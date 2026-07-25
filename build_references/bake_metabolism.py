"""Bake the frozen MetaNetX atom mapping and directionality into compiled tables.

    mamba run -n fabfos python build_references/bake_metabolism.py [--force]
    mamba run -n fabfos python build_references/bake_metabolism.py --selftest

Reads two immutable, DVC-pinned inputs from ``data/raw/mnxref-4_5/`` and writes three
parquets into ``data/reference/metabolism/``. See ``refs.py`` for the encoding and the
bake-identity contract.

What this is NOT
----------------
It is a re-encoding, not a re-derivation. Nothing here recomputes an atom mapping or a
free energy; every number that comes out is the number that went in, at float32. The
selftest exists to prove exactly that, row by row over all 2,455,235 rows -- because the
one failure mode that matters (a too-narrow rank field merging two distinct atoms onto
one node) RAISES the network's conductance and so reads as an improvement rather than a
bug.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).resolve().parent))
import refs  # noqa: E402


def _say(msg=""):
    print(msg, flush=True)


# =====================================================================
# read + census the inputs
# =====================================================================

def read_sources(*, verify: bool = True):
    for p in (refs.SRC_ATOM_PAIRS, refs.SRC_DIRECTION):
        if not p.exists():
            raise SystemExit(
                f"missing input {p}\n"
                f"  the mnxref bundle is DVC-pinned but may not be checked out here:\n"
                f"  mamba run -n dvc dvc checkout data/raw/mnxref-4_5.dvc")

    t0 = time.perf_counter()
    ap_sha = refs.sha256_file(refs.SRC_ATOM_PAIRS)
    dir_sha = refs.sha256_file(refs.SRC_DIRECTION)
    pairs = pd.read_parquet(refs.SRC_ATOM_PAIRS)
    direction = pd.read_parquet(refs.SRC_DIRECTION,
                                columns=["mnxr", "ratio", "dir_tier"])
    _say(f"  read sources                       {(time.perf_counter()-t0)*1000:7.0f} ms")
    _say(f"    atom_pairs  {len(pairs):>9,} rows  sha {ap_sha[:16]}")
    _say(f"    direction   {len(direction):>9,} rows  sha {dir_sha[:16]}")

    manifest = {}
    if refs.SRC_MANIFEST.exists():
        manifest = json.loads(refs.SRC_MANIFEST.read_text())
        # The manifest counts REACTIONS, which is the direction table's row count; it says
        # nothing about atom-pair rows. Check what it actually claims, not what would be
        # convenient.
        if verify and manifest.get("n_reactions") not in (None, len(direction)):
            raise SystemExit(
                f"MANIFEST.json says n_reactions={manifest['n_reactions']} but "
                f"direction.parquet has {len(direction)} rows -- the bundle is "
                f"inconsistent; refusing to bake against it")

    if verify:
        missing = set(pairs["mnxr"].unique()) - set(direction["mnxr"].unique())
        if missing:
            raise SystemExit(
                f"{len(missing)} reactions have atom pairs but no direction row "
                f"(e.g. {sorted(missing)[:3]}); the union coding assumes the direction "
                f"table is the superset")

    return pairs, direction, ap_sha, dir_sha, manifest


# =====================================================================
# vocabulary
# =====================================================================

def build_vocabulary(pairs: pd.DataFrame, direction: pd.DataFrame):
    elements = sorted(pairs["element"].unique())
    unexpected = set(elements) - set(refs.ELEMENT_ORDER)
    if unexpected:
        raise SystemExit(
            f"element(s) {sorted(unexpected)} are in the source but not in "
            f"refs.ELEMENT_ORDER; adding one renumbers the element codes, so it is a "
            f"deliberate edit to refs.py, not something to infer here")

    mets = sorted(set(pairs["substrate"].unique()) | set(pairs["product"].unique()))

    # The reaction vocabulary is the UNION of both tables. 20,175 reactions carry a
    # direction row and no atom pairs; coding against the atom-pair set alone would drop
    # them, and every dropped reaction then falls back to ratio 1.0 -- fully reversible,
    # i.e. MORE conductance than the evidence supports. Silent, and in the flattering
    # direction.
    rxns = sorted(set(pairs["mnxr"].unique()) | set(direction["mnxr"].unique()))

    # The 'EMPTY' sentinel row of the direction table is kept rather than dropped, so the
    # bake preserves the source census exactly and the row-count check is a real check
    # rather than one with a hand-waved exception. Consumers that filter it by symbol
    # (ecspr_build does) still work, because the symbol survives in the vocabulary.
    spaces = {
        "element": list(refs.ELEMENT_ORDER),
        "met": mets,
        "rxn": rxns,
        "method": sorted(pairs["method"].dropna().unique()),
        "source": sorted(pairs["source"].dropna().unique()),
    }
    vocab = refs.build_vocab(spaces)
    _say(f"    vocabulary  " + "  ".join(
        f"{k} {len(v):,}" for k, v in spaces.items()))
    return vocab, refs.Vocab(vocab)


# =====================================================================
# encode
# =====================================================================

def encode_pairs(pairs: pd.DataFrame, V: refs.Vocab, rank_bits: int) -> pd.DataFrame:
    """Int-code every id column, keeping the (metabolite, rank) split.

    Fusing the two into one node code here would cost 7 MB: it makes every value
    distinct and parquet loses the dictionary encoding it gets on a 34k-symbol
    metabolite column. Consumers pack at load in 4 ms -- see ``refs.pack_pairs``.
    """
    met = V.codes("met")
    enc = pd.DataFrame({
        "element": pd.Series(V.encode("element", pairs["element"].to_numpy()),
                             dtype=np.uint8),
        "rxn": pd.Series(V.encode("rxn", pairs["mnxr"].to_numpy()), dtype=np.uint32),
        "tail_met": np.fromiter((met[s] for s in pairs["substrate"]), np.uint16,
                                count=len(pairs)),
        "tail_rank": pairs["sub_idx"].to_numpy().astype(np.uint16),
        "head_met": np.fromiter((met[s] for s in pairs["product"]), np.uint16,
                                count=len(pairs)),
        "head_rank": pairs["prod_idx"].to_numpy().astype(np.uint16),
        "pair_w": pairs["pair_w"].to_numpy().astype(np.float32),
        "method": pd.Series(V.encode("method", pairs["method"].to_numpy()), dtype=np.uint8),
        "source": pd.Series(V.encode("source", pairs["source"].to_numpy()), dtype=np.uint8),
        "confidence": pairs["confidence"].to_numpy().astype(np.float32),
    })
    return enc.sort_values(list(refs.ATOM_PAIRS_SORT),
                           kind="mergesort").reset_index(drop=True)


def encode_direction(direction: pd.DataFrame, V: refs.Vocab) -> pd.DataFrame:
    """``ratio`` stays float64. The consumer's flip test is a threshold at exactly 1.0,
    and four reactions sit within float32 epsilon of it -- narrowing reorients those four
    edges, which is a topology change no re-encoding is allowed to make."""
    enc = pd.DataFrame({
        "rxn": pd.Series(V.encode("rxn", direction["mnxr"].to_numpy()), dtype=np.uint32),
        "ratio": direction["ratio"].to_numpy().astype(np.float64),
        "dir_tier": direction["dir_tier"].to_numpy().astype(np.uint8),
    })
    return enc.sort_values("rxn", kind="mergesort").reset_index(drop=True)


def write_pairs(enc: pd.DataFrame, identity: dict) -> Path:
    """One row group per element, so the element filter prunes row groups rather than
    reading everything and masking. Written slice by slice because parquet's
    ``row_group_size`` is a uniform cap and would not land on element boundaries."""
    table = pa.Table.from_pandas(enc, preserve_index=False)
    md = dict(table.schema.metadata or {})
    md[refs.BAKE_KEY] = json.dumps(identity, sort_keys=True).encode()
    schema = table.schema.with_metadata(md)
    refs.ATOM_PAIRS.parent.mkdir(parents=True, exist_ok=True)
    codes = enc["element"].to_numpy()
    with pq.ParquetWriter(refs.ATOM_PAIRS, schema, compression="zstd",
                          compression_level=refs.ZSTD_LEVEL) as w:
        for code in np.unique(codes):
            w.write_table(table.filter(pa.array(codes == code)))
    return refs.ATOM_PAIRS


# =====================================================================
# bake
# =====================================================================

def bake(force: bool = False) -> dict:
    existing = [p for p in refs.BAKED if p.exists()]
    if existing and not force:
        raise SystemExit(
            "refusing to overwrite an existing bake: "
            + ", ".join(p.name for p in existing) + "\n  pass --force")

    _say("== bake metabolism ==")
    pairs, direction, ap_sha, dir_sha, manifest = read_sources()

    t0 = time.perf_counter()
    vocab, V = build_vocabulary(pairs, direction)
    max_rank = int(max(pairs["sub_idx"].max(), pairs["prod_idx"].max()))
    met_bits, rank_bits = refs.bit_widths(V.size("met"), max_rank)
    _say(f"    widths      met {met_bits} + rank {rank_bits} = node {met_bits+rank_bits} bits, "
         f"edge {2*(met_bits+rank_bits)} of {refs.NODE_KEY_BUDGET}  "
         f"(max atom rank {max_rank}, headroom to {(1<<rank_bits)-1})")
    _say(f"  vocabulary                         {(time.perf_counter()-t0)*1000:7.0f} ms")

    t0 = time.perf_counter()
    enc_pairs = encode_pairs(pairs, V, rank_bits)
    enc_dir = encode_direction(direction, V)
    _say(f"  encode + sort                      {(time.perf_counter()-t0)*1000:7.0f} ms")

    identity = {
        "bake_version": refs.BAKE_VERSION,
        "mnxref_version": manifest.get("mnxref_version"),
        "mnxref_date": manifest.get("mnxref_date"),
        "src_atom_pairs_sha256": ap_sha,
        "src_atom_pairs_rows": int(len(pairs)),
        "src_direction_sha256": dir_sha,
        "src_direction_rows": int(len(direction)),
        "vocab_sha256": refs.vocab_sha256(vocab),
        "n_element": V.size("element"),
        "n_met": V.size("met"),
        "n_rxn": V.size("rxn"),
        "n_method": V.size("method"),
        "n_source": V.size("source"),
        "met_bits": met_bits,
        "rank_bits": rank_bits,
        "max_atom_rank": max_rank,
        "element_order": list(refs.ELEMENT_ORDER),
        "orientation": refs.ORIENTATION,
        "atom_pairs_rows": int(len(enc_pairs)),
        "direction_rows": int(len(enc_dir)),
    }

    t0 = time.perf_counter()
    refs.write_with_identity(pa.Table.from_pandas(vocab, preserve_index=False),
                             refs.VOCAB, identity)
    write_pairs(enc_pairs, identity)
    refs.write_with_identity(pa.Table.from_pandas(enc_dir, preserve_index=False),
                             refs.DIRECTION, identity)
    _say(f"  write                              {(time.perf_counter()-t0)*1000:7.0f} ms")

    total = sum(p.stat().st_size for p in refs.BAKED)
    src = refs.SRC_ATOM_PAIRS.resolve().stat().st_size + \
        refs.SRC_DIRECTION.resolve().stat().st_size
    for p in refs.BAKED:
        _say(f"    {p.name:22s} {p.stat().st_size/1e6:7.2f} MB")
    _say(f"    {'total':22s} {total/1e6:7.2f} MB   (sources {src/1e6:.2f} MB)")
    _say(f"    row groups in atom_pairs: {pq.ParquetFile(refs.ATOM_PAIRS).num_row_groups} "
         f"(one per element)")
    _say(f"    bake {identity['vocab_sha256'][:16]}")
    return identity


# =====================================================================
# selftest
# =====================================================================

def selftest() -> int:
    failures: list[str] = []

    def check(label, ok, detail=""):
        print(f"[{'PASS' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""),
              flush=True)
        if not ok:
            failures.append(label)

    _say("== selftest: the bake is a re-encoding and nothing more ==")
    ident = refs.assert_same_bake()
    check("bake identity agrees across all three files", True,
          f"bake {ident['vocab_sha256'][:16]}")

    pairs = pd.read_parquet(refs.SRC_ATOM_PAIRS)
    direction = pd.read_parquet(refs.SRC_DIRECTION, columns=["mnxr", "ratio", "dir_tier"])
    V = refs.load_vocab()
    enc = refs.load_atom_pairs()
    enc_dir = refs.load_direction()

    check("source atom_pairs unchanged since the bake",
          refs.sha256_file(refs.SRC_ATOM_PAIRS) == ident["src_atom_pairs_sha256"])
    check("source direction unchanged since the bake",
          refs.sha256_file(refs.SRC_DIRECTION) == ident["src_direction_sha256"])
    check("vocabulary hash matches the identity block",
          refs.vocab_sha256(V.df) == ident["vocab_sha256"])

    check("atom_pairs row census preserved", len(enc) == len(pairs),
          f"{len(enc):,} vs {len(pairs):,}")
    check("direction row census preserved", len(enc_dir) == len(direction),
          f"{len(enc_dir):,} vs {len(direction):,}")
    # The 1,280 pair_w == 0 rows are kept. graph_from_pairs drops them itself; removing
    # them at bake time would silently move the census pin.
    check("zero-weight pair rows preserved",
          int((enc["pair_w"].to_numpy() == 0).sum()) == int((pairs["pair_w"] == 0).sum()),
          f"{int((enc['pair_w'].to_numpy()==0).sum()):,} rows")

    # -- the round trip, over every row, not a sample --------------------------------
    t0 = time.perf_counter()
    rank_bits = ident["rank_bits"]
    src = pairs.sort_values(["element", "mnxr", "substrate", "product",
                             "sub_idx", "prod_idx"], kind="mergesort").reset_index(drop=True)
    met_sym = V.symbols("met")
    rxn_sym = V.symbols("rxn")
    el_sym = V.symbols("element")
    # Decode through the PACKED form, not straight off the split columns, so the round
    # trip actually exercises pack_node/unpack_node -- the step that silently merges two
    # atoms onto one node if rank_bits is a bit too narrow.
    tail_node, head_node = refs.pack_pairs(enc, rank_bits)
    tm, ti = refs.unpack_node(tail_node, rank_bits)
    hm, hi = refs.unpack_node(head_node, rank_bits)
    dec = pd.DataFrame({
        "element": el_sym[enc["element"].to_numpy()],
        "mnxr": rxn_sym[enc["rxn"].to_numpy()],
        "substrate": met_sym[tm], "product": met_sym[hm],
        "sub_idx": ti, "prod_idx": hi,
        "pair_w": enc["pair_w"].to_numpy().astype(float),
    }).sort_values(["element", "mnxr", "substrate", "product", "sub_idx", "prod_idx"],
                   kind="mergesort").reset_index(drop=True)

    for col in ("element", "mnxr", "substrate", "product"):
        check(f"round trip exact: {col}",
              bool((dec[col].to_numpy() == src[col].to_numpy()).all()))
    for col in ("sub_idx", "prod_idx"):
        check(f"round trip exact: {col}",
              bool((dec[col].to_numpy() == src[col].to_numpy()).all()))

    # float32 is a deliberate narrowing: ~6e-8 relative, which is four orders below the
    # spread of the evidence weights it multiplies. Asserted with an explicit tolerance
    # rather than equality, and the max observed error is printed so a future widening of
    # the source's dynamic range shows up here instead of in a conductance.
    dw = np.abs(dec["pair_w"].to_numpy() - src["pair_w"].to_numpy())
    rel = dw / np.maximum(np.abs(src["pair_w"].to_numpy()), 1.0)
    check("round trip within float32: pair_w", float(rel.max()) <= 1e-6,
          f"max relative error {float(rel.max()):.2e}")
    _say(f"       ({len(enc):,} rows verified in {(time.perf_counter()-t0)*1000:.0f} ms)")

    # Per-source-atom weight sums are PRESERVED, not conserved. The source table does not
    # in fact carry sum == 1.0 per (reaction, element, substrate atom): 1,360,927 groups
    # do, but 763,983 do not -- 0.5 for a lone-member correspondence, up to 4.0 observed.
    # So the invariant a re-encoding can honestly assert is that it changed none of them.
    key = ["element", "mnxr", "substrate", "sub_idx"]
    gs = src.groupby(key, observed=True)["pair_w"].sum()
    gd = dec.groupby(key, observed=True)["pair_w"].sum()
    gd = gd.reindex(gs.index)
    check("per-source-atom weight sums preserved",
          bool((np.abs(gd.to_numpy() - gs.to_numpy()) <= 1e-6 * np.maximum(gs.to_numpy(), 1.0)).all()),
          f"{len(gs):,} source atoms, max drift "
          f"{float(np.abs(gd.to_numpy()-gs.to_numpy()).max()):.2e}")

    # -- direction --------------------------------------------------------------------
    dsrc = direction.set_index("mnxr")
    ddec = pd.DataFrame({
        "mnxr": rxn_sym[enc_dir["rxn"].to_numpy()],
        "ratio": enc_dir["ratio"].to_numpy().astype(float),
        "dir_tier": enc_dir["dir_tier"].to_numpy(),
    }).set_index("mnxr").reindex(dsrc.index)
    check("round trip exact: dir_tier",
          bool((ddec["dir_tier"].to_numpy() == dsrc["dir_tier"].to_numpy()).all()))
    check("round trip exact: ratio",
          bool((ddec["ratio"].to_numpy() == dsrc["ratio"].to_numpy()).all()),
          f"float64 preserved; ratios span {dsrc['ratio'].min():.1e} to "
          f"{dsrc['ratio'].max():.1e}")
    # ratio > 1 is not a rare edge case, and it FLIPS the edge rather than amplifying it.
    # Pinning the count here is what caught the float32 ratio: four reactions sit within
    # float32 epsilon of 1.0, so narrowing the column silently reoriented four edges.
    n_flip_src = int((dsrc["ratio"] > 1.0).sum())
    check("edge-flipping ratio count preserved",
          int((ddec["ratio"] > 1.0).sum()) == n_flip_src,
          f"{n_flip_src:,} of {len(dsrc):,} ratios exceed 1")

    # -- widths -----------------------------------------------------------------------
    max_rank = int(max(pairs["sub_idx"].max(), pairs["prod_idx"].max()))
    check("recorded max atom rank matches the source", ident["max_atom_rank"] == max_rank,
          f"{max_rank}")
    check("rank field is wide enough for the observed maximum",
          max_rank < (1 << ident["rank_bits"]),
          f"{max_rank} < {1 << ident['rank_bits']}")
    check("metabolite field is wide enough for the vocabulary",
          ident["n_met"] <= (1 << ident["met_bits"]),
          f"{ident['n_met']:,} <= {1 << ident['met_bits']:,}")
    check("edge key fits in int64",
          2 * (ident["met_bits"] + ident["rank_bits"]) <= refs.NODE_KEY_BUDGET)

    _say()
    if failures:
        _say(f"{len(failures)} FAILED: " + "; ".join(failures))
        return 1
    _say("all checks passed")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="overwrite an existing bake")
    ap.add_argument("--selftest", action="store_true",
                    help="verify an existing bake against its sources; does not write")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    bake(force=a.force)
    _say()
    return selftest()


if __name__ == "__main__":
    raise SystemExit(main())
