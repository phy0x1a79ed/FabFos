"""Curated genome-scale model -> per-host gene->reaction table.

    mamba run -n scadc-metabolic-model python build_references/build_host_gem_gpr.py \
        [--host e_coli_k12] [--force]

Writes ``data/reference/hosts/<host>/gpr_gem.parquet`` in the shared schema from
``hosts.py``. cobra is needed for the crosswalk (it normalises the annotation
dictionaries ``crosswalk_gem`` reads), which is why this step runs in a different
environment from the bake.

Three decisions worth knowing about
-----------------------------------
**The GPR is not evaluated here.** A reference table has no perturbation, so evaluating
the boolean rule at bake time would bake in one condition. The table carries the full
crosswalk plus each reaction's ``gpr_rule``, and the knockout evaluation happens in the
consumer -- which must remember that ``cobra.GPR.eval`` takes KNOCKOUTS, not the active
set. Handing it the active set inverts the question, and a previous run of that logic
reported more reactions live after a knockout than before.

**``raw_score`` is uniform 1.0.** A curated model asserts that a reaction is PRESENT,
not how much evidence there is for it. Weighting it by anything would be inventing a
quantity. The evidence-weighted lane is ``gpr_denovo.parquet``.

**Reactions with no gene get a row too.** Exchanges, diffusion and spontaneous chemistry
have no gene to attribute, but they are live in every condition -- ``keep_ruleless`` in
the reference builder exists precisely because dropping them "would make every gene set
look like a starvation". They are emitted with ``feature_kind='ruleless'`` and a null
``feature_id``, so the table carries the model's whole reaction set and a consumer never
has to reopen the model to find them.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src" / "metasmith_libraries" / "resources" / "lib"))
import hosts  # noqa: E402
import refs  # noqa: E402

REAC_XREF = refs.REPO / "data" / "reference" / "metabolism" / "metanetx" / "4.5" / \
    "reac_xref.tsv"
CHANNEL = "gem_gpr"


def _say(msg=""):
    print(msg, flush=True)


def atom_universe_from_bake() -> set:
    """Every MNXR carrying at least one atom transfer, read off the compiled table.

    "In universe" means the reaction has atom-pair coverage, so an edge can exist for it.
    Reading it from the bake rather than the source table is the point of the bake, and
    the identity assertion below is what stops a stale vocab decoding it to the wrong
    reactions.
    """
    refs.assert_same_bake()
    V = refs.load_vocab()
    codes = pd.read_parquet(refs.ATOM_PAIRS, columns=["rxn"])["rxn"].to_numpy()
    return set(V.symbols("rxn")[np.unique(codes)].tolist())


def build_host(host: str, reac_xref, universe: set, *, force: bool) -> dict:
    from ecspr_build import crosswalk_gem, load_model

    spec = hosts.HOSTS[host]
    out = hosts.out_dir(host) / "gpr_gem.parquet"
    if out.exists() and not force:
        raise SystemExit(f"{out} exists; pass --force")

    model_path = hosts.gem_json(host)
    if not model_path.exists():
        raise SystemExit(f"missing model {model_path}")

    t0 = time.perf_counter()
    model = load_model(model_path)
    xw, stats = crosswalk_gem(model, reac_xref, universe)
    _say(f"  {host}: {spec['gem_id']}  {stats['n_reactions']:,} reactions -> "
         f"{stats['n_resolved']:,} resolved ({xw['mnxr'].nunique():,} distinct MNXR), "
         f"{stats['n_unresolved']:,} unresolved, "
         f"{stats['n_aam_gap']:,} outside the atom universe "
         f"[{(time.perf_counter()-t0)*1000:.0f} ms]")

    by_mnxr = dict(zip(xw["rxn_id"], zip(xw["mnxr"], xw["source"], xw["in_universe"])))
    build_id = f"gem_{spec['gem_id']}"
    rows = []
    n_ruleless = 0
    for r in model.reactions:
        hit = by_mnxr.get(r.id)
        if hit is None:
            # no candidate MNXR at all: it can never carry an edge, and emitting a row
            # with a null reaction would put an unjoinable key in the table
            continue
        mnxr, via, in_universe = hit
        rule = (r.gene_reaction_rule or "").strip()
        genes = list(r.genes)
        common = dict(
            build_id=build_id, host=host, unit_id=spec["gem_id"],
            mnxr=mnxr, channel=CHANNEL,
            evidence_id=r.id, evidence_name=r.name or None,
            raw_score=1.0, projection_via=via,
            in_atom_universe=bool(in_universe), gpr_rule=rule or None,
        )
        if not genes:
            n_ruleless += 1
            rows.append(dict(common, feature_id=None, feature_kind="ruleless",
                             feature_name=None))
            continue
        for g in genes:
            rows.append(dict(common, feature_id=g.id, feature_kind="gem_gene",
                             feature_name=g.name or None))

    df = pd.DataFrame(rows, columns=list(hosts.GPR_COLS))
    df["raw_score"] = df["raw_score"].astype(np.float32)
    df["in_atom_universe"] = df["in_atom_universe"].astype(bool)
    df = df.sort_values(["feature_kind", "feature_id", "mnxr", "evidence_id"],
                        kind="mergesort", na_position="last").reset_index(drop=True)

    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False, compression="zstd")

    gene_rows = df[df["feature_kind"] == "gem_gene"]
    summary = dict(
        host=host, model=spec["gem_id"], build_id=build_id,
        rows=len(df), gene_rows=len(gene_rows), ruleless_rows=n_ruleless,
        genes=int(gene_rows["feature_id"].nunique()),
        mnxr=int(df["mnxr"].nunique()),
        distinct_gene_mnxr=int(gene_rows.drop_duplicates(["feature_id", "mnxr"]).shape[0]),
        in_universe_rows=int(df["in_atom_universe"].sum()),
        aam_gap_mnxr=int(df.loc[~df["in_atom_universe"], "mnxr"].nunique()),
        model_genes=len(model.genes),
        bytes=out.stat().st_size,
        **{f"crosswalk_{k}": v for k, v in stats.items()},
    )
    _say(f"      {len(df):>7,} rows  {summary['genes']:>5,} genes  "
         f"{summary['mnxr']:>6,} MNXR  "
         f"{summary['distinct_gene_mnxr']:>6,} distinct (gene, MNXR)  "
         f"{n_ruleless:>5,} ruleless  "
         f"{summary['aam_gap_mnxr']:>5,} MNXR with no atom pairs  "
         f"{out.stat().st_size/1e3:.0f} KB")
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host", action="append", choices=sorted(hosts.HOSTS),
                    help="repeatable; default is every host")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    targets = a.host or sorted(hosts.HOSTS)

    _say("== GEM GPR tables ==")
    if not REAC_XREF.exists():
        raise SystemExit(f"missing {REAC_XREF}")
    t0 = time.perf_counter()
    universe = atom_universe_from_bake()
    _say(f"  atom universe: {len(universe):,} reactions with atom-pair coverage "
         f"[{(time.perf_counter()-t0)*1000:.0f} ms]")

    summaries = [build_host(h, REAC_XREF, universe, force=a.force) for h in targets]

    # EPI300 and DH10B share a model, so their GEM tables must agree everywhere except
    # the host tag. Asserting it here means a future divergence surfaces as a failure
    # rather than as two silently different tables.
    shared = [s for s in summaries if s["model"] == hosts.HOSTS["e_coli_epi300"]["gem_id"]]
    if len({s["host"] for s in shared}) > 1:
        cols = [c for c in hosts.GPR_COLS if c != "host"]
        frames = {s["host"]: pd.read_parquet(hosts.out_dir(s["host"]) / "gpr_gem.parquet")
                  for s in shared}
        ref_host, ref_df = next(iter(frames.items()))
        for h, d in frames.items():
            if h == ref_host:
                continue
            same = d[cols].equals(ref_df[cols])
            _say(f"  {h} vs {ref_host}: identical apart from the host tag -- {same}")
            if not same:
                raise SystemExit(
                    f"{h} and {ref_host} share model "
                    f"{hosts.HOSTS[h]['gem_id']} but produced different tables")

    (hosts.OUT_HOSTS / "gpr_gem_summary.json").write_text(
        json.dumps(summaries, indent=2, sort_keys=True) + "\n")
    _say(f"\nwrote {(hosts.OUT_HOSTS / 'gpr_gem_summary.json').relative_to(refs.REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
