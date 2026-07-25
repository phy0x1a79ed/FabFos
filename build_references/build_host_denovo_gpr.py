"""Per-host annotation lanes -> the de-novo gene->reaction table.

    mamba run -n fabfos python build_references/build_host_denovo_gpr.py \
        [--host e_coli_k12] [--force]

Writes ``data/reference/hosts/<host>/gpr_denovo.parquet`` in the shared schema from
``hosts.py``. Every (host, lane) pair and its reader are declared in ``lanes.yml``; this
module holds the readers and nothing about which files exist.

Readers are per-shape on purpose
--------------------------------
There is no tolerant reader here, and that is the design. A reader pointed at the wrong
shape does not raise -- it produces a frame with the wrong columns, an empty inner join,
and a table that is missing a lane rather than an error. The three shapes that actually
occur in this repo (HMMER detail text vs KOfam CSV; 13-column BLAST6 vs 14 with BSR) each
get their own function, and ``lanes.yml`` says which host gets which.

The projection logic is reused from ``lib/fabfos_evidence.py`` -- the stitle cleaner, the
BLAST6 column names, the best-hit-per-ORF rule and the bridge joins -- so this table and
the deployed evidence tables describe the same evidence. What is added here is the
per-host reader dispatch and the host/unit/feature columns the GPR schema carries.

``raw_score`` is the lane's own number on the lane's own scale and is NOT comparable
across channels. That is the whole reason the downstream weighting normalises within
``(orf, channel)`` before combining.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src" / "metasmith_libraries" / "resources" / "lib"))
import hosts  # noqa: E402
import refs  # noqa: E402
import fabfos_evidence as fe  # noqa: E402

LANES_YML = HERE / "lanes.yml"
BUILD_ID = "denovo_4lane"


def _say(msg=""):
    print(msg, flush=True)


def load_config() -> dict:
    return yaml.safe_load(LANES_YML.read_text())


# =====================================================================
# readers -- one per shape, dispatched by lanes.yml
# =====================================================================

def kofam_csv(path: Path, bridge: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """``contig,orf,ko,hmm_threshold,score,evalue,description``.

    Straight to ``fabfos_evidence.read_kofam``, which rejoins ``contig`` and ``orf``
    into the ORF id. That rejoin is load-bearing for K-12, whose writer split
    ``NP_414542.1`` on its underscore into ``contig=NP, orf=414542.1``.
    """
    df = fe.read_kofam(path, "", bridge)
    return df.drop(columns=["source"])


def kofam_hmmer_detail(path: Path, bridge: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """KOfamScan ``--format detail`` text.

        ``* ECDH10B_0001        K08278   23.70   38.5   1.6e-10 thr operon leader peptide``

    Only ``*``-marked rows are at or above the family's own adaptive threshold, which is
    the same predicate ``read_kofam`` applies to the CSV as ``score >= hmm_threshold``.
    Rows without the mark are reported hits below threshold and must not be kept.
    """
    rows = []
    with open(path) as fh:
        for ln in fh:
            if not ln.startswith("*"):
                continue
            parts = ln[1:].strip().split(None, 5)
            if len(parts) < 5:
                continue
            orf, ko, thr, score, evalue = parts[:5]
            desc = parts[5].strip() if len(parts) > 5 else ""
            rows.append((orf, ko, float(score), float(thr), desc))
    df = pd.DataFrame(rows, columns=["orf", "ko", "raw_score", "thr",
                                     "intermediate_name"])
    df = df[df["raw_score"] >= df["thr"]].drop(columns=["thr"])
    df = df.merge(bridge, on="ko", how="inner")
    df["channel"] = "kofam"
    df["intermediate_id"] = df["ko"]
    df["projection_via"] = "kegg.reaction"
    return df[[c for c in fe.SCHEMA_COLS if c != "source"]]


def clean_tsv(path: Path, bridge: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """CLEAN: ``Query ID``, ``Predicted EC number``, ``clean_score``.

    Not ``read_dl_ec``: that reader expects an EZpred parquet with a ``head_kind``
    column, and no host here has one. Same projection, different input shape.

    The K-12 file carries a known writer defect on 350 of its lines: a protein
    description leaked into the EC field and got quoted, so the quote spans into the
    next line and leaves the real prediction's ORF id reading ``"NP_414593.1``. Left
    alone that is not a parse error but an IDENTITY SPLIT -- the deployed table ships
    3,879 gene ids for 3,792 real ORFs because of it. So the file is read with quoting
    disabled, the stray quotes are stripped off the ORF id, the description rows fall
    out on the EC pattern, and the count of repaired rows is returned for reporting
    rather than swallowed.
    """
    import csv as _csv

    floor = float(spec.get("score_floor", 0.0))
    df = pd.read_csv(path, sep="\t", quoting=_csv.QUOTE_NONE, dtype=str)
    df = df.rename(columns={"Query ID": "orf", "Predicted EC number": "ec",
                            "clean_score": "raw_score"})
    quoted = df["orf"].str.contains('"', na=False)
    df["orf"] = df["orf"].str.strip('"')
    df["raw_score"] = pd.to_numeric(df["raw_score"], errors="coerce")
    is_ec = df["ec"].astype(str).str.match(r"^\d+\.\d+\.\d+\.\d+$", na=False)
    clean_tsv.last_repair = dict(
        quote_marked_rows=int(quoted.sum()),
        non_ec_rows_dropped=int((~is_ec).sum()),
        orfs_before_strip=int(df.loc[is_ec, "orf"].nunique()),
    )
    df = df[is_ec & df["raw_score"].notna() & (df["raw_score"] >= floor)]
    df = df.merge(bridge, on="ec", how="inner")
    df["channel"] = "clean_ec"
    df["intermediate_id"] = df["ec"]
    df["intermediate_name"] = ""
    df["projection_via"] = "ec"
    return df[[c for c in fe.SCHEMA_COLS if c != "source"]]


clean_tsv.last_repair = {}


def _blast6(path: Path, bridge: pd.DataFrame, *, has_bsr: bool) -> pd.DataFrame:
    names = list(fe._BLAST6_BSR_COLS)
    if not has_bsr:
        names = names[:-1]
    df = pd.read_csv(path, sep="\t", header=None, names=names, dtype=str)
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
    score = "bsr" if has_bsr else "bitscore"
    if has_bsr:
        df["bsr"] = pd.to_numeric(df["bsr"], errors="coerce")
    # best hit per ORF, exactly as the deployed reader does
    df = (df.sort_values(["qseqid", "evalue", "bitscore"], ascending=[True, True, False])
            .drop_duplicates(subset=["qseqid"], keep="first"))
    df["uniprot_accession"] = df["sseqid"].str.replace(r"^UniRef50_", "", regex=True)
    df["intermediate_name"] = df["stitle"].apply(fe._clean_stitle)
    out = df[["qseqid", "uniprot_accession", "intermediate_name", score]].rename(
        columns={"qseqid": "orf", score: "raw_score"})
    out["orf"] = out["orf"].str.replace(r"-(\d+)$", r"_\1", regex=True)
    j = out.merge(bridge[["uniprot_accession", "dr_source", "mnxr"]],
                  on="uniprot_accession", how="inner")
    j["channel"] = "uniref50_dr"
    j["intermediate_id"] = j["uniprot_accession"]
    j["projection_via"] = j["dr_source"]
    return j[[c for c in fe.SCHEMA_COLS if c != "source"]]


def blast6_bsr(path: Path, bridge: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """14-column DIAMOND BLAST6 + stitle + BSR. ``raw_score`` is the BSR."""
    return _blast6(path, bridge, has_bsr=True)


def blast6_no_bsr(path: Path, bridge: pd.DataFrame, spec: dict) -> pd.DataFrame:
    """13-column DIAMOND BLAST6 + stitle, no BSR column.

    ``raw_score`` carries the raw bitscore instead, and ``lanes.yml`` declares the
    difference rather than papering over it. It does not change any weight: this reader
    keeps one hit per ORF, so each (orf, channel) group holds a single nomination and
    the belief-conservation weight is 1.0 regardless of the number. Manufacturing a
    pseudo-BSR here would put two different quantities in one column under one name.
    """
    return _blast6(path, bridge, has_bsr=False)


def embed_transfer(path: Path, bridge, spec: dict) -> pd.DataFrame:
    raise NotImplementedError  # never reached; the lane is gated on its bridge


READERS = {
    "kofam_csv": kofam_csv,
    "kofam_hmmer_detail": kofam_hmmer_detail,
    "clean_tsv": clean_tsv,
    "blast6_bsr": blast6_bsr,
    "blast6_no_bsr": blast6_no_bsr,
    "embed_transfer": embed_transfer,
}


# =====================================================================
# bridges
# =====================================================================

def load_bridges(cfg: dict) -> tuple[dict, dict]:
    """``{name: frame}`` for what is available, ``{name: reason}`` for what is not."""
    avail, missing = {}, {}
    for name, spec in cfg["bridges"].items():
        if spec.get("available") is False or not spec.get("path"):
            missing[name] = " ".join((spec.get("reason") or "not available").split())
            continue
        p = refs.REPO / spec["path"]
        if not p.exists():
            missing[name] = f"declared at {spec['path']} but not on disk"
            continue
        t0 = time.perf_counter()
        if name == "ko_to_mnxr":
            avail[name] = fe.load_ko_to_mnxr(p)
        elif name == "ec_to_mnxr":
            avail[name] = pd.read_csv(p, sep="\t")[["ec", "mnxr"]].drop_duplicates()
        elif name == "uniprot_to_mnxr":
            avail[name] = pd.read_parquet(
                p, columns=["uniprot_accession", "dr_source", "mnxr"])
        else:
            raise SystemExit(f"no loader for bridge {name!r}")
        _say(f"  bridge {name:18s} {len(avail[name]):>11,} rows "
             f"[{time.perf_counter()-t0:.1f} s]")
    return avail, missing


# =====================================================================
# build
# =====================================================================

def build_host(host: str, cfg: dict, bridges: dict, missing: dict, *,
               force: bool) -> dict:
    hspec = cfg["hosts"][host]
    out = hosts.out_dir(host) / "gpr_denovo.parquet"
    if out.exists() and not force:
        raise SystemExit(f"{out} exists; pass --force")

    root = hosts.OUT_HOSTS / host
    frames, built, skipped, repairs = [], [], [], {}
    for lane, lspec in hspec["lanes"].items():
        gspec = dict(cfg["lanes"][lane])
        gspec.update(lspec)
        need = gspec.get("bridge")
        if need in missing:
            skipped.append((lane, f"bridge {need!r}: {missing[need]}"))
            continue
        path = root / lspec["file"]
        if not path.exists():
            skipped.append((lane, f"declared file {lspec['file']} is not on disk"))
            continue
        t0 = time.perf_counter()
        df = READERS[lspec["reader"]](path, bridges[need], gspec)
        if df.empty:
            # An empty lane after a successful read is the failure mode the per-shape
            # readers exist to prevent, so it is loud rather than a zero-row concat.
            raise SystemExit(
                f"{host}/{lane}: reader {lspec['reader']!r} returned 0 rows from "
                f"{lspec['file']}. That is the wrong-shape signature -- check the "
                f"reader assignment in lanes.yml before relaxing anything.")
        df["lane"] = lane
        frames.append(df)
        built.append(lane)
        note = ""
        if lspec["reader"] == "clean_tsv" and clean_tsv.last_repair.get("quote_marked_rows"):
            r = clean_tsv.last_repair
            repairs[lane] = dict(r)
            note = (f"  [repaired {r['quote_marked_rows']} quote-split ids, "
                    f"dropped {r['non_ec_rows_dropped']} non-EC rows]")
        _say(f"    {lane:15s} {len(df):>8,} rows  {df['orf'].nunique():>5,} ORFs  "
             f"{df['mnxr'].nunique():>6,} MNXR  [{(time.perf_counter()-t0):.1f} s]{note}")
    for lane, why in skipped:
        _say(f"    {lane:15s} SKIPPED -- {why}")

    if not frames:
        raise SystemExit(f"{host}: no lane produced rows; nothing to write")

    ev = pd.concat(frames, ignore_index=True)
    universe = _atom_universe()
    df = pd.DataFrame({
        "build_id": f"{BUILD_ID}_{'+'.join(sorted(built))}",
        "host": host,
        "unit_id": hspec["unit_id"],
        "feature_id": ev["orf"].astype(str),
        "feature_kind": "orf",
        "feature_name": None,
        "mnxr": ev["mnxr"].astype(str),
        "channel": ev["channel"].astype(str),
        "evidence_id": ev["intermediate_id"].astype(str),
        "evidence_name": ev["intermediate_name"].fillna("").astype(str),
        "raw_score": pd.to_numeric(ev["raw_score"], errors="coerce").astype(np.float32),
        "projection_via": ev["projection_via"].astype(str),
        "in_atom_universe": ev["mnxr"].isin(universe),
        # a de-novo call has no boolean rule; the column exists so both lanes share one
        # schema, and leaving it null is what says "this lane asserts no rule"
        "gpr_rule": None,
    })[list(hosts.GPR_COLS)]
    df = df.drop_duplicates(["feature_id", "channel", "evidence_id", "mnxr"])
    df = df.sort_values(["channel", "feature_id", "mnxr"],
                        kind="mergesort").reset_index(drop=True)

    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False, compression="zstd")

    per_lane = {c: int(n) for c, n in df["channel"].value_counts().items()}
    summary = dict(
        host=host, build_id=df["build_id"].iloc[0],
        lanes_built=built, lanes_skipped={l: w for l, w in skipped},
        rows=len(df), orfs=int(df["feature_id"].nunique()),
        mnxr=int(df["mnxr"].nunique()),
        rows_per_lane=per_lane,
        input_repairs=repairs,
        in_universe_rows=int(df["in_atom_universe"].sum()),
        aam_gap_mnxr=int(df.loc[~df["in_atom_universe"], "mnxr"].nunique()),
        bytes=out.stat().st_size,
    )
    _say(f"      {len(df):>8,} rows  {summary['orfs']:>5,} ORFs  "
         f"{summary['mnxr']:>6,} MNXR  "
         f"{summary['aam_gap_mnxr']:>6,} MNXR with no atom pairs  "
         f"{out.stat().st_size/1e6:.1f} MB")
    return summary


def _atom_universe() -> set:
    refs.assert_same_bake()
    V = refs.load_vocab()
    codes = pd.read_parquet(refs.ATOM_PAIRS, columns=["rxn"])["rxn"].to_numpy()
    return set(V.symbols("rxn")[np.unique(codes)].tolist())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host", action="append", choices=sorted(hosts.HOSTS))
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    cfg = load_config()
    targets = a.host or sorted(cfg["hosts"])
    _say("== de-novo GPR tables ==")
    bridges, missing = load_bridges(cfg)
    for name, why in missing.items():
        _say(f"  bridge {name:18s} UNAVAILABLE -- {why}")

    summaries = []
    for h in targets:
        _say(f"  {h}:")
        summaries.append(build_host(h, cfg, bridges, missing, force=a.force))

    (hosts.OUT_HOSTS / "gpr_denovo_summary.json").write_text(
        json.dumps(summaries, indent=2, sort_keys=True) + "\n")
    _say(f"\nwrote {(hosts.OUT_HOSTS / 'gpr_denovo_summary.json').relative_to(refs.REPO)}")

    lane_sets = {s["host"]: sorted(s["lanes_built"]) for s in summaries}
    if len(set(map(tuple, lane_sets.values()))) > 1:
        _say("\nlane sets differ between hosts -- that is part of each table's identity, "
             "not a defect to average away:")
        for h, ls in lane_sets.items():
            _say(f"  {h:16s} {', '.join(ls)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
