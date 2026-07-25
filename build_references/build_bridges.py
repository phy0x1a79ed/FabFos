"""Build the intermediate -> MNXR bridge tables the annotation lanes project through.

    mamba run -n fabfos python build_references/build_bridges.py [--force]

Writes into ``data/reference/functional_annotation/bridges/``::

    ec_to_mnxr.tsv          EC number  -> MNXR   (MetaNetX reac_prop classifs)
    ko_to_mnxr.tsv          KEGG KO    -> MNXR   (KOfam ko_list EC -> the above)
    uniprot_to_mnxr.parquet UniProt    -> MNXR   (Rhea DR -> MetaNetX reac_xref)

Two of the three already have builders in ``lib/fabfos_evidence.py`` and are called
here rather than reimplemented. ``ko_to_mnxr`` did not -- ``fabfos_evidence`` calls it
"a reused reference table (staged input), not built here" -- so its producer is ported
here from the method that actually built the deployed lane,
``scadc .../betweenness/05_build_ko_to_mnxr.py``:

    KO -> KEGG reaction -> MNXR

where ``KO -> KEGG reaction`` comes from the ``REACTION`` block of each KO's KEGG
flat-file, and ``KEGG reaction -> MNXR`` from the ``kegg.reaction:`` rows of MetaNetX
``reac_xref``. Same regexes, same join.

Routing KO through its ``[EC:...]`` tag instead would have been easy and wrong: measured
on these three hosts it takes the kofam lane from 646 to 8,548 reactions and makes 91%
of them reactions ``clean_ec`` already reaches, so the two lanes stop being independent
evidence -- which is the entire reason for having four of them.

The KO flat-files are read from a repo-local table, ``ko_to_kegg_r.tsv``, built once
from the scadc KEGG request cache (read-only) plus a live REST fetch for the remainder,
and thereafter the only input. That keeps the sibling project's cache a recorded
provenance source rather than a runtime dependency, which is the point of baking these
tables at all.

The fourth thing the lanes want -- ``reference_label_pool``, the labelled reference
embeddings the kNN transfer lane votes against -- is NOT buildable from anything here.
It needs a reference proteome embedded with the same model as the query. See
``build_host_denovo_gpr.py``, which names it as missing rather than skipping the lane
silently.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sqlite3
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src" / "metasmith_libraries" / "resources" / "lib"))
import refs  # noqa: E402
import fabfos_evidence as fe  # noqa: E402

REF = refs.REPO / "data" / "reference"
REAC_PROP = REF / "metabolism" / "metanetx" / "4.5" / "reac_prop.tsv"
REAC_XREF = REF / "metabolism" / "metanetx" / "4.5" / "reac_xref.tsv"
KO_LIST = REF / "functional_annotation" / "kofam" / "ko_list"
RHEA_SWISS = REF / "functional_annotation" / "rhea" / "rhea2uniprot.tsv"
RHEA_TREMBL = REF / "functional_annotation" / "rhea" / "rhea2uniprot_trembl.tsv.gz"

OUT = REF / "functional_annotation" / "bridges"
EC_TO_MNXR = OUT / "ec_to_mnxr.tsv"
KO_TO_MNXR = OUT / "ko_to_mnxr.tsv"
KO_TO_KEGG_R = OUT / "ko_to_kegg_r.tsv"
UNIPROT_TO_MNXR = OUT / "uniprot_to_mnxr.parquet"
BRIDGES = (EC_TO_MNXR, KO_TO_MNXR, UNIPROT_TO_MNXR)

# One-time provenance source for the KO flat-files, read-only. Not a runtime dependency:
# once ko_to_kegg_r.tsv exists in this repo it is the only input.
SCADC_KEGG_CACHE = Path("/home/tony/agentic_workspace/data/scadc/cache/kegg_requests.db")
KEGG_URL = "https://rest.kegg.jp/get/{ko}"
KEGG_PAUSE = 0.4  # polite; KEGG asks for no more than ~3 requests/second

# A complete EC. Partial ECs ("3.6.3.-", "1.2.1") are kept in ec_to_mnxr verbatim
# because reac_prop writes them that way.
_EC4 = re.compile(r"^\d+\.\d+\.\d+\.\d+$")

# The REACTION block of a KEGG K-entry flat-file, and the R-numbers on it. Both lifted
# verbatim from 05_build_ko_to_mnxr.py so this bridge cannot drift from the deployed one.
_REACTION_BLOCK = re.compile(r"^REACTION\s+(.+?)(?=\n[A-Z]+\s|\Z)", re.S | re.M)
_R_ID = re.compile(r"\bR\d{5}\b")


def _say(msg=""):
    print(msg, flush=True)


def build_ec(force: bool) -> pd.DataFrame:
    if EC_TO_MNXR.exists() and not force:
        return pd.read_csv(EC_TO_MNXR, sep="\t")
    t0 = time.perf_counter()
    df = fe.load_ec_to_mnxr(REAC_PROP)
    OUT.mkdir(parents=True, exist_ok=True)
    df.sort_values(["ec", "mnxr"]).to_csv(EC_TO_MNXR, sep="\t", index=False)
    _say(f"  ec_to_mnxr       {len(df):>9,} rows  {df['ec'].nunique():>7,} ECs  "
         f"{df['mnxr'].nunique():>7,} MNXR  [{(time.perf_counter()-t0)*1000:.0f} ms]")
    n4 = int(df["ec"].str.match(_EC4).sum())
    _say(f"                   {n4:>9,} rows carry a complete 4-level EC "
         f"({len(df)-n4:,} partial)")
    return df


def _reactions_of(text) -> set:
    """R-numbers on a KO flat-file's REACTION block. Verbatim from the deployed builder."""
    if not text:
        return set()
    if isinstance(text, bytes):
        text = (gzip.decompress(text) if text[:2] == b"\x1f\x8b" else text).decode()
    if text.startswith('"') and text.endswith('"'):
        try:
            text = json.loads(text)
        except Exception:
            pass
    m = _REACTION_BLOCK.search(text)
    return set(_R_ID.findall(m.group(1))) if m else set()


def declared_host_kos() -> set:
    """Every KO the hosts declared in lanes.yml actually call.

    Read from the declared kofam files, using the same two shapes their readers use, so
    the bridge's scope is derived from the hosts rather than guessed.
    """
    import yaml
    cfg = yaml.safe_load((HERE / "lanes.yml").read_text())
    out: set[str] = set()
    for host, hspec in cfg["hosts"].items():
        lane = hspec["lanes"].get("kofam")
        if not lane:
            continue
        p = refs.REPO / "data" / "reference" / "hosts" / host / lane["file"]
        if not p.exists():
            continue
        if lane["reader"] == "kofam_csv":
            out |= set(pd.read_csv(p)["ko"].dropna())
        elif lane["reader"] == "kofam_hmmer_detail":
            with open(p) as fh:
                for ln in fh:
                    if ln.startswith("*"):
                        parts = ln[1:].strip().split(None, 2)
                        if len(parts) >= 2:
                            out.add(parts[1])
        else:
            raise SystemExit(f"no KO extractor for reader {lane['reader']!r}")
    return out


def build_ko_to_kegg_r(force: bool, *, allow_fetch: bool = True) -> pd.DataFrame:
    """``ko -> kegg_r``, from the scadc KEGG cache plus a live fetch for the remainder.

    Written into this repo so it is built once and then owned here.

    SCOPE: every KO already in the cache, plus every KO the declared hosts call. Not all
    28,277 KOs in ``ko_list`` -- that would be ~24,000 live REST requests for KOs nothing
    here annotates. The covered set is recorded in the sidecar, and adding a host whose
    KOs fall outside it re-runs this step rather than quietly projecting through a
    narrower bridge. The deployed builder was scoped the same way (to the fosmid KOs).
    """
    if KO_TO_KEGG_R.exists() and not force:
        df = pd.read_csv(KO_TO_KEGG_R, sep="\t")
        _say(f"  ko_to_kegg_r     {len(df):>9,} rows  {df['ko'].nunique():>7,} KOs  "
             f"(reusing the repo-local table)")
        return df

    rows, n_cached, n_live, n_no_reaction, n_failed = [], 0, 0, 0, 0

    cached: dict[str, object] = {}
    if SCADC_KEGG_CACHE.exists():
        con = sqlite3.connect(f"file:{SCADC_KEGG_CACHE}?mode=ro", uri=True)
        try:
            for url, data in con.execute(
                    "SELECT id, data FROM json_cache WHERE id LIKE ?",
                    ("%rest.kegg.jp/get/K%",)):
                cached[url.rsplit("/", 1)[1]] = data
        finally:
            con.close()
        _say(f"  KEGG flat-files from {SCADC_KEGG_CACHE.name}: {len(cached):,} cached")
    else:
        _say(f"  {SCADC_KEGG_CACHE} not readable; every KO must be fetched live")

    host_kos = declared_host_kos()
    known = set(pd.read_csv(KO_LIST, sep="\t", dtype=str)["knum"].dropna())
    wanted = sorted((set(cached) | host_kos) & known)
    _say(f"  scope: {len(host_kos):,} KOs called by the declared hosts, "
         f"{len(host_kos - set(cached)):,} of them uncached; "
         f"resolving {len(wanted):,} in total")

    todo = [k for k in wanted if k not in cached]
    if todo and not allow_fetch:
        raise SystemExit(
            f"{len(todo):,} of {len(wanted):,} KOs are not in the cache and "
            f"--no-fetch was given. A bridge built from the cached subset would be "
            f"silently narrower than the deployed one; refusing.")
    if todo:
        _say(f"  fetching {len(todo):,} uncached KOs from KEGG REST "
             f"(~{len(todo)*KEGG_PAUSE/60:.1f} min at {KEGG_PAUSE}s spacing)")

    for i, ko in enumerate(wanted, 1):
        text = cached.get(ko)
        if text is not None:
            n_cached += 1
        else:
            try:
                with urllib.request.urlopen(KEGG_URL.format(ko=ko), timeout=20) as r:
                    text = r.read().decode()
                n_live += 1
                time.sleep(KEGG_PAUSE)
            except Exception as e:
                n_failed += 1
                _say(f"    [warn] {ko}: {type(e).__name__}: {e}")
                continue
        rs = _reactions_of(text)
        if not rs:
            n_no_reaction += 1
            continue
        rows.extend((ko, r) for r in sorted(rs))
        if n_live and i % 500 == 0:
            _say(f"    {i:,}/{len(wanted):,} (cached {n_cached:,}, live {n_live:,})")

    df = pd.DataFrame(rows, columns=["ko", "kegg_r"]).drop_duplicates()
    OUT.mkdir(parents=True, exist_ok=True)
    df.sort_values(["ko", "kegg_r"]).to_csv(KO_TO_KEGG_R, sep="\t", index=False)
    _say(f"  ko_to_kegg_r     {len(df):>9,} rows  {df['ko'].nunique():>7,} KOs  "
         f"(cached {n_cached:,}, fetched {n_live:,}, "
         f"{n_no_reaction:,} KOs carry no REACTION block, {n_failed:,} failed)")
    return df


def _kegg_r_to_mnxr() -> pd.DataFrame:
    """``kegg.reaction:R#####`` rows of reac_xref. Both prefixes, as the deployed
    builder accepts both."""
    rows = []
    with open(REAC_XREF) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            src, mnxr = parts[0], parts[1]
            if not mnxr.startswith("MNXR") or mnxr == "EMPTY":
                continue
            if src.startswith("kegg.reaction:") or src.startswith("keggR:"):
                rows.append((src.split(":", 1)[1], mnxr))
    return pd.DataFrame(rows, columns=["kegg_r", "mnxr"]).drop_duplicates()


def build_ko(force: bool, *, allow_fetch: bool = True) -> pd.DataFrame:
    if KO_TO_MNXR.exists() and not force:
        return pd.read_csv(KO_TO_MNXR, sep="\t")
    t0 = time.perf_counter()
    ko_r = build_ko_to_kegg_r(force, allow_fetch=allow_fetch)
    r_mnxr = _kegg_r_to_mnxr()
    df = ko_r.merge(r_mnxr, on="kegg_r", how="inner")[["ko", "mnxr"]].drop_duplicates()
    OUT.mkdir(parents=True, exist_ok=True)
    df.sort_values(["ko", "mnxr"]).to_csv(KO_TO_MNXR, sep="\t", index=False)
    n_unmapped_r = int((~ko_r["kegg_r"].isin(set(r_mnxr["kegg_r"]))).sum())
    _say(f"  ko_to_mnxr       {len(df):>9,} rows  {df['ko'].nunique():>7,} KOs  "
         f"{df['mnxr'].nunique():>7,} MNXR  [{time.perf_counter()-t0:.1f} s]")
    _say(f"                   reac_xref carries {r_mnxr['kegg_r'].nunique():,} KEGG "
         f"reactions; {n_unmapped_r:,} of the KO links point at one it does not have")
    _say(f"                   MNXR per KO: median "
         f"{df.groupby('ko')['mnxr'].nunique().median():.0f}, "
         f"max {df.groupby('ko')['mnxr'].nunique().max()}")
    return df


def build_uniprot(force: bool) -> Path:
    if UNIPROT_TO_MNXR.exists() and not force:
        return UNIPROT_TO_MNXR
    t0 = time.perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    # TrEMBL is ~10x SwissProt and is what gives the lane its reach on non-model ORFs;
    # it is passed when present and the builder degrades to reviewed-only when not.
    fe.build_uniprot_bridge(REAC_XREF, RHEA_SWISS,
                            RHEA_TREMBL if RHEA_TREMBL.exists() else None,
                            UNIPROT_TO_MNXR)
    d = pd.read_parquet(UNIPROT_TO_MNXR, columns=["uniprot_accession", "mnxr",
                                                  "evidence_quality"])
    _say(f"  uniprot_to_mnxr  {len(d):>9,} rows  "
         f"{d['uniprot_accession'].nunique():>7,} accessions  "
         f"{d['mnxr'].nunique():>7,} MNXR  [{time.perf_counter()-t0:.1f} s]")
    _say(f"                   {d['evidence_quality'].value_counts().to_dict()}")
    return UNIPROT_TO_MNXR


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--no-fetch", action="store_true",
                    help="fail rather than fetch uncached KO flat-files from KEGG REST")
    a = ap.parse_args()

    missing = [p for p in (REAC_PROP, REAC_XREF, KO_LIST, RHEA_SWISS) if not p.exists()]
    if missing:
        raise SystemExit("missing bridge inputs:\n  "
                         + "\n  ".join(str(p.relative_to(refs.REPO)) for p in missing))

    _say("== intermediate -> MNXR bridges ==")
    build_ec(a.force)
    build_ko(a.force, allow_fetch=not a.no_fetch)
    build_uniprot(a.force)
    _say("")
    for p in BRIDGES:
        _say(f"  {p.relative_to(refs.REPO)}  {p.stat().st_size/1e6:.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
