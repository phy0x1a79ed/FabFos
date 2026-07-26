"""Read the four benchmark cohorts out of the acquired literature products.

Ported from the sibling project's ``benchmark_v4/_v4common.py`` cohort loaders
(``load_gof_pairs``, ``load_lof_genes``, ``load_eydallin``) and step 20's heterologous
half. One departure, and it is the thing to know about this module:

    THE COHORT SOURCES ARE THE EXTRACTED TABLES, NOT THE PAPERS.

``_v4common`` read ``gof_observations.tsv`` and ``lof_observations.tsv`` -- tabular
extractions of the LASER records and the Keio supplement that a human made once. This
module reads the same extracted tables, looked for inside each literature product. It
does NOT re-derive them from the LASER checkout's own tables or from Keio's .xls, because
extraction from a paper supplement is human judgement, and a re-extraction that disagreed
with the deployed one would silently move the benchmark rather than fail.

Under the tier rule that makes these extractions ``curated/`` artifacts rather than
``raw/`` ones, and ``build_references/REFERENCES.md`` flags exactly that as an open
decision ("`benchmark/v3/ground_truth/gof_*.tsv` ... are either a `curated/` artifact or
they get re-extracted by a transform"). This module implements the first reading and
REFUSES BY NAME when a table is absent, so the choice stays visible instead of shipping a
cohort that is quietly smaller.

Every loader returns rows in one schema, so the caller never branches on cohort:

    cohort  condition_id  gene_label  uniprot  action  host_hint  source_organism

``uniprot`` is set only for heterologous proteins. A native gene is named by symbol and
its edges are the host's own, resolved against the per-host GPR tables at network
construction time -- which is why a null there is a routing fact, not a missing value.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pandas as pd

COLUMNS = ("cohort", "condition_id", "gene_label", "uniprot", "action",
           "host_hint", "source_organism")

# The extracted table each product is expected to carry. Named here so a refusal can say
# which paper's extraction is missing rather than which path did not exist.
EXPECTED = {
    "laser":    ("gof_observations.tsv",
                 "LASER engineering records, extracted to one row per observation with a "
                 "genes_json column"),
    "keio":     ("lof_observations.tsv",
                 "Keio single-gene knockouts (Baba 2006), extracted to one row per "
                 "knockout with its b-number"),
    "eydallin": ("eydallin2010_hits.tsv",
                 "Eydallin 2010 screen hits, extracted from the PDF supplement"),
    "het":      ("heterologous_uniprot.tsv",
                 "the heterologous / metagenomic screening cohort's protein list, one row "
                 "per protein with its UniProt accession"),
}

# LASER writes a literal "None" where a field is absent (v3 learned this the hard way);
# treat those as native rather than as an unnamed organism.
_NULLISH = {"", "none", "n/a", "na", "null", "not specified", "-"}


def nullish(s: str) -> bool:
    return (s or "").strip().lower() in _NULLISH


def is_native(source: str) -> bool:
    """A LASER GeneSource of '' / 'None' / anything *coli* means the host genome."""
    return nullish(source) or "coli" in (source or "").lower()


def _find(product: Path, key: str) -> Path:
    filename, what = EXPECTED[key]
    product = Path(product)
    hits = [product / filename] + sorted(product.rglob(filename))
    for h in hits:
        if h.exists():
            return h
    raise SystemExit(
        f"cohort source missing: {filename} is not under {product}.\n"
        f"  It is {what}.\n"
        f"  This is a NAMED refusal. Building without it would ship a conditions table "
        f"that looks complete and is missing an arm.")


def load_gof(laser: Path) -> pd.DataFrame:
    """Distinct (gene, source) pairs from every LASER observation's genes_json.

    Splits into the two GOF arms on the SAME test the deployed loader uses: a gene whose
    source is blank / 'None' / *coli* is the host's own (``gof_native``); anything else is
    heterologous (``gof_het``). Getting that test wrong moves proteins between the arm
    that resolves against a host GPR and the arm that projects through the bridge.
    """
    rows = []
    with open(_find(laser, "laser")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gj = (row.get("genes_json") or "").strip()
            if not gj:
                continue
            for g in json.loads(gj):
                gene = (g.get("gene") or "").strip()
                src = (g.get("source") or "").strip()
                if nullish(gene):
                    continue
                native = is_native(src)
                rows.append(dict(
                    cohort="gof_native" if native else "gof_het",
                    condition_id=row.get("obs_id") or "",
                    gene_label=gene,
                    uniprot=None,          # filled from the het table below, where known
                    action=",".join(sorted(g.get("actions") or [])),
                    host_hint=row.get("host", ""),
                    source_organism="" if native else src,
                ))
    return pd.DataFrame(rows, columns=list(COLUMNS)).drop_duplicates()


def load_lof(keio: Path) -> pd.DataFrame:
    rows = []
    with open(_find(keio, "keio")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = (row.get("gene_set") or "").strip()
            if not gene:
                continue
            rows.append(dict(
                cohort="lof",
                condition_id=row.get("obs_id") or f"KEIO:{gene}",
                gene_label=gene,
                uniprot=None,
                # A Keio strain is a deletion; the action is the cohort's meaning, not a
                # per-row field, but it is written out so every arm reads the same way.
                action="del",
                host_hint=row.get("host", ""),
                source_organism="",
            ))
    return pd.DataFrame(rows, columns=list(COLUMNS)).drop_duplicates()


def load_eydallin(eydallin: Path) -> pd.DataFrame:
    rows = []
    with open(_find(eydallin, "eydallin")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            gene = (r.get("gene_norm") or r.get("gene") or "").strip()
            if not gene:
                continue
            rows.append(dict(
                cohort="eydallin",
                condition_id=f"EYDALLIN:{gene}",
                gene_label=gene,
                uniprot=None,
                action=(r.get("phenotype") or "").strip(),
                host_hint="",
                source_organism="",
            ))
    return pd.DataFrame(rows, columns=list(COLUMNS)).drop_duplicates()


def load_het(het: Path) -> pd.DataFrame:
    """The heterologous protein list -- the only arm that carries UniProt accessions."""
    rows = []
    with open(_find(het, "het")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            acc = (r.get("uniprot") or "").strip()
            gene = (r.get("gene") or "").strip()
            if not acc:
                # Named in the screen but never resolved to a sequence. Kept as a row with
                # no accession so the attrition is visible downstream; dropping it here
                # would make the arm look smaller than the screen actually was.
                acc = None
            rows.append(dict(
                cohort="gof_het",
                condition_id=(r.get("obs_id") or f"HET:{gene or acc}"),
                gene_label=gene or (acc or ""),
                uniprot=acc,
                action="add",
                host_hint="",
                source_organism=(r.get("source") or "").strip(),
            ))
    return pd.DataFrame(rows, columns=list(COLUMNS)).drop_duplicates()


def load_all_cohorts(*, laser, keio, eydallin, het) -> pd.DataFrame:
    """All four arms in one frame, with LASER's heterologous entries given accessions.

    LASER names a heterologous gene by label and source organism; the het table is what
    carries the accession. They are joined on (gene, source) -- the same key the deployed
    step 20 uses -- rather than on the label alone, because two screens can nominate the
    same gene symbol from different organisms and merging them would attribute one
    organism's protein to the other's condition.
    """
    gof = load_gof(Path(laser))
    lof = load_lof(Path(keio))
    eyd = load_eydallin(Path(eydallin))
    hets = load_het(Path(het))

    key = hets.dropna(subset=["uniprot"]).drop_duplicates(["gene_label", "source_organism"])
    lut = dict(zip(zip(key["gene_label"], key["source_organism"]), key["uniprot"]))
    mask = gof["cohort"] == "gof_het"
    gof.loc[mask, "uniprot"] = [
        lut.get((g, s)) for g, s in zip(gof.loc[mask, "gene_label"],
                                        gof.loc[mask, "source_organism"])
    ]

    # The het table's own rows are kept alongside LASER's: the screening cohort includes
    # proteins LASER never recorded, and they are what the arm is for.
    out = pd.concat([gof, lof, eyd, hets], ignore_index=True)
    return out.drop_duplicates(subset=["cohort", "condition_id", "gene_label",
                                       "uniprot"]).reset_index(drop=True)
