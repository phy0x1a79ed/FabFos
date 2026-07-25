"""Host registry, genome readers and the shared GPR schema.

One place that knows which organisms we build reference tables for, where each one's
genome, proteome and curated model live, and what a row of a GPR table looks like.

The GPR schema
--------------
Long format, one row per (feature x channel x MNXR). Both lanes -- the curated-model
lane and the annotation lane -- write the SAME columns, so they can be concatenated and
compared rather than conflated. What differs between them is the values, and the
``channel`` column says which lane a row came from.

Fourteen columns, and every one past the obvious seven exists to stop one specific
silent failure:

``build_id``          what produced this table; carries the model or lane-set identity
``host``              the organism. Two hosts can share a build_id (see EPI300)
``unit_id``           the physical or logical carrier of the feature -- a model, a
                      contig, a fosmid. Currently guessed downstream by
                      ``orf.rsplit("_", 1)[0]``, which is a guess this column retires
``feature_id``        the gene or ORF
``feature_kind``      ``gem_gene`` or ``orf``. A GEM gene and an ORF are NOT the same
                      identifier; the incumbent table ships a known 3,879-vs-3,792
                      mismatch from conflating them
``feature_name``      gene symbol where one exists. Without it every consumer reopens
                      the model or the GFF to label a row
``mnxr``              the reaction, as a STRING and never a vocabulary code -- 23.6% of
                      the incumbent's reactions have no atom pairs, and a code cannot
                      represent a reaction outside the atom universe
``channel``           the evidence lane
``evidence_id``       the intermediate the mapping went through: a model reaction id, a
                      KO, an EC, a UniProt accession
``evidence_name``     that intermediate's label
``raw_score``         the lane's own score, on the lane's own scale. NOT comparable
                      across channels; that is what the belief-conservation weighting
                      downstream is for
``projection_via``    how evidence_id reached mnxr (``bigg``, ``kegg``, ``embedded``,
                      a bridge name)
``in_atom_universe``  whether this MNXR has any atom-pair coverage. Rows are KEPT when
                      false; an inner join here would report perfect coverage
``gpr_rule``          the model's boolean rule, so a knockout does not require
                      reopening the model
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
RAW_HOSTS = REPO / "data" / "raw" / "hosts"
OUT_HOSTS = REPO / "data" / "reference" / "hosts"

GPR_COLS = (
    "build_id", "host", "unit_id", "feature_id", "feature_kind", "feature_name",
    "mnxr", "channel", "evidence_id", "evidence_name", "raw_score",
    "projection_via", "in_atom_universe", "gpr_rule",
)

# EPI300 carries no model of its own and needs none: measured against the genome
# announcement (doi 10.1128/mra.01124-25) and both genomes, the EPI300<->DH10B delta is
# three genes lost (yfdH, rhsA, ypjCB) and three gained (araC, dfrA1, trfA), and NONE of
# the six is a gene in iECDH10B_1368 -- so the GEM-side edit list is empty. Writing a
# copied iEPI300_derived.json would mint a distinct hash asserting a difference that
# does not exist. check_epi300_identity.py re-measures this on every build.
HOSTS: dict[str, dict] = {
    "e_coli_dh10b": dict(
        label="E. coli K-12 DH10B",
        accession="NC_010473.1",
        gem_host="e_coli_dh10b",
        gem_file="GEM/iECDH10B_1368.json",
        gem_id="iECDH10B_1368",
    ),
    "e_coli_k12": dict(
        label="E. coli K-12 MG1655",
        accession="NC_000913.3",
        gem_host="e_coli_k12",
        gem_file="GEM/iML1515.json",
        gem_id="iML1515",
    ),
    "e_coli_epi300": dict(
        label="E. coli EPI300",
        accession="CP189566.1",
        gem_host="e_coli_dh10b",
        gem_file="GEM/iECDH10B_1368.json",
        gem_id="iECDH10B_1368",
        gem_note=("no model of its own; the EPI300<->DH10B edit list is empty at the "
                  "model level -- see gem_identity_report.md"),
    ),
}

# P. putida DOT-T1E has a genome and a proteome but no annotation lanes and no curated
# model, so neither table can be built for it until the annotation pipeline has run.
OUT_OF_SCOPE = {"p_putida_dot_t1e": "no annotation lanes and no curated model"}


def host_dir(host: str) -> Path:
    return RAW_HOSTS / host


def genome_gbk(host: str) -> Path:
    return RAW_HOSTS / host / "genome" / f"{HOSTS[host]['accession']}.gbk"


def proteome_faa(host: str) -> Path:
    return RAW_HOSTS / host / "genome" / f"{HOSTS[host]['accession']}.faa"


def gem_json(host: str) -> Path:
    h = HOSTS[host]
    return RAW_HOSTS / h["gem_host"] / h["gem_file"]


def out_dir(host: str) -> Path:
    return OUT_HOSTS / host


# =====================================================================
# genome readers
# =====================================================================

def _open(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


# A GenBank feature block starts at column 5 and its key at column 6. Splitting on that
# is what keeps a /note containing the word "gene" from being read as a new feature.
_FEATURE_SPLIT = re.compile(r"\n     (?=\S)")


def parse_genbank_features(path, kinds=("gene",)) -> pd.DataFrame:
    """``locus_tag, old_locus_tag, gene, kind, pseudo`` for the requested feature kinds.

    Deliberately a regex pass rather than Biopython: this needs four qualifiers off a
    12 MB flatfile, and a full SeqIO parse costs seconds and a dependency for nothing.
    """
    with _open(path) as f:
        txt = f.read()
    rows = []
    for block in _FEATURE_SPLIT.split(txt):
        head = block.lstrip().split(None, 1)
        if not head or head[0] not in kinds:
            continue
        lt = re.search(r'/locus_tag="([^"]+)"', block)
        ol = re.search(r'/old_locus_tag="([^"]+)"', block)
        gn = re.search(r'/gene="([^"]+)"', block)
        rows.append(dict(
            locus_tag=lt.group(1) if lt else None,
            old_locus_tag=ol.group(1) if ol else None,
            gene=gn.group(1) if gn else None,
            kind=head[0],
            pseudo="/pseudo" in block,
        ))
    return pd.DataFrame(rows, columns=["locus_tag", "old_locus_tag", "gene", "kind",
                                       "pseudo"])


_FAA_TAG = re.compile(r"\[locus_tag=([^\]]+)\]")
_FAA_GENE = re.compile(r"\[gene=([^\]]+)\]")
_FAA_PROT = re.compile(r"\[protein=([^\]]+)\]")


def parse_faa(path) -> pd.DataFrame:
    """``locus_tag, gene, product, seq`` from an NCBI-style protein FASTA."""
    heads, seqs, cur = [], [], []
    with _open(path) as f:
        for ln in f:
            if ln.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
                heads.append(ln[1:].strip())
            else:
                cur.append(ln.strip())
    if cur:
        seqs.append("".join(cur))
    if len(heads) != len(seqs):
        raise ValueError(f"{path}: {len(heads)} headers but {len(seqs)} sequences")
    return pd.DataFrame({
        "locus_tag": [(_FAA_TAG.search(h).group(1) if _FAA_TAG.search(h) else h)
                      for h in heads],
        "gene": [(_FAA_GENE.search(h).group(1) if _FAA_GENE.search(h) else None)
                 for h in heads],
        "product": [(_FAA_PROT.search(h).group(1) if _FAA_PROT.search(h) else None)
                    for h in heads],
        "seq": seqs,
    })


def locus_tag_aliases(gbk_path) -> dict:
    """``old_locus_tag -> locus_tag``, for models keyed on the pre-RefSeq tags.

    iECDH10B's gene ids are ``ECDH10B_1296``-style old tags while the current annotation
    uses ``ECDH10B_RS06655``. Without this map every join on the model's gene ids returns
    nothing, and an empty join is indistinguishable from a strain with no metabolism.
    """
    df = parse_genbank_features(gbk_path)
    df = df[df["old_locus_tag"].notna() & df["locus_tag"].notna()]
    return dict(zip(df["old_locus_tag"], df["locus_tag"]))


# =====================================================================
# curated model
# =====================================================================

def load_gem_json(path) -> dict:
    """The cobra JSON as plain data.

    Reading it as JSON rather than through cobra is deliberate for the checks that only
    need gene ids and boolean rules: it keeps them runnable in the same environment as
    the rest of the bake. ``build_host_gem_gpr.py`` does load it through cobra, because
    the crosswalk needs the annotation dictionaries cobra normalises.
    """
    with _open(path) as f:
        return json.load(f)


_RULE_TOKEN = re.compile(r"[A-Za-z_][A-Za-z0-9_.\-]*")


def rule_genes(rule: str | None, known: set) -> list:
    """Gene ids mentioned in a ``gene_reaction_rule``, restricted to ids the model has.

    ``and``/``or`` fall out because they are not model gene ids -- no operator list to
    keep in sync.
    """
    if not rule:
        return []
    return [t for t in dict.fromkeys(_RULE_TOKEN.findall(rule)) if t in known]


def gem_gene_index(gem: dict, gbk_path=None) -> pd.DataFrame:
    """``gem_gene_id -> (name, locus_tag)``, resolving the model's MIXED id space.

    iECDH10B keys 1,294 of its 1,327 genes on ``ECDH10B_*`` old locus tags and the other
    33 on bare gene symbols (``acpS``, ``adiA``, ``cyaA``, ``frdB``, ``pabB``, ...), and
    33 entries carry an empty ``name``. A join on either key alone silently drops rows,
    so this resolves by old tag first, then by symbol, and reports what it could not
    place rather than dropping it.
    """
    genes = pd.DataFrame([{"gem_gene_id": g["id"], "name": g.get("name") or None}
                          for g in gem["genes"]])
    genes["locus_tag"] = None
    if gbk_path is not None:
        feat = parse_genbank_features(gbk_path)
        by_old = dict(zip(feat["old_locus_tag"].dropna(),
                          feat.loc[feat["old_locus_tag"].notna(), "locus_tag"]))
        by_tag = set(feat["locus_tag"].dropna())
        by_name = {}
        for g, lt in zip(feat["gene"], feat["locus_tag"]):
            if g is not None and g not in by_name:
                by_name[g] = lt
        def resolve(row):
            gid = row["gem_gene_id"]
            if gid in by_old:
                return by_old[gid]
            if gid in by_tag:
                return gid
            for key in (gid, row["name"]):
                if key and key in by_name:
                    return by_name[key]
            return None
        genes["locus_tag"] = genes.apply(resolve, axis=1)
    return genes
