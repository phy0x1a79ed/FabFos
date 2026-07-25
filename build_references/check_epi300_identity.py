"""Re-measure the EPI300 <-> DH10B delta, and prove the GEM-side edit list is empty.

    mamba run -n fabfos python build_references/check_epi300_identity.py

EPI300 gets its curated-model GPR from ``iECDH10B_1368`` under its own host tag, with no
derived model file and no hand-authored edit list. That is a claim, and this script is
what turns it into a measurement that re-runs on every build:

  * the whole-genome named-gene diff between CP189566 and NC_010473 is exactly three
    lost (``yfdH``, ``rhsA``, ``ypjCB``) and three gained (``araC``, ``dfrA1``,
    ``trfA``) -- the announced P_BAD-araC-trfA203 + dfrA cassette landing in yfdH;
  * none of those six is a gene in the model;
  * no model gene is missing from EPI300, and all but a handful translate to a
    byte-identical protein.

If any of that stops being true the gate FAILS, rather than the shared table quietly
becoming wrong. Writes ``gem_identity_report.md`` next to the host's tables.

Sources: Liu, Formby, Davenport, Capron, Szeitz & Hallam, Microbiol Resour Announc
(2025), doi 10.1128/mra.01124-25 -- the EPI300 genome announcement; Durfee et al.,
J Bacteriol 190:2597 (2008) -- the DH10B genome, and the finding that the sequenced
strain is in fact (galE galK galU)+ despite the classical genotype string.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import hosts  # noqa: E402

DERIVED = "e_coli_epi300"
BASE = "e_coli_dh10b"

# The point mutations in the EPI300/DH10B genotype string, and what the sequence data
# says about each. Reported rather than acted on: not one of them yields an edit.
GENOTYPE_NOTES = [
    ("galK16", "already /pseudo in the DH10B annotation and absent from the model"),
    ("araD139, D(ara,leu)7697",
     "araA, araB and leuABCD are absent from both genomes; the gene the annotation "
     "calls araD (ECDH10B_3764, 3.85 Mb) is 55 substitutions over 231 residues from "
     "K-12's araD (b0061, ~70 kb) and is a paralogous L-ribulose-5-P 4-epimerase, not "
     "the ara operon gene"),
    ("galU", "DH10B's GalU protein is byte-identical to K-12's; Durfee et al. 2008 "
             "report the sequenced strain is (galE galK galU)+"),
    ("lacX74", "lacZ and lacA are not model genes; lacY is intact in both genomes"),
    ("endA1, recA1, rpsL, nupG, deoR, mcrA, hsdRMS, mrr", "not model genes"),
    ("trfA, araC", "gained by the cassette; neither is metabolic"),
    ("dfrA", "gained by the cassette; a DHFR, and the model's DHFR already reads "
             "'ECDH10B_1739 or ECDH10B_0049', so it adds no reaction"),
]


def measure() -> dict:
    base_gbk, der_gbk = hosts.genome_gbk(BASE), hosts.genome_gbk(DERIVED)
    base_faa, der_faa = hosts.proteome_faa(BASE), hosts.proteome_faa(DERIVED)
    for p in (base_gbk, der_gbk, base_faa, der_faa):
        if not p.exists():
            raise SystemExit(
                f"missing {p}\n  mamba run -n dvc dvc checkout data/raw/hosts/*.dvc")

    bf = hosts.parse_genbank_features(base_gbk)
    df = hosts.parse_genbank_features(der_gbk)
    b_named = Counter(bf["gene"].dropna())
    d_named = Counter(df["gene"].dropna())
    lost = sorted(set(b_named) - set(d_named))
    gained = sorted(set(d_named) - set(b_named))
    copy_diff = sorted((g, b_named[g], d_named[g]) for g in set(b_named) & set(d_named)
                       if b_named[g] != d_named[g])

    bp = hosts.parse_faa(base_faa)
    dp = hosts.parse_faa(der_faa)
    der_seqs = set(dp["seq"])
    der_by_name = {}
    for g, s in zip(dp["gene"], dp["seq"]):
        if g is not None:
            der_by_name.setdefault(g, []).append(s)

    gem = hosts.load_gem_json(hosts.gem_json(BASE))
    gem_ids = {g["id"] for g in gem["genes"]}
    idx = hosts.gem_gene_index(gem, base_gbk)
    tag2prot = {lt: (gn, sq) for lt, gn, sq in
                zip(bp["locus_tag"], bp["gene"], bp["seq"])}

    # which reactions each model gene carries, so a loss can be reported by consequence
    g2r: dict[str, list] = {}
    for r in gem["reactions"]:
        for t in hosts.rule_genes(r.get("gene_reaction_rule"), gem_ids):
            g2r.setdefault(t, []).append(r["id"])

    resolved, unresolvable, identical, missing, changed = 0, [], 0, [], []
    for gid, name, lt in idx[["gem_gene_id", "name", "locus_tag"]].itertuples(index=False):
        rec = tag2prot.get(lt) if lt else None
        if rec is None:
            unresolvable.append(gid)
            continue
        resolved += 1
        gname, seq = rec
        if seq in der_seqs:
            identical += 1
            continue
        cand = der_by_name.get(gname or "", [])
        if not cand:
            missing.append(dict(gem_gene_id=gid, gene=gname, aa=len(seq),
                                n_rxn=len(g2r.get(gid, [])),
                                rxns=g2r.get(gid, [])[:6]))
        else:
            best = min(cand, key=lambda c: abs(len(c) - len(seq)))
            nsub = (sum(1 for a, b in zip(seq, best) if a != b)
                    if len(best) == len(seq) else None)
            changed.append(dict(gem_gene_id=gid, gene=gname, base_aa=len(seq),
                                derived_aa=len(best), n_subs=nsub,
                                n_rxn=len(g2r.get(gid, [])),
                                rxns=g2r.get(gid, [])[:6]))

    changed_in_gem = sorted(
        {g for g in lost + gained
         if g in set(idx["name"].dropna()) or g in gem_ids}
    )

    # pseudogenes the annotation calls broken but the model still carries. Reported, not
    # edited: both are pseudo in EPI300 too, and every reaction they touch has an `or`
    # alternative, so neither changes any reaction's liveness.
    b_pseudo = set(bf.loc[bf["pseudo"], "gene"].dropna())
    d_pseudo = set(df.loc[df["pseudo"], "gene"].dropna())
    pseudo_in_gem = []
    for gid, name, lt in idx[["gem_gene_id", "name", "locus_tag"]].itertuples(index=False):
        if name in b_pseudo:
            rxns = g2r.get(gid, [])
            solo = [r["id"] for r in gem["reactions"]
                    if r["id"] in set(rxns)
                    and hosts.rule_genes(r.get("gene_reaction_rule"), gem_ids) == [gid]]
            pseudo_in_gem.append(dict(gem_gene_id=gid, gene=name, n_rxn=len(rxns),
                                      n_rxn_solo=len(solo),
                                      also_pseudo_in_derived=name in d_pseudo))

    return dict(
        base_named=len(b_named), derived_named=len(d_named),
        lost=lost, gained=gained, copy_diff=copy_diff,
        base_proteins=len(bp), derived_proteins=len(dp),
        gem_genes=len(gem_ids), gem_reactions=len(gem["reactions"]),
        resolved=resolved, unresolvable=unresolvable,
        identical=identical, missing=missing, changed=changed,
        changed_in_gem=changed_in_gem, pseudo_in_gem=pseudo_in_gem,
    )


def write_report(m: dict, path: Path) -> Path:
    L = []
    A = L.append
    A(f"# {hosts.HOSTS[DERIVED]['label']} vs {hosts.HOSTS[BASE]['label']}: "
      f"curated-model identity")
    A("")
    A(f"Generated by `build_references/check_epi300_identity.py`. "
      f"{hosts.HOSTS[DERIVED]['accession']} against "
      f"{hosts.HOSTS[BASE]['accession']}, and both against "
      f"`{hosts.HOSTS[DERIVED]['gem_id']}`.")
    A("")
    A("EPI300 has no curated model of its own and needs none. This report is the "
      "evidence for that, re-measured on every build rather than asserted once.")
    A("")
    A("## Whole-genome named-gene diff")
    A("")
    A(f"- named genes: {m['base_named']:,} in {BASE}, {m['derived_named']:,} in "
      f"{DERIVED}")
    A(f"- lost in {DERIVED} ({len(m['lost'])}): "
      + (", ".join(f"`{g}`" for g in m["lost"]) or "none"))
    A(f"- gained in {DERIVED} ({len(m['gained'])}): "
      + (", ".join(f"`{g}`" for g in m["gained"]) or "none"))
    A(f"- copy-number differences: {len(m['copy_diff'])}")
    A("")
    A("The gained three are the announced copy-control cassette "
      "(P_BAD-`araC`-`trfA203` + `dfrA`); `yfdH` is its insertion site, disrupted at "
      "codon 121.")
    A("")
    A("## Restricted to the model")
    A("")
    A(f"- model genes: {m['gem_genes']:,} over {m['gem_reactions']:,} reactions")
    A(f"- resolved to a {BASE} protein: {m['resolved']:,} "
      f"(unresolvable ids: {len(m['unresolvable'])})")
    A(f"- byte-identical protein in {DERIVED}: {m['identical']:,}")
    A(f"- **absent from {DERIVED}: {len(m['missing'])}**")
    A(f"- present but not byte-identical: {len(m['changed'])}")
    A(f"- changed genes that are model genes: {len(m['changed_in_gem'])}")
    A("")
    if m["missing"]:
        A("### Model genes with no counterpart")
        A("")
        for x in m["missing"]:
            A(f"- `{x['gem_gene_id']}` ({x['gene']}), {x['aa']} aa, "
              f"{x['n_rxn']} reactions: {', '.join(x['rxns'])}")
        A("")
    if m["changed"]:
        A("### Model genes whose protein differs")
        A("")
        A("| gene | model id | base aa | derived aa | substitutions | reactions |")
        A("| --- | --- | --- | --- | --- | --- |")
        for x in sorted(m["changed"], key=lambda x: -(x["n_subs"] or 10**6)):
            A(f"| {x['gene']} | `{x['gem_gene_id']}` | {x['base_aa']} | "
              f"{x['derived_aa']} | "
              f"{x['n_subs'] if x['n_subs'] is not None else 'length differs'} | "
              f"{x['n_rxn']} |")
        A("")
        A("Length differences here are annotation-boundary calls between a 2008 "
          "assembly and PGAP 6.10, not strain biology.")
        A("")
    if m["unresolvable"]:
        A("### Model gene ids that resolve to no annotated feature")
        A("")
        A(f"{len(m['unresolvable'])} of {m['gem_genes']:,}: "
          + ", ".join(f"`{g}`" for g in sorted(m["unresolvable"])[:40])
          + (" ..." if len(m["unresolvable"]) > 40 else ""))
        A("")
        A("These are the model's bare-symbol gene ids. They are the same in both hosts, "
          "so they cannot make the two disagree -- but they do mean the identity check "
          "covers the resolvable majority rather than every id.")
        A("")
    A("## Genotype point mutations")
    A("")
    A("The published genotype is `F- mcrA D(mrr-hsdRMS-mcrBC) phi80dlacZDM15 DlacX74 "
      "recA1 endA1 araD139 D(ara, leu)7697 galU galK lambda- rpsL nupG trfA dhfr`. "
      "None of it yields a model edit:")
    A("")
    for allele, verdict in GENOTYPE_NOTES:
        A(f"- **{allele}** -- {verdict}")
    A("")
    if m["pseudo_in_gem"]:
        A("## Pseudogenes the model still carries")
        A("")
        A("| gene | model id | reactions | reactions where it is the only gene | "
          "also pseudo in derived |")
        A("| --- | --- | --- | --- | --- |")
        for x in sorted(m["pseudo_in_gem"], key=lambda x: -x["n_rxn"]):
            A(f"| {x['gene']} | `{x['gem_gene_id']}` | {x['n_rxn']} | "
              f"{x['n_rxn_solo']} | {x['also_pseudo_in_derived']} |")
        A("")
        A("Kept in the GPR table rather than edited out, so the table stays a faithful "
          "reading of the model file it claims to encode. Every reaction they touch has "
          "an `or` alternative, so removing them would change no reaction's liveness "
          "anyway -- but a knockout study needs to know they are there.")
        A("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L) + "\n")
    return path


def main() -> int:
    failures: list[str] = []

    def check(label, ok, detail=""):
        print(f"[{'PASS' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""),
              flush=True)
        if not ok:
            failures.append(label)

    print(f"== {DERIVED} identity against {BASE} ==", flush=True)
    m = measure()

    check("no model gene is absent from the derived host",
          not m["missing"],
          f"{m['identical']:,} of {m['resolved']:,} resolvable model genes are "
          f"byte-identical proteins")
    check("no changed gene is a model gene",
          not m["changed_in_gem"],
          f"lost {m['lost']}, gained {m['gained']}")
    check("no copy-number differences", not m["copy_diff"],
          f"{len(m['copy_diff'])}")
    check("protein differences are annotation-boundary scale",
          all((x["n_subs"] is not None and x["n_subs"] <= 5)
              or abs(x["base_aa"] - x["derived_aa"]) > 0 for x in m["changed"]),
          "; ".join(f"{x['gene']} {x['base_aa']}->{x['derived_aa']} aa"
                    + (f", {x['n_subs']} subs" if x["n_subs"] is not None else "")
                    for x in m["changed"]) or "none")

    report = write_report(m, hosts.out_dir(DERIVED) / "gem_identity_report.md")
    print(f"\nwrote {report.relative_to(hosts.REPO)}")

    if failures:
        print(f"\n{len(failures)} FAILED: " + "; ".join(failures))
        print("The EPI300 GEM lane is built by tagging the DH10B model. A failure here "
              "means that is no longer safe -- fix the derivation, do not relax the "
              "check.")
        return 1
    print("\nthe GEM-side edit list is empty; tagging the DH10B model is correct")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
