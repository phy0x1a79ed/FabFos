"""Resolve LASER GOF conditions -> the MNXR reaction set each one ADDS.

Join: conditions_manifest.units (';'-separated gene tokens) -> phenotype_table.gene ->
phenotype_table.condition, whose leading token is the MNXR id. Genes that do not resolve are
REPORTED, never dropped silently: an unresolved gene means the condition is measured with a
smaller insert than the paper describes, which changes its effect size mechanically.
"""
from __future__ import annotations
import csv, os, re
from pathlib import Path

REF = Path(os.environ.get(
    "ECSPR_REF", "/home/tony/agentic_workspace/projects/fabfos/scadc/.awm/data/ref"))
DN = REF / "validation/dual_network"
GT = DN / "ground_truth"
_MNXR = re.compile(r"^(MNXR\d+)")


def _rows(p):
    with open(p, newline="") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def gene_to_rxns() -> dict:
    """gene -> {MNXR, ...} from the GOF arm of the phenotype table."""
    out = {}
    for r in _rows(GT / "phenotype_table.tsv"):
        if r.get("arm") != "GOF":
            continue
        m = _MNXR.match((r.get("condition") or "").strip())
        if m:
            out.setdefault((r.get("gene") or "").strip(), set()).add(m.group(1))
    return out


def gof_conditions() -> tuple:
    """``(conditions, report)``; each condition is a dict with its resolved MNXR set."""
    g2r = gene_to_rxns()
    conds, unresolved = [], {}
    for r in _rows(GT / "conditions_manifest.tsv"):
        if r.get("arm") != "gof":
            continue
        units = [u.strip() for u in (r.get("units") or "").split(";") if u.strip()]
        rxns, miss = set(), []
        for u in units:
            hit = g2r.get(u)
            if hit:
                rxns |= hit
            else:
                miss.append(u)
        if miss:
            unresolved[r["condition_id"]] = miss
        conds.append(dict(condition_id=r["condition_id"], host=r.get("host", ""),
                          units=units, n_units=len(units), rxns=sorted(rxns),
                          n_rxns=len(rxns), n_unresolved=len(miss),
                          expected=(r.get("expected") or "").strip(),
                          is_control=(r.get("is_control") or "").strip(),
                          note=(r.get("note") or "").strip()))
    return conds, unresolved


if __name__ == "__main__":
    conds, unres = gof_conditions()
    ok = [c for c in conds if c["n_rxns"]]
    print(f"GOF conditions in manifest      : {len(conds)}")
    print(f"  with >=1 resolved reaction    : {len(ok)}")
    print(f"  fully unresolved (no reaction): {len(conds) - len(ok)}")
    print(f"  partially unresolved          : {sum(1 for c in conds if c['n_unresolved'] and c['n_rxns'])}")
    import collections
    print(f"  insert size (n_rxns) distribution: "
          f"{sorted(collections.Counter(c['n_rxns'] for c in ok).items())}")
    print(f"  hosts: {collections.Counter(c['host'] for c in ok).most_common()}")
    print(f"  controls: {collections.Counter(c['is_control'] for c in conds).most_common()}")
    print("\n  first 8 resolved conditions:")
    for c in ok[:8]:
        print(f"    {c['condition_id']:26s} n_rxn={c['n_rxns']:2d} {','.join(c['rxns'][:3]):36s} "
              f"| {c['expected'][:60]}")
    if unres:
        print(f"\n  UNRESOLVED GENES in {len(unres)} conditions (first 8):")
        for k, v in list(unres.items())[:8]:
            print(f"    {k:26s} {v}")
