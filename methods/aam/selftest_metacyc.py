"""Self-test for the MetaCyc curated AAM member (T3 acceptance), mirroring
``direction/selftest.py``.

Three checks, each a LOUD FLOOR -- a failure means the wiring broke, not that the
chemistry is interesting:

  1. Crosswalk floor. MetaCyc atom-mapping lines must join to an MNXR at >= the
     ``metacyc_member.JOIN_FLOOR`` rate. A curated member that joins almost nothing is a
     broken reac_xref prefix, not a data property (the analogue of ``curated.py``'s
     ``decided > 0.5`` floor).

  2. Map-number collision floor. Within a reaction, a single atom-map number must not span
     two DIFFERENT C/N/S/P elements -- that would let the engine pair a carbon source to a
     nitrogen product. MetaCyc reuses map numbers almost exclusively across equivalent
     oxygens (excluded from ``ELEMENTS``); the C/N/S/P collision rate is ~0 and must stay
     negligible.

  3. Spot-check: MetaCyc agrees with a neural mapper. On a deterministic sample of
     reactions BOTH MetaCyc and RXNMapper cover and both map, the curated member's argmax
     product-rank must match RXNMapper's on the great majority of shared C/N/S/P atoms --
     the "does it argue" check at the atom level, floored so a namespace regression trips it.

Runs under ``scadc-metabolic-model`` (rdkit + the engine import):

    mamba run -n scadc-metabolic-model python selftest_metacyc.py
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
ENGINE_LIB = Path("/home/tony/agentic_workspace/projects/metasmith-libraries/fabfos/resources/lib")
sys.path.insert(0, str(ENGINE_LIB))

from ecspr_atom_pairs import (ELEMENTS, load_equations, parse_equation,       # noqa: E402
                              load_mnxm_smiles, canon_smiles, pairs_from_mapped)
import metacyc_member as mcm                                                  # noqa: E402
from combine_aam import load_mapper, correspondence, _argmax                  # noqa: E402

# Floors. Deliberately loose -- they catch a broken wiring, not a few odd reactions. The
# real collision rate is ~0 and the real agreement is well above these.
COLLISION_FLOOR = 0.005          # fraction of reactions with a cross-element mapnum collision
SPOT_AGREE_FLOOR = 0.80          # metacyc-vs-rxnmapper argmax agreement on the sample
SPOT_SAMPLE = 400                # deterministic sample size for the agreement spot-check

DATA = Path("/home/tony/agentic_workspace/data/scadc")
SMILES_DAT = DATA / "references/metacyc26_flatfiles/atom-mappings-smiles.dat"
REAC_XREF = DATA / "references/metanetx/reac_xref.tsv"
REAC_PROP = DATA / "references/metanetx/reac_prop.tsv"
CHEM_PROP = DATA / "references/metanetx/chem_prop.tsv"
TRY1 = Path("/home/tony/agentic_workspace/projects/scadc/metabolic-modelling/main/"
            "metabolic-modelling/_reference_try1/betweenness/cache")
RXNMAPPER = TRY1 / "aam_rxnmapper_universe.tsv"


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    return ok


def collision_rate(smiles_dat: Path):
    """Fraction of MetaCyc reactions where one atom-map number spans >1 distinct C/N/S/P
    element -- the corruption the engine's ``sel != el`` guard would silently drop."""
    elems = set(ELEMENTS)
    n = coll = 0
    with open(smiles_dat) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or "\t" not in line:
                continue
            _mc, smi = line.split("\t", 1)
            try:
                rxn = AllChem.ReactionFromSmarts(smi, useSmiles=True)
            except Exception:
                rxn = None
            if rxn is None:
                continue
            n += 1
            seen = defaultdict(set)          # mapnum -> {C/N/S/P element symbols}
            for tmpls in (
                [rxn.GetReactantTemplate(i) for i in range(rxn.GetNumReactantTemplates())],
                [rxn.GetProductTemplate(j) for j in range(rxn.GetNumProductTemplates())],
            ):
                for mol in tmpls:
                    for a in mol.GetAtoms():
                        k = a.GetAtomMapNum()
                        if k > 0 and a.GetSymbol() in elems:
                            seen[k].add(a.GetSymbol())
            if any(len(v) > 1 for v in seen.values()):
                coll += 1
    return coll, n


def spot_agreement(sample):
    """metacyc-vs-rxnmapper argmax agreement over a deterministic MNXR sample both cover."""
    rxnm = load_mapper(RXNMAPPER)
    mc = mcm.load_member(SMILES_DAT, REAC_XREF)
    shared = sorted(set(rxnm) & set(mc))[:sample]
    eqs = load_equations(REAC_PROP, set(shared))
    parsed, want = {}, set()
    for r, eq in eqs.items():
        pe = parse_equation(eq)
        if pe:
            parsed[r] = pe
            want |= set(pe[0]) | set(pe[1])
    raw = load_mnxm_smiles(CHEM_PROP, want)
    canon = {m: cs for m, s in raw.items() if (cs := canon_smiles(s))}
    n_shared = n_agree = n_rxn = 0
    for mnxr in shared:
        pe = parsed.get(mnxr)
        if pe is None:
            continue
        subs, prods = pe
        pr, _ = pairs_from_mapped(rxnm[mnxr][0], subs, prods, canon)
        pm, _ = pairs_from_mapped(mc[mnxr][0], subs, prods, canon)
        if not pr or not pm:
            continue
        cr, cm = correspondence(pr), correspondence(pm)
        keys = set(cr) & set(cm)
        if keys:
            n_rxn += 1
        for k in keys:
            n_shared += 1
            if _argmax(cr[k]) == _argmax(cm[k]):
                n_agree += 1
    return n_agree, n_shared, n_rxn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=SPOT_SAMPLE)
    a = ap.parse_args()
    allok = True

    # 1. crosswalk floor (build() already asserts it loudly; here we surface the number)
    per_mnxr, rep = mcm.build(SMILES_DAT, REAC_XREF)
    allok &= check("MetaCyc->MNXR join clears floor",
                   rep["join_frac"] >= mcm.JOIN_FLOOR,
                   f"{rep['join_frac']:.2%} joined (floor {mcm.JOIN_FLOOR:.0%}); "
                   f"{rep['n_mnxr']:,} MNXR, {rep['n_collapsed']:,} collapsed")

    # 2. map-number collision floor
    coll, ntot = collision_rate(SMILES_DAT)
    rate = coll / ntot if ntot else 1.0
    allok &= check("C/N/S/P map-number collision rate negligible",
                   rate < COLLISION_FLOOR,
                   f"{coll:,}/{ntot:,} = {rate:.4%} (floor {COLLISION_FLOOR:.2%})")

    # 3. spot-check agreement with the neural mapper
    ag, sh, nr = spot_agreement(a.sample)
    frac = ag / sh if sh else 0.0
    allok &= check("MetaCyc argmax agrees with RXNMapper on sample",
                   sh > 0 and frac >= SPOT_AGREE_FLOOR,
                   f"{ag:,}/{sh:,} = {frac:.1%} over {nr:,} co-mapped reactions "
                   f"(floor {SPOT_AGREE_FLOOR:.0%})")

    print("\nRESULT:", "ALL PASS" if allok else "FAILURES ABOVE")
    return 0 if allok else 1


if __name__ == "__main__":
    raise SystemExit(main())
