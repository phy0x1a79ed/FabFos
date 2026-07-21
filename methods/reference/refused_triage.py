"""Triage the refused (unmapped) reactions by RECOVERABILITY -- the reclaim worklist.

The closure ledger adjudicates every reaction, but ~25k carry a REFUSED verdict (unmapped).
This asks: of those, how many are refused for a reason a RULE could honestly convert, vs a
reason that is terminal by construction? The pivotal distinction is the *reason the reaction
is unmapped*:

  * CARRIER-VETOED -- the reaction has real, structured chemistry (>=1 substrate WITH a SMILES
    and >=1 product WITH a SMILES), and is refused only because SOME participant (a generic
    Acceptor / e- / ACP / [protein] / Enzyme-ligand complex) is structureless. Today ANY
    structureless participant vetoes the whole reaction. It should not: a participant that
    carries no scored atom is a conserved carrier the atom lane can EXCISE (a super-node),
    mapping the structured remainder. The honest gate is per-element MASS BALANCE: if the
    structured substrates and products already balance in element X, the structureless
    participants carry no X -> excising them is exact for X (no fabrication). If they don't
    balance, X is hiding on a structureless participant -> abstain on X (diluted), never guess.

  * TERMINAL -- no structured substrate or no structured product (a whole side is
    structureless), or a pseudo/transport reaction, or a fully-generic R-group/class. There is
    no concrete chemistry to map without inventing atoms; these stay refused.

Output: REF/sanity/refused_recoverability.tsv (one row per refused reaction: reason,
reachability, structured-both-sides?, which elements balance on the structured side, tier)
plus a printed summary. This is a WORKLIST, not a mutation -- it changes no verdict; it names
what a balance-gated carrier-excision rule could recover and bounds the honest ceiling.

Env: p312.  Read-only.
"""
from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import reference as ref  # noqa: E402
sys.path.insert(0, "/home/tony/agentic_workspace/projects/metasmith-libraries/fabfos/resources/lib")
from ecspr_aam_rescue import count_el  # noqa: E402

CNSP = ["C", "N", "S", "P"]
_TOK = re.compile(r"(\d+(?:\.\d+)?)?\s*(MNXM\d+|MNXM[A-Za-z0-9]+|WATER|BIOMASS)")


def ce(f: str, X: str) -> int:
    return count_el(f, X) or 0


def sided(reac_prop: Path, want: set[str]) -> dict[str, tuple[list, list]]:
    out = {}
    with open(reac_prop) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or p[0] not in want or "=" not in p[1]:
                continue
            L, R = p[1].split("=", 1)
            f = lambda s: [(float(c) if c else 1.0, m) for c, m in _TOK.findall(s)]
            out[p[0]] = (f(L), f(R))
    return out


def main() -> int:
    d = pd.read_parquet(ref.LEDGER)
    refused = d[d.aam_state == "refused"]
    reason_of = dict(zip(refused.mnxr.astype(str), refused.aam_reason))
    reach_of = dict(zip(d.mnxr.astype(str), d.in_evidence))
    refset = set(refused.mnxr.astype(str))

    side = sided(ref.REAC_PROP, refset)
    allmet = {m for L, R in side.values() for _, m in L + R}
    formula, smi = {}, set()
    with open(ref.CHEM_PROP) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if p[0] in allmet:
                formula[p[0]] = p[3] if len(p) > 3 else ""
                if len(p) > 8 and p[8].strip():
                    smi.add(p[0])
                if len(formula) == len(allmet):
                    break

    rows = []
    for r in refset:
        L, R = side.get(r, ([], []))
        Ss = [(c, m) for c, m in L if m in smi and formula.get(m, "").strip()]
        Ps = [(c, m) for c, m in R if m in smi and formula.get(m, "").strip()]
        both = bool(Ss and Ps)
        bal = []
        if both:
            for X in CNSP:
                ls = sum(c * ce(formula[m], X) for c, m in Ss)
                ps = sum(c * ce(formula[m], X) for c, m in Ps)
                if ls > 0 and ls == ps:
                    bal.append(X)
        tier = ("recoverable-balanced" if bal else
                "carrier-vetoed-unbalanced" if both else "terminal-structureless")
        rows.append(dict(mnxr=r, reason=reason_of[r], reachable=bool(reach_of.get(r, False)),
                         structured_both_sides=both, balances=",".join(bal), tier=tier))

    df = pd.DataFrame(rows)
    out = ref.REF_ROOT / "sanity" / "refused_recoverability.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)

    print(f"refused reactions: {len(df):,}  (reachable {int(df.reachable.sum()):,})\n")
    print("tier                         total   reachable")
    for tier in ("recoverable-balanced", "carrier-vetoed-unbalanced", "terminal-structureless"):
        s = df[df.tier == tier]
        print(f"  {tier:26s} {len(s):6,} {int(s.reachable.sum()):11,}")
    print("\nrecoverable-balanced, by element (a reaction can balance in several):")
    for X in CNSP:
        m = df.balances.str.contains(X, na=False)
        print(f"    {X}: {int(m.sum()):6,} / {int((m & df.reachable).sum()):5,} reachable")
    print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
