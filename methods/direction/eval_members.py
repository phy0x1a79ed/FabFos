"""Evaluate one thermo member (eQuilibrator or dGbyG) over a list of MNXRs.

Emits a parquet (mnxr, dg, sigma, flag, reason) -- flag is uses_gc for eQ (True
= group-contribution arm) and wildcard for dGbyG (True = abstained on an R-group).
dg is null when the member is silent on that reaction. This is the per-member
evaluation the combiner joins; keeping it a separate step mirrors the lane
structure (each member staged independently, the combiner fuses).

Run once per member, each under its own env:
  eq    : mamba run -n equilibrator python eval_members.py --member eq ...
  dgbyg : (LD_LIBRARY_PATH=$CONDA_PREFIX/lib) mamba run -n dgbyg python eval_members.py --member dgbyg ...
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from refdata import load_mnxr_stoich, load_mnxm_props


def run(member, mnxrs, stoich, props):
    rows = []
    for i, mnxr in enumerate(mnxrs):
        s = stoich.get(mnxr)
        if s is None:
            rows.append(dict(mnxr=mnxr, dg=None, sigma=None, flag=None, reason="no_stoich"))
            continue
        st, is_bal, is_tr = s
        dg, sig, flag, reason = member.dgr(st, props)
        rows.append(dict(mnxr=mnxr, dg=dg, sigma=sig, flag=flag, reason=reason))
        if (i + 1) % 500 == 0:
            print(f"[eval] {i+1}/{len(mnxrs)}...", flush=True)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--member", required=True, choices=["eq", "dgbyg"])
    ap.add_argument("--mnxrs", required=True, help="json list of MNXR ids")
    ap.add_argument("--reac-prop", required=True)
    ap.add_argument("--chem-prop", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    mnxrs = json.load(open(a.mnxrs))
    stoich = load_mnxr_stoich(Path(a.reac_prop))
    props = load_mnxm_props(Path(a.chem_prop))
    if a.member == "eq":
        from thermo_eq import EquilibratorMember
        member = EquilibratorMember()
    else:
        from thermo_dgbyg import DgbygMember
        member = DgbygMember()

    df = run(member, mnxrs, stoich, props)
    df.to_parquet(a.out, index=False)
    print(f"[eval:{a.member}] {len(df)} reactions -> {a.out}")
    print(df["reason"].value_counts().to_string())
    ok = df[df["dg"].notna()]
    print(f"[eval:{a.member}] answered {len(ok)} ({len(ok)/len(df):.1%})")


if __name__ == "__main__":
    main()
