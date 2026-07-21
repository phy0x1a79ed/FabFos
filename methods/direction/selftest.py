"""Self-test / verification for the direction annotator (T4 acceptance).

Runs the plan's verification section against a produced table:
  * no-op parity: every no-evidence reaction is ratio 1.0 EXACTLY, so the
    annotator reproduces the undirected model where it is silent;
  * ratio<->dG' sign coherence (ratio = exp(dG'/RT), strictly positive);
  * the canon assert accepts the real table and rejects corrupted copies
    (missing row, duplicate MNXR, null ratio).

Runs in p312. `python selftest.py --table <parquet> --base-mnxrs <json>`.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # methods/ on path, for siblings
from fabfos import canon


def check(name, ok):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--base-mnxrs", required=True)
    a = ap.parse_args()
    df = pd.read_parquet(a.table)
    base = set(json.load(open(a.base_mnxrs)))
    allok = True

    # 1. one row per base-graph MNXR, no duplicates
    allok &= check("one row per base-graph MNXR",
                   len(df) == len(base) and set(df["mnxr"]) == base
                   and not df["mnxr"].duplicated().any())

    # 2. no-op parity: no-evidence rows are ratio 1.0 exactly
    noev = df[df["dir_tier"] == 0]
    allok &= check(f"no-evidence rows ({len(noev)}) all ratio == 1.0 exactly",
                   bool((noev["ratio"] == 1.0).all()))

    # 3. ratio is never null and strictly positive
    allok &= check("ratio never null and > 0",
                   bool(df["ratio"].notna().all() and (df["ratio"] > 0).all()))

    # 4. ratio <-> dG' sign coherence: ratio = exp(dG'/RT)
    expected = np.exp(df["dG_prime"].to_numpy() / canon.DIR_RT)
    allok &= check("ratio == exp(dG'/RT) within 1e-9",
                   bool(np.allclose(expected, df["ratio"].to_numpy(), rtol=1e-9, atol=1e-12)))

    # 5. canon assert accepts the real table
    try:
        canon.assert_canonical_direction_table(df, n_reactions=len(base))
        allok &= check("canon.assert_canonical_direction_table accepts real table", True)
    except canon.CanonError as e:
        allok &= check(f"canon assert accepts real table -- RAISED: {e}", False)

    # 6. canon assert rejects corrupted copies
    for label, corrupt in (
        ("missing row", df.iloc[1:].copy()),
        ("duplicate MNXR", pd.concat([df, df.iloc[[0]]], ignore_index=True)),
        ("null ratio", df.assign(ratio=df["ratio"].mask(df.index == 0))),
    ):
        try:
            canon.assert_canonical_direction_table(corrupt, n_reactions=len(base))
            allok &= check(f"canon assert REJECTS {label}", False)
        except canon.CanonError:
            allok &= check(f"canon assert rejects {label}", True)

    print("\nRESULT:", "ALL PASS" if allok else "FAILURES ABOVE")
    sys.exit(0 if allok else 1)


if __name__ == "__main__":
    main()
