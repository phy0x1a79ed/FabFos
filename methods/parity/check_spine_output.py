#!/usr/bin/env python
"""T3 gate: the SPINE's output must reproduce the incumbent.

    python check_spine_output.py [--experiment scadc_main]

The T2 gate proves the SCORER reproduces the incumbent when driven by hand. This
one proves the whole path does: staging, planning, the container, the transform
wiring, and the scorer together. A regression in any of them shows up here as a
number, not as a plausible new table.

Outputs are content-addressed, so the lane is in the DIRECTORY name and the
filename is a hash. Resolve by directory; never by filename.

Note this gate FAILS when it finds nothing to compare. That is not defensive
padding -- the first draft of this check "passed" by comparing zero lanes, which
is exactly the class of vacuous green this whole exercise exists to eliminate.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fabfos import canon  # noqa: E402

HERE = Path(__file__).resolve().parent
RUNS = HERE.parent / ".runs"


def latest_results(experiment: str) -> Path:
    runs = sorted((RUNS / experiment / "agent_home" / "runs").glob("*/results"),
                  key=lambda p: p.stat().st_mtime)
    if not runs:
        raise SystemExit(f"no completed runs for [{experiment}] under {RUNS}")
    return runs[-1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--experiment", default="scadc_main")
    a = ap.parse_args()

    results = latest_results(a.experiment)
    print(f"results: {results}")

    worst, compared = 0.0, 0
    for lane in canon.LANES:
        d = results / f"ecspr-{lane}_significance"
        tsvs = sorted(d.glob("*.tsv")) if d.is_dir() else []
        if not tsvs:
            print(f"  {lane}: NO OUTPUT at {d}")
            continue
        # the spine's output must survive the same assertion a figure generator
        # applies -- if it would not, the run is not usable regardless of parity
        mine = canon.assert_canonical_significance(tsvs[0])
        ref = pd.read_csv(canon.incumbent_sig_table(lane), sep="\t")
        keys = ["fosmid", "element", "axis_id", "null"]
        m = mine.merge(ref, on=keys, suffixes=("_mine", "_ref"))
        if not len(m):
            print(f"  {lane}: joined 0 rows against the incumbent")
            continue
        dp = float(np.abs(m["p_mine"] - m["p_ref"]).max())
        dpe = float(np.abs(m["p_emp_mine"] - m["p_emp_ref"]).max())
        af = int((m["at_floor_mine"] != m["at_floor_ref"]).sum())
        unjoined = len(ref) - len(m)
        worst = max(worst, dp, dpe)
        compared += 1
        print(f"  {lane}: rows={len(mine)} ref={len(ref)} joined={len(m)} "
              f"unjoined={unjoined} | max|dp|={dp:.3e} max|dp_emp|={dpe:.3e} "
              f"| at_floor disagreements={af}")
        if af or unjoined:
            worst = float("inf")

    print()
    if compared != len(canon.LANES):
        print(f"FAIL: compared {compared}/{len(canon.LANES)} lanes. "
              f"A gate that compares nothing is not a gate.")
        return 1
    print(f"worst |d| = {worst:.3e}   (canon.PARITY_TOL = {canon.PARITY_TOL:.0e})")
    if worst >= canon.PARITY_TOL:
        print("T3 GATE: FAIL -- the spine does not reproduce the incumbent.")
        return 1
    print("T3 GATE: PASS -- the spine reproduces the incumbent end to end.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
