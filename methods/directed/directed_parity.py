#!/usr/bin/env python
"""T4 symmetric-limit parity gate: the DIRECTED solve at ratio->1 reproduces the
UNDIRECTED reference axes_report to solver tolerance.

This is the replacement solver gate that lands with the canon repoint to the directed
model. The undirected reference-solve parity (parity/run_solver_parity.py) proves the
engine's undirected solve reproduces the frozen reference report; this proves the NEW
directed path degrades to that same undirected answer exactly where direction is unknown
(``--force-noop`` forces every ratio to 1.0). Together they mean the directed model is
never silently different from the undirected one it generalises -- any divergence in a
real directed run (ratios != 1) is attributable to genuine forward/backward asymmetry,
never to the wiring.

What it exercises that the engine self-tests do not: the real reference base graphs, the
MetaNetX role parse, the LCC extraction, the addition-side ``build_ar2m_directed``
orientation, and the reff/ieff report schema -- the whole ``cmd_solve_directed`` path on
production inputs, not a toy net.

The engine primitive's symmetric-limit parity is separately gated exhaustively inside
``ecspr_directed.py`` (``_selftest_directed_parity`` ~5e-17, ``_selftest_orient_and_solve``,
``_selftest_build_ar2m_directed`` 2.9e-15). This is the integration counterpart.

Env: ml (or any env with numpy/pandas/scipy/networkx + the engine on sys.path).

    python directed_parity.py            # 3 axes/element (fast wiring gate ~3 min)
    python directed_parity.py --full     # every testable axis (~15-20 min)
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from fabfos import canon                                   # noqa: E402

ENGINE = canon.ENGINE_LIB / "resources" / "lib" / "ecspr_network.py"
KEY = ["element", "axis_id", "fosmid"]
LANE_COLS = {"reff": ["r_base", "r_aug", "delta_reff"],
             "ieff": ["g_base", "g_aug", "delta_ieff"]}
# The base/aug potentials are the direct solver quantity -- these must match to ~machine
# epsilon. delta is a difference of two near-equal large numbers, so its ABSOLUTE error
# is what is meaningful (its relative error blows up wherever the true delta -> 0); gate
# delta on abs, base/aug on both.
PARITY_ABS = 1e-9


def run_forcenoop(out_reff: Path, out_ieff: Path, *, full: bool) -> None:
    cmd = [
        sys.executable, str(ENGINE), "solve-directed", "--force-noop",
        "--addition", str(canon.ADDITION_WEIGHTS),
        "--base-dir", str(canon.SOLVE_BASE_DIR),
        "--bipartite-dir", str(canon.BIPARTITE_DIR),
        "--axes", str(canon.AXES_JSON),
        "--testable", str(canon.AXES_TESTABLE_JSON),
        "--roles", str(canon.DIR_REAC_PROP),
        "--elements", *canon.ELEMENTS,
        "--out-reff", str(out_reff), "--out-ieff", str(out_ieff),
    ]
    if not full:
        cmd.append("--smoke")
    print("[gate]", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def compare(lane: str, directed: Path) -> bool:
    dz = pd.read_csv(directed, sep="\t")
    # Compare against the UNDIRECTED reference explicitly: this gate proves the directed
    # solve at ratio->1 reduces to the undirected 1/R_eff. canon.CANONICAL_ORIENTATION is
    # now "directed", so a bare reference_axes_report(lane) would dispatch to the directed
    # table (the very thing under test / pending production) -- always request undirected.
    ref = pd.read_csv(canon.reference_axes_report(lane, orientation="undirected"), sep="\t")
    j = dz.merge(ref, on=KEY, suffixes=("_d", "_r"))
    if not len(j):
        raise SystemExit(f"[{lane}] no cells joined directed vs reference -- gate cannot run")
    ok = True
    print(f"[{lane}] {len(j)} cells (directed force-noop vs undirected reference)")
    for c in LANE_COLS[lane]:
        a, b = j[f"{c}_d"].values, j[f"{c}_r"].values
        max_abs = float(np.abs(a - b).max())
        status = "ok" if max_abs < PARITY_ABS else "FAIL"
        ok = ok and max_abs < PARITY_ABS
        print(f"   {c:12s} max|abs|={max_abs:.3e}   [{status}]")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--full", action="store_true",
                    help="every testable axis (default: 3 axes/element)")
    a = ap.parse_args()

    with tempfile.TemporaryDirectory() as td:
        tr, ti = Path(td) / "reff.tsv", Path(td) / "ieff.tsv"
        run_forcenoop(tr, ti, full=a.full)
        ok_r = compare("reff", tr)
        ok_i = compare("ieff", ti)
    ok = ok_r and ok_i
    print(f"\n[directed-parity] {'PASS' if ok else 'FAIL'} "
          f"(symmetric limit reproduces the undirected reference to < {PARITY_ABS:.0e})")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
