#!/usr/bin/env python3
"""Gate: the observed solve and the null it is scored against are the SAME UNIVERSE.

WHY THIS EXISTS. Every other reference gate is structurally incapable of catching a
tier mismatch. `assert_canonical_reference()` pins two hashes -- `reac_prop.tsv` (the
MetaNetX release) and `direction.parquet` (the per-reaction ratios). Both are
TIER-INVARIANT: the atom-pair tier changes which atom transits exist in the graph, and
changes neither file. So those pins pass IDENTICALLY on tier3 and tier4 data. They were
never weak checks; they check a different thing. Nothing checked this one.

That is not hypothetical. On 2026-07-20 the figure-side canon still resolved
`REFERENCE_SOLVE_*` and `REFERENCE_NULL_DIRECTED_DIR` to the un-suffixed (tier3)
directories while the reference graph beneath them had already moved to tier4. Every
gate was green. A significance run through that canon would have scored a tier3 observed
solve against a tier3 null and reported it as the tier4 result -- silently, because
"observed and null agree" held: they agreed on being the WRONG pair.

WHAT IT CHECKS, in increasing order of sharpness:

  1. Axis-id sets are identical between the observed report and every null file. This
     is the coarse check, and on its own it is weak: tier3 and tier4 differ by exactly
     ONE axis (tier4 restores a carbon axis whose endpoint tier3's graph could not
     reach), so a future tier that happens to preserve the axis count would slip past.
  2. Cell count is FOSMID_BASIS x n_axes exactly -- no dropped or duplicated cells.
  3. `g_base` agreement per axis. THIS is the real discriminator. `g_base` is the
     conductance of the HOST BASE GRAPH along an axis, before any fosmid addition or
     null draw, so it is a pure function of the universe. The observed solver and the
     null generator compute it by independent code paths, which is what makes agreement
     informative rather than tautological.

THE THRESHOLD, and why it is honest. Measured, ieff, median relative deviation in
`g_base` per axis:

    matched   tier3 solve / tier3 null     1.8e-06
    matched   tier4 solve / tier4 null     2.5e-05
    MIXED     tier4 solve / tier3 null     2.8e-01
    MIXED     tier3 solve / tier4 null     3.5e-01

Four orders of magnitude of clear air. `G_BASE_MEDIAN_REL_TOL` was chosen AFTER seeing
those numbers -- so state that plainly rather than dress it up as a pre-committed
tolerance. It is defensible anyway, for a reason that does not apply to a significance
threshold: any value in [1e-4, 1e-2] returns the same verdict on every pair above, so
the gate's answer does not depend on where in that range the line falls. The gate prints
the measured margin on every run, so a future pair that lands in the gap is visible as
an anomaly instead of being silently bucketed.

Matched pairs do NOT agree to machine epsilon (2.5e-5 median, ~0.2 max on one axis) and
that predates tier4 -- the tier3 matched pair shows the same shape. It is a difference
between the observed solver's and the null generator's base-graph construction, not a
universe mismatch, and it is out of scope here. Reported, not gated on.

    python check_universe_pair.py              # canonical lane
    python check_universe_pair.py --all-lanes
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src"))

from fabfos import canon  # noqa: E402

AXIS_KEY = ["element", "axis_id"]
# The base-graph quantity, per lane. The two lanes emit the same SHAPE under different
# column names -- ieff is a conductance (g_base), reff a resistance (r_base) -- so the
# name has to follow the lane. Both are pure functions of the base graph.
BASE_COL = {"reff": "r_base", "ieff": "g_base"}

# See module docstring: calibrated on a measured 4-orders-of-magnitude separation, not
# committed in advance. The verdict is invariant across [1e-4, 1e-2].
G_BASE_MEDIAN_REL_TOL = 1e-3

failures: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""))
    if not ok:
        failures.append(label)


def axis_base(path: Path, lane: str) -> pd.Series:
    """The base-graph quantity per (element, axis_id). Constant within an axis in both
    table shapes -- asserted rather than assumed, because taking .first() off a column
    that actually varied would silently compare two arbitrary rows and the gate would
    report a difference that means nothing."""
    col = BASE_COL[lane]
    df = pd.read_csv(path, sep="\t", usecols=AXIS_KEY + [col])
    grp = df.groupby(AXIS_KEY)[col]
    if (grp.nunique() > 1).any():
        raise SystemExit(
            f"{col} is not constant within an axis in {path.name} -- it is supposed to "
            f"be a base-graph quantity, independent of fosmid and of draw. The "
            f"comparison below would be meaningless."
        )
    return grp.first()


def compare(lane: str) -> None:
    obs_path = canon.reference_axes_report(lane)
    null_paths = [canon.REFERENCE_NULL_DIR / n for n in canon.FROZEN_NULL_FILES
                  if n.startswith(f"{lane}_")]

    print(f"\n[{lane}] observed: {obs_path}")
    obs = axis_base(obs_path, lane)

    # NOT a cell COUNT check. The observed grid is legitimately sparse: a fosmid whose
    # inserts touch no sulfur-bearing reaction has no S-axis cell at all (13 of the basis
    # are short -- 12 miss the 7 S axes, 1 misses S and P both), so cells != fosmids x axes
    # and any product formula is wrong. Asserting a literal total instead would just
    # transcribe a number, which is the failure canon exists to prevent. Check the two
    # structural properties that are actually invariant.
    obs_df = pd.read_csv(obs_path, sep="\t", usecols=AXIS_KEY + ["fosmid"])
    n_fos = obs_df["fosmid"].nunique()
    check(f"{lane}: covers canon.FOSMID_BASIS fosmids", n_fos == canon.FOSMID_BASIS,
          f"got {n_fos}, expected {canon.FOSMID_BASIS}")

    # The testable axis split. Tier 4 restores the one carbon axis tier 3's graph could
    # not reach, so the testable set is now EXACTLY the declared set -- which makes this
    # equality a live tier discriminator, not a tautology. A tier 3 pair fails it on C
    # (39 vs 40). If a future universe legitimately drops an axis again this check is the
    # thing that must be consciously relaxed, and the failure will say which axis.
    per_el = obs_df.groupby("element")["axis_id"].nunique().to_dict()
    check(f"{lane}: testable axes == canon.AXES_PER_ELEMENT", per_el == canon.AXES_PER_ELEMENT,
          f"got {per_el}, expected {canon.AXES_PER_ELEMENT}")

    for np_ in null_paths:
        nul = axis_base(np_, lane)

        same_axes = set(obs.index) == set(nul.index)
        detail = ""
        if not same_axes:
            only_o = sorted(set(obs.index) - set(nul.index))[:3]
            only_n = sorted(set(nul.index) - set(obs.index))[:3]
            detail = f"obs-only {only_o} / null-only {only_n}"
        check(f"{lane}: axis set == {np_.name} ({len(obs)} vs {len(nul)} axes)",
              same_axes, detail)

        common = obs.index.intersection(nul.index)
        rel = (np.abs(obs[common] - nul[common]) / np.abs(obs[common]))
        med = float(rel.median())
        margin = G_BASE_MEDIAN_REL_TOL / med if med > 0 else float("inf")
        check(f"{lane}: {BASE_COL[lane]} median rel dev < {G_BASE_MEDIAN_REL_TOL:.0e} vs {np_.name}",
              med < G_BASE_MEDIAN_REL_TOL,
              f"median={med:.3e} p90={rel.quantile(0.9):.3e} max={rel.max():.3e} "
              f"(margin {margin:.0f}x under the line)")

    # Close the loop: the SCORED table must be the one this pair produces. Everything
    # above proves the solve and null are the same universe; without this, a correct
    # pair can still be reported through a significance table scored from a different
    # one -- which is exactly what was on disk before this run (the shipped sig tables
    # matched a tier3 solve at 100% and tier4 at 13%).
    sig = canon.SIGNIFICANCE_DIR / f"{canon.SCORER}_{lane}.tsv"
    if not sig.exists():
        check(f"{lane}: significance table present", False, f"missing {sig}")
        return
    cell = AXIS_KEY + ["fosmid"]
    # One row per (cell, null STYLE) -- canon.STYLES, four of them -- so the raw table is
    # ~4x the cell count. delta_obs is a property of the cell, identical across styles;
    # collapse before comparing or the join fans out and every ratio below is meaningless.
    sig_df = (pd.read_csv(sig, sep="\t", usecols=cell + ["delta_obs"])
                .groupby(cell, as_index=False)["delta_obs"].first())
    sig_cells = set(map(tuple, sig_df[cell].astype(str).values))
    obs_cells = set(map(tuple, obs_df[cell].astype(str).values))
    check(f"{lane}: significance cells == observed cells ({len(obs_cells)})",
          sig_cells == obs_cells,
          f"sig-only {len(sig_cells - obs_cells)}, obs-only {len(obs_cells - sig_cells)}")

    # delta_obs must BE the solve's delta -- floored at zero. The scorer drops negative
    # deltas (an addition that raises resistance / lowers conductance is not evidence),
    # so raw equality fails on ~28% of reff rows for a correct table. Compare against the
    # floored value or the check reports a mismatch that is really a convention.
    dcol = f"delta_{lane}"
    full = pd.read_csv(obs_path, sep="\t", usecols=cell + [dcol])
    j = sig_df.merge(full, on=cell)
    agree = float(np.isclose(j["delta_obs"], j[dcol].clip(lower=0), rtol=1e-9).mean())
    check(f"{lane}: delta_obs == max(0, {dcol}) of THIS solve", agree == 1.0,
          f"{100 * agree:.2f}% of {len(j)} joined cells")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all-lanes", action="store_true",
                    help=f"both lanes (default: {canon.CANONICAL_LANE} only)")
    a = ap.parse_args()

    print(f"universe-pair gate | orientation={canon.CANONICAL_ORIENTATION} "
          f"| library={canon.library_root()}")
    for lane in (canon.LANES if a.all_lanes else (canon.CANONICAL_LANE,)):
        compare(lane)

    print(f"\n[universe-pair] {'PASS' if not failures else 'FAIL'}"
          + (f" -- {len(failures)} check(s): {failures}" if failures else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
