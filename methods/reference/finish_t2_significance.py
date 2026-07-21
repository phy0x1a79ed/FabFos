#!/usr/bin/env python
"""T2 finish: score canon.AXIS_SET significance on the HONEST reference graph, compare.

The reference-graph null (``build_reference_null.py``) and the reference-graph
observed solve (``REF_ROOT/solve/{reff,ieff}_axes_report.tsv``) are now on the SAME
graph -- the invariant the whole promotion turns on. This driver runs the canonical
mixture scorer on that matched pair and reports how the promoted result compares to
the incumbent (a different graph), so the honesty promotion is delivered as a number,
not a claim.

Steps, per lane:
  1. Build a CURATED nulls dir: symlinks to exactly ``canon.FROZEN_NULL_FILES`` from
     ``canon.REFERENCE_NULL_DIR`` -- never a glob (the cache carries retired sizes and
     the publish tree carries a ``_work`` copy; either would silently redefine the
     draw-size basis the scorer reads).
  2. Run the engine scorer ``ecspr_significance.py run`` (mixture SF) with the reference
     observed report + the curated reference null.
  3. ``canon.assert_canonical_significance`` on the output (schema + basis + splits).
  4. Compare to ``canon.incumbent_sig_table(lane)`` BY EFFECT SIZE -- the way ECSPr
     significance is actually read (the mixture ``survives`` q<0.05 flag saturates on the
     heavy tail and is a diagnostic, not the ranking). Primary metric: axis-level
     effect-size rank correlation (max ``delta_obs`` per axis) + top-K axis overlap;
     secondary: per-cell effect-size Spearman. High agreement is the expected/honest
     outcome -- the reference refuses fabricated transits but preserves the ranking the
     paper reads. The survivor-flag Jaccard is printed only as a labelled diagnostic (it
     reshuffles near-threshold even when the ranking is preserved).

Outputs land in ``REF_ROOT/significance/sig_mix_{lane}.tsv`` (a NEW dir; the frozen
``solve/`` artifacts are hash-pinned and never touched).

Env: ml.

    python finish_t2_significance.py            # both lanes, full compare
    python finish_t2_significance.py --lane reff
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
from fabfos import canon                                   # noqa: E402

SIG_SCRIPT = canon.ENGINE_LIB / "resources" / "lib" / "ecspr_significance.py"
CURATED = canon.REFERENCE_NULL_DIR / "_curated"
# Tier-suffixed, from canon -- NOT `REFERENCE_ROOT / "significance"`. That literal wrote
# into the tier3 directory, so a tier4 run would have overwritten the only tables that
# record what the tier change moved.
SIG_DIR = canon.REFERENCE_SIGNIFICANCE_DIR
CELL_KEY = ["fosmid", "element", "axis_id"]


def build_curated() -> None:
    """Symlink exactly canon.FROZEN_NULL_FILES into a clean dir the scorer can glob."""
    CURATED.mkdir(parents=True, exist_ok=True)
    # clear any stragglers so the derived draw-size basis is exactly the curated list
    for old in CURATED.glob("*.tsv"):
        old.unlink()
    for name in canon.FROZEN_NULL_FILES:
        src = canon.REFERENCE_NULL_DIR / name
        if not src.exists():
            raise SystemExit(f"reference null missing: {src} -- null-gen incomplete")
        # The scorer (ecspr_significance._NULL_RE) globs the fixed `{lane}_null_canonical_N{n}`
        # convention and reads draw sizes only from names matching it. The curated dir is this
        # module's private staging for that scorer, so normalise the orientation stem to what it
        # reads -- the DIRECTED null is the canonical null now; its content is unchanged.
        curated_name = name.replace("_null_directed_", "_null_canonical_")
        (CURATED / curated_name).symlink_to(src)


def run_scorer(lane: str) -> Path:
    SIG_DIR.mkdir(parents=True, exist_ok=True)
    out = SIG_DIR / f"{canon.SCORER}_{lane}.tsv"
    cmd = [
        sys.executable, str(SIG_SCRIPT), "run",
        "--lane", lane,
        "--nulls-dir", str(CURATED),
        "--report", str(canon.reference_axes_report(lane)),
        "--faa", str(canon.ORFS_FAA),
        "--out", str(out),
    ]
    print("[scorer]", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    return out


def compare(lane: str, new_path: Path) -> dict:
    """Compare reference vs incumbent by EFFECT SIZE, which is how ECSPr significance
    is read (the mixture q<0.05 `survives` flag saturates on the heavy tail and is a
    diagnostic, not the ranking -- see canon/MIGRATION). The primary agreement metric
    is therefore axis-level effect-size rank correlation + top-K axis overlap; the
    survivor-set Jaccard is kept only as a labelled secondary diagnostic.
    """
    new = canon.assert_canonical_significance(new_path, lane=lane)
    inc = pd.read_csv(canon.incumbent_sig_table(lane), sep="\t")

    # --- PRIMARY: axis-level effect size (max delta_obs across fosmids per axis) ---
    na = new.groupby("axis_id")["delta_obs"].max()
    ia = inc.groupby("axis_id")["delta_obs"].max()
    common = na.index.intersection(ia.index)
    axis_rho, _ = spearmanr(na[common], ia[common]) if len(common) > 2 else (float("nan"), None)
    def topk_axes(k):
        nt, it = set(na.nlargest(k).index), set(ia.nlargest(k).index)
        return len(nt & it)
    top10, top20 = topk_axes(10), topk_axes(20)

    # --- SECONDARY: per-cell effect-size rank corr over jointly-present cells ---
    nc = new.groupby(CELL_KEY, as_index=False)["delta_obs"].max()
    ic = inc.groupby(CELL_KEY, as_index=False)["delta_obs"].max()
    j = nc.merge(ic, on=CELL_KEY, suffixes=("_new", "_inc"))
    cell_rho, _ = spearmanr(j["delta_obs_new"], j["delta_obs_inc"]) if len(j) else (float("nan"), None)

    # --- DIAGNOSTIC ONLY: saturating survivor-flag Jaccard (not the ranking) ---
    def surv_set(df):
        s = df[df["survives"].astype(bool)]
        return set(map(tuple, s[CELL_KEY].astype(str).values.tolist()))
    a, b = surv_set(new), surv_set(inc)
    jacc = len(a & b) / len(a | b) if (a | b) else float("nan")

    r = dict(lane=lane, n_axes=len(common),
             axis_spearman=float(axis_rho), top10_axes=top10, top20_axes=top20,
             cell_spearman=float(cell_rho), n_cells=len(j),
             surv_jaccard_diag=jacc, n_new_surv=len(a), n_inc_surv=len(b))
    print(f"[compare:{lane}] EFFECT-SIZE (primary): axis Spearman={axis_rho:.3f} over "
          f"{len(common)} axes | top-10 axes {top10}/10 | top-20 {top20}/20 | "
          f"cell Spearman={cell_rho:.3f}", flush=True)
    print(f"[compare:{lane}] survivor-flag Jaccard (DIAGNOSTIC, saturating q<0.05, "
          f"not the ranking) = {jacc:.3f}  [ref {len(a)} / inc {len(b)} survivors]",
          flush=True)
    return r


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lane", choices=list(canon.LANES), default=None,
                    help="one lane; default both")
    a = ap.parse_args()

    lanes = [a.lane] if a.lane else list(canon.LANES)
    build_curated()
    rows = []
    for lane in lanes:
        out = run_scorer(lane)
        rows.append(compare(lane, out))
    print("\n[T2 summary]")
    print(pd.DataFrame(rows).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
