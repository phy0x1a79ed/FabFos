"""Head-to-head: score the directed observations and tabulate directional-specificity
AUC against the undirected incumbent baseline, carbon, same facet, same cells.

Runs the benchmark's own reference scorer (30_score.py) on the directed observations,
then joins the resulting specificity table with the incumbent baseline's. Rank-only
metric, so we report AUCs and CIs, never raw magnitudes.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
BENCH_DIR = Path("/home/tony/agentic_workspace/projects/scadc/fig-model/main/ecspr/benchmark")
BASELINE = Path("/home/tony/agentic_workspace/data/scadc/ecspr/benchmark/v1/baseline")
FACET_MAP = {"A_iECDH10B": "netA_iECDH10B", "A_iML1515": "netA_iML1515",
             "B_iECDH10B": "netB_iECDH10B", "B_iML1515": "netB_iML1515"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--facet", default="A_iECDH10B")
    ap.add_argument("--obs", type=Path, default=None)
    ap.add_argument("--out", type=Path, default=HERE / "bench_out")
    args = ap.parse_args()
    bench_facet = FACET_MAP[args.facet]
    obs = args.obs or (args.out / f"observations_directed_{args.facet}.tsv")
    if not obs.exists():
        sys.exit(f"observations not found: {obs} (has the --run finished?)")
    scored = args.out / f"scored_directed_{args.facet}"
    scored.mkdir(parents=True, exist_ok=True)

    # 1. score the directed arm with the benchmark's own scorer
    cmd = [sys.executable, str(BENCH_DIR / "30_score.py"), "--obs", str(obs),
           "--mode", "signed", "--out", str(scored)]
    print("scoring:", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)

    # 2. join directed vs incumbent baseline specificity, carbon
    dsp = pd.read_csv(scored / "specificity_signed.tsv", sep="\t")
    bsp = pd.read_csv(BASELINE / "specificity.tsv", sep="\t")
    d = dsp[(dsp.facet == bench_facet) & (dsp.element == "C")].iloc[0]
    b = bsp[(bsp.facet == bench_facet) & (bsp.element == "C")].iloc[0]

    print("\n" + "=" * 68)
    print(f"DIRECTIONAL-SPECIFICITY AUC (carbon, {bench_facet}) — same cells")
    print("=" * 68)
    print(f"{'arm':<26}{'AUC':>8}{'CI_lo':>9}{'CI_hi':>9}{'n_cond':>8}{'passes':>8}")
    print("-" * 68)
    print(f"{'undirected (incumbent)':<26}{b.auc:>8.4f}{b.ci_lo:>9.4f}{b.ci_hi:>9.4f}"
          f"{int(b.n_conditions):>8}{str(bool(b.passes)):>8}")
    print(f"{'directed (real ratios)':<26}{d.auc:>8.4f}{d.ci_lo:>9.4f}{d.ci_hi:>9.4f}"
          f"{int(d.n_conditions):>8}{str(bool(d.passes)):>8}")
    delta = d.auc - b.auc
    print("-" * 68)
    print(f"{'Δ AUC (directed − undirected)':<26}{delta:>+8.4f}")
    overlap = not (d.ci_lo > b.ci_hi or b.ci_lo > d.ci_hi)
    print(f"\nCIs {'OVERLAP (no significant difference)' if overlap else 'DISJOINT (significant)'}")
    print(f"direction {'improves' if delta > 0 else 'does not improve'} biological alignment "
          f"by ΔAUC={delta:+.4f} on this facet.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
