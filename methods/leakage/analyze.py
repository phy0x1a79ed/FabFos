"""Are the LASER GOF measurements SENSIBLE and DISTINGUISHABLE, and does the ground matter?

Three questions, each answered by a number rather than an impression:

1. DISTINGUISHABLE -- can the instrument tell conditions apart? A probe that returns the
   same signature for every insert measures nothing, however well-conditioned it looks.
   Read: the spread of |delta_total| across conditions, how many clear a noise floor set by
   the conditions that CANNOT move (unreachable inserts -- a built-in negative control), and
   the pairwise cosine similarity of per-precursor delta signatures. Near-identical
   signatures mean one degree of freedom, not 22.

2. SENSIBLE -- does a condition move the precursors its stated expectation names? The
   manifest's `expected` field names a near-target metabolite per condition; a sensible
   measurement puts its largest response at or near that target rather than uniformly.

3. DOES THE GROUND MATTER -- Spearman of the per-condition effect ranking between the three
   ground modes. If hard and leak_norm rank conditions identically, the ground construction
   is not where the answer lives, and that is the finding.

Unreachable conditions are the load-bearing control: their edges land outside the source's
component, so NO ground can move them. They are what a true zero looks like on this
instrument, which is what makes "above the floor" mean something.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _sig_matrix(d: pd.DataFrame):
    """conditions x precursors matrix of delta_draw, plus the condition index."""
    piv = d.pivot_table(index="condition_id", columns="metabolite", values="delta_draw",
                        aggfunc="first").fillna(0.0)
    return piv.index.to_numpy(), piv.to_numpy()


def _cos(M):
    n = np.linalg.norm(M, axis=1, keepdims=True)
    keep = (n[:, 0] > 0)
    if keep.sum() < 2:
        return np.array([]), keep.sum()
    U = M[keep] / n[keep]
    C = U @ U.T
    iu = np.triu_indices(C.shape[0], k=1)
    return C[iu], int(keep.sum())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True, nargs="+")
    ap.add_argument("--out", default=None)
    ap.add_argument("--reach", default=None,
                    help="corrected reachability table to join on (lane, element, condition_id)")
    a = ap.parse_args(argv)

    df = pd.concat([pd.read_csv(f, sep="\t") for f in a.tsv], ignore_index=True)
    if a.reach:
        fx = pd.read_csv(a.reach, sep="\t")
        df = df.drop(columns=[c for c in ("reachable", "n_edges_added",
                                          "n_edges_touching_source_component")
                              if c in df.columns])
        df = df.merge(fx, on=["lane", "element", "condition_id"], how="left")
        # The HONEST negative control: a graph with no edges appended is byte-identical to
        # the base, so its delta must be exactly zero. Anything else is a solver artefact.
        z = df[df.n_edges_added == 0]
        if len(z):
            print(f"negative control (n_edges_added==0): {z.condition_id.nunique()} conditions, "
                  f"max |rel_total| = {z.rel_total.abs().max():.3e}  "
                  f"max |delta_draw| = {z.delta_draw.abs().max():.3e}")
    print(f"rows {len(df):,}   conditions {df.condition_id.nunique()}   "
          f"lanes {sorted(df.lane.unique())}   elements {sorted(df.element.unique())}   "
          f"modes {sorted(df['mode'].unique())}")

    # per (lane, mode, element, condition): one total-level row
    per = (df.groupby(["lane", "mode", "element", "condition_id"])
             .agg(delta_total=("delta_total", "first"), total_base=("total_base", "first"),
                  rel_total=("rel_total", "first"), reachable=("reachable", "first"),
                  n_edges_added=("n_edges_added", "first"),
                  is_novel=("is_novel", "first"),
                  prec_share_base=("prec_share_base", "first"))
             .reset_index())

    print("\n" + "=" * 96)
    print("REACHABILITY -- the built-in negative control")
    print("=" * 96)
    rr = (per[per["mode"] == "leak_norm"].groupby(["lane", "element"])
          .agg(n=("condition_id", "nunique"), reachable=("reachable", "sum"),
               novel=("is_novel", "sum")).reset_index())
    print(rr.to_string(index=False))

    print("\n" + "=" * 96)
    print("Q1 DISTINGUISHABLE -- effect spread, and the floor set by unreachable inserts")
    print("=" * 96)
    print(f"{'lane':6s} {'mode':10s} {'el':3s} {'base G':>10s} {'n':>4s} {'reach':>6s} "
          f"{'floor |dG/G|':>13s} {'med |dG/G|':>11s} {'max |dG/G|':>11s} {'>10x floor':>11s}")
    rows = []
    for (lane, mode, el), d in per.groupby(["lane", "mode", "element"]):
        if d.total_base.iloc[0] == 0:
            continue
        rel = d.rel_total.abs()
        unreach = rel[d.reachable == 0]
        reach = rel[d.reachable == 1]
        floor = float(np.nanmax(unreach)) if len(unreach) else 0.0
        nabove = int((reach > max(floor, 0) * 10).sum()) if len(reach) else 0
        rows.append(dict(lane=lane, mode=mode, element=el, floor=floor,
                         med=float(np.nanmedian(reach)) if len(reach) else np.nan,
                         mx=float(np.nanmax(reach)) if len(reach) else np.nan,
                         n_above=nabove, n_reach=len(reach)))
        print(f"{lane:6s} {mode:10s} {el:3s} {d.total_base.iloc[0]:10.4g} "
              f"{len(d):4d} {len(reach):6d} {floor:13.3e} "
              f"{rows[-1]['med']:11.3e} {rows[-1]['mx']:11.3e} {nabove:11d}")

    print("\n" + "=" * 96)
    print("Q1b SIGNATURE DIVERSITY -- pairwise cosine of per-precursor delta vectors")
    print("   (1.000 everywhere = every insert moves the network the same way = 1 dof)")
    print("=" * 96)
    print(f"{'lane':6s} {'mode':10s} {'el':3s} {'n_nonzero':>10s} {'median cos':>11s} "
          f"{'p05':>8s} {'p95':>8s} {'frac>0.99':>10s}")
    for (lane, mode, el), d in df[df.reachable == 1].groupby(["lane", "mode", "element"]):
        _, M = _sig_matrix(d)
        c, nz = _cos(M)
        if not len(c):
            continue
        print(f"{lane:6s} {mode:10s} {el:3s} {nz:10d} {np.median(c):11.4f} "
              f"{np.percentile(c, 5):8.4f} {np.percentile(c, 95):8.4f} "
              f"{float((c > 0.99).mean()):10.3f}")

    print("\n" + "=" * 96)
    print("Q3 DOES THE GROUND MATTER -- Spearman of per-condition |dG/G| ranking, mode vs mode")
    print("=" * 96)
    from scipy.stats import spearmanr
    print(f"{'lane':6s} {'el':3s} {'pair':26s} {'spearman':>9s} {'n':>5s}")
    for (lane, el), d in per.groupby(["lane", "element"]):
        w = d.pivot_table(index="condition_id", columns="mode", values="rel_total")
        w = w.dropna()
        modes = [m for m in ("hard", "leak_raw", "leak_norm") if m in w.columns]
        for i in range(len(modes)):
            for j in range(i + 1, len(modes)):
                x, y = w[modes[i]].abs(), w[modes[j]].abs()
                if x.nunique() < 3:
                    continue
                r = spearmanr(x, y).statistic
                print(f"{lane:6s} {el:3s} {modes[i]+' vs '+modes[j]:26s} {r:9.4f} {len(x):5d}")

    print("\n" + "=" * 96)
    print("Q2 SENSIBLE -- top responders per condition vs the stated expectation (netB/C)")
    print("=" * 96)
    sub = df[(df.lane == "netB") & (df.element == "C") & (df["mode"] == "leak_norm")
             & (df.reachable == 1)]
    top = (sub.assign(ad=sub.delta_draw.abs())
              .sort_values("ad", ascending=False)
              .groupby("condition_id")
              .head(2))
    order = (sub.groupby("condition_id").delta_total.first().abs()
                .sort_values(ascending=False).head(15).index)
    for cid in order:
        t = top[top.condition_id == cid]
        exp = t.expected.iloc[0] if len(t) else ""
        hits = ", ".join(f"{r.metabolite}({r.delta_draw:+.2e})" for r in t.itertuples())
        print(f"  {cid:26s} dG/G={t.rel_total.iloc[0]:+.3e}")
        print(f"      top precursors : {hits}")
        print(f"      expected       : {exp[:95]}")

    if a.out:
        per.to_csv(a.out, sep="\t", index=False)
        print(f"\nper-condition table -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
