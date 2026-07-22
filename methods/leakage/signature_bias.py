"""Does adding a universal-leakage ground BIAS the per-precursor conductance signature?

The hard ground reports only set4 precursors -- outside the contracted sink, net current is
identically zero by KCL, so there is nothing else it *could* report. A leakage ground lets
every atom node drain into Omega, so the precursors now compete with a background sink. The
worry this answers: does that background pull the per-precursor attribution around, i.e. is
the leak's readout a biased version of the hard readout rather than the same signature?

Two levels, both restricted to the set4 precursor set so the comparison is like-for-like:

BASELINE SIGNATURE
    The unperturbed share vector over precursors, ``draw_base / sum(draw_base)``. One vector
    per (lane, element, mode). Compared hard-vs-leak by cosine, Spearman, and the largest
    per-precursor share difference. A bias shows up here as a systematic reweighting.

DELTA SIGNATURE
    Per condition, the vector of ``delta_draw`` over precursors -- WHERE the perturbation
    moved the carbon/nitrogen, which is the thing the instrument is actually read for.
    Compared hard-vs-leak per condition; reported as the distribution over conditions.

Also reported, because a bias would most plausibly enter through metabolite SIZE: the
correlation between a precursor's share difference (leak - hard) and its atom count. That is
the confound ``leak_norm`` was built to remove, so hard-vs-``leak_raw`` and
hard-vs-``leak_norm`` are both shown.

Reads the run's TSV shards. No solving. Env: pandas + scipy.
"""
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def _cos(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    return float(a @ b / (na * nb)) if na > 0 and nb > 0 else float("nan")


def _spear(a, b):
    if len(a) < 3 or np.allclose(a, a[0]) or np.allclose(b, b[0]):
        return float("nan")
    return float(spearmanr(a, b).statistic)


def load(out_dir: Path) -> pd.DataFrame:
    fs = sorted(glob.glob(str(out_dir / "gof_*.tsv")))
    if not fs:
        raise SystemExit(f"no gof_*.tsv under {out_dir}")
    df = pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True)
    # the in-file reachability columns are known-wrong (endpoint-key diagnostic); this
    # analysis never uses them, but drop them so nothing downstream can pick them up
    return df.drop(columns=[c for c in ("reachable", "n_edges_added",
                                        "n_edges_touching_source_component")
                            if c in df.columns])


def baseline(df: pd.DataFrame) -> pd.DataFrame:
    """One share vector per (lane, element, mode) -> pairwise hard-vs-leak agreement."""
    b = (df.groupby(["lane", "element", "mode", "metabolite"], as_index=False)
           .draw_base.first())
    rows = []
    for (lane, el), sub in b.groupby(["lane", "element"]):
        piv = sub.pivot(index="metabolite", columns="mode", values="draw_base")
        if "hard" not in piv.columns:
            continue
        piv = piv.fillna(0.0)
        tot = piv.sum(axis=0)
        if float(tot.get("hard", 0.0)) <= 0:
            rows.append(dict(lane=lane, element=el, note="hard total is zero (unreachable)",
                             n_prec=len(piv)))
            continue
        share = piv / tot.replace(0.0, np.nan)
        h = share["hard"].to_numpy(float)
        for mode in ("leak_raw", "leak_norm"):
            if mode not in share.columns:
                continue
            k = share[mode].to_numpy(float)
            ok = np.isfinite(h) & np.isfinite(k)
            rows.append(dict(
                lane=lane, element=el, mode=mode, n_prec=int(ok.sum()),
                cosine=_cos(h[ok], k[ok]), spearman=_spear(h[ok], k[ok]),
                max_abs_share_diff=float(np.max(np.abs(k[ok] - h[ok]))) if ok.any() else np.nan,
                l1_share_diff=float(np.sum(np.abs(k[ok] - h[ok]))) if ok.any() else np.nan,
                top_prec_hard=share["hard"].idxmax(), top_prec_leak=share[mode].idxmax()))
    return pd.DataFrame(rows)


def deltas(df: pd.DataFrame) -> pd.DataFrame:
    """Per-condition delta_draw vectors, hard vs each leak mode.

    Carries the condition's own effect size (``rel_total`` under ``hard``) so that a poor
    cosine can be attributed. A near-zero delta vector has no well-defined direction: its
    cosine is dominated by the last bits of the solve, and calling that "the ground biased
    the signature" would be reading noise as a finding.
    """
    rows = []
    for (lane, el, cond), sub in df.groupby(["lane", "element", "condition_id"]):
        piv = sub.pivot_table(index="metabolite", columns="mode", values="delta_draw",
                              aggfunc="first")
        if "hard" not in piv.columns:
            continue
        piv = piv.fillna(0.0)
        h = piv["hard"].to_numpy(float)
        if not np.any(h):                       # negative control: exact zero, no direction
            continue
        hard_rows = sub[sub["mode"] == "hard"]
        rel = float(hard_rows["rel_total"].abs().max()) if len(hard_rows) else float("nan")
        for mode in ("leak_raw", "leak_norm"):
            if mode not in piv.columns:
                continue
            k = piv[mode].to_numpy(float)
            if not np.any(k):
                continue
            rows.append(dict(lane=lane, element=el, condition_id=cond, mode=mode,
                             abs_rel_total_hard=rel,
                             cosine=_cos(h, k), spearman=_spear(h, k),
                             ratio_l2=float(np.linalg.norm(k) / np.linalg.norm(h)),
                             argmax_agrees=int(np.argmax(np.abs(h)) == np.argmax(np.abs(k)))))
    return pd.DataFrame(rows)


def by_effect_size(dl: pd.DataFrame) -> pd.DataFrame:
    """Delta-signature agreement bucketed by the condition's own effect size."""
    d = dl.dropna(subset=["abs_rel_total_hard"]).copy()
    if d.empty:
        return d
    d["bucket"] = pd.cut(d.abs_rel_total_hard,
                         [-np.inf, 1e-9, 1e-7, 1e-6, 1e-5, np.inf],
                         labels=["<1e-9", "1e-9..1e-7", "1e-7..1e-6",
                                 "1e-6..1e-5", ">1e-5"])
    return (d.groupby(["lane", "element", "mode", "bucket"], observed=True)
             .agg(n=("cosine", "size"), cos_median=("cosine", "median"),
                  frac_cos_gt_099=("cosine", lambda s: float((s > 0.99).mean())),
                  frac_cos_lt_0=("cosine", lambda s: float((s < 0).mean())))
             .reset_index())


def per_precursor_baseline(df: pd.DataFrame, lane: str, el: str) -> pd.DataFrame:
    """The actual share numbers for one cell -- a summary statistic can hide a big move
    on a small precursor, so the cell with the largest L1 gets shown in full."""
    b = (df[(df.lane == lane) & (df.element == el)]
         .groupby(["mode", "metabolite"], as_index=False).draw_base.first())
    piv = b.pivot(index="metabolite", columns="mode", values="draw_base").fillna(0.0)
    share = piv / piv.sum(axis=0).replace(0.0, np.nan)
    share = share.sort_values("hard", ascending=False)
    for m in ("leak_raw", "leak_norm"):
        if m in share.columns:
            share[f"d_{m}"] = share[m] - share["hard"]
    return share


def size_confound(df: pd.DataFrame, natoms: dict | None) -> pd.DataFrame:
    if not natoms:
        return pd.DataFrame()
    b = (df.groupby(["lane", "element", "mode", "metabolite"], as_index=False)
           .draw_base.first())
    rows = []
    for (lane, el), sub in b.groupby(["lane", "element"]):
        piv = sub.pivot(index="metabolite", columns="mode", values="draw_base").fillna(0.0)
        if "hard" not in piv.columns or piv["hard"].sum() <= 0:
            continue
        share = piv / piv.sum(axis=0).replace(0.0, np.nan)
        n = np.array([natoms.get((el, m), np.nan) for m in share.index], float)
        for mode in ("leak_raw", "leak_norm"):
            if mode not in share.columns:
                continue
            d = (share[mode] - share["hard"]).to_numpy(float)
            ok = np.isfinite(d) & np.isfinite(n)
            if ok.sum() < 3:
                continue
            rows.append(dict(lane=lane, element=el, mode=mode, n=int(ok.sum()),
                             spearman_shareDiff_vs_natoms=_spear(d[ok], n[ok])))
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="run dir holding out/gof_*.tsv")
    ap.add_argument("--natoms", default=None,
                    help="optional TSV: element<TAB>metabolite<TAB>n_atoms")
    a = ap.parse_args(argv)
    run = Path(a.run)
    df = load(run / "out")

    natoms = None
    if a.natoms:
        t = pd.read_csv(a.natoms, sep="\t")
        natoms = {(r.element, r.metabolite): r.n_atoms for r in t.itertuples()}

    pd.set_option("display.width", 200, "display.max_columns", 40,
                  "display.float_format", lambda v: f"{v:.4f}")

    print("=" * 78)
    print("BASELINE SIGNATURE -- unperturbed per-precursor share, hard vs leak")
    print("=" * 78)
    bl = baseline(df)
    print(bl.to_string(index=False))

    print()
    print("=" * 78)
    print("DELTA SIGNATURE -- per-condition delta_draw vector, hard vs leak")
    print("=" * 78)
    dl = deltas(df)
    if dl.empty:
        print("  no conditions with a nonzero hard delta")
    else:
        g = (dl.groupby(["lane", "element", "mode"])
               .agg(n=("cosine", "size"),
                    cos_median=("cosine", "median"),
                    cos_p05=("cosine", lambda s: float(np.nanpercentile(s, 5))),
                    cos_min=("cosine", "min"),
                    frac_cos_gt_099=("cosine", lambda s: float((s > 0.99).mean())),
                    spear_median=("spearman", "median"),
                    argmax_agree=("argmax_agrees", "mean"),
                    l2_ratio_median=("ratio_l2", "median")).reset_index())
        print(g.to_string(index=False))
        worst = dl.sort_values("cosine").head(8)
        print("\n  lowest-cosine conditions:")
        print(worst.to_string(index=False))
        print("\n  agreement bucketed by the condition's OWN effect size "
              "(abs rel_total under hard):")
        print(by_effect_size(dl).to_string(index=False))

    print()
    print("=" * 78)
    print("PER-PRECURSOR BASELINE SHARE -- the cell with the largest L1 difference")
    print("=" * 78)
    bl2 = bl.dropna(subset=["l1_share_diff"])
    if not bl2.empty:
        r = bl2.loc[bl2.l1_share_diff.idxmax()]
        print(f"  {r.lane}/{r.element}  (L1 = {r.l1_share_diff:.4f} under {r['mode']})")
        print(per_precursor_baseline(df, r.lane, r.element).to_string())

    sc = size_confound(df, natoms)
    if not sc.empty:
        print()
        print("=" * 78)
        print("SIZE CONFOUND -- (leak share - hard share) vs precursor atom count")
        print("=" * 78)
        print(sc.to_string(index=False))

    bl.to_csv(run / "signature_bias_baseline.tsv", sep="\t", index=False)
    if not dl.empty:
        dl.to_csv(run / "signature_bias_delta.tsv", sep="\t", index=False)
    print(f"\nwrote {run}/signature_bias_baseline.tsv"
          f"{' and signature_bias_delta.tsv' if not dl.empty else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
