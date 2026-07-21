# Directed prototype vs the validation benchmark — head-to-head

The decisive experiment the assembly exists to make: does feeding **real per-reaction
direction** into the effective-conductance model change its **biological alignment** on the
X/Y validation benchmark, versus the undirected incumbent — and what does the benchmark see
that the pulse-chase behavioural gate cannot?

## Setup

- Facet `netA_iECDH10B`, carbon, the same cells the incumbent validation panel scored.
- Undirected arm = the incumbent's `s = log(G_pert/G_base)`. Directed arm = the same panel
  re-solved with the rectified "diode" solver, forward `gp = w`, backward `gm = ratio·w`,
  `ratio = exp(ΔG'/RT)` from the directionality ensemble (`direction_annotation.parquet`).
- Pipeline validated first by the symmetric limit: with every ratio forced to 1.0 the directed
  solver reproduces the incumbent's `s` to **5.2e-12 over 6955 cells**. So any difference in
  the directed arm is direction, not reimplementation.
- Wall-clock: the directed solve is ~0.35 s/cell (semismooth Newton, ~30 iters on the
  2620-node carbon graph — active-set oscillation, intrinsic). The full facet is many hours,
  so the head-to-head runs on a **logged subsample**: all 202 conditions kept, all target
  cells kept, off-target cells capped at 12/condition (seed 0) → 4,275 of 142,675 cells (3.0%).
  Directional specificity is target-vs-off-target *within* a condition, so this keeps the
  metric valid with a stated, smaller negative set.

## Result — directional-specificity AUC (same cells, same conditions)

| arm | AUC | 95% CI | n_cond |
|---|---|---|---|
| undirected (incumbent) | 0.7800 | [0.7272, 0.8353] | 88 |
| directed (real ratios) | 0.7751 | [0.7215, 0.8319] | 88 |

**ΔAUC = −0.005.** The CIs overlap almost entirely. Direction neither improves nor degrades
biological alignment on this facet.

(The incumbent's headline carbon AUC on the *full* cell set is 0.833; the ~0.78 here is the
subsample, not direction — both arms sit at ~0.78 on the same 4,275 cells. The comparison that
isolates direction is same-cells-vs-same-cells above, not against the full-set number.)

## Why ΔAUC ≈ 0 while direction is a real refinement

Direction is **not** a global no-op — but its effect is concentrated and rank-neutral:

- Per-cell change `|effect_directed − effect_undirected|`: **median 0.000**, mean 0.035; **17.5%**
  of cells move by >1e-2 and **5.0%** by >0.1 (max 6.1). Most routes run through symmetric or
  no-evidence edges (ratio 1.0, the proven no-op); a sizeable minority genuinely shift.
- Within-condition rank correlation between the two arms: **Spearman ρ ≈ 0.75** (median over 200
  conditions) — direction *does* reorder cells within a condition; it is not a monotone rescale.
- But that reordering is roughly orthogonal to the target/off-target split, so the AUC — which
  reads only that ordering — does not move.

In one line: **direction changes the conductances (17.5% of cells, substantially for 5%) and
reorders them within conditions, but does not change which cells this benchmark counts as
lit — so its directional-specificity ranking is unchanged.**

## What the pulse-chase gate structurally could not see

The behavioural gate (`pulsechase/`) proves the method is *self-consistent*: atoms conserved,
all-paths measured, and now (GROUP IV) that an irreversible reaction throttles in reverse and
that the symmetric limit reproduces the undirected answer. **None of that can tell you whether
direction helps the biology.** That is a different instrument — the benchmark AUC — and it
returns a specific, falsifiable answer the gate cannot: on this facet, real direction is
biologically **net-neutral** for the directional-specificity ranking. Neither a bug nor a win;
a measurement the behavioural gate is blind to by construction. This is exactly the
underspecification flagged at the outset: direction is now *required* by the acceptance
criteria (GROUP IV) and *measured* by the benchmark, where before neither happened.

## Caveats and next cuts

- One facet, carbon, 3% off-target subsample. The point estimate is stable (both arms move
  together under the subsample), but a full-coverage run and the other three facets would
  tighten the CI and test whether netB (58% of edges directed, vs 36% here) shifts the verdict.
- Potency (needs the enumerated null pool) was not computed for the directed arm; the
  essentiality diagnostic moved trivially (0.558 → 0.567, carbon).
- The directed solve wants the one principled speedup the prior work identified (symbolic-
  factorisation reuse) before a full-panel or HPC run is worthwhile.

## Reproduce

```
python benchmark_bridge.py --facet A_iECDH10B --validate 5           # pipeline parity gate
python benchmark_bridge.py --facet A_iECDH10B --run --offtarget 12 --tol 1e-6
python same_cells_incumbent.py                                       # incumbent on the same cells
python compare_auc.py --facet A_iECDH10B                             # score + tabulate
```
