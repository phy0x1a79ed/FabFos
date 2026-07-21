# What moved between the tier-3 and tier-4 results

Measured 2026-07-20, directed solve, full-mixture SF scorer, both lanes. Baseline is the
**tier-3 directed** significance (`significance/sig_mix_{lane}_directed.tsv`), not the
star-graph incumbent — the incumbent is two generations back, and comparing to it would
attribute the graph promotion's movement to this step.

## Read this first: the comparison is CONFOUNDED

The two solves differ by **two** changes, not one, and the smaller one is the tier.

**1. The atom-pair tier (tier 3 → tier 4).** Changes edge *weights* only. Carbon universe
graph:

| | nodes | edges |
|---|---|---|
| tier 3 | 126,636 | 362,482 |
| tier 4 | 126,636 | 362,482 |

Identical topology. The curated MetaCyc increment revises `w_C` on edges that already
existed. It adds no reachability.

**2. The host-base evidence-weights fix.** Much larger. The base was being induced from
weights hand-carried from the incumbent cache while every other input had moved to the
199-fosmid CLEAN basis; the mismatch silently zeroed ~4,688 CLEAN-nominated reactions.
Carbon host base:

| | nodes | edges |
|---|---|---|
| tier 3 solve | 7,199 | 12,127 |
| tier 4 solve | 13,524 | 25,049 |

It roughly doubles.

**Everything below is the JOINT effect of both.** Nothing here isolates the tier. A clean
tier-3-vs-tier-4 measurement needs both solves built on the same fixed-weights base, and
that run does not exist. Given a doubled base against an unchanged topology, expect the
weights fix to dominate.

## Coverage

|  | reff | ieff |
|---|---|---|
| cells, tier 3 | 15,422 | 15,422 |
| cells, tier 4 | 15,621 | 15,621 |
| new | 199 | 199 |
| dropped | **0** | **0** |

The 199 new cells are one carbon axis × the whole fosmid basis: the testable set goes
78 → 79, i.e. becomes the full declared set (`AXES_PER_ELEMENT` C=40 N=23 S=7 P=9).

**Attribute this to the weights fix, not the tier.** fig-model tested it directly:
rebuilding the base from weights coherent with `EVIDENCE_TABLE` restores the axis on
*both* the old and the tier-4 universe. It could not have come from the tier — the tier
does not change topology, and reachability is topology.

## Effect sizes move a lot

| | reff | ieff |
|---|---|---|
| cells with a tier-3 effect | 11,122 | 11,293 |
| moved >10% relative | **95.0%** | **89.2%** |
| median relative move | **93.9%** | **98.5%** |
| per-cell effect-size Spearman | 0.647 | 0.540 |

The typical cell's effect size roughly doubles or halves. Two things are in play and both
point the same way: the base doubled, and Δ is a difference of two near-equal large
numbers, so its *absolute* error is what is bounded while its *relative* value is free to
move wherever the true Δ approaches zero.

**Any per-cell number quoted off tier 3 has to be re-read off tier 4.**

## Rankings hold

| | reff | ieff |
|---|---|---|
| axis-level effect-size Spearman | **0.954** | **0.965** |
| top-10 axis overlap | 8/10 | 9/10 |
| top-20 axis overlap | 17/20 | 18/20 |

This is the level the paper reads, and it survives both changes. Axis-level conclusions
stand; the per-cell numbers behind them do not.

## Survivors reshuffle — a diagnostic, not the result

| | reff | ieff |
|---|---|---|
| survivors, tier 3 | 10 | 26 |
| survivors, tier 4 | 21 | 37 |
| kept / lost / gained | 4 / 6 / 17 | 15 / 11 / 22 |
| Jaccard | 0.148 | 0.312 |

The `survives` flag is `q < 0.05` on a scorer that **saturates**, so significance is not
the discriminating quantity here — effect size is. A survivor set that reshuffles near
threshold while the ranking is preserved is the expected behaviour of a saturating flag
under perturbation, not evidence the answer changed.

Survivor counts roughly double. That tracks the base doubling — a coverage effect, not a
sensitivity claim.

## What this does NOT license

The tier-4 benchmark AUC moved by ~0.20× its own noise floor, with tier 4 inside the star
CI in 40/40 slices. **Do not cite tier 4, or the atom-mapping work behind it, as evidence
that the model got better.** The case for tier 4 is coverage and curation provenance —
16,890 expert-assigned MetaCyc atom maps that nothing in the tree had ever read — not
score. And do not credit tier 4 with the restored axis or the survivor gain; those are
the weights fix.
