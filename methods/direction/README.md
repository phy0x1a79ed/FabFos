# Per-reaction direction annotator

Gives every base-graph reaction a directional conductance ratio
`g_reverse / g_forward = exp(ΔG'/RT)`, so a downstream edge-builder can split an
undirected edge into a forward/reverse pair. Per the standing "ratios, never
one-way gates" ruling, direction is a continuous ratio, never a hard gate; a
reaction with no evidence gets ΔG'=0 → ratio 1.0 → reversible, so the annotator
is a **provable no-op against today's undirected model where it is silent**.

Scope is the **table only** — not the directed-edge rebuild, not the solver.

## The three members (method families, not databases)

| # | member | family | routing | file |
|---|--------|--------|---------|------|
| 1 | eQuilibrator | measured (TECRDB) + group contribution | InChIKey/InChI | `thermo_eq.py` |
| 2 | dGbyG | learned GNN, 100-head deep ensemble | SMILES | `thermo_dgbyg.py` |
| 3 | BioCyc `REACTION-DIRECTION` | curated physiology (TECRDB-independent) | MNXR crosswalk | `curated.py` |

Members 1 and 2 are both TECRDB-fitted (correlated); member 3 is independent.

## Two load-bearing correctness facts

- **Orientation flip (~60%).** `REACTION-DIRECTION` is stated in MetaCyc's equation
  orientation, but MNXref re-canonicalises orientation on import, so a naive
  metacyc→MNXR join inverts the curated direction on **59.9%** of reactions
  (measured; cross-checked against the reac_xref col-3 equation at 99.88%).
  `curated.py` re-expresses every curated call in MNXR orientation by comparing
  compound sets (`chem_xref` MNXM ↔ `reac_prop` MNXR sides). Undecidable cases
  (transport ties, unmappable compounds) are recorded, never guessed.
- **dGbyG wildcard accept.** dGbyG silently returns a confident number for R-group
  wildcards (RDKit gives `*` atomic number 0 → a valid feature vector). The member
  **abstains** on any reaction with a wildcard atom — confirmed on 25.3% of the
  base graph.

## The combiner (`combine.py`)

A weighted vote in ΔG' space, the continuous analogue of the functional-lane
belief pattern:

1. **Thermo vote.** A measured eQ value takes precedence (dGbyG is a lossy readback
   of the same TECRDB number). Two *predictions* (eQ group-contribution arm and/or
   dGbyG) are fused and floored by `TAU_SHARED` (the common-mode TECRDB error their
   spread cannot see) so two correlated predictors never vote as two independent.
2. **Curated prior.** The category is mapped to an empirical ΔG' via the T3
   calibration — never scored as a classifier (that is what hit the MCC≈0.4 ceiling).
3. **Fuse** thermo + curated by inverse-variance precision (they are independent).
4. **Shrink** toward 0 by evidence: `λ = σ₀²/(σ₀²+s²)`; `μ_eff = λ·μ`. No evidence →
   λ=0 → ratio 1.0 (reversible **as a limit**, not an if-branch).
5. **Ratio** = `exp(μ_eff/RT)` (median transform; never `E[ratio]`), clamped to a
   finite two-way range.

`dir_tier`/`dir_method`/`dir_confidence` are a **provenance record, not a selection**
— tiers 1 (measured) / 2 (predicted) / 3 (curated) / 0 (none). **`dir_tier > 0`
means "carries direction", NOT "usable": every row is usable and the edge-builder
must read all of them.**

## Calibration (T3) findings

- REVERSIBLE bin median = **0.00** on the measured arm — the reversible default is
  empirically anchored, not just stipulated.
- PHYSIOL antisymmetry holds (−12.3 vs +14.0); plain LEFT-TO-RIGHT/RIGHT-TO-LEFT sit
  near 0 (the names mislead — binning beats ranking-by-name).
- Measured arm is only ~15% of eQ answers, so calibration uses all MetaCyc∩eQ.

## Constants — all in `canon.py`

`DIR_RT`, `DIR_DECADE`, `DIR_TAU_SHARED`, `DIR_TAU_CUR_FLOOR`, `DIR_SIGMA_CEILING`,
`DIR_SIGMA_0` (=9.505, frozen from the calibration run; band [5,40]), `DIR_DG_CLAMP`,
`DIR_CATEGORIES`, `DIR_COLUMNS`, plus artifact paths (`DIR_TABLE`, …). The assert
`canon.assert_canonical_direction_table` refuses a table with a missing/duplicate
MNXR or a null ratio.

## Run / verify

```
bash run.sh            # regenerate end to end across the three envs
```

Envs: `envs/{equilibrator,dgbyg}.yml` (+ `p312`). dGbyG is vendored at
`lib/vendor/dGbyG` (weights bundled). Outputs → `data/scadc/direction/`
(`direction_annotation.parquet` is the deliverable). `selftest.py` checks the
no-op parity, ratio well-formedness, and the canon assert.
