# ECSPr against benchmark V3, on the tier-4 universe

Run: sockeye job `12319671` (solve) + `12321473` (merge), 2026-07-20.
Agent home `/scratch/st-shallam-1/txyliu/ecspr_bench_v3_1784548824`, run key `cPieN6LE`.
17.4 min wall; 16 shards (4 facets x 4 elements), 32 workers, one node.
Output 5,690,248 rows — exactly the previous run's count.

Method `0.3.1+37a1d54`. **The `+hash` is valid as of this run and no longer resolves
once the frozen null is regenerated**, because the data library's index is one of the
seven hashed components and re-staging moves it. The component values are therefore
recorded here rather than the id alone:

| component | value |
|---|---|
| canon | unchanged from 0.3.0 (`sha256` of `canon.py`; the repoint touched no code) |
| transform_library | `data.epi300.1-58-g508cb77`, dirty=false |
| metasmith | 0.18.8 |
| containers | ecspr `sha256:97f7afad99f55fe66a213851f48f6293f6b328ce7db4fe035a93de6f24f582db` |
| data_library | index moved (tier-4 graph, universe, base graphs, answer key) |
| type_contract | unchanged |
| domains | unchanged |

## What changed, and what it was worth

The previous result was scored against the **incumbent star graph** — dated 13 June,
weighted from `aam_unified_expanded.tsv`, the graph `build_atom_graph.py` indicts for
fabricating transits for structureless stubs and giving atom-transit edges to
non-molecules that carry no atom. It was never any tier of the honest reference lane.

This run replaces it with tier 4. Measured on the carbon graph:

| pair | shared edges | Jaccard | turnover |
|---|---|---|---|
| star vs tier 4 | 172,471 | 0.808 | **19.2%** |
| tier 3 vs tier 4 | — | 0.997 | 0.3% |

So the change is roughly 40x more "stop using the star" than "adopt tier 4". Anyone
reading a shift here as the value of the MetaCyc curated increment would be
misattributing it.

## The headline: the benchmark cannot see the difference

Directional specificity, 925 scored conditions. Star universe -> tier 4:

| facet | mode | star ALL | tier 4 ALL | delta |
|---|---|---|---|---|
| netA_iML1515 | signed | 0.8307 | 0.8268 | −0.0039 |
| netA_iECDH10B | signed | 0.8289 | 0.8251 | −0.0038 |
| netB_iML1515 | signed | 0.8492 | 0.8455 | −0.0037 |
| netB_iECDH10B | signed | 0.8452 | 0.8419 | −0.0033 |
| netA_iML1515 | unsigned | 0.8724 | 0.8709 | −0.0015 |
| netA_iECDH10B | unsigned | 0.8705 | 0.8692 | −0.0013 |
| netB_iML1515 | unsigned | 0.8919 | 0.8930 | +0.0012 |
| netB_iECDH10B | unsigned | 0.8904 | 0.8904 | +0.0000 |

Across all 40 (facet x element x mode) slices:

- mean |delta| = **0.0040**, against a mean CI half-width of **0.0222**
- ratio |delta| / half-CI = **0.20**
- the tier-4 AUC falls inside the star run's 95% CI in **40 of 40 slices**

**Replacing 19% of the carbon edge set moved the benchmark by one fifth of its own
noise floor.** That is the result, and it is a negative one worth more than a
confirmation would have been.

Two readings, both true and neither optional:

1. **The published 0.83–0.92 is robust.** It was not an artifact of the star's
   fabricated transits. A reader who feared the headline number rested on a
   discredited graph can stop fearing it.
2. **This benchmark does not validate universe quality.** It cannot discriminate a
   graph carrying ~24k fabricated carbon edges from one without them. Whatever it is
   measuring — broadly, that conductance to a condition's own target outranks its
   off-targets — survives a fifth of the topology being replaced. So it must not be
   cited as evidence that the atom-mapping work improved the model. It is evidence
   that the model's ranking behaviour is insensitive to it.

The signed deltas are consistently, very slightly negative (19 of 20 slices down,
mean −0.0043). Tier 4 is not measurably better on this benchmark and is arguably a
hair worse. That does not make the star better science — the star fabricates transits
that do not exist, and a benchmark that cannot see the difference is not a licence to
keep them. It does mean the case for tier 4 rests on the chemistry, not on this AUC.

## Essentiality diagnostic

Carbon, 131 essential vs 32 non-essential genes. Effectively unmoved:

| facet | star | tier 4 |
|---|---|---|
| netA_iML1515 | 0.7560 | 0.7512 |
| netA_iECDH10B | 0.7615 | 0.7609 |
| netB_iML1515 | 0.7765 | 0.7767 |
| netB_iECDH10B | 0.7653 | 0.7655 |

netB still beats netA on every element in both modes, as before.

## Why the comparison is trustworthy

Re-freezing the answer key under a new universe could have made this delta
meaningless: the panel's anchors are re-gated live against the base graphs, so a
universe swap can silently change which conditions are scoreable. Measured instead of
assumed — star key vs tier-4 key:

- anchors 49 -> 49, **zero liveness flips**, all live on all four facets
- panel edges 3,865 -> 3,865, identical `edge_id` set
- conditions 974 -> 974, identical `condition_id` set
- expectations 24,229 -> 24,229
- coverage 925/974 both, and the **covered sets are identical**

The scored condition set does not move, so the AUC delta is attributable to the
universe and not to the panel shifting underneath it. The identical 5,690,248 output
row count is independent corroboration.

`50_audit_v3.py` passes, including GATE 4 (answer-key precedence — the key freezes
before any solve output exists). `check_canon.py` and `check_provenance.py` green.

## The scorer is ~60x faster, and provably unchanged

Scoring took ~15 min per mode. It now takes **~15 seconds**.

The cost was not the data volume — reading the 423 MB table takes 2.0 s. It was
`cluster_bootstrap`: 2,000 resamples x ~900 conditions x a sort-based
`mannwhitney_u` each, repeated for all 20 (facet, element) slices — roughly 36M U
computations per mode. But AUC_S is a ratio of sums over strata, and the bootstrap
resamples *which conditions are included*, never the data inside a stratum. So `U_j`
and `|P_j||N_j|` are invariant across all 2,000 iterations and were being recomputed
every time. Precomputing them once per condition makes each iteration two sums.

A second fix replaced a per-cell Python loop in `directional_specificity` (~5.7M dict
lookups per slice) with a vectorised merge. That one was real but an order of
magnitude smaller — it is recorded because the first diagnosis was wrong, and finding
the per-cell loop first is the mistake a reader is likely to repeat.

Equivalence is verified, not asserted:

- the rewritten `cluster_bootstrap` reproduces the reference implementation's `(lo, hi)`
  **exactly** — same float64 bits — on both a small and a realistic (902-condition)
  case, including degenerate strata (empty side, multi-stratum condition). The single
  `(n_boot, n)` draw consumes the same RNG stream in the same C-order as `n_boot`
  successive draws, so the resamples are identical rather than merely equivalent in
  distribution.
- re-scoring the same tier-4 observations with the optimised code produced
  `specificity_signed.tsv`, `essentiality_signed.tsv` and `potency_signed.tsv`
  **byte-identical** to the original-code run.
- `--selftest` passes: rank-invariance under x^3, 5x, arcsinh, x|x| (signed) and
  exp|x|, log|x|+1 (unsigned); perfect-oracle AUC 1.0; the no-null path returns
  NaN/undefined; the planted-signal permutation test recovers.

Both universes in this report were scored with the *same* optimised code, so the
comparison is internally consistent regardless.

Caveat: `_panelstats.py` documents itself as a VERBATIM copy of
`fig-validate/main/ecspr/validation/_panelstats.py`. It is no longer verbatim. The
sibling should take this change or the claim should be dropped.

## Potency: still undefined, by decision

Every potency cell reports `status=undefined:no_null_pool`. V3 ships no enumerated
null, so a vs-random claim has no reference distribution and the scorer says so
instead of emitting a number. This is also why the benchmark could run before the
frozen null was regenerated — it consumes no null.

## Caveats a reader needs

- **The frozen null is currently mismatched with this universe.** It was built on the
  old graph. The benchmark does not consume it, so this run is unaffected, but
  significance scoring is not licensed until the null is regenerated on tier 4.
  **No gate refuses this pair.** This caveat previously said
  `canon.assert_canonical_reference()` "refuses the pair by design" — corrected
  2026-07-20, it does not. That function pins `reac_prop.tsv` and
  `direction.parquet` by sha256, and a tier move changes neither (it moves
  `atom_pairs.parquet` and the graph weights), so a tier-4 solve scored against
  an old-graph null passes it. `check_canon`'s frozen-null check is PRESENCE-only.
  Solve/null agreement holds only **by construction** — both built from the same
  base graphs — and has to be checked by hand.
- **The 341-row target-resolution table still has no human review**, unchanged from
  the previous run. 124 component rows resolved to a metabolite, 217 carry no target
  cell, 5 ambiguous strings were refused rather than guessed. A review that moves rows
  moves these numbers.
- **The answer key is set-stable but row-order unstable** between builds, so "frozen"
  means frozen by content, not by bytes. A sort key would fix it.
- **Image reference**: this run executed against
  `quay.io/hallamlab/external_ecspr:2026.07.14`. The rename to `hallamlab/ecspr` had
  left a PRIVATE repo, so the first attempt refused at pre-pull with a 401 — correctly,
  since a pull scheduled on a compute node with no outbound internet would have become
  a silent green run with empty outputs under `errorStrategy='ignore'`. `hallamlab/ecspr`
  was made public mid-session and both references now resolve to the identical digest
  `sha256:97f7afad...`, verified against the registry.
- The star-universe baseline is archived at
  `transforms/build/benchmark/v3_build/run_star_universe_baseline/` and its answer key
  at `reports/benchmark_v3_star_universe_key/`, both independent of the live tree, so
  this comparison is reproducible without re-running sockeye.
