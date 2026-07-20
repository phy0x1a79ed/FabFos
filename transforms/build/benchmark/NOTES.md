# benchmark v3 answer-key build — working notes

State and decisions for the v3 answer key. Written as the build proceeds so the
reasoning survives the session that produced it.

## Where v3 lives

One self-contained tree, `validation/benchmark/v3/` in the data library
(`.awm/data/ref/`), 55 files / 15 MB, reachable through canon:

| canon symbol | holds |
|---|---|
| `BENCH_V3_OBSERVATIONS` | the 548 measured perturbations + curation prose |
| `BENCH_V3_DECISIONS` | host lineage, gene→reaction overrides, target resolution |
| `BENCH_V3_X` | the frozen X — four network builds, carried from v1 unchanged |
| `BENCH_V3_CONTRACT` | v1's Y, staged as the contract SHAPE |
| `BENCH_V3_GROUND_TRUTH` | the upstream tables the key is derived FROM |
| `BENCH_V3_BASELINE` | v1's scored reference numbers |
| `BENCH_V3_V1_PROVENANCE` | v1's manifest/audit/README, nested under `v1/` |
| `BENCH_V3_Y` | the answer key — **raises CanonError until built** |

## Measured facts (this build, not inherited)

- **548 observations**: 382 GOF + 166 LOF, across 58 distinct hosts.
- GOF perturbation split: 127 add-only, 177 add-and-delete, 78 delete-only.
- **Direction comes from `measured`, not `expected`.** `expected` is populated on
  only 10 of 382 GOF rows. `measured` is **322 up / 3 down / 57 unknown**, so
  **491 of 548** observations carry a usable direction (325 GOF + 166 LOF).
- Nearly all sign contrast is the 166 uniformly-negative LOF conditions plus 3
  GOF rows. The GOF arm is therefore close to a pure "does this raise
  conductance toward its own target more than elsewhere" test. This belongs in
  the result, not in a footnote.

## Why the mechanical path carries the GOF arm

A condition's target is derived from the reactions it perturbs, not from the
free-text `target` string:

- 766 distinct MNXR are referenced; **697 (91%) are present in X**.
- **380 of 382** GOF rows have at least one perturbed reaction in X.
- 299 rows have an `add_mnxr` in X, so the "direct non-currency product of the
  inserted reaction" rule (v1's rule) applies directly.
- LOF: 138/166 rows have a `del_mnxr` in X; all 166 resolve through
  `pheno_edges.tsv`, which is authoritative for that arm and is what v1 used.

## Why free-text target coverage is low, and why that is honest

`target_resolution.tsv` resolves 274 raw strings → 341 components, of which 124
land on a metabolite in X. The residual is **not** mostly a name-matching
failure:

- **27** components resolve in MetaNetX but are **outside X** — catechin,
  pinocembrin, lycopene, 6-deoxyerythronolide B, GFP. These are heterologous
  products of engineered strains. ECSPr cannot score conductance toward a
  metabolite the network does not contain, so they correctly carry no target
  cell.
- A further group are **compound classes, not metabolites** — "Carotenoids",
  "Flavonoids", "Biodiesel", "Carboxylic acids" — and one is a protein (GFP).
  There is no single id to target and inventing one would be a fabricated
  answer.
- The remainder is a tail of spelling and naming-convention variants
  ("4-courmaric acid"; `1_2-Propanediol`, where the source escaped a comma as an
  underscore and MetaNetX spells it `propane-1,2-diol`).

This is a real property of a LASER-derived GOF set read against an *E. coli*
network, not a defect in the resolver. It is reported as coverage, never
silently dropped, and the mechanical path above is what keeps the arm scorable.

Ambiguity is refused rather than guessed: a folded name matching more than one
X metabolite (e.g. `L-lactate`) resolves to nothing and records why. Picking
the first would be a wrong target that scores as though it were right.

## Ordering constraint discovered during the build

The panel's anchor **liveness gate** re-measures every anchor in each facet's
largest connected component, which requires the per-facet base graphs. Those
are built from X by the `ecsprBenchmark` domain (T4). So the order is

    X → base graphs (T4) → panel → Y → freeze → solve

Base-graph construction is *network construction*, not ECSPr output, so using it
here does not violate answer-key independence — the audit's forbidden-token list
already draws that line correctly, and it must stay drawn.
