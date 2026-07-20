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

## T3 outcome (built, gated, frozen)

`20_build_panel_v3.py` -> `21_build_Y_v3.py` -> `50_audit_v3.py`, all writing into
`transforms/build/benchmark/v3_build/` in the REPO and hardlinked into the
library by the staging script. Never written into the library directly: the
library is built from declarations, and a file written straight in is one no
record describes and no hash covers.

- panel: 49 anchors (all live on all 4 facets), 325 products, **3,865 edges**
  (v1's 175 products carried + 150 v3 adds)
- key: **974 conditions**, **24,229** sparse expectation cells (0.64% of the
  full cross; the rest are the declared default role=off_target dir=0)
- direction: 805 `+` / 169 `-`; 57 GOF observations carry no usable direction
  and are reported as coverage, never scored as null results
- one row per (observation, ELEMENT) -- the panel is element-resolved, so a
  perturbation whose targets carry both C and N is two claims

Four gates pass. Gate 4 is new in v3: **answer-key precedence**, which refuses
to freeze if solve output already exists in the tree. That is the one ordering
error no care inside the builders can catch.

### Two defects this found

1. **The LOF join matched 0 of 166.** v3 writes the gene token as `argA:del`;
   `pheno_edges.tsv` keys on `argA`. Every row fell through to the free-text
   fallback and ~half resolved by luck -- 88 of 166, which reads exactly like
   ordinary coverage loss. Gate 3 now asserts the full 166.
2. **canon resolved against an EMPTY manifest.** The fallback reader asked
   `index.yml` for a `manifest:` wrapper it does not have, so every symbol
   failed as "not in the manifest" -- pointing at the declaration instead of at
   the reader. An empty manifest is now a refusal that names the real fault.

## Base graphs are NOT derivable from X

X withholds the atom mapping deliberately (`X/_provenance.json`: "aam: WITHHELD
-- the atom mapping is under test"). Both base-graph builders are pure
inductions of an atom-mapped universe: the edge weight IS the atom-transit
count, and `load_bipartite` deletes every `w<=0` edge, so the AAM fixes
connectivity and LCC membership, not just weights. Formulas plus stoichiometry
are not a substitute.

This is the benchmark working as designed -- X withholds the AAM so an
independent implementation derives its own; the incumbent's own AAM-derived
graphs are its answer to that, and are therefore an INPUT to scoring the
incumbent. Declared as `benchmark.v3.base_graphs` (16 pkl) and
`benchmark.v3.universe` (4 pkl, 57.6 MB), the latter because
`base_plus_reactions` re-reads the universe at SOLVE time -- a tree without it
fails on the first of the 382 GOF conditions.

## T4/T5 sizing (measured, not estimated)

`resources/lib/ecspr_benchmark.py`, subcommands `solve` and `merge`.

- unit of work is a CONDITION: one graph rebuild + one factorization, then a
  cheap solve per panel edge against that shared factorization
- emits **signed delta conductance**, not a log ratio: a GOF insertion that
  creates a route has `g_base = 0`, where a log ratio is undefined and any
  sentinel is an arbitrary rank injection
- silent conditions emit ZERO, never dropped -- that is the redundancy signal
  the essentiality diagnostic reads
- parallelism is fork PROCESSES with BLAS pinned to 1 before numpy loads;
  SuperLU is serial so intra-process threading buys nothing
- **parity gate holds**: 4 workers byte-identical to `--workers 1`

Measured on the worst case (netB / C, 2,508 panel edges): **~40 s per condition
serial**. So 16 shards (4 facets x 4 elements) x 32 worker processes = 512-way
fan-out, wall clock ~10-15 min. That is the "batches of 32" shape.

fir is reachable as user `phyberos` (checked via the ssh domain's status verb,
NOT by opening a new connection -- never loop-retry connect, never delete the
ControlMaster socket; a prior run here was halted by a Duo lockout that way).
