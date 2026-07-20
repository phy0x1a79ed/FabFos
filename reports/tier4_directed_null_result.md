# The directed null, regenerated on the tier-4 universe

Closes the reference chain that the tier-4 benchmark run left half-moved. As of that
run, `mnxref.graph` sat on `graph_tier4` but everything downstream of it —
`solve`, `solve_directed`, `null_directed` — still lived on the pre-tier-4 graph. A
null built on different edges than the observed solve is not a valid reference, so
significance scoring was **not licensed** until this work landed.

Method `0.3.1+7b3976c` (the `+hash` moved from the benchmark run's `+37a1d54` exactly as
that report predicted: the data-library index is one of the seven hashed components, and
re-staging the null moved it).

## What was rebuilt, and on what

The whole chain was regenerated on one graph — `graph_tier4` — and on the **coherent**
evidence/weights pair, not the frozen mismatched pair the 17-Jul solve reproduced from
(see the weights-fix note below). Three new sibling directories under
`ecspr_reference/mnxref-4_5/`, staged by hardlink into the library, leaving the old
`solve`/`solve_directed`/`null_directed` untouched:

| item | new src | contents |
|---|---|---|
| `mnxref.solve` | `solve_tier4/` | `base_{C,N,S,P}.pkl` + `axes_testable.json` + undirected `{ieff,reff}_axes_report.tsv` (the parity gates' symmetric-limit referent) |
| `mnxref.solve_directed` | `solve_directed_tier4/` | directed `{ieff,reff}_axes_report.tsv`, 15,621 rows each |
| `mnxref.null_directed` | `null_directed_tier4/` | 10 files, `{ieff,reff}_null_directed_N{14,25,30,35,43}.tsv`, 316,000 rows each |

## The observed directed solve

Local, `p312`, softplus-diode directed solve, ~82 min serial. Verified by row count and
per-element census, not exit status:

| element | cells | wall | s/cell | LCC | directed edges | fosmids w/ additions |
|---|---|---|---|---|---|---|
| C | 7,960 | 3,011.7 s | 0.378 | 12,884 | 11,152/24,565 | 199/199 |
| N | 4,577 | 1,445.7 s | 0.316 | 9,085 | 7,610/16,210 | 199/199 |
| S | 1,302 | 27.2 s | 0.021 | 2,838 | 2,079/4,198 | **186**/199 |
| P | 1,782 | 489.7 s | 0.275 | 7,430 | 6,897/13,980 | **198**/199 |

**S scores over 186 of 199 fosmids and P over 198**, where C and N are 199/199. Not a
defect — not every fosmid carries sulfur- or phosphorus-bearing reactions — but it means
S's and P's per-axis denominators differ from C/N's, and any per-fosmid summary must say so.

## The null sweep

Sockeye, apptainer (`ecspr:2026.07.14`), job `12325885`, `--array=0-19`, 32 workers/shard,
`--mem=96G`. **All 20 shards COMPLETED** — 47/47 sacct entries, zero FAILED/OOM/TIMEOUT.

Per-shard row count = axes(X) × 4,000 draws, audited on the cluster and again locally after
retrieval, every shard exact:

| element | axes | rows/shard | × 5 sizes |
|---|---|---|---|
| C | 40 | 160,000 | 800,000 |
| N | 23 | 92,000 | 460,000 |
| S | 7 | 28,000 | 140,000 |
| P | 9 | 36,000 | 180,000 |
| **total** | **79** | | **1,580,000 per lane** |

Cluster and local totals both **3,160,040 lines = 3,160,000 rows + 40 headers = 2 lanes ×
1,580,000**.

**N was the slow lane, and it is not a fault.** N ran at 104–113 ms/cell against C's
25–46 ms, despite a smaller graph (LCC 9,085 vs 12,884) — where LCC-scaling predicts N
should be *faster*. Chased to the compute node rather than assumed: `ssh se231` showed the
fork-pool workers pinned at 99.5% CPU, i.e. computing, not deadlocked. N's directed Newton
solve is simply harder per cell (more iterations on N's structure). The 24 h walltime made
the slowness harmless.

## The load-bearing check: 316,000, not 312,000

The merge from 20 per-element shards to the 10 canonical per-size files
(`transforms/run/null_directed/merge_null.py`) **asserts** the row count rather than
reporting it, because three separate mechanisms on this path narrow the null basis silently
(`--resume` is row-count-only; `discover_draw_sizes` enumerates rather than checks; the
scorer's `.get()` + `continue` drops missing keys). All 10 merged files landed at exactly
**316,000 = 79 axes × 4,000**.

That number is itself the test. The **previous** canonical null held 312,000 = **78** axes ×
4,000. The +4,000 is the carbon axis the stale pre-CLEAN evidence weights had been silently
zeroing (C went 39/40 → 40/40 testable when the weights were repointed at the CLEAN
199-fosmid basis). So 316,000 is a positive confirmation that the coherent weights took — the
same conclusion the base-graph rebuild reached, arriving here by a completely independent
route (the sweep never touches the base rebuild). A merged file at 312,000 would have meant
the weights fix did not propagate and the whole chain was suspect. It did not happen.

## Gates

- `check_canon.py` — **green**. Includes "10 frozen null files present — missing: []" and
  the canon↔provenance hash agreements.
- `check_provenance.py` — **green**. reac_prop sha256 agrees across the record and
  `MANIFEST.json`; all 53 present items hashed.
- All 10 `canon.FROZEN_NULL_FILES` resolve through the library to files at 316,000 rows;
  both `ieff` and `reff` carry the full grid `{14,25,30,35,43}`.

## Solve/null agreement holds BY CONSTRUCTION, not by a gate

The observed solve and the null were built from the **same** base graphs in one pass, which
is the only thing that makes them a valid pair. Nothing mechanical enforces it:
`canon.assert_canonical_reference()` pins the *reference* inputs (`reac_prop.tsv`,
`direction.parquet`) by sha256, and a tier move changes neither — it moves
`atom_pairs.parquet` and the graph weights. `check_canon`'s frozen-null check is
presence-only. So a mismatched solve/null pair would pass both gates cleanly; agreement here
rests on the build discipline, not on being stopped. (This corrects the earlier false
"refuses the pair by design" claim; committed `07aba73`.)

## Caveats a reader needs

- **The frozen manifest's `aam_source` still names the tier-1 atom-pairs table**, and the
  closure ledger's refused-bucket census is ~28% wrong under tier 4 (7,119 of 25,436 refused
  reactions now have atom pairs). `assert_reference` checks that table's columns, never its
  hash or census, so the gate certifies a stale manifest. Re-running the ledger is out of this
  work's scope — filed, not fixed.
- **The `benchmark.v3.observations` source has accumulated a stray `bench_run_1784548824/`
  subdir** (the tier-4 benchmark run's staged X inputs, 53 files / 36 MB, leaked into the
  declared observations dataset). Left the committed record describing the clean 14-file
  dataset; the leak is a benchmark-run hygiene item, filed separately.
- **`validation.dual_network`'s declared `src` was stale** (`ecspr_validation`; the dataset
  had moved under `ecspr_benchmark/`). Repointed to `ecspr_benchmark/ecspr_validation` so the
  gate resolves — a pre-existing bug discovered while running the T5 gates, not part of the
  null regeneration.
- **K stays at 1000.** K=10000 buys resolution in the null tail, not correctness, and the
  benchmark finding says this instrument is insensitive to universe quality; raising K also
  needs new draws written to a distinctly-named path (the draws filename carries no K suffix,
  a known silent-overwrite trap).

## Provenance

Sweep job `12325885`; observed solve local. New source dirs under
`/home/tony/agentic_workspace/data/scadc/ecspr_reference/mnxref-4_5/{solve,solve_directed,null_directed}_tier4`.
Merged null and retrieved shards under the job scratch. The old `solve`/`solve_directed`/
`null_directed` dirs are preserved as the historical pre-tier-4 record.
