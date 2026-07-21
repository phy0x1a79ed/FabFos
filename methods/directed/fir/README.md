# Directed regeneration — fir CPU deployment

The path-parametrised, canon-free drivers that run the directed observed solve + directed
null on fir CPU nodes. They are the cluster-deployment siblings of the canon-integrated
`../build_directed_observed.py` and `../build_directed_null.py`: same engine primitives
(`orient_and_weight` / `build_ar2m_directed` / `directed_solve_element`), same fork-per-core
parallelisation and output schema, but every path is a `--stage-dir` argument and the inputs
are pre-resolved to pickles (draws→weights, direction→ratios) so the venv needs no pyarrow.
This canon-free shape is also the metasmith-transform precursor: the transform receives these
same inputs as staged parents.

## Staging (on fir, `/scratch/$USER/scadc_dir/`)
- `lib/` — the engine (`ecspr_directed.py` **must be the post-O1 build**, MMD ordering),
  `ecspr_network.py`, `ecspr_solver.py`.
- `ref/solve/base_{C,N,S,P}.pkl`, `ref/graph/mnx_bipartite_{X}.pkl` — the reference graph.
- `draws_resolved.pkl` (per draw-size: `keys` + `draw_w`), `ratios.pkl`, `reac_prop.tsv`,
  `axes.json`, `axes_testable.json`, `fosmid_addition_weights.pkl`.

## Run
- `sbatch run_null.sbatch` — array 0-19 (`manifest.txt` = element × draw-size), one 192-core
  node per shard, `fir_null.py --workers $SLURM_CPUS_PER_TASK --resume`. Output:
  `out/{reff,ieff}_{X}_N{n}.tsv`, concatenated per draw-size across elements into the
  canonical `{lane}_null_directed_N{n}.tsv`.
- `sbatch run_obs.sbatch` — one 192-core node, `fir_observed.py` (all elements) →
  `obs/{reff,ieff}_axes_report.tsv` (the directed `reference_axes_report`).

## Notes
- No `--partition` (fir's submit wrapper auto-routes; a full-node core count lands on the
  bynode partitions). Batch context exposes the full node cpuset (`nproc=192`).
- The venv has no `sksparse`, so every solve takes the engine's `_reg_spsolve` path — which is
  exactly where the O1 MMD ordering applies, so the speedup is realised on every solve here.
- `--resume` skips a shard whose row count already matches; delete `out/` to force a clean
  rebuild (required whenever the engine changes, e.g. the pre-O1 → post-O1 re-run).
