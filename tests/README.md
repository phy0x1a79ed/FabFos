# FabFos tests

Compile-checks for the FabFos pipeline as wired in `src/metasmith_libraries`.
These are **planning-only**: they ask the metasmith planner to resolve a
workflow and assert its shape. Nothing is staged, containerised, or executed, so
they run in a couple of seconds and need no test data.

## `test_assembly_workflow.py`

Verifies the **assembly stage** compiles: from a raw pooled-clone short-read
sample the planner must resolve a workflow that produces *both* the megahit and
the spades assembly, routing through read QC (`seqkit_reads`) and read cleaning
(`bbduk`). Asserts the plan is valid and contains every expected transform, then
renders the resolved DAG to `tests/artifacts/assembly_dag.svg` for inspection.

## `test_fosmids_workflow.py`

Verifies the **full fosmid-recovery pipeline** compiles: from raw pooled-clone
reads (two pools, so the cross-pool aggregation is exercised) through host
filtering, both assemblers, contig clustering, pool coverage, chimera splitting
and coverage trimming to ORF prediction on the final inserts. Renders
`tests/artifacts/fosmids_dag.svg`.

## `test_gpr_workflow.py`

Verifies the **annotation → GPR-table stage** compiles: from one ORF FASTA the
planner must resolve all seven functional-annotation run tools (kofamscan, clean,
deepec, ezpred, diamond_uniref50, proteinbert, esm_c) feeding both GPR mappers —
`gpr_4lane` (kofam + CLEAN + uniref50 + ProteinBERT) and `gpr_7lane` (all seven).
EZpred consumes the ESM-C embeddings rather than re-embedding, so the plan orders
`esm_c → ezpred`. Both mapper targets are demanded on purpose: their outputs are
distinct subtypes of a shared `annotation::gpr_table` parent, so asking for only
one would let the planner drop the other mapper. The MNXR bridges, reference
label pool and model bundles are staged type-only stand-ins — their build
transforms are deferred, and planning never reads the bytes. Renders
`tests/artifacts/gpr_dag.svg`.

> Note: editing a `data_types/*.yml` or transform in `src/metasmith_libraries`
> requires regenerating the per-library `_metadata/` snapshots before these tests
> see the change — `./dev.sh -b` only bundles. See the submodule's build step
> (`python -m metasmith build all --types … --uniques … --transforms …`).

## Running

Needs an environment with the `metasmith` package **and** graphviz's `dot` on
`PATH`. The `msm` conda env has both (prepend its `bin/` so `dot` is found):

```bash
PATH="/home/tony/lib/miniforge3/envs/msm/bin:$PATH" python -m pytest tests/ -v
```

Each test file also runs standalone to just print the plan and (re)write its SVG:

```bash
PATH="/home/tony/lib/miniforge3/envs/msm/bin:$PATH" python tests/test_assembly_workflow.py
```
