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
