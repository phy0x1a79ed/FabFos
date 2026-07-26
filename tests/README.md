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

## `build_references_stage{1,2}_*.py` — not tests

These two are **drivers**, not compile-checks, and they are the exception to
everything above: they can execute, they read real data, and they write to
`data/`. They live here because they are the other half of
`examples/build_references_dag.py` and share its target list, and because
plan-only is their default — run either with no flags and it prints the plan,
renders a DAG to `artifacts/` and stops, exactly like the tests around it.

The reference build splits at the tier boundary:

| | | |
|---|---|---|
| `build_references_stage1_acquire.py` | fills `data/raw/` | network-bound; 11 acquisitions, 6 on a machine already holding the DVC-pinned bulk |
| `build_references_stage2_curate.py`  | fills `data/reference/`, `data/benchmark/` | compute-bound; 17 steps, no network |

Stage 2 also loads `acquire/` and then asserts none of it is scheduled. That
looks backwards — the usual way to enforce a split is to withhold the other
half's transforms — and it is the point: withholding them turns an unmet input
into an unresolvable plan whose error names a *type*, leaving you to guess which
fetch should have covered it. Loading them names the acquisition instead.

`--standins` on stage 2 fills anything stage 1 has not produced with empty files,
so its curation DAG renders before stage 1 has ever run. Planning never reads
bytes; running would, so that flag refuses to combine with `--run`.

Both need the pinned engine — `src/metasmith` on `feat/fabfos`, which carries the
mamba executor. They put it ahead of whatever is installed in the env and fail
with an explanation if `Runtime.MAMBA` is missing, because three of this build's
envs are conda-only and metasmith picks `conda:` vs `container:` off one global
runtime.

### The `msm-fabfos` agent env

`Agent.container` is read as a *conda environment name* under MAMBA, so the agent
needs an env holding the pinned metasmith — and each step resolves its own env the
same way, from the `conda:` key of its `env::*.env` declaration. Nothing creates
these on demand. Build all of them with:

```bash
bash envs/setup_agent_env.sh
```

It is idempotent; re-run it after a `git pull` to repoint the hook at the pinned
engine. Two parts of it are not obvious and are why this is a script rather than
two commands:

- **`base.yml` never installs metasmith.** It is metasmith's *dependency* env and
  carries a standing comment that it can hold no pip dependencies, so the env has
  no `metasmith` console script. The MAMBA launcher runs `mamba run -n <env>
  metasmith api …`, so every step's bootstrap would die on `metasmith: not found`.
  The script installs `metasmith`/`msm` shims (`exec python -m metasmith "$@"`)
  rather than `pip install -e`, which would write `*.egg-info` into the pinned
  submodule and dirty it.
- **The activation hook is not optional and not a `.pth`.** This machine exports
  `PYTHONPATH=~/lib/locals`, whose `metasmith` symlink points into the metasmith
  `dev` worktree — the one *without* the mamba executor — and PYTHONPATH is placed
  ahead of site-packages, so only an activation-time prepend wins. The script
  verifies afterwards that `import metasmith` resolves inside `src/metasmith`.

### Staging under MAMBA — fixed

This used to be a hard blocker: `Deploy` worked and `StageWorkflow` did not, because
the staged-agent API read `AgentPaths.HOME_ROOT` (`/msm_home`) and `WORK_ROOT`
(`/ws`) unconditionally, and both are *container bind* paths that the relay-free
MAMBA launcher never creates.

Fixed upstream in the engine (`metasmith` `feat/mamba-executor`, merged here):
those two roots now resolve from `METASMITH_HOME_ROOT` / `METASMITH_WORK_ROOT`,
which the relay-free `msm` / `msm_bootstrap` scripts export with the real host
paths. Container runtimes are byte-identical — the defaults are unchanged, and
separate immutable `CONTAINER_*` literals are used for shell text that runs
*inside* a container and must not follow a host override.

Run directories go under `data/scratch/`, which is gitignored and safe to delete.

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
