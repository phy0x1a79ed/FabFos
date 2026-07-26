"""Compile-check for the full FabFos fosmid-recovery pipeline.

Plans the whole intended chain from raw pooled-clone reads through to ORFs:

    reads (per pool) -> host filter -> assembly (megahit + spades)
      -> junction split (per pool) -> cluster contigs (all pools)
      -> putative inserts -> ORF prediction

and asserts the resolved plan contains every stage. Two read pools are supplied
so both aggregations are exercised: `junction_split` runs per pool
(`group_by=meta`), while `cluster_contigs` groups on the experiment
(`group_by=exp`) and must see *both* pools' pieces in one job.

Two lineage mechanisms are what make this resolve without editing any existing
transform:

* **host-filter coercion** -- `junction_split` pins its consumed `assembly` to a
  `host_filtered_short_reads` ancestor, so the planner is forced to insert
  `background_filter` upstream of the assemblers (which accept host-filtered
  reads because that type extends `clean_short_reads`). This pin lives in
  whichever transform consumes the assemblies first; that used to be
  `cluster_contigs`, and moving it is what keeps host filtering in the plan.
* **ORF reuse** -- `fabfos::putative_inserts` carries the `origin: assembler`
  property, making it a property-superset of `sequences::assembly`, so the
  existing `prodigal` calls ORFs on the final inserts unchanged. The `orfs`
  target is pinned to the inserts so prodigal runs on them, not the raw assembly.

Renders the resolved DAG to `tests/artifacts/fosmids_dag.svg`.

Requirements: an env with `metasmith` + graphviz `dot` (e.g. `msm`):

    PATH="/home/tony/lib/miniforge3/envs/msm/bin:$PATH" python -m pytest tests/test_fosmids_workflow.py -v

Runs standalone too: `python tests/test_fosmids_workflow.py`.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from metasmith.python_api import (
    Agent,
    Runtime,
    Source,
    DataInstanceLibrary,
    TransformInstanceLibrary,
    TargetBuilder,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
MLIB = REPO_ROOT / "src" / "metasmith_libraries"
ARTIFACTS = Path(__file__).resolve().parent / "artifacts"

N_POOLS = 2

# Every stage of the recovery pipeline, keyed by the transform-file stem the
# planner must include to reach ORFs on the final inserts.
EXPECTED_TRANSFORMS = {
    "seqkit_reads",      # read QC (feeds bbduk)
    "bbduk",             # read cleaning
    "background_filter", # host removal (coerced in by junction_split's lineage)
    "megahit",           # assembler 1
    "spades",            # assembler 2
    "junction_split",    # vector-backbone junctions -> cut at them
    "cluster_contigs",   # cross-assembler dedup
    "prodigal",          # ORF prediction on the final inserts
}


def _plan_fosmids(work: Path):
    """Resolve the full recovery pipeline for two read pools. Planning only."""
    inputs = DataInstanceLibrary(work / "inputs.xgdb")
    for tl in ("sequences.yml", "fabfos.yml"):
        inputs.AddTypeLibrary(MLIB / "data_types" / tl)

    # one recovery_experiment groups the whole run; each pool's read_metadata
    # (and thus its reads/assembly) descends from it, so the cross-pool
    # clustering/coverage jobs (group_by=experiment) see all pools at once.
    exp = inputs.AddValue(
        "recovery_experiment.txt", "fabfos_demo", "fabfos::experiment"
    )
    for i in range(N_POOLS):
        meta = inputs.AddValue(
            f"read_metadata_{i}.json",
            {"parity": "paired", "length_class": "short"},
            "sequences::read_metadata",
            parents={exp},
        )
        reads = work / f"pool_{i}.fq.gz"
        reads.touch()
        inputs.AddItem(reads, "sequences::short_reads_pe", parents={meta})

    # per-run resource: the host/vector reference for background_filter
    host = work / "host.fna"
    host.touch()
    inputs.AddItem(host, "sequences::background_genome")

    # per-run reference: the pCC1 backbone junction_split blasts the contigs
    # against. Pinned to the experiment so a plan cannot silently reach for a
    # different run's vector -- the same reason the host pin hangs off `exp`.
    vector = work / "pCC1.fna"
    vector.touch()
    inputs.AddItem(vector, "fabfos::vector_backbone", parents={exp})
    inputs.Save()

    resources = [
        DataInstanceLibrary.Load(MLIB / "resources" / "env"),
        DataInstanceLibrary.Load(MLIB / "resources" / "lib"),
        inputs,
    ]
    transforms = [
        TransformInstanceLibrary.Load(MLIB / "transforms" / d)
        for d in ("assembly", "fabfos", "metagenomics")
    ]

    # targets: final inserts + report, and ORFs pinned to the inserts so
    # prodigal runs on them (not the raw assembly).
    targets = TargetBuilder()
    ins = targets.Add("fabfos::putative_inserts")
    targets.Add("fabfos::putative_insert_report")
    targets.Add("sequences::orfs", parents={ins})

    agent = Agent(home=Source.FromLocal(work / "agent_home"), runtime=Runtime.APPTAINER)
    task = agent.GenerateWorkflow(
        samples=list(inputs.AsSamples("fabfos::experiment")),
        resources=resources,
        transforms=transforms,
        targets=targets,
    )
    return task


def _step_transform_names(task) -> set[str]:
    return {Path(step.transform._path).stem for step in task.plan.steps}


def _render_dag(task, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    task.plan.RenderDAG(out, blacklist_namespaces={"lib", "env"})
    return out


def test_fosmids_pipeline_compiles(tmp_path):
    """The full recovery pipeline plans and wires in every expected stage."""
    task = _plan_fosmids(tmp_path)

    assert task.ok, f"fosmids workflow failed to plan: {task.plan}"

    names = _step_transform_names(task)
    missing = EXPECTED_TRANSFORMS - names
    assert not missing, (
        f"pipeline plan is missing expected stage(s): {sorted(missing)}; "
        f"plan used: {sorted(names)}"
    )

    svg = _render_dag(task, ARTIFACTS / "fosmids_dag.svg")
    assert svg.exists() and svg.stat().st_size > 0, f"DAG SVG not written: {svg}"


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        task = _plan_fosmids(Path(td))
        if not task.ok:
            raise SystemExit(f"FAILED to plan: {task.plan}")
        names = _step_transform_names(task)
        print(f"Plan OK -- {len(task.plan.steps)} steps")
        for step in sorted(task.plan.steps, key=lambda s: s.order):
            prods = [i.dtype_name for g in step.produces for i in g]
            print(f"  step {step.order}: {Path(step.transform._path).stem} -> {prods}")
        missing = EXPECTED_TRANSFORMS - names
        print(f"missing stages: {sorted(missing) if missing else 'none'}")
        svg = _render_dag(task, ARTIFACTS / "fosmids_dag.svg")
        print(f"DAG written to: {svg}")
