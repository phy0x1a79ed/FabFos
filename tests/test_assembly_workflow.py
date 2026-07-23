"""Compile-check for the FabFos assembly stage.

The intended FabFos pipeline assembles each pooled-clone read set with **both**
megahit and spades before the contigs are clustered downstream. This test asks
the metasmith planner to resolve a workflow that produces both assemblies from
raw short reads and asserts:

1. the plan is valid (``task.ok``), and
2. every transform the assembly stage is expected to use is present in the
   resolved plan — read QC (``seqkit_reads``), read cleaning (``bbduk``), and
   both assemblers (``megahit`` and ``spades``).

As a side effect it renders the resolved DAG to ``tests/artifacts/assembly_dag.svg``
so the wiring can be eyeballed.

Requirements: an environment with the ``metasmith`` package **and** graphviz's
``dot`` on ``PATH`` — e.g. the ``msm_fresh`` conda env:

    /home/tony/lib/miniforge3/envs/msm_fresh/bin/python -m pytest tests/ -v

The file also runs standalone (``python tests/test_assembly_workflow.py``),
which just prints the plan and writes the SVG.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from metasmith.python_api import (
    Agent,
    ContainerRuntime,
    Source,
    DataInstanceLibrary,
    TransformInstanceLibrary,
    TargetBuilder,
)

# repo_root/tests/this_file  ->  repo_root
REPO_ROOT = Path(__file__).resolve().parent.parent
MLIB = REPO_ROOT / "src" / "metasmith_libraries"
ARTIFACTS = Path(__file__).resolve().parent / "artifacts"

# The assembly stage is expected to route raw short reads through read QC and
# cleaning into *both* assemblers. Each name is the stem of a transform file
# under transforms/assembly/.
EXPECTED_TRANSFORMS = {"seqkit_reads", "bbduk", "megahit", "spades"}


def _plan_assembly(work: Path):
    """Resolve a plan that produces both the megahit and spades assemblies.

    Planning only — nothing is staged or executed, so mock (empty) read files
    are sufficient.
    """
    # ── inputs: one pooled-clone short-read sample ────────────────────────────
    inputs = DataInstanceLibrary(work / "inputs.xgdb")
    inputs.AddTypeLibrary(MLIB / "data_types" / "sequences.yml")

    reads = work / "reads.fq.gz"
    reads.touch()
    meta = inputs.AddValue(
        "reads_metadata.json",
        {"parity": "paired", "length_class": "short"},
        "sequences::read_metadata",
    )
    inputs.AddItem(reads, "sequences::short_reads_pe", parents={meta})
    inputs.Save()

    # ── resources: tool envs + the sample library itself ──────────────────────
    # (the env-migration renamed resources/containers -> resources/env and the
    #  containers:: namespace -> env::)
    containers = DataInstanceLibrary.Load(MLIB / "resources" / "env")

    # ── transforms: the assembly domain ───────────────────────────────────────
    transforms = [TransformInstanceLibrary.Load(MLIB / "transforms" / "assembly")]

    # ── targets: force both assemblers to be resolved ─────────────────────────
    targets = TargetBuilder()
    targets.Add("sequences::megahit_assembly")
    targets.Add("sequences::spades_assembly")

    agent = Agent(home=Source.FromLocal(work / "agent_home"), runtime=ContainerRuntime.APPTAINER)
    task = agent.GenerateWorkflow(
        samples=list(inputs.AsSamples("sequences::read_metadata")),
        resources=[containers, inputs],
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


def test_assembly_stage_compiles(tmp_path):
    """The assembly stage plans cleanly and wires in every expected transform."""
    task = _plan_assembly(tmp_path)

    assert task.ok, f"assembly workflow failed to plan: {task.plan}"

    names = _step_transform_names(task)
    missing = EXPECTED_TRANSFORMS - names
    assert not missing, (
        f"assembly plan is missing expected transform(s): {sorted(missing)}; "
        f"plan used: {sorted(names)}"
    )

    # Render for inspection; the SVG is a deliverable of this test, so a failure
    # to render (e.g. graphviz `dot` missing) should fail the test loudly.
    svg = _render_dag(task, ARTIFACTS / "assembly_dag.svg")
    assert svg.exists() and svg.stat().st_size > 0, f"DAG SVG not written: {svg}"


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        task = _plan_assembly(Path(td))
        if not task.ok:
            raise SystemExit(f"FAILED to plan: {task.plan}")
        names = _step_transform_names(task)
        print(f"Plan OK — {len(task.plan.steps)} steps")
        for step in sorted(task.plan.steps, key=lambda s: s.order):
            prods = [i.dtype_name for g in step.produces for i in g]
            print(f"  step {step.order}: {Path(step.transform._path).stem} -> {prods}")
        print(f"expected transforms present: {EXPECTED_TRANSFORMS <= names}")
        svg = _render_dag(task, ARTIFACTS / "assembly_dag.svg")
        print(f"DAG written to: {svg}")
