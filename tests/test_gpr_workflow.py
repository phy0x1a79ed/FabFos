"""Compile-check for the fabfos annotation -> GPR-table stage.

From a single ORF FASTA the intended pipeline runs all seven functional-annotation
tools and folds their outputs into a gene-attributed **GPR table** two ways:

    orfs -> {kofamscan, clean, deepec, ezpred, diamond_uniref50, proteinbert, esm_c}
         -> gpr_4lane  (kofam + CLEAN + uniref50 + ProteinBERT) -> gpr_table_4lane
         -> gpr_7lane  (all seven lanes)                        -> gpr_table_7lane

This test asks the metasmith planner to resolve a workflow that produces *both*
GPR tables and asserts:

1. the plan is valid (``task.ok``), and
2. every transform the stage is expected to use is present — the seven run tools
   plus the two mappers.

Both mapper targets are added on purpose: ``gpr_table_4lane`` and
``gpr_table_7lane`` are distinct subtypes of a shared ``annotation::gpr_table``
parent, so unless *both* are demanded the planner would resolve only one producer
and drop a mapper from the DAG.

Planning is type-driven: nothing is staged, containerised, or executed, so an
empty ORF FASTA and empty type-only stand-ins for every staged ``ref::*`` (the
KO/EC/UniProt->MNXR bridges, the reference label pool, the tool databases and
model weights) are sufficient. The bridges and pool have no producer here on
purpose -- their build transforms are deferred out of scope.

As a side effect it renders the resolved DAG to ``tests/artifacts/gpr_dag.svg``.

Requirements: an env with ``metasmith`` + graphviz ``dot`` (e.g. ``msm``):

    PATH="/home/tony/lib/miniforge3/envs/msm/bin:$PATH" python -m pytest tests/test_gpr_workflow.py -v

Runs standalone too: ``python tests/test_gpr_workflow.py``.
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

REPO_ROOT = Path(__file__).resolve().parent.parent
MLIB = REPO_ROOT / "src" / "metasmith_libraries"
ARTIFACTS = Path(__file__).resolve().parent / "artifacts"

# The seven run tools + the two GPR mappers. Each name is the stem of a transform
# file under transforms/functionalAnnotation/.
EXPECTED_TRANSFORMS = {
    "kofamscan",         # KO (HMM bitscore)
    "clean",             # CLEAN contrastive EC
    "deepec",            # DeepEC EC (score-less)
    "ezpred",            # EZpred DL-only EC head
    "diamond_uniref50",  # UniRef50 homology + analytic BSR
    "proteinbert",       # ProteinBERT embeddings
    "esm_c",             # ESM-C embeddings
    "gpr_4lane",         # chosen-4 mapper -> gpr_table_4lane
    "gpr_7lane",         # full-7 mapper   -> gpr_table_7lane
}

# Staged ``ref::*`` inputs the run tools + mappers consume. Type-only stand-ins:
# planning only needs the type, not the bytes. The bridges / pool have no
# producer transform (deferred), so they must be supplied as inputs.
STAGED_REFS = [
    "ref::kofamscan_profiles",
    "ref::kofamscan_ko_list",
    "ref::uniref50_diamond_db",
    "ref::esm_c_600m_weights",
    "ref::ezpred_model",
    "ref::ko_to_mnxr",
    "ref::ec_to_mnxr",
    "ref::uniprot_to_mnxr",
    "ref::reference_label_pool",
]


def _plan_gpr(work: Path):
    """Resolve a plan that produces both the 4-lane and 7-lane GPR tables.

    Planning only — an empty ORF FASTA + empty ref stand-ins are enough.
    """
    inputs = DataInstanceLibrary(work / "inputs.xgdb")
    for tl in ("sequences.yml", "annotation.yml", "ref.yml", "lib.yml"):
        inputs.AddTypeLibrary(MLIB / "data_types" / tl)

    # the one real input: an ORF FASTA (empty; type-driven planning)
    orfs = work / "orfs.faa"
    orfs.touch()
    inputs.AddItem(orfs, "sequences::orfs")

    # type-only stand-ins for every staged reference the stage consumes
    for i, dtype in enumerate(STAGED_REFS):
        stub = work / f"ref_{i}.dat"
        stub.touch()
        inputs.AddItem(stub, dtype)
    inputs.Save()

    resources = [
        DataInstanceLibrary.Load(MLIB / "resources" / "env"),
        DataInstanceLibrary.Load(MLIB / "resources" / "lib"),
        inputs,
    ]
    transforms = [
        TransformInstanceLibrary.Load(MLIB / "transforms" / "functionalAnnotation"),
    ]

    # both GPR tables must be demanded so both mappers survive planning
    targets = TargetBuilder()
    targets.Add("annotation::gpr_table_4lane")
    targets.Add("annotation::gpr_table_7lane")

    agent = Agent(home=Source.FromLocal(work / "agent_home"), runtime=ContainerRuntime.APPTAINER)
    task = agent.GenerateWorkflow(
        samples=list(inputs.AsSamples("sequences::orfs")),
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


def test_gpr_pipeline_compiles(tmp_path):
    """The annotation -> GPR-table stage plans and wires in all nine transforms."""
    task = _plan_gpr(tmp_path)

    assert task.ok, f"gpr workflow failed to plan: {task.plan}"

    names = _step_transform_names(task)
    missing = EXPECTED_TRANSFORMS - names
    assert not missing, (
        f"gpr plan is missing expected transform(s): {sorted(missing)}; "
        f"plan used: {sorted(names)}"
    )

    svg = _render_dag(task, ARTIFACTS / "gpr_dag.svg")
    assert svg.exists() and svg.stat().st_size > 0, f"DAG SVG not written: {svg}"


if __name__ == "__main__":
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        task = _plan_gpr(Path(td))
        if not task.ok:
            raise SystemExit(f"FAILED to plan: {task.plan}")
        names = _step_transform_names(task)
        print(f"Plan OK -- {len(task.plan.steps)} steps")
        for step in sorted(task.plan.steps, key=lambda s: s.order):
            prods = [i.dtype_name for g in step.produces for i in g]
            print(f"  step {step.order}: {Path(step.transform._path).stem} -> {prods}")
        missing = EXPECTED_TRANSFORMS - names
        print(f"missing transforms: {sorted(missing) if missing else 'none'}")
        svg = _render_dag(task, ARTIFACTS / "gpr_dag.svg")
        print(f"DAG written to: {svg}")
