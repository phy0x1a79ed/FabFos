"""Resolve and render the DAG that builds every compiled reference.

    PATH="/home/tony/lib/miniforge3/envs/msm/bin:$PATH" python examples/build_references_dag.py

THE DELIVERABLE IS THE DAG. None of the build transforms are implemented -- every
`protocol` raises NotImplementedError. What is asserted here is that their
CONTRACTS compose: that each compiled reference has exactly one producer, that
producer's requirements are themselves produced, and that the whole thing closes
from ONE given.

That one given is the MetaCyc drop-in, and it is given only because it is licensed
and cannot be fetched. Everything else -- including which three hosts the
references are built for -- is produced by a transform in the library. The host set
in particular is declared by acquire/host_accessions.py rather than passed in here,
so "the reference build" means one thing rather than whatever accessions a caller
happened to supply.

Planning is type-driven -- nothing is staged, containerised or executed -- so empty
stand-ins are sufficient and this renders on a machine holding none of the 41 GB.

What the graph should show, reading from the leaves:

    metacyc_flatfiles (GIVEN) -> metacyc_licensed -> the two curated members
    host_accessions -> accession -> host_genome -> orfs
                                      -> {kofamscan, clean, diamond_uniref50,
                                          proteinbert} -> gpr_4lane
                                      -> host_gpr_denovo   (B2, x3 hosts)
                                -> host_gem -> host_gpr_gem  (B1, x3 hosts)
    metanetx + metacyc  -> aam_ensemble       -\\
    metanetx + metacyc + equilibrator                       -> bake_metabolism (R6)
                        -> direction_ensemble -/
    metanetx + kegg + rhea -> mnxr_lookup                  (R5)
                           -> reference_label_pool     (R7, with uniref50)
    uniref50 -> uniref50_dmnd                              (R4)
    host_genome -> isolate_assembly                        (R1, the background itself)
    kofam    -> kofam_ref                                  (R3)
    literature -> condition_gpr -> conditions              (B3, B4)

Three things this is really checking, because all three are silent failures rather
than loud ones:

0. Exactly ONE node in the graph is a given. If a second appears, something that
   should be fetched or derived is instead being handed in, and the build stops
   being reproducible from the library alone.

1. `mnxr_lookup` is a SINGLE producer feeding the mapper. The three former bridge
   types are gone; if the mapper still asked for them the plan would simply not
   close, which is the honest outcome. (Measured before consolidating: the three id
   spaces share no ids, and the dedup to distinct (id, mnxr) is behaviour-preserving.)

2. The metabolism trio comes out of ONE transform. Three separate producers would
   plan just as happily and then hand a reader an atom_pairs baked against a
   different vocab, which decodes every node to the wrong metabolite without raising.
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

from metasmith.python_api import (
    Agent,
    ContainerRuntime,
    Source,
    DataInstanceLibrary,
    TransformInstanceLibrary,
    TargetBuilder,
)

REPO = Path(__file__).resolve().parent.parent
MLIB = REPO / "src" / "metasmith_libraries"
BREF = REPO / "build_references"
# Beside the other rendered DAGs (assembly, fosmids, gpr). Gitignored -- a DAG is
# regenerated from the transforms in seconds, so committing one only creates a
# second source of truth that can disagree with the library.
ARTIFACTS = REPO / "tests" / "artifacts"

# Every compiled reference, by artifact id in build_references/REFERENCES.md.
TARGETS = [
    ("R1", "sequences::isolate_assembly"),
    ("R3", "ref::kofamscan_profiles"),
    ("R3", "ref::kofamscan_ko_list"),
    ("R4", "ref::uniref50_diamond_db"),
    ("R5", "ref::mnxr_lookup"),
    ("R6", "ref::atom_pairs"),
    ("R6", "ref::metabolism_vocab"),
    ("R6", "ref::direction_ratios"),
    ("R7", "ref::reference_label_pool"),
    ("B1", "ref::gpr_table_gem"),
    ("B2", "ref::gpr_table_denovo"),
    ("B3", "bench::condition_gpr"),
    ("B4", "bench::conditions"),
]

# Every transform that must appear, so a plan that quietly drops a branch fails
# rather than rendering a smaller graph.
EXPECTED = {
    # acquire
    "host_accessions", "host_genome", "host_gem", "metanetx", "rhea", "kofam",
    "uniref50",
    "kegg_ko_reactions", "metacyc_licensed", "equilibrator_cache",
    "literature",
    # compile
    "kofam_ref", "uniref50_dmnd", "mnxr_lookup",
    "aam_ensemble", "direction_ensemble", "bake_metabolism",
    "reference_label_pool",
    # the shipped annotation lanes + mapper, on the host proteomes
    "kofamscan", "clean", "diamond_uniref50", "proteinbert", "gpr_4lane",
    # benchmark
    "host_gpr_gem", "host_gpr_denovo", "condition_gpr", "conditions",
}


def plan(work: Path):
    inputs = DataInstanceLibrary(work / "inputs.xgdb")
    for tl in ("ncbi.yml", "sequences.yml", "annotation.yml", "ref.yml",
               "lib.yml"):
        inputs.AddTypeLibrary(MLIB / "data_types" / tl)
    for tl in ("raw.yml", "interm.yml", "bench.yml"):
        inputs.AddTypeLibrary(BREF / "data_types" / tl)

    # THE ONE GIVEN: the licensed MetaCyc drop-in. Type-only stand-in, since
    # planning never opens it. Everything else in the graph is produced -- which is
    # the claim being tested.
    metacyc = work / "metacyc"
    metacyc.mkdir()
    inputs.AddItem(metacyc, "raw::metacyc_flatfiles")
    inputs.Save()

    resources = [
        DataInstanceLibrary.Load(MLIB / "resources" / "env"),
        DataInstanceLibrary.Load(MLIB / "resources" / "lib"),
        inputs,
    ]
    # functionalAnnotation only, from the shipped library. Its logistics/ sibling
    # carries downloaders that produce the same ref:: types compile/ does, and two
    # producers for one reference is a tiebreak deciding provenance.
    transforms = [
        TransformInstanceLibrary.Load(MLIB / "transforms" / "functionalAnnotation"),
        TransformInstanceLibrary.Load(BREF / "transforms" / "acquire"),
        TransformInstanceLibrary.Load(BREF / "transforms" / "compile"),
        TransformInstanceLibrary.Load(BREF / "transforms" / "benchmark"),
    ]

    targets = TargetBuilder()
    for _id, dtype in TARGETS:
        targets.Add(dtype)

    agent = Agent(home=Source.FromLocal(work / "agent_home"),
                  runtime=ContainerRuntime.APPTAINER)
    return agent.GenerateWorkflow(
        samples=list(inputs.AsSamples("raw::metacyc_flatfiles")),
        resources=resources,
        transforms=transforms,
        targets=targets,
    )


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        task = plan(Path(td))
        if not task.ok:
            print(f"FAILED to plan:\n{task.plan}")
            return 1

        used = {Path(s.transform._path).stem for s in task.plan.steps}
        print(f"Plan OK -- {len(task.plan.steps)} steps, {len(used)} distinct transforms\n")
        for step in sorted(task.plan.steps, key=lambda s: s.order):
            prods = [i.dtype_name for g in step.produces for i in g]
            print(f"  {step.order:>3}  {Path(step.transform._path).stem:<28} -> {prods}")

        missing = EXPECTED - used
        extra = used - EXPECTED
        print()
        if missing:
            print(f"MISSING expected transforms: {sorted(missing)}")
        if extra:
            print(f"note -- transforms used but not in EXPECTED: {sorted(extra)}")
        if not missing:
            print("every expected transform is in the plan")

        ARTIFACTS.mkdir(parents=True, exist_ok=True)
        svg = ARTIFACTS / "build_references_dag.svg"
        task.plan.RenderDAG(svg, blacklist_namespaces={"lib", "env"})
        print(f"\nDAG -> {svg}")
        return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
