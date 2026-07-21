"""SCADC main: score the canonical basis for ECSPr significance.

Stages the canonical solve and asks for the significance targets. Every path and
every number comes from `canon`; this file transcribes neither.

WHY THE SOLVE IS STAGED RATHER THAN RECOMPUTED
----------------------------------------------
The observed solve is the expensive, already-verified artifact, and re-deriving it
here would change nothing except the risk. What this spec exists to exercise is the
SPINE and the SCORER -- and staging the solve makes that a real gate: the output must
reproduce the incumbent table, so a regression anywhere in staging, planning, the
container, or the scorer shows up as a numeric difference rather than as a plausible
new table.

WHY THERE IS NO ANNOTATION LANE HERE
------------------------------------
Deliberate, and the single most important line in this file. The driver this replaces
re-ran annotation from gated databases on every invocation. Fresh annotation changes
the evidence, which changes the numbers, which destroys the parity signal that makes
any of this trustworthy -- so the pipeline could never be used to check itself. A
canonical driver stages evidence as an INPUT. Running annotation fresh is a separate,
later, deliberate question, and it should be asked by an experiment whose name says so.

Only the `ecspr` domain is loaded. A loaded domain can be selected structurally by the
planner, so leaving the annotation domain in is how you get a lane you did not ask for.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fabfos import canon
from spec import DirInput, ExperimentSpec, Value


def assert_axes():
    """Assert the axis set by COUNT, not by filename -- the filename is what drifted.

    Assert the graph-INDEPENDENT axis definition (AXES_JSON, exactly canon.AXES_N of
    the canon.AXIS_SET axes), not the graph-derived testable subset
    (AXES_TESTABLE_JSON), which on the honest reference graph is a strict subset (one
    phospholipid C axis lost its LCC endpoint).
    """
    canon.assert_canonical_axes(canon.AXES_JSON)


def assert_basis():
    """Assert the ORF FASTA carries the expected basis, including the split contigs.

    Cheap, and it fails here rather than after the planner has staged everything.
    The split contigs are the ones a `\\w`-based ORF regex drops silently.
    """
    contigs = set()
    with open(canon.ORFS_FAA) as fh:
        for line in fh:
            if line.startswith(">"):
                contigs.add(line[1:].split()[0].rsplit("_", 1)[0])
    missing = [c for c in canon.SPLIT_CONTIGS if c not in contigs]
    if missing:
        raise SystemExit(f"basis FASTA is missing split contigs: {missing}")
    if len(contigs) != canon.FOSMID_BASIS:
        raise SystemExit(
            f"basis FASTA has {len(contigs)} contigs; expected canon.FOSMID_BASIS")


def note_production_status():
    """Surface the canonical method's PRODUCTION status so a pending state reads as
    'produce these outputs', not 'this experiment is misconfigured'. Non-fatal: the
    input-existence check that follows is what stops a run; this only explains it. The
    canonical orientation is the declared standard regardless of whether its frozen
    outputs exist yet (see canon.CANONICAL_ORIENTATION)."""
    st = canon.production_status()
    if st["pending"]:
        print(
            f"  PENDING PRODUCTION: canonical orientation is '{st['orientation']}'; its "
            f"frozen outputs are not on disk yet.\n"
            f"    This experiment STAGES them (it does not produce them) -- run the "
            f"directed reference solve + null\n"
            f"    with the fixed solver to produce them, then re-run. The standard is "
            f"declared; this is validation work, not a misconfiguration.")
    else:
        print(f"  ok: canonical '{st['orientation']}' outputs produced")


SPEC = ExperimentSpec(
    name="scadc_main",
    inputs={
        # the staged, canonical solve -- the frozen reference solve on the honest
        # reference graph (delta_obs). ORIENTATION-DISPATCHED: reference_axes_report()
        # returns the undirected OR the directed (diode) frozen table per
        # canon.CANONICAL_ORIENTATION. This spec never re-solves -- it stages the frozen
        # report -- so it needs NO roles/direction inputs (those are the solve_directed
        # transform's, not the scorer's); flipping the orientation moves this spec onto
        # the directed reference table with no edit here, only the one-line canon flip.
        "ecspr::reff_axes_report": canon.reference_axes_report("reff"),
        "ecspr::ieff_axes_report": canon.reference_axes_report("ieff"),
        # ORF counts for null size matching
        "sequences::open_reading_frames": canon.ORFS_FAA,
        # the null -- the reference-GRAPH null, the matched pair to the solve above:
        # delta_obs and its null are on the SAME graph AND the SAME orientation.
        # reference_null_files() dispatches on canon.CANONICAL_ORIENTATION in lockstep
        # with the solve above, so a flip can never leave delta_obs and its null on
        # different orientations. Curated from an EXPLICIT list (the scorer derives its
        # draw sizes from whatever is here, so this list IS the basis -- a glob would
        # silently redefine it).
        "ecspr::frozen_null": DirInput(canon.reference_null_files()),
        # device/dtype as a hashed input, so a GPU run and a CPU run of the same
        # step cannot collide on one cache entry
        "ecspr::compute_profile": Value("compute_profile.yml", canon.COMPUTE_CPU),
    },
    per_experiment=frozenset({
        "ecspr::reff_axes_report",
        "ecspr::ieff_axes_report",
        "sequences::open_reading_frames",
    }),
    targets=(
        "ecspr::reff_significance",
        "ecspr::ieff_significance",
    ),
    domains=("ecspr",),
    namespaces=("sequences", "fosmids", "ecspr"),
    preflight=(assert_axes, assert_basis, note_production_status),
)
