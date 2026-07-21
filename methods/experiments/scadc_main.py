"""SCADC main: score the canonical basis for ECSPr ground-probe significance.

Stages the reference ingredients and asks for `ecspr::ground_significance`. Every path
and every number comes from `canon`; this file transcribes neither.

WHAT CHANGED HERE (2026-07-20)
------------------------------
This spec used to STAGE a frozen two-terminal solve (`ecspr::{reff,ieff}_axes_report`)
and ask only for the scorer, so that the output had to reproduce the incumbent table and
any regression in staging / planning / the container / the scorer showed up as a numeric
difference. That gate is gone with the lane it gated.

The star topology joined every metabolite to reaction-node HUBS, and eliminating a
reaction node -- which is exactly what the Woodbury update does -- left a CLIQUE over its
participants, giving two participants that share no atom a conductance between them
(measured on MNXR106432, carbon: a zero-carbon channel 21x a real one). Reproducing that
table would be reproducing the artifact. So this spec now stages the INGREDIENTS -- the
atom-pair table, the direction ensemble, the evidence weights, the addition maps -- and
the atom lane builds, probes and scores from them.

The replacement gate is not parity against a frozen table; it is CONSERVATION plus the
SYMMETRIC LIMIT, both asserted in the engine's own self-tests (`ecspr_graph.py`), plus
the pulse-chase suite. Those gate the measurement rather than its agreement with an
earlier run of itself.

WHY THERE IS NO ANNOTATION LANE HERE
------------------------------------
Deliberate, and still the single most important line in this file. The driver this
replaces re-ran annotation from gated databases on every invocation. Fresh annotation
changes the evidence, which changes the numbers, which destroys any hope of attributing a
difference to the thing you changed -- so the pipeline could never be used to check
itself. A canonical driver stages evidence as an INPUT. Running annotation fresh is a
separate, later, deliberate question, and it should be asked by an experiment whose name
says so.

WHICH DOMAINS ARE LOADED, AND WHY EXACTLY THESE
-----------------------------------------------
A loaded domain can be selected structurally by the planner, so leaving one in is how you
get a lane you did not ask for. THREE domains produce `ecspr::atom_graph` -- `ecsprAtomB`
(evidence-weighted), `ecsprAtomA` (curated GEM), `ecsprAtomGPR` (GEM + active gene set).
Exactly one is loaded here, and which one is a claim about the method.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fabfos import canon
from spec import DirInput, ExperimentSpec


def assert_axes():
    """Assert the axis set by COUNT, not by filename -- the filename is what drifted.

    The ground probe derives its media / central / precursor sets from THIS file's
    contents, so an unnoticed substitution would silently redefine what "biomass" means
    for every share in the output.
    """
    canon.assert_canonical_axes(canon.AXES_JSON)


def assert_basis():
    """Assert the ORF FASTA carries the expected basis, including the split contigs.

    Cheap, and it fails here rather than after the planner has staged everything. The
    split contigs are the ones a `\\w`-based ORF regex drops silently.
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


def note_null_status():
    """Surface whether the ground null has been produced yet.

    Non-fatal: the input-existence check that follows is what stops a run; this only
    explains it. `ecspr::ground_null` is BUILD-side and expensive, and unlike its
    predecessor it HAS a producing transform (`ecsprGround/null.py`) -- so a pending
    state reads as "produce this artifact", not "this experiment is misconfigured".
    """
    missing = [p.name for p in canon.ground_null_paths() if not p.exists()]
    if missing:
        print(f"  PENDING PRODUCTION: {len(missing)} ground-null draw file(s) are not on "
              f"disk yet: {missing}\n"
              f"    Produce them with the ecsprGround/null.py transform, then re-run. The "
              f"draw-size grid comes from canon.DRAW_SIZES (the run's own ORF\n"
              f"    percentiles); a size that is missing must be GENERATED, never "
              f"interpolated across.")
    else:
        print(f"  ok: all {len(canon.DRAW_SIZES)} ground-null draw sizes present")


SPEC = ExperimentSpec(
    name="scadc_main",
    inputs={
        # The atom-transfer universe and its directionality. Both are network-AGNOSTIC
        # static functions of the MNXR ids, so they are shared staged inputs.
        "ecspr::atom_pairs": canon.REFERENCE_ATOM_PAIRS,
        "ecspr::direction_ratios": canon.REFERENCE_DIRECTION,
        # The host's per-reaction conductances and the per-fosmid addition maps -- the
        # matched pair. Staged rather than recomputed so that two runs differ only in
        # what they were asked to differ in.
        "ecspr::evidence_weights": canon.EVIDENCE_WEIGHTS,
        "ecspr::addition_weights": canon.ADDITION_WEIGHTS,
        # The biomass DAG. The probe derives media / central / precursor from it.
        "ecspr::biomass_axes": canon.AXES_JSON,
        # ORF counts for null size matching.
        "sequences::open_reading_frames": canon.ORFS_FAA,
        # The null. Curated from an EXPLICIT list: the scorer derives its draw sizes by
        # listing this directory, so this list IS the basis -- a glob would silently
        # redefine it, and a cache holds retired sizes beside canonical ones.
        "ecspr::ground_null": DirInput(canon.ground_null_paths()),
    },
    per_experiment=frozenset({
        "ecspr::evidence_weights",
        "ecspr::addition_weights",
        "sequences::open_reading_frames",
    }),
    targets=(
        "ecspr::ground_significance",
    ),
    # ecsprAtomB is the evidence-weighted graph builder. It is ONE of three producers of
    # ecspr::atom_graph and the others must not be loaded beside it.
    domains=("ecspr", "ecsprAtomB", "ecsprGround"),
    namespaces=("sequences", "fosmids", "ecspr"),
    preflight=(assert_axes, assert_basis, note_null_status),
)
