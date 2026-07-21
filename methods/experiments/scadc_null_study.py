"""SCADC null-family study: the same transform, planned into the DAG unchanged.

    FABFOS_NULL_FAMILY=gpd       python run_experiment.py scadc_null_study ...
    FABFOS_NULL_FAMILY=empirical python run_experiment.py scadc_null_study ...

The point of this spec is that `significance_study.py` -- the file iterated solo
via `dev_null_family.py` with no nextflow -- is planned here UNCHANGED. There is no
separate production variant to drift from the one that was tested.

The family is a staged, content-hashed input, so two families produce different
cache keys and cannot overwrite one another. That is the structural fix for the bug
that motivated all of this: a draws file whose name did not record K, silently
overwritten in place by a run at a different K, mid-run, with the file's mtime moving
BACKWARDS so nothing looked wrong. It was recovered only because the seed happened to
be deterministic. Here the overwrite is not mitigated, it is unrepresentable --
outputs are addressed by the hash of their inputs.

This produces STUDY tables. It cannot produce the canonical significance type; that
type has exactly one producer. See the engine's ecspr_null_study.py for why.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fabfos import canon
from spec import DirInput, ExperimentSpec, Value

FAMILY = os.environ.get("FABFOS_NULL_FAMILY", "gpd")


def assert_axes():
    # AXES_JSON = the graph-independent canon.AXIS_SET definition; AXES_TESTABLE_JSON is
    # now a strict subset on the reference graph (see scadc_main.assert_axes).
    canon.assert_canonical_axes(canon.AXES_JSON)


SPEC = ExperimentSpec(
    # Deliberately CONSTANT across families, so every family shares one staging
    # dir and one nextflow cache. That is the honest test: it recreates exactly
    # the situation the overwrite bug lived in -- same run, same paths, one knob
    # different -- and forces the engine to distinguish them on content alone.
    # Putting the family in the name would make them different by construction
    # and prove nothing.
    name="scadc_null_study",
    inputs={
        "ecspr::reff_axes_report": canon.incumbent_axes_report("reff"),
        "ecspr::ieff_axes_report": canon.incumbent_axes_report("ieff"),
        "sequences::open_reading_frames": canon.ORFS_FAA,
        "ecspr::frozen_null": DirInput(
            [canon.INCUMBENT_K1000 / f for f in canon.FROZEN_NULL_FILES]),
        # the family: staged, hashed, and therefore part of the address
        "ecspr::null_model_spec": Value("null_model_spec.yml", f"family: {FAMILY}\n"),
    },
    per_experiment=frozenset({
        "ecspr::reff_axes_report",
        "ecspr::ieff_axes_report",
        "sequences::open_reading_frames",
    }),
    targets=(
        "ecspr::reff_significance_study",
        "ecspr::ieff_significance_study",
    ),
    domains=("ecspr",),
    namespaces=("sequences", "fosmids", "ecspr"),
    preflight=(assert_axes,),
)
