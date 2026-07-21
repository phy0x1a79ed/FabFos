#!/usr/bin/env python
"""Iterate a null-scoring family SOLO -- no solver, no nextflow.

    python dev_null_family.py --family gpd
    python dev_null_family.py --family empirical

This is the answer to "can transforms be run solo without invoking nextflow?":
yes, via `RunTransform`. It skips the workflow solver and nextflow entirely, but
routes through the SAME `ExecuteStep` the DAG uses -- it builds the same
ExecutionContext and runs the same protocol. So what you iterate on here is what
runs in the pipeline: no separate dev harness to drift from production, which is
the usual way a "tested" transform still breaks once planned.

The loop:
  1. add a family to FAMILIES in <engine>/resources/lib/ecspr_null_study.py
  2. run it here until it does what you meant
  3. plan the SAME file into the DAG unchanged -- see experiments/scadc_null_study.py

Because the family arrives as a staged, content-hashed input, two families cannot
collide on one cache entry. The bug that motivated this -- a draws file whose name
did not record K, silently overwritten in place by a run at a different K, with the
mtime moving backwards -- is not mitigated here; it is unrepresentable.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from fabfos import canon  # noqa: E402
from metasmith.python_api import RunTransform  # noqa: E402

WORK = HERE / ".dev" / "null_family"


def build_curated_nulls(dest: Path) -> Path:
    """Explicit list from canon -- never a glob. See MIGRATION.md."""
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    for name in canon.FROZEN_NULL_FILES:
        (dest / name).symlink_to(canon.INCUMBENT_K1000 / name)
    return dest


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", default="gpd",
                    help="a key of FAMILIES in the engine's ecspr_null_study.py")
    a = ap.parse_args()

    work = WORK / a.family
    work.mkdir(parents=True, exist_ok=True)
    nulls = build_curated_nulls(work / "nulls")

    spec = work / "null_model_spec.yml"
    spec.write_text(f"family: {a.family}\n")

    exp = work / "experiment.txt"
    exp.write_text("dev_null_family")

    lib = canon.ENGINE_LIB
    inputs = [
        ("fosmids::recovery_experiment", exp),
        ("ecspr::reff_axes_report", canon.incumbent_axes_report("reff")),
        ("ecspr::ieff_axes_report", canon.incumbent_axes_report("ieff")),
        ("sequences::open_reading_frames", canon.ORFS_FAA),
        ("ecspr::frozen_null", nulls),
        ("ecspr::null_model_spec", spec),
        ("containers::ecspr.oci", lib / "resources/containers/ecspr.oci"),
        ("lib::ecspr_null_study.py", lib / "resources/lib/ecspr_null_study.py"),
        ("lib::ecspr_significance.py", lib / "resources/lib/ecspr_significance.py"),
    ]

    print(f"=== solo run: family={a.family} (no solver, no nextflow) ===")
    for t, p in inputs:
        print(f"  {t:36s} {p.name}")
    res = RunTransform(
        transform_lib=lib / "transforms" / "ecspr",
        transform="significance_study.py",
        inputs=inputs,
        work_dir=work / "out",
    )
    print(f"\nsuccess: {getattr(res, 'success', None)}")
    for f in sorted((work / "out").rglob("*.tsv")):
        print(f"  -> {f}")
    return 0 if getattr(res, "success", False) else 1


if __name__ == "__main__":
    sys.exit(main())
