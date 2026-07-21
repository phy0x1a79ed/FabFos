#!/usr/bin/env python
"""The parity gate: prove the engine reproduces the incumbent, or fail.

This is the task that earns everything else. Until it is green, `canon.STATUS` stays
PROVISIONAL and the incumbent tables remain authoritative -- the FabFos path is a
claim, and there are already four of those.

    python run_parity.py                # canonical lane
    python run_parity.py --lane reff    # the other one
    python run_parity.py --all

What it does: builds a CURATED nulls directory from `canon.FROZEN_NULL_FILES` -- an
explicit list, never a glob of a cache -- and hands it to the engine's `selftest`,
which joins against `canon.incumbent_sig_table()` and exits non-zero past
`canon.PARITY_TOL`.

Two things it checks that a numeric diff alone would miss:

  * the split contigs (`canon.SPLIT_CONTIGS`) appear in the output. A `\\w`-based
    ORF-id regex drops them SILENTLY -- no error, the fosmids simply vanish -- so
    their presence has to be asserted positively rather than noticed.
  * `at_floor` agrees EXACTLY. It marks which rows rest on the fitted tail rather
    than on counted draws; a row that quietly changes side is a row whose provenance
    claim changed.

If parity cannot be reached: STOP. Do not tune the port to match and do not relax
`canon.PARITY_TOL`. A real divergence is a finding about the incumbent, and it needs
a human.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fabfos import canon  # noqa: E402

HERE = Path(__file__).resolve().parent
STAGED_NULLS = HERE / "_staged_nulls"
ENGINE_SIG = canon.ENGINE_LIB / "resources" / "lib" / "ecspr_significance.py"


def build_curated_nulls() -> Path:
    """Symlink exactly the canonical nulls into a clean directory.

    Explicit list only. The live cache holds retired sizes and an un-suffixed draws
    file left by the overwrite incident; even the byte-copy backup holds retired
    sizes. Globbing either would silently widen the basis.
    """
    if STAGED_NULLS.exists():
        for f in STAGED_NULLS.iterdir():
            f.unlink()
    STAGED_NULLS.mkdir(parents=True, exist_ok=True)
    for name in canon.FROZEN_NULL_FILES:
        src = canon.INCUMBENT_K1000 / name
        if not src.exists():
            raise SystemExit(f"curated null missing from the frozen backup: {src}")
        (STAGED_NULLS / name).symlink_to(src)
    return STAGED_NULLS


def assert_fitter_deps(python: str) -> None:
    """Refuse to run in an interpreter that cannot fit the mixture.

    Without sklearn the fitter does not raise -- it falls back to a single
    log-normal, flags the row `no-sklearn`, and returns finite numbers. Those
    numbers are then diffed against an incumbent that WAS fitted with sklearn, so
    the gate fails at ~3e-1 and reports FAIL directly beneath a banner saying a real
    divergence is a finding about the incumbent. It is not: it is a missing package.

    This gate exists to make a wrong number impossible to mistake for a right one,
    so it must not itself hand you a wrong diagnosis. Fail here, naming the cause.
    `sys.executable` is the default `--python`, and the env that runs this spine
    (metasmith) is not the env that can fit the mixture -- so the DOCUMENTED
    invocation is exactly the one that trips this.
    """
    probe = "import sklearn, numpy, pandas, scipy"
    r = subprocess.run([python, "-c", probe], capture_output=True, text=True)
    if r.returncode != 0:
        missing = r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "?"
        raise SystemExit(
            f"parity cannot run under {python}\n"
            f"  {missing}\n\n"
            f"The mixture fitter DEGRADES rather than failing when sklearn is absent "
            f"(single log-normal, flag `no-sklearn`), so this would report a FAIL that "
            f"looks like a port divergence and is not one.\n"
            f"Pass an interpreter that has scikit-learn:\n"
            f"  python parity/run_parity.py --all --python <env>/bin/python"
        )


def run_lane(lane: str, python: str) -> bool:
    ref = canon.incumbent_sig_table(lane)
    cmd = [
        python, "-u", str(ENGINE_SIG), "selftest",
        "--lane", lane,
        "--nulls-dir", str(STAGED_NULLS),
        "--report-dir", str(canon.INCUMBENT_CACHE),
        "--reference", str(ref),
        "--faa", str(canon.ORFS_FAA),
        "--tol", repr(canon.PARITY_TOL),
        "--kmax", str(canon.KMAX),
        "--seed", str(canon.FIT_SEED),
        "--split-contigs", *canon.SPLIT_CONTIGS,
    ]
    print(f"\n=== parity: {lane} vs {ref.name} ===", flush=True)
    return subprocess.run(cmd).returncode == 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lane", choices=list(canon.LANES), default=canon.CANONICAL_LANE)
    ap.add_argument("--all", action="store_true", help="both lanes")
    ap.add_argument("--python", default=sys.executable)
    a = ap.parse_args()

    assert_fitter_deps(a.python)
    build_curated_nulls()
    print(f"curated {len(canon.FROZEN_NULL_FILES)} null files -> {STAGED_NULLS}")

    lanes = list(canon.LANES) if a.all else [a.lane]
    results = {ln: run_lane(ln, a.python) for ln in lanes}

    print("\n=== parity summary ===")
    for ln, ok in results.items():
        print(f"  {ln}: {'PASS' if ok else 'FAIL'}")
    if not all(results.values()):
        print("\nSTOP. Do not tune the port to match and do not relax canon.PARITY_TOL.")
        print("A real divergence is a finding about the incumbent.")
        return 1
    print(f"\nGate green. canon.STATUS may now be flipped (currently {canon.STATUS}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
