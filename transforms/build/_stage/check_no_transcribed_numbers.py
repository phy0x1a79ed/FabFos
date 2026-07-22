#!/usr/bin/env python
"""Tripwire: fail if a canonical basis value has been transcribed into prose.

`canon.py` is the only place the basis exists. That rule is worth exactly as much
as the thing that checks it -- which is this file. Run it in CI, in the parity
gate, and before any commit to the method or its consumers.

    python check_no_transcribed_numbers.py
    python check_no_transcribed_numbers.py --strict-all   # advisory roots fail too

What is forbidden, and where:

  * CANONICAL values -- the live basis -- may not appear as literals anywhere in a
    scanned root except canon.py. If you need one, import it. This is the rule that
    stops a document from drifting away from the method.

  * RETIRED values may additionally appear in MIGRATION.md. A tombstone has to
    name its dead, and a dead value cannot drift by definition -- that is what
    makes it safe there and nowhere else.

This is a tripwire, not a proof. It catches transcription, which is the failure
that actually happened. It cannot catch a value that is merely described
("roughly two hundred inserts"), and it is not trying to.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]

# WHAT IS SCANNED, and why it is a LIST now. This check used to live at the root of the
# figure scope's `main/fabfos/` and root itself at `Path(__file__).parent` -- so it saw
# exactly the directory it sat in, and NOTHING else. That is how three canon.DRAW_SIZES
# literals drifted into `main/ecspr/`: the tripwire could not see the tree that was
# drifting. Rooting it at the code it governs, explicitly, is the fix.
#
# The figure scopes are included by path because they are the CONSUMERS the rule protects
# and they live in sibling checkouts; a missing one is skipped rather than fatal, so this
# runs on a machine that has only the fabfos project.
_FIGURE_SCOPES = Path("/home/tony/agentic_workspace/projects/scadc")
# (root, strict). STRICT roots fail the gate. ADVISORY roots are reported and counted
# but do not fail it -- a deliberate, and temporary, concession.
#
# WHY ADVISORY EXISTS AT ALL. Widening the scan from `main/fabfos/` to the consumers it
# governs surfaced a backlog of ~190 pre-existing hits in `main/ecspr/`, accumulated over
# the whole period the tripwire could not see that tree. Most are the token `set4`, and a
# large share of those are not transcription at all -- `set4` is also a ROLE LABEL in
# these dataframes ("registry" vs "set4"), i.e. a value-collision, which is what the
# `# canon-ok:` marker is for. Triaging them is real work and it is not this change.
#
# The alternative was to leave the gate red on arrival, and a gate that is red on arrival
# is a gate everyone learns to ignore -- which is how the original blind spot persisted.
# So: the rule now SEES that tree and reports it every run, and the number is the backlog.
# Promote it to strict once triaged; `--strict-all` does it today.
SCAN_ROOTS = [
    (REPO / "src" / "fabfos", True),
    (REPO / "methods", True),
    (REPO / "transforms", True),
    (REPO / "examples", True),
    (_FIGURE_SCOPES / "figures-main" / "main" / "ecspr", False),
]

SELF = {"canon.py", Path(__file__).name}
TOMBSTONE = "MIGRATION.md"

# (label, pattern). Patterns are matched against the raw file text.
CANONICAL = [
    ("fosmid basis (canon.FOSMID_BASIS)", r"\b199\b"),
    ("axis count (canon.AXES_N)", r"\b79\b"),
    ("draw count K (canon.DRAW_K)", r"\b1000\b"),
    ("draw seed (canon.DRAW_SEED)", r"\bseed[\s=:]+42\b"),
    ("axis set (canon.AXIS_SET)", r"\bset4\b"),
    ("scorer name (canon.SCORER)", r"\bsig_mix\b"),
    ("CLEAN floor (canon.CLEAN_FLOOR)", r"(?<![\d.])0\.01(?![\d])"),
    # CORRECTED 2026-07-20 (was the RETIRED grid 14,28,34,42,51). Listing the retired
    # grid here INVERTED the gate: it flagged files transcribing the OLD values as
    # canonical transcriptions, and never detected a transcription of the LIVE grid -- so
    # a file hardcoding the real draw sizes passed silently, which is precisely what this
    # checker exists to catch. The old pattern moved to RETIRED, where it belongs.
    ("draw sizes (canon.DRAW_SIZES)", r"14\s*,\s*25\s*,\s*30\s*,\s*35\s*,\s*43"),
    ("axis split (canon.AXES_PER_ELEMENT)", r"\b40\s*/\s*23\s*/\s*[79]\s*/\s*[79]\b"),
]

RETIRED = [
    ("retired figure basis (canon.FIGURE_BASIS_OPEN)", r"\b132\b"),
    ("retired axis set (canon.RETIRED_AXIS_SETS)", r"\bset2cat\b"),
    ("retired scorer (canon.RETIRED_SCORERS)", r"\bsig_emp\b"),
    ("retired scorer (canon.RETIRED_SCORERS)", r"\bsig_negbin\b"),
    ("retired draw sizes (canon.RETIRED_DRAW_SIZES)", r"\{?\s*21\s*,\s*34\s*,\s*56\s*\}?"),
    # Moved down from CANONICAL 2026-07-20. The pre-2026-07-18 grid; fir's committed
    # manifest.txt and the pre-T1 e2e mirror still enumerate it, so it is worth detecting
    # as a RETIRED transcription.
    ("retired draw sizes (canon.RETIRED_DRAW_SIZES)", r"14\s*,\s*28\s*,\s*34\s*,\s*42\s*,\s*51"),
    ("retired column (canon.SIG_STALE_COLUMNS)", r"\bmatched_N\b"),
]

SCAN_SUFFIXES = {".py", ".md", ".sh", ".yml", ".yaml", ".json", ".toml", ".cfg"}

# A line carrying this marker (with a mandatory reason after the colon) is exempt: it is a
# value-collision, not a transcription. Audit every use -- the reason is the justification.
EXEMPT = re.compile(r"#\s*canon-ok:\s*\S")

# The rule governs SOURCE WE AUTHOR here, not artifacts that land in the tree at
# runtime. A run stages the engine's own library files under .runs/ -- those
# legitimately carry these values (that is what the engine is), and they are
# gitignored. Scanning them reports the engine's constants as our transcription,
# which is both wrong and noisy enough to get the whole check ignored.
#
# `v3_build` is the same class and arrived with the benchmark-v3 merge: it is the
# FROZEN benchmark build tree -- MANIFEST.tsv, the answer key under Y/, the scored
# tables -- hashed and pinned in provenance. Its provenance JSON records which axis
# set the panel was built from, which is a RECORD of a value, not a restatement of
# one; and it is frozen, so the only way to satisfy the gate by editing it would be
# to break the hash the freeze exists to protect.
#
# `_library` is the .runs/ case again, one directory up: `dev.sh -b` bundles
# src/metasmith_libraries/ into src/fabfos/_library/ for shipping, and that path is
# gitignored generated output, not source we author. Scanning it charged us 66 of 75
# strict hits -- every one of them the ENGINE restating its own constants, in files
# whose originals live in a submodule this gate deliberately does not govern. It also
# double-counted: the same line is scanned once in the bundle and never in its source.
SKIP_DIRS = {"__pycache__", ".runs", "_staged_nulls", ".git", "v3_build", "_library"}


def _skipped(rel: Path) -> bool:
    # any dot-directory, plus the named ephemerals
    return any(p in SKIP_DIRS or (p.startswith(".") and p not in (".",))
               for p in rel.parts[:-1])


def scan(strict_all: bool = False) -> tuple[list[str], list[str]]:
    """Returns (blocking hits, advisory hits)."""
    blocking: list[str] = []
    advisory: list[str] = []
    for root, strict in SCAN_ROOTS:
        if not root.exists():
            # a sibling figure scope that is not checked out on this machine
            print(f"  (skipping absent root {root})")
            continue
        found = _scan_root(root)
        (blocking if (strict or strict_all) else advisory).extend(found)
    return blocking, advisory


def _scan_root(root: Path) -> list[str]:
    hits: list[str] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix not in SCAN_SUFFIXES:
            continue
        rel = path.relative_to(root.parent)
        if _skipped(rel):
            continue
        if path.name in SELF:
            continue
        rules = list(CANONICAL)
        if path.name != TOMBSTONE:
            rules += RETIRED
        try:
            text = path.read_text()
        except UnicodeDecodeError:
            continue
        for lineno, line in enumerate(text.splitlines(), 1):
            # A line may opt out ONLY with an explicit, reasoned marker -- the escape hatch for
            # a literal that collides by VALUE with a basis but is a different quantity (e.g. an
            # Indigo timeout in ms that happens to equal DRAW_K). The reason is mandatory so the
            # exemption is auditable, never a silent suppression.
            if EXEMPT.search(line):
                continue
            for label, pat in rules:
                if re.search(pat, line):
                    hits.append(f"{rel}:{lineno}: transcribed {label}\n      {line.strip()}")
    return hits


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strict-all", action="store_true",
                    help="fail on advisory roots too (use once the backlog is triaged)")
    ap.add_argument("--show-advisory", action="store_true",
                    help="list the advisory hits instead of only counting them")
    a = ap.parse_args()

    blocking, advisory = scan(strict_all=a.strict_all)

    if advisory:
        print(f"\nADVISORY: {len(advisory)} transcribed value(s) in non-strict root(s) "
              f"-- reported, not failing. Re-run with --show-advisory to list them, "
              f"--strict-all to enforce.")
        if a.show_advisory:
            for h in advisory:
                print(f"  {h}")

    if blocking:
        print(f"\nFAIL: {len(blocking)} transcribed basis value(s) in strict root(s)\n")
        for h in blocking:
            print(f"  {h}")
        print("\nImport the name from canon.py instead of restating its value.")
        return 1
    print("\nPASS: no transcribed basis values in strict roots "
          f"({', '.join(r.name for r, s in SCAN_ROOTS if s)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
