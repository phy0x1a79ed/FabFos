#!/usr/bin/env python3
"""Merge the sweep's per-(element, draw-size) shards into the canonical null files.

The sweep emits 20 shards per lane -- `{ieff,reff}_{X}_N{N}.tsv`, one per
(element x draw-size). `mnxref.null_directed` selects something different: FIVE
files per lane, `{ieff,reff}_null_directed_N{N}.tsv`, each holding all four
elements. This closes that gap.

It is a script rather than a shell one-liner for the same reason `stage_null.py`
is: the fir harness's staging was only ever done by hand, which is why it was
never reproducible. A hand-run `cat` here would recreate exactly that gap at the
last step of the chain.

WHY THE ROW COUNT IS ASSERTED, not merely reported. Three separate mechanisms on
this path narrow the null basis SILENTLY rather than failing:

  * `null_shard.py --resume` checks expected-vs-actual ROW COUNT only, and the
    count barely moves under a universe change.
  * `discover_draw_sizes` ENUMERATES what is present; a missing size is simply
    never discovered.
  * the canonical scorer looks its keys up with a bare `.get()` + `continue`, so
    missing rows shrink the significance table instead of raising.

A short merge would therefore sail through every downstream check and produce a
smaller, wrong answer that looks finished. So the count is a hard gate here.

THE EXPECTED COUNT IS ALSO A REAL TEST, not bookkeeping. Expected rows per merged
file = (total testable axes) x (draws per size) = 79 x 4,000 = 316,000. The
PREVIOUS canonical null holds 312,000 = 78 x 4,000. The +4,000 difference is the
carbon axis that the stale evidence weights had been silently zeroing (C went
39/40 -> 40/40 testable when the weights were repointed at the CLEAN 199-fosmid
basis). So a merged file landing at 312,000 does not mean "close enough" -- it
means the coherent weights did not take, and the whole chain is suspect.

Usage:
    python merge_null.py --shard-dir <sc>/out --out-dir <dir> \\
        --testable <base-dir>/axes_testable.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ELEMENTS = ("C", "N", "S", "P")
LANES = ("ieff", "reff")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shard-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--testable", required=True, type=Path,
                    help="axes_testable.json from the SAME base graphs the sweep ran on")
    ap.add_argument("--draw-sizes", nargs="+", type=int, default=None,
                    help="default: canon.DRAW_SIZES (never a literal here)")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(repo / "src"))
    from fabfos import canon  # noqa: E402

    draw_sizes = args.draw_sizes or list(canon.DRAW_SIZES)
    testable = json.loads(args.testable.read_text())
    n_axes = sum(len(testable.get(X, [])) for X in ELEMENTS)
    per_el = {X: len(testable.get(X, [])) for X in ELEMENTS}

    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[merge] {n_axes} testable axes {per_el}; grid {tuple(draw_sizes)}")

    failures: list[str] = []
    for lane in LANES:
        for N in draw_sizes:
            shards = [args.shard_dir / f"{lane}_{X}_N{N}.tsv" for X in ELEMENTS]
            missing = [p.name for p in shards if not p.exists()]
            if missing:
                failures.append(f"{lane} N={N}: missing shards {missing}")
                continue

            dest = args.out_dir / f"{lane}_null_directed_N{N}.tsv"
            header = None
            n_rows = 0
            with open(dest, "w") as out:
                for p, X in zip(shards, ELEMENTS):
                    with open(p) as fh:
                        h = fh.readline()
                        if header is None:
                            header = h
                            out.write(h)
                        elif h != header:
                            failures.append(f"{lane} N={N}: header mismatch in {p.name}")
                            break
                        rows = 0
                        for line in fh:
                            out.write(line)
                            rows += 1
                        exp_el = per_el[X] * 4000
                        if rows != exp_el:
                            failures.append(
                                f"{lane} N={N} element {X}: {rows:,} rows, expected "
                                f"{exp_el:,} ({per_el[X]} axes x 4,000)")
                        n_rows += rows

            expected = n_axes * 4000
            status = "OK" if n_rows == expected else "SHORT"
            print(f"[merge] {dest.name}: {n_rows:,} rows (expected {expected:,}) {status}")
            if n_rows != expected:
                failures.append(f"{dest.name}: {n_rows:,} rows, expected {expected:,}")

    if failures:
        print("\nMERGE FAILED -- the null basis would be silently narrowed:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1

    print(f"\nmerged {len(LANES) * len(draw_sizes)} files to {args.out_dir}")
    print(f"each {n_axes * 4000:,} rows; {n_axes} axes x 4,000 draws")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
