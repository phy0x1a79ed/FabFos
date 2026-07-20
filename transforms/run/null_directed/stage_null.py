#!/usr/bin/env python3
"""Build the self-contained stage directory the directed-null shards run from.

The fir harness (`fig-model/main/fabfos/directed/fir/`) was written to be
canon-free and `--stage-dir` parametrised precisely so it could be staged onto
any cluster -- but the staging step itself was only ever done by hand, so the
harness has never been reproducible. `grep -rl draws_resolved` over the whole
tree hits only `fir/README.md` and `fir/fir_null.py`: the producer of the thing
they both read did not exist. This is that producer.

Two resolutions happen HERE, on the workstation, rather than on the cluster:

  * the frozen draws parquet -> `draws_resolved.pkl` (`{N: {"keys": [...],
    "draw_w": {(style, iter): {mnxr: E}}}}`), via `_null_canon.draw_to_weights`
  * the direction ensemble parquet -> `ratios.pkl` (`{mnxr: g_rev/g_fwd}`)

so the cluster venv needs neither pyarrow nor the incumbent tree on sys.path.
That is the whole reason the fir shards are as thin as they are, and it is worth
preserving: the incumbent tree is not something we want to have to replicate on
a login node.

WHERE THE DRAWS COME FROM, and why it is not what canon says. `canon.py` warns
that the live `04_reaction_network/cache` "is never a safe source" and names
`INCUMBENT_K1000` instead. That guidance is stale and, for the live grid,
inverted: K1000 holds only the RETIRED sizes (14/21/28/34/42/51/56, all stamped
2026-07-14) and has none of 25/30/35/43, while the live cache holds the whole
live grid written in one run on 2026-07-18. So the live cache is both the only
source and the internally-consistent one. The single divergence between them --
`N14`, 1,000 of 4,000 rows -- is confined entirely to null style E (A/B/D are
byte-identical) and is old-E vs new-E, not corruption.

Style E is staged, deliberately. It was retired FIG-MODEL-SIDE ONLY (that lane's
`_data.STYLES = ["A","B","D"]`); on the canonical path `canon.STYLES` still
carries it and `resources/lib/ecspr_significance.py:75` still scores it. Dropping
E would cut the sweep ~25% and would ALSO silently shrink the canonical
significance table, because the scorer looks its keys up with a bare
`cell_p.get(...)` + `continue`. A quiet 25% smaller answer is not a saving.

Usage:
    python stage_null.py --stage-dir /path/to/stage [--draw-sizes 14 25 30 35 43]

Then rsync the stage dir to the cluster and submit `run_null.sbatch`.
"""
from __future__ import annotations

import argparse
import json
import pickle
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
ENGINE_LIB = REPO / "src" / "metasmith_libraries" / "resources" / "lib"

# The frozen ORF samples. Graph- AND solver-independent (they are ORF identity
# samples), so they are reused verbatim across universes -- a tier change does
# not invalidate them. See the module docstring for why this is the live cache.
INCUMBENT_RN = Path(
    "/home/tony/agentic_workspace/projects/scadc/metabolic-modelling"
    "/main/metabolic-modelling/04_reaction_network"
)
FROZEN_CACHE = INCUMBENT_RN / "cache"

# The three engine modules a shard imports. Listed explicitly rather than
# globbed: the lib directory holds the whole ECSPr surface, and staging a module
# we do not mean to run is how a shard ends up executing something that was
# never gated.
ENGINE_MODULES = ("ecspr_directed.py", "ecspr_network.py", "ecspr_solver.py")

ELEMENTS = ("C", "N", "S", "P")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage-dir", required=True, type=Path)
    ap.add_argument("--draw-sizes", nargs="+", type=int, default=None,
                    help="default: canon.DRAW_SIZES (never a literal here)")
    ap.add_argument("--base-dir", type=Path, required=True,
                    help="directory holding base_{X}.pkl + axes_testable.json")
    ap.add_argument("--elements", nargs="+", default=list(ELEMENTS))
    args = ap.parse_args()

    sys.path.insert(0, str(REPO / "src"))
    from fabfos import canon  # noqa: E402

    # The grid is READ from canon, never transcribed. Three stale copies of this
    # tuple already exist in the trees (fir/manifest.txt, the incumbent
    # check_no_transcribed_numbers.py, and the pre-T1 run_fabfos_e2e mirror), and
    # every one of them narrows the null basis SILENTLY rather than raising --
    # `discover_draw_sizes` enumerates what is present, so a missing size is
    # simply never discovered. This file will not become the fourth.
    draw_sizes = args.draw_sizes or list(canon.DRAW_SIZES)
    if args.draw_sizes:
        retired = sorted(set(draw_sizes) & set(canon.RETIRED_DRAW_SIZES))
        if retired:
            print(f"REFUSING: {retired} are in canon.RETIRED_DRAW_SIZES. "
                  f"The live grid is {tuple(canon.DRAW_SIZES)}.", file=sys.stderr)
            return 2

    sd: Path = args.stage_dir
    (sd / "lib").mkdir(parents=True, exist_ok=True)
    (sd / "ref" / "solve").mkdir(parents=True, exist_ok=True)
    (sd / "ref" / "graph").mkdir(parents=True, exist_ok=True)
    (sd / "out").mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(ENGINE_LIB))
    from ecspr_network import load_direction_ratios  # noqa: E402

    # ---- engine ----
    for m in ENGINE_MODULES:
        src = ENGINE_LIB / m
        if not src.exists():
            print(f"missing engine module {src}", file=sys.stderr)
            return 1
        shutil.copy2(src, sd / "lib" / m)
    print(f"[stage] lib/: {len(ENGINE_MODULES)} engine modules")

    # ---- static inputs ----
    shutil.copy2(canon.REAC_PROP, sd / "reac_prop.tsv")
    shutil.copy2(canon.AXES_JSON, sd / "axes.json")
    shutil.copy2(Path(args.base_dir) / "axes_testable.json", sd / "axes_testable.json")
    n_ax = sum(len(v) for v in json.loads((sd / "axes_testable.json").read_text()).values())
    print(f"[stage] reac_prop.tsv, axes.json, axes_testable.json ({n_ax} testable axes)")

    # ---- direction ratios, pre-resolved ----
    ratios = load_direction_ratios(canon.REFERENCE_DIRECTION)
    with open(sd / "ratios.pkl", "wb") as fh:
        pickle.dump(ratios, fh, protocol=4)
    print(f"[stage] ratios.pkl: {len(ratios):,} reactions")

    # ---- base graphs + bipartite ----
    for X in args.elements:
        shutil.copy2(Path(args.base_dir) / f"base_{X}.pkl", sd / "ref" / "solve" / f"base_{X}.pkl")
        shutil.copy2(canon.BIPARTITE_DIR / f"mnx_bipartite_{X}.pkl",
                     sd / "ref" / "graph" / f"mnx_bipartite_{X}.pkl")
    print(f"[stage] ref/solve + ref/graph: {len(args.elements)} elements")

    # ---- draws, pre-resolved ----
    sys.path.insert(0, str(INCUMBENT_RN))
    import pandas as pd  # noqa: E402
    from _null_canon import draw_to_weights  # noqa: E402

    mg_w = pickle.load(open(FROZEN_CACHE / "metag_orf_weights.pkl", "rb"))
    resolved: dict = {}
    for N in draw_sizes:
        dp = FROZEN_CACHE / f"null_canonical_N{N}_draws.parquet"
        if not dp.exists():
            print(f"missing frozen draws {dp}", file=sys.stderr)
            return 1
        draws = pd.read_parquet(dp)
        keys = [(r.null, int(r.iter)) for r in draws.itertuples()]
        draw_w = {(r.null, int(r.iter)): draw_to_weights(r.orf_list, mg_w)
                  for r in draws.itertuples()}
        n_ne = sum(1 for w in draw_w.values() if w)
        styles = sorted({k[0] for k in keys})
        resolved[N] = {"keys": keys, "draw_w": draw_w}
        print(f"[stage] draws N={N}: {len(keys):,} ({n_ne:,} non-empty) styles={styles}")
    with open(sd / "draws_resolved.pkl", "wb") as fh:
        pickle.dump(resolved, fh, protocol=4)

    # ---- the shard manifest, generated from the live grid ----
    # fir's committed manifest.txt enumerates the RETIRED grid, so an unedited
    # sbatch there reproduces retired sizes. Generating it here from
    # canon.DRAW_SIZES is what stops that recurring.
    lines = [f"{X} {N}" for X in args.elements for N in draw_sizes]
    (sd / "manifest.txt").write_text("\n".join(lines) + "\n")
    print(f"[stage] manifest.txt: {len(lines)} shards -> --array=0-{len(lines)-1}")

    print(f"\nstaged to {sd}")
    print(f"grid {tuple(draw_sizes)}; elements {tuple(args.elements)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
