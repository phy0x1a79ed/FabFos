#!/usr/bin/env python
"""Regenerate the canonical null on the HONEST reference graph (T2).

The observed reference solve (``REF_ROOT/solve/{reff,ieff}_axes_report.tsv``) already
runs on the reference graph. Its null must run on the SAME graph, or ``delta_obs`` and
the null it is scored against sit on different graphs -- the exact invariant the whole
promotion turns on.

METHOD -- deliberately identical to how the incumbent null was made, only the graph
changed:

  * The DRAWS are graph-independent ORF lists (``null_canonical_N{n}_draws.parquet``);
    reused VERBATIM from the frozen cache -- not regenerated, so the sampling is
    byte-for-byte the incumbent's.
  * ``metag_orf_weights.pkl`` (ORF -> {mnxr: E}) is graph-independent; reused verbatim.
  * The reff null SMW pass (frozen ``17r_reff_null.py``) is re-run with two redirects and
    NOTHING else: (1) ``load_bipartite`` -> the reference universe graph, (2) the base
    read -> the reference-fed base. This is the same float32-GPU ``reff_batch`` the
    incumbent null used, so the reference null is the incumbent null with the graph
    swapped -- attributable.
  * The ieff null is the SAME column transform of the reff null the observed solve uses
    (engine ``derive_ieff``: dI = 1/r_aug - 1/r_base), NOT a second GPU pass.

Nothing in the frozen tree is edited on disk; the redirect is injected at runtime, the
way ``main/ecspr/validation/_ref_patch.py`` does. Outputs land in ``REF_ROOT/null/``.

Env: ml (torch+cuda for the SMW batch; numpy-2.x pickles).

    python build_reference_null.py --smoke          # 2 axes, tiny, validates the redirect
    python build_reference_null.py                  # full: canon.DRAW_SIZES, all elements
"""
from __future__ import annotations

import argparse
import importlib.util
import pickle
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
from fabfos import canon                                   # noqa: E402
import reference as ref                        # noqa: E402  (sibling reference module)

# The frozen betweenness/null tree the redirect drives.
RN = Path("/home/tony/agentic_workspace/projects/scadc/metabolic-modelling"
          "/main/metabolic-modelling/04_reaction_network")
FROZEN_CACHE = RN / "cache"

# Engine transform for the ieff derivation -- the SAME one cmd_solve uses for observed.
sys.path.insert(0, str(canon.ENGINE_LIB / "resources" / "lib"))
from ecspr_solver import derive_ieff           # noqa: E402

OUT_DIR = canon.REFERENCE_ROOT / "null"
WORK = OUT_DIR / "_work"                        # 17r's redirected CACHE (inputs + reff outputs)

ELEMENTS = list(canon.ELEMENTS)


# =====================================================================
# The reference-graph loader (mirrors _ref_patch.ref_load_bipartite)
# =====================================================================
def ref_load_bipartite(resolution: str):
    """Replacement for `_graphs.load_bipartite`: the honest reference skeleton, w_X<=0
    edges dropped, exactly as the frozen element branch does."""
    if resolution not in ELEMENTS:
        raise SystemExit(f"reference null asked for non-element resolution {resolution!r}")
    with open(canon.REFERENCE_GRAPH_DIR / f"mnx_bipartite_{resolution}.pkl", "rb") as fh:
        G = pickle.load(fh)
    wk = f"w_{resolution}"
    drop = [(u, v) for u, v, d in G.edges(data=True) if d.get(wk, 0) <= 0]
    G.remove_edges_from(drop)
    return G


# =====================================================================
# Stage the redirected CACHE: reference bases + reused graph-independent inputs
# =====================================================================
def stage_work_dir():
    """Populate WORK with what 17r reads from CACHE, redirected where the graph matters.

    base_{X}.pkl        -> the reference-fed base (graph-DEPENDENT: reference)
    *_draws.parquet     -> frozen draws           (graph-INDEPENDENT: reused verbatim)
    metag_orf_weights   -> frozen weights         (graph-INDEPENDENT: reused verbatim)
    biomass_dag_axes    -> canon.AXES_JSON        (graph-INDEPENDENT axis definitions)
    axes_testable_set4  -> reference testable      (graph-DEPENDENT: reference, C is a
                                                    strict subset of the incumbent's)
    """
    WORK.mkdir(parents=True, exist_ok=True)
    def link(src: Path, name: str):
        dst = WORK / name
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(src)

    for X in ELEMENTS:
        link(canon.SOLVE_BASE_DIR / f"base_{X}.pkl", f"base_{X}.pkl")
    for n in canon.DRAW_SIZES:
        link(FROZEN_CACHE / f"null_canonical_N{n}_draws.parquet",
             f"null_canonical_N{n}_draws.parquet")
    link(FROZEN_CACHE / "metag_orf_weights.pkl", "metag_orf_weights.pkl")
    link(canon.AXES_JSON, "biomass_dag_axes_set4.json")
    link(canon.AXES_TESTABLE_JSON, "axes_testable_set4.json")


# =====================================================================
# Load the frozen 17r under the redirect
# =====================================================================
def load_frozen_17r():
    """Import 17r_reff_null.py, then override its module-global CACHE + load_bipartite.

    17r's helpers (draws_parquet/null_tsv/load_axes) and run_smw all read the module
    global CACHE and the imported name load_bipartite; rebinding them in the module
    namespace redirects every read without editing the frozen file.
    """
    for p in (str(RN), str(RN / "reff")):
        if p not in sys.path:
            sys.path.insert(0, p)
    spec = importlib.util.spec_from_file_location("_frozen_17r", RN / "17r_reff_null.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.CACHE = WORK
    mod.load_bipartite = ref_load_bipartite
    return mod


# =====================================================================
# ieff null <- reff null (the observed solve's own column transform)
# =====================================================================
IEFF_COLS = ["element", "axis_id", "source_hub", "sink_hub", "null", "iter",
             "n_draw", "delta_ieff", "g_base", "g_aug", "n_added_rxns"]


def derive_ieff_null(reff_tsv: Path, ieff_tsv: Path):
    df = pd.read_csv(reff_tsv, sep="\t")
    di, gb, ga = zip(*[derive_ieff(rb, ra)
                       for rb, ra in zip(df["r_base"], df["r_aug"])])
    out = df.copy()
    out["delta_ieff"] = di
    out["g_base"] = gb
    out["g_aug"] = ga
    out[IEFF_COLS].to_csv(ieff_tsv, sep="\t", index=False)
    print(f"[ieff] {ieff_tsv.name}: {len(out):,} rows <- {reff_tsv.name}", flush=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--draw-sizes", nargs="+", type=int, default=list(canon.DRAW_SIZES))
    ap.add_argument("--elements", nargs="+", default=ELEMENTS)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--batch-size", type=int, default=512)
    ap.add_argument("--c-batch-size", type=int, default=128)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stage_work_dir()
    mod = load_frozen_17r()

    # Mirror 17r's arg surface for run_smw.
    smw_args = argparse.Namespace(
        draw_sizes=a.draw_sizes, elements=a.elements, device=a.device,
        batch_size=a.batch_size, c_batch_size=a.c_batch_size, smoke=a.smoke)

    for N in a.draw_sizes:
        reff_out = mod.null_tsv(N)              # WORK/reff_null_canonical_N{N}.tsv
        if reff_out.exists():
            reff_out.unlink()                   # 17r appends; a fresh file per run
        print(f"\n########## reff null smw  N={N}  (reference graph) ##########", flush=True)
        mod.run_smw(smw_args, N)
        # publish reff + derive ieff into OUT_DIR
        final_reff = OUT_DIR / f"reff_null_canonical_N{N}.tsv"
        final_reff.write_bytes(reff_out.read_bytes())
        derive_ieff_null(final_reff, OUT_DIR / f"ieff_null_canonical_N{N}.tsv")

    print(f"\n[done] reference nulls -> {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
