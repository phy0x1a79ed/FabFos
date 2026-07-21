#!/usr/bin/env python
"""Directed (diode) null on the reference graph -- the lockstep companion to the directed
observed solve (engine ``cmd_solve_directed``).

The frozen DRAWS (``null_canonical_N{n}_draws.parquet``) and ``metag_orf_weights.pkl`` are
graph- AND solver-independent ORF samples; they are reused VERBATIM from the frozen cache,
exactly as the undirected reference null does. The ONLY change from that null is the solver:
each draw is resolved to ``{mnxr: E}`` (``_null_canon.draw_to_weights``, byte-for-byte the
observed scheme restricted to one draw's ORFs), oriented per reaction roles + the
directionality ensemble's ``g_rev/g_fwd`` ratio, and solved as a rectified-network ``C_eff``
through the engine's ``directed_solve_element`` -- the SAME primitive the observed directed
solve calls. So ``delta_obs`` and this null share one model by construction; that is the
lockstep the significance scorer relies on.

Why not 17r's GPU ``reff_batch``: the diode has a per-cell active set and does not batch.
Instead each solve warm-starts from the undirected grounded potential (``OrientedNet``,
convex energy so the seed changes only the iteration count, never the answer) and the work
shards over (element x draw-size). ``--force-noop`` forces every ratio to 1.0 -> this must
reproduce the UNDIRECTED reference null (the null-side symmetric-limit lockstep gate).

Schema is the undirected null's COLS exactly, so ``ecspr_significance.py`` consumes it
unchanged. The ieff null is the same ``derive_ieff`` column transform the observed solve
uses (done inside ``directed_solve_element``), not a second pass.

Env: ml (numpy/pandas/scipy/networkx + the engine on sys.path). Sharded on Sockeye for the
full run; a single (N, element) shard is a self-contained invocation.

    # one shard (Sockeye array task): size 14, element C
    python build_directed_null.py --draw-sizes 14 --elements C \
        --out-reff out/reff_null_directed_N14_C.tsv --out-ieff out/ieff_null_directed_N14_C.tsv
    # local smoke: 2 axes, first 40 draws, one (N, element)
    python build_directed_null.py --draw-sizes 14 --elements S --smoke \
        --out-reff /tmp/r.tsv --out-ieff /tmp/i.tsv
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import pickle
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

# Single-threaded BLAS: the null fans single-threaded solves across FORK processes (one per
# core), so per-process BLAS threads would oversubscribe every core -- measured >10x slower
# unpinned. Set before numpy/scipy load (via pandas) so forked workers inherit it too.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
from fabfos import canon                                          # noqa: E402

sys.path.insert(0, str(canon.ENGINE_LIB / "resources" / "lib"))
from ecspr_directed import (orient_and_weight, directed_solve_element,      # noqa: E402
                            build_ar2m_directed)
from ecspr_network import (load_reaction_roles, load_direction_ratios,      # noqa: E402
                           load_bipartite, _base_lcc)

# Frozen draws + weights live beside the incumbent null (graph/solver-independent).
RN = Path("/home/tony/agentic_workspace/projects/scadc/metabolic-modelling"
          "/main/metabolic-modelling/04_reaction_network")
FROZEN_CACHE = RN / "cache"
for _p in (str(RN), str(RN / "reff")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from _null_canon import draw_to_weights                                     # noqa: E402

REFF_COLS = ["element", "axis_id", "source_hub", "sink_hub", "null", "iter",
             "n_draw", "delta_reff", "r_base", "r_aug", "n_added_rxns"]
IEFF_COLS = ["element", "axis_id", "source_hub", "sink_hub", "null", "iter",
             "n_draw", "delta_ieff", "g_base", "g_aug", "n_added_rxns"]


def _open(path: Path, cols):
    path.parent.mkdir(parents=True, exist_ok=True)
    fh = open(path, "w", buffering=1)
    fh.write("\t".join(cols) + "\n")
    return fh


def _row_tuples(row, X, N, axes):
    """The (reff, ieff) TSV field tuples for one directed_solve_element row -- the single
    place row schema is materialised, so the serial and worker paths emit identical bytes."""
    ax = axes[row["axis_id"]]
    kind, it = row["unit"]
    src, snk = ax["source"][0], ax["sink"][0]
    reff = (X, row["axis_id"], src, snk, kind, it, N,
            row["delta_reff"], row["r_base"], row["r_aug"], row["n_added_rxns"])
    ieff = (X, row["axis_id"], src, snk, kind, it, N,
            row["delta_ieff"], row["g_base"], row["g_aug"], row["n_added_rxns"])
    return reff, ieff


def _chunks(keys, workers, per_worker=4):
    """Contiguous, in-order draw-key slices (~``per_worker`` per worker for load balance).
    Contiguity + ProcessPoolExecutor.map's ordered results keep output byte-equivalent to the
    serial (draw-major) path; more chunks than workers lets a fast worker steal the tail."""
    n = len(keys)
    nchunks = min(n, max(1, workers * per_worker)) or 1
    size = (n + nchunks - 1) // nchunks
    return [keys[i:i + size] for i in range(0, n, size)]


# Set by the parent per (N, element) BEFORE the fork pool is created; workers read it via
# fork copy-on-write (never pickled). Holds the oriented net + the full additions dict + axes,
# so a worker receives only a small slice of draw KEYS -- not the graph or the edge lists.
_WORK: dict = {}


def _solve_key_slice(slice_keys):
    """Worker: solve one contiguous slice of draws on the fork-inherited graph. Reads the
    graph/additions from ``_WORK`` (CoW), builds this slice's sub-additions locally, and
    returns the materialised (reff, ieff) row tuples in draw-major order."""
    w = _WORK
    onet, axes_list, additions = w["onet"], w["axes_list"], w["additions"]
    axes, X, N, tol = w["axes"], w["X"], w["N"], w["tol"]
    sub = {k: additions[k] for k in slice_keys}
    reff_rows, ieff_rows = [], []
    for row in directed_solve_element(onet, axes_list, sub, warm=True, tol=tol):
        r, i = _row_tuples(row, X, N, axes)
        reff_rows.append(r)
        ieff_rows.append(i)
    return reff_rows, ieff_rows


def _write_rows(fr, fi, reff_rows, ieff_rows):
    for r in reff_rows:
        fr.write("\t".join(str(v) for v in r) + "\n")
    for i in ieff_rows:
        fi.write("\t".join(str(v) for v in i) + "\n")


def run(args) -> int:
    roles = load_reaction_roles(args.roles)
    ratios = {} if args.force_noop else load_direction_ratios(args.direction)
    axes = json.loads(Path(args.axes).read_text())
    testable = json.loads(Path(args.testable).read_text())
    mg_w = pickle.load(open(FROZEN_CACHE / "metag_orf_weights.pkl", "rb"))

    workers = args.workers if getattr(args, "workers", None) else (os.cpu_count() or 1)

    fr = _open(Path(args.out_reff), REFF_COLS)
    fi = _open(Path(args.out_ieff), IEFF_COLS)
    n_rows = 0
    try:
        for N in args.draw_sizes:
            dp = FROZEN_CACHE / f"null_canonical_N{N}_draws.parquet"
            if not dp.exists():
                raise SystemExit(f"need frozen draws {dp}")
            draws = pd.read_parquet(dp)
            if args.smoke:
                draws = draws.head(args.smoke_draws)
            keys = [(r.null, int(r.iter)) for r in draws.itertuples()]
            draw_w = {(r.null, int(r.iter)): draw_to_weights(r.orf_list, mg_w)
                      for r in draws.itertuples()}
            n_ne = sum(1 for w in draw_w.values() if w)
            print(f"[null-directed] N={N}: {len(keys):,} draws "
                  f"({n_ne:,} non-empty)", flush=True)

            for X in args.elements:
                wk = f"w_{X}"
                with open(Path(args.base_dir) / f"base_{X}.pkl", "rb") as fh:
                    base = pickle.load(fh)
                t = time.time()
                lcc = _base_lcc(base, wk)
                onet = orient_and_weight(lcc, X, roles, ratios,
                                         force_noop=args.force_noop)
                base_nodes = set(onet.nodes)
                G_X = load_bipartite(Path(args.bipartite_dir) / f"mnx_bipartite_{X}.pkl", X)
                # EVERY draw becomes a unit, INCLUDING empty ones: the undirected 17r
                # null writes a row per (axis, draw) with delta=0 for a draw that touches
                # no base reaction, so the null distribution carries those zeros. An empty
                # edge list -> directed_solve_element's base-only branch -> delta 0, n=0,
                # matching exactly. (The OBSERVED solve, by contrast, skips empty fosmids,
                # matching the undirected cmd_solve -- the two paths differ on purpose.)
                additions = {k: build_ar2m_directed(draw_w[k], G_X, base_nodes, wk,
                                                     roles, ratios, reinforce=True,
                                                     force_noop=args.force_noop)
                             for k in keys}
                n_ne = sum(1 for e in additions.values() if e)
                ax_ids = [a for a in testable.get(X, []) if a in axes]
                if args.smoke:
                    ax_ids = ax_ids[:2]
                axes_list = []
                for ax_id in ax_ids:
                    ax = axes[ax_id]
                    src, snk = ("met", ax["source"][0]), ("met", ax["sink"][0])
                    if src in onet.idx and snk in onet.idx and src != snk:
                        axes_list.append((ax_id, src, snk))
                print(f"=== [N{N}/{X}] n_LCC={onet.n} axes {len(axes_list)} "
                      f"draws-with-additions {n_ne}/{len(keys)} "
                      f"orient {time.time()-t:.1f}s ===", flush=True)

                t = time.time()
                nb = 0
                if workers <= 1:
                    # Serial path: byte-for-byte the pre-parallel behaviour (the reference
                    # the fork path is validated against).
                    for row in directed_solve_element(onet, axes_list, additions,
                                                      warm=True, tol=args.tol):
                        r_row, i_row = _row_tuples(row, X, N, axes)
                        _write_rows(fr, fi, [r_row], [i_row])
                        nb += 1
                else:
                    # Fan the draws across fork workers. The parent has already built onet +
                    # additions (once); publish them to the module global so the forked
                    # workers inherit them copy-on-write and receive only draw-KEY slices --
                    # the graph/edge-lists are never pickled. Contiguous slices + ordered
                    # map results keep output byte-equivalent to the serial path.
                    global _WORK
                    _WORK = dict(onet=onet, axes_list=axes_list, additions=additions,
                                 axes=axes, X=X, N=N, tol=args.tol)
                    slices = _chunks(keys, workers)
                    ctx = mp.get_context("fork")
                    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as ex:
                        for reff_rows, ieff_rows in ex.map(_solve_key_slice, slices):
                            _write_rows(fr, fi, reff_rows, ieff_rows)
                            nb += len(reff_rows)
                    _WORK = {}
                n_rows += nb
                print(f"  [N{N}/{X}] {nb:,} rows in {time.time()-t:.1f}s "
                      f"(workers={workers})", flush=True)
    finally:
        fr.close()
        fi.close()
    print(f"[null-directed] wrote {args.out_reff} + {args.out_ieff} ({n_rows:,} rows each)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--draw-sizes", nargs="+", type=int, default=list(canon.DRAW_SIZES))
    ap.add_argument("--elements", nargs="+", default=list(canon.ELEMENTS))
    ap.add_argument("--base-dir", default=str(canon.SOLVE_BASE_DIR))
    ap.add_argument("--bipartite-dir", default=str(canon.BIPARTITE_DIR))
    ap.add_argument("--axes", default=str(canon.AXES_JSON))
    ap.add_argument("--testable", default=str(canon.AXES_TESTABLE_JSON))
    ap.add_argument("--direction", default=str(canon.REFERENCE_ROOT / "direction.parquet"),
                    help="frozen directionality ensemble parquet (mnxr, ratio)")
    ap.add_argument("--roles", default=str(canon.DIR_REAC_PROP))
    ap.add_argument("--out-reff", required=True)
    ap.add_argument("--out-ieff", required=True)
    ap.add_argument("--force-noop", action="store_true",
                    help="ratios->1.0: must reproduce the undirected reference null")
    ap.add_argument("--tol", type=float, default=None)
    ap.add_argument("--workers", type=int, default=None,
                    help="fork workers to fan the draws across (default: os.cpu_count(); "
                         "1 = serial path, byte-equivalent reference)")
    ap.add_argument("--smoke", action="store_true", help="2 axes/element")
    ap.add_argument("--smoke-draws", type=int, default=40)
    a = ap.parse_args()
    return run(a)


if __name__ == "__main__":
    raise SystemExit(main())
