#!/usr/bin/env python
"""Directed (diode) OBSERVED solve on the reference graph -- the delta_obs the significance
scorer joins against the directed null (build_directed_null.py). The lockstep twin of that
null: both call the SAME engine primitive (``directed_solve_element``) on the SAME oriented
reference graph, so the observed signal and the null it is scored against are one model by
construction -- the whole point of the directed promotion.

This is the directed counterpart of the undirected reference solve
(``canon.reference_axes_report``); it emits that report's schema EXACTLY (``fosmid`` unit,
``null="obs"``, ``iter=0``) so ``ecspr_significance.py`` and the canon accessors consume it
unchanged. Each fosmid's evidence-network additions (``canon.ADDITION_WEIGHTS``) reinforce
the host base graph; every canonical axis is solved as a rectified-network ``C_eff`` through the
diode. ``--force-noop`` forces every ratio to 1.0 -> must reproduce the UNDIRECTED reference
report (the observed-side symmetric-limit lockstep, gated by directed_parity.py).

Mirrors build_directed_null.py's fork parallelisation: the parent builds the oriented net +
per-fosmid additions once, workers are forked (inheriting them copy-on-write) and receive
only fosmid-KEY slices -- the graph is never pickled. Contiguous slices + ordered map results
keep the output byte-equivalent to the serial path.

This module is the transform-ready precursor to the engine's ``ecsprDirected`` /
``solve_directed`` metasmith transform: its core (orient -> per-unit additions ->
directed_solve_element) is exactly the transform body, with canon supplying the pinned inputs
here that the transform would receive as staged parents.

Env: ml (numpy/pandas/scipy/networkx + the engine on sys.path). The full cluster run uses the
path-parametrised sibling directed/fir/fir_observed.py; this is the canon-integrated twin.

    # all elements, canon inputs
    python build_directed_observed.py --out-reff reff_axes_report.tsv --out-ieff ieff_axes_report.tsv
    # symmetric-limit lockstep (must reproduce the undirected reference report)
    python build_directed_observed.py --force-noop --out-reff /tmp/r.tsv --out-ieff /tmp/i.tsv
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

# Single-threaded BLAS before numpy/scipy load (inherited by forks): we fan single-threaded
# solves across cores, so per-process BLAS threads would oversubscribe. Measured >10x slower.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = Path(__file__).resolve().parent
from fabfos import canon                                          # noqa: E402

sys.path.insert(0, str(canon.ENGINE_LIB / "resources" / "lib"))
from ecspr_directed import (orient_and_weight, directed_solve_element,      # noqa: E402
                            build_ar2m_directed)
from ecspr_network import (load_reaction_roles, load_direction_ratios,      # noqa: E402
                           load_bipartite, _base_lcc)

# The undirected reference report's schema, verbatim (canon.reference_axes_report columns):
# the directed observed solve is a drop-in replacement for it.
REFF_COLS = ["element", "axis_id", "source_hub", "sink_hub", "fosmid", "null", "iter",
             "delta_reff", "r_base", "r_aug", "n_added_rxns"]
IEFF_COLS = ["element", "axis_id", "source_hub", "sink_hub", "fosmid", "null", "iter",
             "delta_ieff", "g_base", "g_aug", "n_added_rxns"]


def _open(path: Path, cols):
    path.parent.mkdir(parents=True, exist_ok=True)
    fh = open(path, "w", buffering=1)
    fh.write("\t".join(cols) + "\n")
    return fh


def _row_tuples(row, X, axes):
    """The (reff, ieff) TSV field tuples for one directed_solve_element row (observed schema:
    the unit IS the fosmid, null='obs', iter=0). Single place the schema is materialised, so
    the serial and worker paths emit identical bytes."""
    ax = axes[row["axis_id"]]
    fos = row["unit"]
    src, snk = ax["source"][0], ax["sink"][0]
    reff = (X, row["axis_id"], src, snk, fos, "obs", 0,
            row["delta_reff"], row["r_base"], row["r_aug"], row["n_added_rxns"])
    ieff = (X, row["axis_id"], src, snk, fos, "obs", 0,
            row["delta_ieff"], row["g_base"], row["g_aug"], row["n_added_rxns"])
    return reff, ieff


def _chunks(keys, workers, per_worker=4):
    """Contiguous, in-order fosmid-key slices (~per_worker per worker). Contiguity + ordered
    map results keep output byte-equivalent to the serial path."""
    n = len(keys)
    nchunks = min(n, max(1, workers * per_worker)) or 1
    size = (n + nchunks - 1) // nchunks
    return [keys[i:i + size] for i in range(0, n, size)] or [[]]


# Set by the parent per element BEFORE the fork pool is created; workers read via CoW.
_WORK: dict = {}


def _solve_key_slice(slice_keys):
    """Worker: solve one contiguous slice of fosmids on the fork-inherited graph."""
    w = _WORK
    onet, axes_list, additions = w["onet"], w["axes_list"], w["additions"]
    axes, X, tol = w["axes"], w["X"], w["tol"]
    sub = {k: additions[k] for k in slice_keys}
    reff_rows, ieff_rows = [], []
    for row in directed_solve_element(onet, axes_list, sub, warm=True, tol=tol):
        r, i = _row_tuples(row, X, axes)
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
    fos_w = pickle.load(open(args.addition, "rb"))

    workers = args.workers if getattr(args, "workers", None) else (os.cpu_count() or 1)

    fr = _open(Path(args.out_reff), REFF_COLS)
    fi = _open(Path(args.out_ieff), IEFF_COLS)
    n_rows = 0
    try:
        for X in args.elements:
            wk = f"w_{X}"
            with open(Path(args.base_dir) / f"base_{X}.pkl", "rb") as fh:
                base = pickle.load(fh)
            t = time.time()
            lcc = _base_lcc(base, wk)
            onet = orient_and_weight(lcc, X, roles, ratios, force_noop=args.force_noop)
            base_nodes = set(onet.nodes)
            G_X = load_bipartite(Path(args.bipartite_dir) / f"mnx_bipartite_{X}.pkl", X)
            # The OBSERVED solve skips empty fosmids (matching the undirected cmd_solve): a
            # fosmid touching no base reaction contributes no row. (The null, by contrast,
            # keeps empty draws as delta-0 rows -- the two paths differ on purpose.)
            additions = {}
            for c, w in fos_w.items():
                e = build_ar2m_directed(w, G_X, base_nodes, wk, roles, ratios,
                                        reinforce=True, force_noop=args.force_noop)
                if e:
                    additions[c] = e
            keys = list(additions.keys())
            ax_ids = [a for a in testable.get(X, []) if a in axes]
            axes_list = []
            for ax_id in ax_ids:
                ax = axes[ax_id]
                src, snk = ("met", ax["source"][0]), ("met", ax["sink"][0])
                if src in onet.idx and snk in onet.idx and src != snk:
                    axes_list.append((ax_id, src, snk))
            print(f"=== [{X}] n_LCC={onet.n} axes {len(axes_list)} "
                  f"fosmids-with-additions {len(keys)}/{len(fos_w)} "
                  f"orient {time.time()-t:.1f}s ===", flush=True)

            t = time.time()
            nb = 0
            if workers <= 1 or not keys:
                for row in directed_solve_element(onet, axes_list, additions,
                                                  warm=True, tol=args.tol):
                    r_row, i_row = _row_tuples(row, X, axes)
                    _write_rows(fr, fi, [r_row], [i_row])
                    nb += 1
            else:
                global _WORK
                _WORK = dict(onet=onet, axes_list=axes_list, additions=additions,
                             axes=axes, X=X, tol=args.tol)
                slices = _chunks(keys, workers)
                ctx = mp.get_context("fork")
                with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as ex:
                    for reff_rows, ieff_rows in ex.map(_solve_key_slice, slices):
                        _write_rows(fr, fi, reff_rows, ieff_rows)
                        nb += len(reff_rows)
                _WORK = {}
            n_rows += nb
            print(f"  [{X}] {nb:,} rows in {time.time()-t:.1f}s (workers={workers})",
                  flush=True)
    finally:
        fr.close()
        fi.close()
    print(f"[obs-directed] wrote {args.out_reff} + {args.out_ieff} ({n_rows:,} rows each)")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", nargs="+", default=list(canon.ELEMENTS))
    ap.add_argument("--base-dir", default=str(canon.SOLVE_BASE_DIR))
    ap.add_argument("--bipartite-dir", default=str(canon.BIPARTITE_DIR))
    ap.add_argument("--axes", default=str(canon.AXES_JSON))
    ap.add_argument("--testable", default=str(canon.AXES_TESTABLE_JSON))
    ap.add_argument("--addition", default=str(canon.ADDITION_WEIGHTS),
                    help="per-fosmid evidence-network addition weights pkl")
    ap.add_argument("--direction", default=str(canon.REFERENCE_ROOT / "direction.parquet"),
                    help="frozen directionality ensemble parquet (mnxr, ratio)")
    ap.add_argument("--roles", default=str(canon.DIR_REAC_PROP))
    ap.add_argument("--out-reff", required=True)
    ap.add_argument("--out-ieff", required=True)
    ap.add_argument("--force-noop", action="store_true",
                    help="ratios->1.0: must reproduce the undirected reference report")
    ap.add_argument("--tol", type=float, default=None)
    ap.add_argument("--workers", type=int, default=None,
                    help="fork workers to fan the fosmids across (default: os.cpu_count(); "
                         "1 = serial path, byte-equivalent reference)")
    a = ap.parse_args()
    return run(a)


if __name__ == "__main__":
    raise SystemExit(main())
