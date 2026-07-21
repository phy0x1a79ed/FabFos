"""Self-contained directed-null driver for fir CPU nodes.

Reuses the ENGINE primitives (orient_and_weight / build_ar2m_directed /
directed_solve_element) exactly as main/fabfos/directed/build_directed_null.py does, but
with every path parametrised and a fork-based multiprocessing pool so one node's cores are
all used. Inputs are pre-resolved to pickles (draws->weights, direction->ratios) so no
pyarrow is needed. One invocation = one (element, draw-size) shard; output is a per-shard
TSV pair with the exact undirected-null schema, resumable (a complete shard is skipped).

    python fir_null.py --stage-dir $SC --element C --draw-size 14 --workers 48 \
        --out-reff $SC/out/reff_C_N14.tsv --out-ieff $SC/out/ieff_C_N14.tsv
"""
from __future__ import annotations
import argparse, json, os, pickle, sys, time
from pathlib import Path
import multiprocessing as mp

REFF_COLS = ["element", "axis_id", "source_hub", "sink_hub", "null", "iter",
             "n_draw", "delta_reff", "r_base", "r_aug", "n_added_rxns"]
IEFF_COLS = ["element", "axis_id", "source_hub", "sink_hub", "null", "iter",
             "n_draw", "delta_ieff", "g_base", "g_aug", "n_added_rxns"]

# module globals shared into forked workers (set before the pool is created)
_G = {}


def _init_engine(stage_dir):
    sys.path.insert(0, str(Path(stage_dir) / "lib"))
    global orient_and_weight, directed_solve_element, build_ar2m_directed
    global load_bipartite, _base_lcc, load_reaction_roles
    from ecspr_directed import (orient_and_weight, directed_solve_element,   # noqa
                                build_ar2m_directed)
    from ecspr_network import load_bipartite, _base_lcc, load_reaction_roles  # noqa


def _solve_unit(unit):
    # unit = (axis_index, draw_chunk_keys): ONE axis over a chunk of draws, so the per-axis
    # base solve is computed once per (axis, chunk) -- keep chunks few so base is not redone
    # many times (base << the chunk's augmented solves).
    ai, chunk_keys = unit
    onet, additions, axes_list, axes, X, N, tol = (
        _G["onet"], _G["additions"], _G["axes_list"], _G["axes"], _G["X"], _G["N"], _G["tol"])
    sub = {k: additions[k] for k in chunk_keys}
    rr, ri = [], []
    for row in directed_solve_element(onet, [axes_list[ai]], sub, warm=True, tol=tol):
        ax = axes[row["axis_id"]]; kind, it = row["unit"]
        src, snk = ax["source"][0], ax["sink"][0]
        rr.append((X, row["axis_id"], src, snk, kind, it, N,
                   row["delta_reff"], row["r_base"], row["r_aug"], row["n_added_rxns"]))
        ri.append((X, row["axis_id"], src, snk, kind, it, N,
                   row["delta_ieff"], row["g_base"], row["g_aug"], row["n_added_rxns"]))
    return rr, ri


def _expected_rows(stage_dir, X, N):
    testable = json.loads((Path(stage_dir) / "axes_testable.json").read_text())
    axes = json.loads((Path(stage_dir) / "axes.json").read_text())
    draws = pickle.load(open(Path(stage_dir) / "draws_resolved.pkl", "rb"))
    n_ax = sum(1 for a in testable.get(X, []) if a in axes)
    return n_ax * len(draws[N]["keys"])


def _complete(path, expected):
    if not Path(path).exists():
        return False
    with open(path) as fh:
        n = sum(1 for _ in fh) - 1  # minus header
    return n >= expected


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage-dir", required=True)
    ap.add_argument("--element", required=True)
    ap.add_argument("--draw-size", type=int, required=True)
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--out-reff", required=True)
    ap.add_argument("--out-ieff", required=True)
    ap.add_argument("--force-noop", action="store_true")
    ap.add_argument("--tol", type=float, default=None)
    ap.add_argument("--smoke-draws", type=int, default=0, help="use only first K draws (calibration)")
    ap.add_argument("--smoke-axes", type=int, default=0, help="use only first K axes (calibration)")
    ap.add_argument("--resume", action="store_true")
    a = ap.parse_args()
    SC = Path(a.stage_dir); X = a.element; N = a.draw_size
    t0 = time.time()

    if a.resume and not a.smoke_draws and _complete(a.out_reff, _expected_rows(a.stage_dir, X, N)) \
            and _complete(a.out_ieff, _expected_rows(a.stage_dir, X, N)):
        print(f"[{X}/N{N}] already complete, skip", flush=True); return 0

    _init_engine(a.stage_dir)
    roles = load_reaction_roles(str(SC / "reac_prop.tsv"))
    ratios = {} if a.force_noop else pickle.load(open(SC / "ratios.pkl", "rb"))
    axes = json.loads((SC / "axes.json").read_text())
    testable = json.loads((SC / "axes_testable.json").read_text())
    resolved = pickle.load(open(SC / "draws_resolved.pkl", "rb"))[N]
    keys, draw_w = resolved["keys"], resolved["draw_w"]
    if a.smoke_draws:
        keys = keys[:a.smoke_draws]

    wk = f"w_{X}"
    base = pickle.load(open(SC / "ref" / "solve" / f"base_{X}.pkl", "rb"))
    lcc = _base_lcc(base, wk)
    onet = orient_and_weight(lcc, X, roles, ratios, force_noop=a.force_noop)
    base_nodes = set(onet.nodes)
    G_X = load_bipartite(SC / "ref" / "graph" / f"mnx_bipartite_{X}.pkl", X)
    t_orient = time.time() - t0

    t = time.time()
    additions = {k: build_ar2m_directed(draw_w[k], G_X, base_nodes, wk, roles, ratios,
                                        reinforce=True, force_noop=a.force_noop) for k in keys}
    t_add = time.time() - t

    ax_ids = [aid for aid in testable.get(X, []) if aid in axes]
    if a.smoke_axes:
        ax_ids = ax_ids[:a.smoke_axes]
    axes_list = []
    for aid in ax_ids:
        ax = axes[aid]; src, snk = ("met", ax["source"][0]), ("met", ax["sink"][0])
        if src in onet.idx and snk in onet.idx and src != snk:
            axes_list.append((aid, src, snk))

    _G.update(onet=onet, additions=additions, axes_list=axes_list, axes=axes,
              X=X, N=N, tol=a.tol)
    W = max(1, a.workers)
    n_ax = len(axes_list)
    ks = list(keys)
    # draw-chunks per axis: enough total units (axis x chunk) to fill 2*W workers, but few
    # enough that the per-axis base solve is redone only a handful of times (base << aug*chunk).
    n_dc = max(1, min(16, -(-(2 * W) // max(n_ax, 1))))
    dcsz = (len(ks) + n_dc - 1) // n_dc
    draw_chunks = [ks[i:i + dcsz] for i in range(0, len(ks), dcsz)]
    units = [(ai, dc) for ai in range(n_ax) for dc in draw_chunks]
    print(f"[{X}/N{N}] n_LCC={onet.n} axes={n_ax} draws={len(ks)} orient={t_orient:.1f}s "
          f"add={t_add:.1f}s workers={W} draw_chunks={len(draw_chunks)} units={len(units)}", flush=True)

    ts = time.time()
    if W == 1:
        results = [_solve_unit(u) for u in units]
    else:
        with mp.get_context("fork").Pool(W) as pool:
            results = pool.map(_solve_unit, units)
    Path(a.out_reff).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out_reff, "w") as fr, open(a.out_ieff, "w") as fi:
        fr.write("\t".join(REFF_COLS) + "\n"); fi.write("\t".join(IEFF_COLS) + "\n")
        nb = 0
        for rr, ri in results:
            for row in rr:
                fr.write("\t".join(str(v) for v in row) + "\n"); nb += 1
            for row in ri:
                fi.write("\t".join(str(v) for v in row) + "\n")
    dt = time.time() - ts
    print(f"[{X}/N{N}] {nb:,} rows solve={dt:.1f}s "
          f"({dt/max(nb,1)*1e3:.2f} ms/cell) total={time.time()-t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
