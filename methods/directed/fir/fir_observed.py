"""Directed OBSERVED solve on fir (new engine): per-fosmid additions to each canonical axis base,
rectified-network C_eff. Mirrors ecspr_network.cmd_solve_directed exactly (fosmid=unit,
null='obs', iter=0) but reads ratios from a pickle (no pyarrow) and parallelises over
(axis, fosmid-chunk). One invocation does all elements -> one reff/ieff pair.

    python fir_observed.py --stage-dir $SC --workers 64 \
        --out-reff $SC/obs/reff_axes_report.tsv --out-ieff $SC/obs/ieff_axes_report.tsv
"""
from __future__ import annotations
import argparse, json, os, pickle, sys, time
from pathlib import Path
import multiprocessing as mp

REFF_COLS = ["element", "axis_id", "source_hub", "sink_hub", "fosmid", "null", "iter",
             "delta_reff", "r_base", "r_aug", "n_added_rxns"]
IEFF_COLS = ["element", "axis_id", "source_hub", "sink_hub", "fosmid", "null", "iter",
             "delta_ieff", "g_base", "g_aug", "n_added_rxns"]
_G = {}


def _solve_unit(unit):
    ai, chunk_keys = unit
    onet, additions, axes_list, axes, X, tol = (
        _G["onet"], _G["additions"], _G["axes_list"], _G["axes"], _G["X"], _G["tol"])
    sub = {k: additions[k] for k in chunk_keys}
    rr, ri = [], []
    for row in directed_solve_element(onet, [axes_list[ai]], sub, warm=True, tol=tol):
        ax = axes[row["axis_id"]]; src, snk = ax["source"][0], ax["sink"][0]
        fos = row["unit"]
        rr.append((X, row["axis_id"], src, snk, fos, "obs", 0,
                   row["delta_reff"], row["r_base"], row["r_aug"], row["n_added_rxns"]))
        ri.append((X, row["axis_id"], src, snk, fos, "obs", 0,
                   row["delta_ieff"], row["g_base"], row["g_aug"], row["n_added_rxns"]))
    return rr, ri


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage-dir", required=True)
    ap.add_argument("--elements", nargs="+", default=["C", "N", "S", "P"])
    ap.add_argument("--workers", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--out-reff", required=True)
    ap.add_argument("--out-ieff", required=True)
    ap.add_argument("--force-noop", action="store_true")
    ap.add_argument("--tol", type=float, default=None)
    a = ap.parse_args()
    SC = Path(a.stage_dir)
    sys.path.insert(0, str(SC / "lib"))
    global orient_and_weight, directed_solve_element, build_ar2m_directed
    global load_bipartite, _base_lcc, load_reaction_roles
    from ecspr_directed import (orient_and_weight, directed_solve_element,  # noqa
                                build_ar2m_directed)
    from ecspr_network import load_bipartite, _base_lcc, load_reaction_roles  # noqa

    roles = load_reaction_roles(str(SC / "reac_prop.tsv"))
    ratios = {} if a.force_noop else pickle.load(open(SC / "ratios.pkl", "rb"))
    axes = json.loads((SC / "axes.json").read_text())
    testable = json.loads((SC / "axes_testable.json").read_text())
    fos_w = pickle.load(open(SC / "fosmid_addition_weights.pkl", "rb"))
    W = max(1, a.workers)

    all_rr, all_ri = [], []
    for X in a.elements:
        t0 = time.time(); wk = f"w_{X}"
        base = pickle.load(open(SC / "ref" / "solve" / f"base_{X}.pkl", "rb"))
        lcc = _base_lcc(base, wk)
        onet = orient_and_weight(lcc, X, roles, ratios, force_noop=a.force_noop)
        base_nodes = set(onet.nodes)
        G_X = load_bipartite(SC / "ref" / "graph" / f"mnx_bipartite_{X}.pkl", X)
        additions = {}
        for c, w in fos_w.items():
            e = build_ar2m_directed(w, G_X, base_nodes, wk, roles, ratios,
                                    reinforce=True, force_noop=a.force_noop)
            if e:
                additions[c] = e
        ax_ids = [aid for aid in testable.get(X, []) if aid in axes]
        axes_list = []
        for aid in ax_ids:
            ax = axes[aid]; src, snk = ("met", ax["source"][0]), ("met", ax["sink"][0])
            if src in onet.idx and snk in onet.idx and src != snk:
                axes_list.append((aid, src, snk))
        _G.update(onet=onet, additions=additions, axes_list=axes_list, axes=axes, X=X, tol=a.tol)
        keys = list(additions.keys())
        n_ax = len(axes_list)
        n_dc = max(1, min(8, -(-(2 * W) // max(n_ax, 1))))
        dcsz = (len(keys) + n_dc - 1) // n_dc
        dchunks = [keys[i:i + dcsz] for i in range(0, len(keys), dcsz)] or [[]]
        units = [(ai, dc) for ai in range(n_ax) for dc in dchunks]
        print(f"[obs {X}] n_LCC={onet.n} fosmids={len(keys)} axes={n_ax} units={len(units)} "
              f"orient={time.time()-t0:.1f}s", flush=True)
        with mp.get_context("fork").Pool(W) as pool:
            for rr, ri in pool.map(_solve_unit, units):
                all_rr.extend(rr); all_ri.extend(ri)

    Path(a.out_reff).parent.mkdir(parents=True, exist_ok=True)
    with open(a.out_reff, "w") as fr, open(a.out_ieff, "w") as fi:
        fr.write("\t".join(REFF_COLS) + "\n"); fi.write("\t".join(IEFF_COLS) + "\n")
        for row in all_rr:
            fr.write("\t".join(str(v) for v in row) + "\n")
        for row in all_ri:
            fi.write("\t".join(str(v) for v in row) + "\n")
    print(f"wrote {len(all_rr)} rows each -> {a.out_reff}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
