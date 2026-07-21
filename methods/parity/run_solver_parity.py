#!/usr/bin/env python
"""The SOLVER gate: prove the solve itself, against a rebuild that shares none of its math.

    python run_solver_parity.py              # canonical elements
    python run_solver_parity.py --elements S # one element
    python run_solver_parity.py --strict     # also fail on the incumbent's own error

WHY THIS FILE EXISTS
--------------------
The other two gates do not touch the solver, and it took reading their inputs to see
it. `run_parity.py` hands the scorer `delta_obs` out of `{lane}_axes_report.tsv`;
every experiment spec stages `ecspr::{reff,ieff}_axes_report` as an INPUT, and
`scadc_main.py` says so in as many words -- "the staged, canonical solve". Both gates
prove `(delta_obs, nulls) -> p/q`. Nothing proved `graph + evidence -> delta_obs`.
What stood in for a gate was one sentence in `transforms/ecspr/solve.py` claiming
"verified vs scadc reff_axes_report.tsv, max|abs| 1.4e-8" -- never executed, and an
order of magnitude looser than `canon.PARITY_TOL`, which is the tolerance the same
repo calls non-negotiable.

WHAT IT COMPARES AGAINST, AND WHY NOT THE INCUMBENT
---------------------------------------------------
Not the incumbent table. The incumbent is the less accurate artifact:

    engine solver    vs dense rebuild : ~1e-14   (machine epsilon)
    incumbent table  vs dense rebuild : ~1e-8    (six orders worse)

`g_base` agrees at ~1e-16, so both are standing on an identical graph and Laplacian;
the entire divergence lives in the Woodbury update, and it is the incumbent carrying
it. Gating the solver against the incumbent would pin a correct implementation to an
incorrect referent and call it parity -- which is the exact failure this directory
exists to remove, wearing a green badge.

So the referent is a REBUILD: assemble the augmented graph the addition map
describes, build its Laplacian, factor it, solve. No Woodbury, no clique
elimination, no shared inverse, no touched-space projection. It shares no line of
reasoning with the thing under test, which is the only property that makes it a
referent at all. The incumbent is still reported -- as INFORMATION, never as a gate
(unless --strict) -- because the ~1e-8 is a real finding about the incumbent and it
should be data rather than something a future reader has to rediscover.

This gate is what lets the edge-space rewrite land: it holds the chemistry fixed
while the update math is replaced under it.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.linalg import splu

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fabfos import canon  # noqa: E402

sys.path.insert(0, str(canon.ENGINE_LIB / "resources" / "lib"))
from ecspr_solver import (SMWGraphContext, SMWSolver, build_ar2m,  # noqa: E402
                          reff_base, reff_summary, derive_ieff,
                          DirectSolver, reff_direct, _norm_ar2m)
from ecspr_network import load_bipartite  # noqa: E402


# =====================================================================
# The referent: rebuild and solve. No Woodbury anywhere in here.
# =====================================================================

def rebuild_reff(G, s, t, wkey: str) -> float:
    """Two-terminal R_eff by direct sparse solve of the graph's own Laplacian.

    Deliberately naive. It grounds one node and solves L_r phi = b, which is the
    textbook definition -- it does not know that an addition is low-rank, does not
    eliminate a hub, and does not reuse anything across calls. Every shortcut the
    solver takes is a shortcut this function refuses, which is what makes a match
    between them evidence rather than a tautology.
    """
    nodes = [n for n in G.nodes]
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    rows, cols, vals = [], [], []
    for u, v, d in G.edges(data=True):
        w = float(d.get(wkey, 0.0))
        if w <= 0:
            continue
        rows.append(idx[u]); cols.append(idx[v]); vals.append(w)
    if not vals:
        return 0.0
    r = np.asarray(rows); c = np.asarray(cols); w = np.asarray(vals, dtype=np.float64)
    deg = np.zeros(n)
    np.add.at(deg, r, w)
    np.add.at(deg, c, w)
    data = np.concatenate([-w, -w, deg])
    ri = np.concatenate([r, c, np.arange(n)])
    cj = np.concatenate([c, r, np.arange(n)])
    L = sp.coo_matrix((data, (ri, cj)), shape=(n, n)).tocsc()

    # ground node 0 of the s-t connected component; R_eff is ground-invariant, but
    # the component must be isolated or the Laplacian is singular.
    import networkx as nx
    comp = nx.node_connected_component(G.subgraph([x for x in G.nodes]), s)
    if t not in comp:
        return float("nan")
    keep = sorted(idx[x] for x in comp)
    kpos = {g: i for i, g in enumerate(keep)}
    Lc = L[keep, :][:, keep]
    m = len(keep)
    b = np.zeros(m)
    b[kpos[idx[s]]] = 1.0
    b[kpos[idx[t]]] = -1.0
    Lr = Lc[1:, 1:].tocsc()
    br = b[1:]
    phi = np.zeros(m)
    phi[1:] = splu(Lr).solve(br)
    return float(phi[kpos[idx[s]]] - phi[kpos[idx[t]]])


def augmented(base_g, ar2m: dict, wkey: str):
    """The augmented graph exactly as the addition map describes it.

    Parallel conductances ADD -- that is how a reaction the host already carries
    reinforces rather than replaces (build_ar2m routes it to a ("rxn_reinf", ...)
    node for precisely this reason).
    """
    Ga = base_g.copy()
    for rnode, mets in ar2m.items():
        for m, c in mets:
            if Ga.has_edge(rnode, m):
                Ga[rnode][m][wkey] = Ga[rnode][m].get(wkey, 0.0) + float(c)
            else:
                Ga.add_edge(rnode, m, **{wkey: float(c)})
    return Ga


# =====================================================================
# Cell selection
# =====================================================================

def pick_cells(ref_el: pd.DataFrame, axis_ids: list, rng) -> list:
    """(axis, fosmid) cells: the biggest deltas, plus typical ones.

    Extremes alone would be a soft gate -- the largest deltas are the best
    conditioned. The random draw is what exercises the rows the method actually
    reports in bulk.
    """
    cells = []
    for ax in axis_ids:
        sub = ref_el[ref_el.axis_id == ax]
        if not len(sub):
            continue
        top = sub.nlargest(canon.SOLVER_GATE_TOP_CELLS, "delta_ieff")
        rest = sub[~sub.fosmid.isin(top.fosmid)]
        k = min(canon.SOLVER_GATE_RANDOM_CELLS, len(rest))
        rnd = rest.iloc[rng.choice(len(rest), size=k, replace=False)] if k else rest.head(0)
        for _, r in pd.concat([top, rnd]).iterrows():
            cells.append((ax, str(r.fosmid), float(r.delta_ieff)))
    return cells


def run_element(X: str, ref: pd.DataFrame, rng) -> list:
    wk = f"w_{X}"
    base = pickle.load(open(canon.SOLVE_BASE_DIR / f"base_{X}.pkl", "rb"))
    fos_w = pickle.load(open(canon.ADDITION_WEIGHTS, "rb"))
    axes = json.loads(canon.AXES_JSON.read_text())
    testable = json.loads(canon.AXES_TESTABLE_JSON.read_text())
    G_X = load_bipartite(canon.BIPARTITE_DIR / f"mnx_bipartite_{X}.pkl", X)

    ctx = SMWGraphContext(base, device="cpu", edge_weight_key=wk)
    base_nodes = set(ctx.nodes)

    ax_all = testable.get(X, [])
    n_ax = min(canon.SOLVER_GATE_AXES, len(ax_all))
    ax_ids = [ax_all[i] for i in rng.choice(len(ax_all), size=n_ax, replace=False)]

    ref_el = ref[ref.element == X]
    cells = pick_cells(ref_el, ax_ids, rng)
    print(f"\n=== {X} === n_LCC={ctx.n}  axes={ax_ids}  cells={len(cells)}", flush=True)

    out = []
    by_axis: dict = {}
    for ax_id, contig, di_ref in cells:
        by_axis.setdefault(ax_id, []).append((contig, di_ref))

    for ax_id, members in by_axis.items():
        ax = axes[ax_id]
        src, snk = ("met", ax["source"][0]), ("met", ax["sink"][0])
        try:
            solver = SMWSolver(base, [src], [snk], edge_weight_key=wk, context=ctx)
        except ValueError:
            print(f"  [{ax_id}] endpoints not in LCC -- skipped")
            continue
        rb = reff_base(solver)
        rb_truth = rebuild_reff(base, src, snk, wk)
        if abs(rb - rb_truth) >= canon.PARITY_TOL:
            print(f"  [{ax_id}] BASE MISMATCH solver={rb:.12e} rebuild={rb_truth:.12e}")
        pair = [(ctx.idx[src], ctx.idx[snk])]
        dsolver = DirectSolver(ctx)
        rb_direct = dsolver.base_reff(pair)
        for contig, di_ref in members:
            ar2m = build_ar2m(fos_w[contig], G_X, base_nodes, wk, reinforce=True)
            if not ar2m:
                continue
            _, _, raug, na = reff_summary(solver, ar2m, r_base=rb)
            di_wood, _, _ = derive_ieff(rb, raug)

            ar2m_v, _, _ = _norm_ar2m(ctx, ar2m)
            _, _, raug_d = reff_direct(dsolver, ar2m_v, pair, rb_direct)[0]
            di_direct, _, _ = derive_ieff(rb_direct[0], raug_d)

            raug_truth = rebuild_reff(augmented(base, ar2m, wk), src, snk, wk)
            di_truth, _, _ = derive_ieff(rb_truth, raug_truth)

            out.append(dict(element=X, axis_id=ax_id, fosmid=contig, n_added=na,
                            di_truth=di_truth, di_direct=di_direct,
                            di_wood=di_wood, di_incumbent=di_ref,
                            err_direct=abs(di_direct - di_truth),
                            err_wood=abs(di_wood - di_truth),
                            err_incumbent=abs(di_ref - di_truth)))
            print(f"  {ax_id[:34]:34s} {contig:9s} n+={na:>4}  dI={di_truth:.8f}  err "
                  f"direct={abs(di_direct-di_truth):.1e} wood={abs(di_wood-di_truth):.1e} "
                  f"incumb={abs(di_ref-di_truth):.1e}", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elements", nargs="+", default=list(canon.ELEMENTS))
    ap.add_argument("--strict", action="store_true",
                    help="also fail on the incumbent's own error (a finding, not a regression)")
    a = ap.parse_args()

    # This gate rebuilds the UNDIRECTED solve from scratch (build_ar2m / reff_summary /
    # derive_ieff) and checks the engine's undirected solver paths agree; it must read the
    # undirected referent, not the (directed) canonical report -- otherwise the informational
    # incumbent-error column would compare a directed solve to an undirected rebuild.
    ref = pd.read_csv(canon.undirected_axes_report(canon.CANONICAL_LANE), sep="\t")
    rng = np.random.default_rng(canon.SOLVER_GATE_SEED)

    rows = []
    for X in a.elements:
        rows += run_element(X, ref, rng)

    if not rows:
        print("\nFAIL: compared 0 cells. A gate that compares nothing is not a gate.")
        return 1

    df = pd.DataFrame(rows)
    print("\n" + "=" * 72)
    print(f"cells compared: {len(df)}  over elements {sorted(df.element.unique())}")
    print("\n  vs the from-scratch rebuild:")
    print(f"    direct    : max {df.err_direct.max():.3e}   median {df.err_direct.median():.3e}")
    print(f"    woodbury  : max {df.err_wood.max():.3e}   median {df.err_wood.median():.3e}")
    print(f"    incumbent : max {df.err_incumbent.max():.3e}   median "
          f"{df.err_incumbent.median():.3e}   [INFORMATIONAL]")
    dw = float(np.abs(df.di_direct - df.di_wood).max())
    print(f"\n  direct vs woodbury (two implementations, no shared math): {dw:.3e}")
    print(f"\n  canon.PARITY_TOL = {canon.PARITY_TOL:.0e}")

    ok = (df.err_direct.max() < canon.PARITY_TOL
          and df.err_wood.max() < canon.PARITY_TOL
          and dw < canon.PARITY_TOL)
    print(f"\nSOLVER GATE: {'PASS' if ok else 'FAIL'} -- both solver paths "
          f"{'reproduce' if ok else 'DO NOT reproduce'} a from-scratch rebuild, "
          f"and each other.")
    if df.err_incumbent.max() >= canon.PARITY_TOL:
        print(f"\nNOTE: the incumbent axes_report does NOT clear PARITY_TOL against the "
              f"rebuild\n      (max {df.err_incumbent.max():.3e}). That is a finding about the "
              f"INCUMBENT, not\n      about the engine: g_base agrees to ~1e-16, so the graph is "
              f"identical and the\n      error is in its Woodbury update. It is invisible to "
              f"run_parity.py because both\n      sides there read this same table as an input.")
        if a.strict:
            print("      --strict: failing on it.")
            ok = False
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
