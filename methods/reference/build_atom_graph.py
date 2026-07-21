"""Rebuild the per-element STAR atom graph FROM THE FROZEN REFERENCE table.

The solve runs on a met-rxn STAR: a reaction node ``("rxn", mnxr)`` joined to each of
its metabolites ``("met", mnxm)``, the edge weight ``w_X`` = the number of element-X atoms
that TRANSIT that reaction through that metabolite (substrate role + product role). That
weight is the only place the atom-atom mapping enters the production solve.

The incumbent graph (``04_reaction_network/cache/mnx_bipartite_{X}.pkl``) was weighted from
``aam_unified_expanded.tsv`` -- a mapper run that fabricated transits for structureless
stubs, gave atom-transit edges to non-molecules (electron/photon acceptors that carry no
atom), and treated pseudo-reactions as chemistry. The closure ledger REFUSES exactly those
(``reference.assert_reference``). This builder re-weights the SAME unweighted skeleton with
``w_X`` taken from the frozen ``atom_pairs.parquet`` instead:

    w_X(rxn, met) = sum of pair_w over reference rows with element==X, substrate != product,
                    and met in {substrate, product}.

``pair_w`` is the margin-conserving ensemble weight (consensus 1.0 -> integer atom count,
matching the incumbent on agreed maps; diluted < 1 -> fractional transit). Self-mapped atoms
(substrate == product: an atom that stays in one metabolite) are dropped, exactly as the
engine ``atom_edges`` drops ``u == v`` -- w_X counts CROSS-metabolite transit only.

Output: ``reference.GRAPH_DIR / mnx_bipartite_{X}.pkl`` -- a NEW location. The incumbent cache
stays frozen so ``run_solver_parity`` / ``run_parity`` keep validating the canonical solve;
this reference-fed graph is a separate, honest measurement, not a replacement.

Env: p312 (the skeleton is numpy-2.x pickled). Reads only frozen inputs; writes only GRAPH_DIR.
"""
from __future__ import annotations

import pickle
import sys
from pathlib import Path

import networkx as nx
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reference as ref  # noqa: E402


def per_element_weights(pairs: pd.DataFrame, element: str) -> dict:
    """{(mnxr, mnxm) -> w_X} = sum of pair_w over cross-metabolite transiting X atoms,
    both substrate and product roles. Vectorised over the 2.16M-row table."""
    p = pairs[(pairs.element == element) & (pairs.substrate != pairs.product)]
    sub = p[["mnxr", "substrate", "pair_w"]].rename(columns={"substrate": "mnxm"})
    prod = p[["mnxr", "product", "pair_w"]].rename(columns={"product": "mnxm"})
    both = pd.concat([sub, prod], ignore_index=True)
    agg = both.groupby(["mnxr", "mnxm"], sort=False)["pair_w"].sum()
    return {(str(r), str(m)): float(w) for (r, m), w in agg.items()}


def frozen_positive(element: str) -> int | None:
    """Count of w_X>0 edges in the FROZEN incumbent graph, for the honesty delta."""
    inc = ref.BASE_GRAPH_DIR / f"mnx_bipartite_{element}.pkl"
    if not inc.exists():
        return None
    with open(inc, "rb") as fh:
        G = pickle.load(fh)
    wk = f"w_{element}"
    return sum(1 for _, _, d in G.edges(data=True) if d.get(wk, 0) > 0)


def main() -> None:
    print(f"[build-atom-graph] loading skeleton {ref.SKELETON}")
    with open(ref.SKELETON, "rb") as fh:
        G = pickle.load(fh)
    print(f"                   nodes={G.number_of_nodes():,} edges={G.number_of_edges():,}")

    pairs = pd.read_parquet(ref.ATOM_PAIRS)
    pairs["mnxr"] = pairs["mnxr"].astype(str)
    print(f"[build-atom-graph] frozen atom_pairs rows={len(pairs):,} "
          f"mnxr={pairs['mnxr'].nunique():,}")

    ref.GRAPH_DIR.mkdir(parents=True, exist_ok=True)
    for X in ref.BASE_GRAPH_ELEMENTS:
        w = per_element_weights(pairs, X)
        H = G.copy()
        wk = f"w_{X}"
        nx.set_edge_attributes(H, 0, wk)
        n_set, n_off_skeleton = 0, 0
        for (mnxr, mnxm), val in w.items():
            u, v = ("rxn", mnxr), ("met", mnxm)
            if H.has_edge(u, v):
                H.edges[u, v][wk] = float(val)
                n_set += 1
            else:
                n_off_skeleton += 1
        n_pos = sum(1 for _, _, d in H.edges(data=True) if d.get(wk, 0) > 0)
        fp = frozen_positive(X)
        delta = f"  frozen={fp:,}  net {n_pos - fp:+,}" if fp is not None else ""
        out = ref.GRAPH_DIR / f"mnx_bipartite_{X}.pkl"
        with open(out, "wb") as fh:
            pickle.dump(H, fh, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"[build-atom-graph] {X}: w_{X}>0 edges = {n_pos:,} / {H.number_of_edges():,}"
              f"  (set {n_set:,}; off-skeleton dropped {n_off_skeleton:,}){delta}  -> {out.name}")

    print(f"[build-atom-graph] wrote {len(ref.BASE_GRAPH_ELEMENTS)} graphs -> {ref.GRAPH_DIR}")


if __name__ == "__main__":
    main()
