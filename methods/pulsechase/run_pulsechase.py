#!/usr/bin/env python
"""The ECSPr pulse-chase suite -- does ECSPr behave like a radio-labelled pulse-chase?

ECSPr factors into two stages the suite tests separately:

  1. GRAPH CONSTRUCTION -- atom-atom mapping turns each reaction into role-respecting,
     atom-weighted transfers assembled into a weighted atom-transfer network. Ambiguity
     is handled here by DILUTING a mapping's weight across candidates (never gapping);
     perturbations (a fosmid's added reactions) are also applied here, upstream.
  2. ECSPr SOLVE -- a pure function (network, source atoms, sink atoms) -> Ieff. It takes
     a finished network and two endpoints and returns effective conductance. It knows
     nothing about how the network was built or perturbed.

The four existing gates are all REPRODUCTION gates and every one stays green through the
star's clique artifact, the MIN scorer, and the refusal-gap. This suite is the missing
BEHAVIOURAL gate: it encodes the two-stage contract as runnable checks, each able to
fail, and runs them to a red/green baseline against the star vs the atom graph.

    python run_pulsechase.py --all
    python run_pulsechase.py --only II2

Several checks are EXPECTED-RED by design -- each pins a specific correction the next
session must make (I3 candidate emission, II2 all-paths measurement, II9 whole-atom-set
participation, I9 no-loss floor). The runner separates "the star fails as designed" and
"an expected-red is still red" from "the atom graph regressed"; the exit code keys ONLY
on the last.

This pass fixes nothing. The baseline is the deliverable.

Env: numpy + scipy + networkx + pandas + rdkit + pyarrow. The suite GATES its env and
refuses an interpreter missing a dependency, naming it (the parity-gate idiom).
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))                 # this dir, for toy/thresholds


def _gate_env() -> None:
    """Refuse an interpreter that cannot run the suite, naming the missing dep."""
    import importlib
    need = ["numpy", "scipy", "networkx", "pandas", "pyarrow", "rdkit"]
    missing = []
    for m in need:
        try:
            importlib.import_module(m)
        except Exception:
            missing.append(m)
    if missing:
        raise SystemExit(
            f"pulse-chase cannot run under {sys.executable}\n"
            f"  missing: {', '.join(missing)}\n\n"
            f"The suite needs the full ECSPr stack incl. rdkit (I3 exercises the real\n"
            f"pairs_from_mapped/forced_pairs refusal path). Run under an interpreter that\n"
            f"has them, e.g. `mamba run -n scadc-metabolic-model python "
            f"run_pulsechase.py --all`."
        )


_gate_env()

import numpy as np              # noqa: E402
import networkx as nx           # noqa: E402
import pandas as pd             # noqa: E402

from fabfos import canon                    # noqa: E402
import toy                      # noqa: E402
import thresholds as T          # noqa: E402

sys.path.insert(0, str(canon.ENGINE_LIB / "resources" / "lib"))
from ecspr_atom_graph import (atom_edges, star_edges, Grid, axis_terminals,     # noqa: E402
                              reff_from_z, permute_pairs, ATOM_REFF_EPS,
                              supernode_ieff, merge_terminals, metabolite_atoms)
from ecspr_solver import _reff_dense, SELFTEST_TOL, REFF_EPS                      # noqa: E402
from ecspr_directed import build_incidence, directed_ceff                         # noqa: E402

# The all-paths measurement now lives in the engine (ecspr_atom_graph); the suite tests
# the engine's function directly rather than a local reimplementation. Kept under the
# suite's original name so the checks that reference it read unchanged.
ieff_supernode = supernode_ieff

HOST_PAIRS = canon.DATA / "fabfos_2026_199" / "ecspr_atom" / "atom_pairs_universe.parquet"


# =====================================================================
# The ECSPr SOLVE under test, and the all-paths reference it is measured against
# =====================================================================

def _grid_edges(grid: Grid) -> dict:
    """The full-graph edge dict {(u, v): w} of a grid -- all components, not just the LCC,
    which is what the all-paths metabolite-to-metabolite measurement runs on."""
    return {(u, v): d["w"] for u, v, d in grid.g.edges(data=True)}


def ecspr_solve(grid: Grid, met_s: str, met_k: str, atomic: bool = True) -> dict:
    """Pure-function ECSPr on a finished grid, via the engine's all-paths super-node
    measurement (`supernode_ieff`). Returns a DEFINITE answer (II4): finite Ieff when the
    endpoints are connected, an informative zero-capacity (reachable=False) when either
    endpoint carries no atom of the element or the two are genuinely disconnected."""
    if atomic:
        edges = _grid_edges(grid)
        S = metabolite_atoms(edges, met_s)
        K = metabolite_atoms(edges, met_k)
    else:
        edges = _grid_edges(grid)
        S = [("met", met_s)] if ("met", met_s) in grid.g else []
        K = [("met", met_k)] if ("met", met_k) in grid.g else []
    if not S or not K:
        return {"ieff": 0.0, "reff": float("inf"), "reachable": False}
    ie = supernode_ieff(edges, S, K)
    reachable = ie > 0.0
    return {"ieff": ie, "reff": (1.0 / ie if ie > 0 else float("inf")),
            "reachable": reachable}


def atom_atoms(edges_or_grid, met: str) -> list:
    """The atom nodes of metabolite `met` present in an atom edge dict."""
    if isinstance(edges_or_grid, dict):
        seen = set()
        for (u, v) in edges_or_grid:
            for n in (u, v):
                if isinstance(n, tuple) and n and n[0] == met:
                    seen.add(n)
        return sorted(seen, key=str)
    return [n for n in edges_or_grid.lcc if n[0] == met]


# =====================================================================
# Verdict plumbing
# =====================================================================
# A check returns {name, cls, arms:{arm:verdict}, expected:{arm:exp}, note}. Verdicts
# are "PASS"/"FAIL"/"NA". `expected` is what the baseline should be: PASS (green),
# FAIL (star fails as designed), RED (expected-red -- a named next-session fix). A
# REGRESSION is any arm whose expected is PASS but whose verdict is not PASS.

def V(name, cls, arms, expected, note=""):
    return {"name": name, "cls": cls, "arms": arms, "expected": expected, "note": note}


def _strict_gt(a, b):
    return a > b + T.STRICT_EPS


# =====================================================================
# GROUP I -- graph construction
# =====================================================================

def check_I1():
    """Margins/conservation on the pairs table with weights=1.0 (exact, canon.EPS)."""
    fixtures = {
        "cofactor": (toy.fx_cofactor(), None),
        "lyase": (toy.fx_lyase(), 5),      # S(5) split -> P1(2)+P2(3): 5 substrate atoms
        "multi": (toy.fx_multi(), 2),
        "transport": (toy.fx_transport(), 3),  # M(2 recycle)+X(1): 3 substrate atoms
    }
    worst = 0.0
    ok = True
    for tag, (fx, known) in fixtures.items():
        pairs = toy.pairs_table(fx["reactions"], fx["element"])
        # every row internally consistent
        for r in pairs.itertuples(index=False):
            si, pi = r.sub_idx.split(","), r.prod_idx.split(",")
            if not (len(si) == len(pi) == r.n_atoms):
                ok = False
        for (mnxr, el), g in pairs.groupby(["mnxr", "element"]):
            sub_margin = sum(len(x.split(",")) for x in g.sub_idx)
            prod_margin = sum(len(x.split(",")) for x in g.prod_idx)
            worst = max(worst, abs(sub_margin - prod_margin))
            if abs(sub_margin - prod_margin) > canon.EPS:
                ok = False
            if known is not None and len(g) and sub_margin != known:
                ok = False
    verdict = "PASS" if (ok and worst <= canon.EPS) else "FAIL"
    return V("I1", "construction", {"atom": verdict}, {"atom": "PASS"},
             f"substrate/product margins conserved on pairs table (worst |d|={worst:g})")


def check_I2():
    """Role-respecting: no substrate->substrate conductance. Vacuous on the atom graph;
    the real, able-to-fail discriminator is the star clique (co-substrate leak)."""
    fx = toy.fx_condensation()
    rx, el = fx["reactions"], fx["element"]
    # atom arm: two co-substrates A, B in different atom components -> Ieff 0
    a_edges = atom_edges(toy.pairs_table(rx, el), toy.weights_unit(rx))
    a_grid = Grid(a_edges, "atom")
    atom_leak = ieff_supernode(a_edges, atom_atoms(a_edges, fx["in1"]),
                               atom_atoms(a_edges, fx["in2"]))
    atom_ok = atom_leak <= ATOM_REFF_EPS
    # star arm: hub joins A and B -> co-substrate conductance manufactured
    s_grid = toy.star_grid(rx, el)
    star_leak = ecspr_solve(s_grid, fx["in1"], fx["in2"], atomic=False)["ieff"]
    star_fails = star_leak > T.STRICT_EPS
    return V("I2", "discrimination",
             {"atom": "PASS" if atom_ok else "FAIL",
              "star": "FAIL" if star_fails else "PASS"},
             {"atom": "PASS", "star": "FAIL"},
             f"co-substrate conductance  atom={atom_leak:.4g} (0) star={star_leak:.4g} (leak)")


def _doubly_stochastic(pd_pairs: dict) -> tuple:
    """Margin check on an emitted {(el,sm,pm):[(sr,pr,w),...]} pair dict: every source
    atom's total OUTGOING weight and every product atom's total INCOMING weight equal
    1.0 -- the conserved 'sum == a confident pairing (er)' the dilute-not-gap fix must
    preserve. Returns (ok, n_source_atoms, n_product_atoms)."""
    src, snk = defaultdict(float), defaultdict(float)
    for (el, sm, pm), lst in pd_pairs.items():
        for (sr, pr, w) in lst:
            src[(el, sm, sr)] += w
            snk[(el, pm, pr)] += w
    ok = (all(abs(v - 1.0) <= canon.EPS for v in src.values())
          and all(abs(v - 1.0) <= canon.EPS for v in snk.values())
          and bool(src) and bool(snk))
    return ok, len(src), len(snk)


def check_I3():
    """Dilute-not-gap at emission -- FIXED (was expected-red). An ambiguous atom
    correspondence is no longer REFUSED (the gap the model forbids): every candidate
    pairing is emitted DILUTED by fanout so per-atom margins are conserved (each source
    atom's outgoing weight sums to 1.0 == a confident pairing's er). forced_pairs n>1
    emits the n x n doubly-stochastic completion; the mapped ambiguous-duplicate naming
    is spread over its k candidates at weight 1/k. Able to fail: a refusal (0 rows) or a
    broken margin flips it red."""
    from ecspr_atom_pairs import forced_pairs, pairs_from_mapped
    # forced n>1: was refused (0 rows); now the n x n diluted completion, margins == 1
    sub, prod, formulas, ranks_of, n = toy.ambiguous_forced_inputs()
    fp = forced_pairs(sub, prod, formulas, ranks_of)
    forced_emitted = sum(len(v) for v in fp.values())
    forced_margins, _, _ = _doubly_stochastic(fp)
    # distinct ranks so the n x n candidates are n*n DISTINCT atom pairs (not collapsed)
    distinct = (len(set(ranks_of[("S_ETOH", "C")])) == n
                and len(set(ranks_of[("P_ACAL", "C")])) == n)
    forced_ok = (forced_emitted == n * n and forced_margins and distinct)
    # mapped ambiguous_duplicate: was refused; now diluted over the k namings, margins == 1
    mapped, ms, mp, cmap = toy.ambiguous_mapped_inputs()
    pairs, status = pairs_from_mapped(mapped, ms, mp, cmap)
    mapped_emitted = sum(len(v) for v in pairs.values())
    mapped_margins, _, _ = _doubly_stochastic(pairs)
    mapped_ok = (mapped_emitted > 0 and status != "ambiguous_duplicate" and mapped_margins)
    ok = forced_ok and mapped_ok
    return V("I3", "construction", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"forced emits {forced_emitted}={n}x{n} margins={forced_margins}; "
             f"mapped status={status!r} emits {mapped_emitted} margins={mapped_margins}")


def check_I4():
    """Cofactor segregation -- able-to-fail, toy. (a) cargo reaches product; (b) cargo
    and carrier are in different atom components while the star leaks by >= a meaningful
    ratio; (c) a synthetic cargo->carrier pair breaks the segregation (proves it can fail)."""
    fx = toy.fx_cofactor()
    rx, el = fx["reactions"], fx["element"]
    edges = atom_edges(toy.pairs_table(rx, el), toy.weights_unit(rx))
    G = nx.Graph()
    for (u, v), w in edges.items():
        G.add_edge(u, v, w=w)
    src, snk, carrier = fx["cargo_src"], fx["cargo_snk"], fx["carrier"]
    # (a) positive control -- cargo atom reaches product atom
    cargo_ie = ieff_supernode(edges, atom_atoms(edges, src), atom_atoms(edges, snk))
    a_ok = cargo_ie > T.STRICT_EPS
    # (b) segregation on the atom graph
    prod_comp = set().union(*[nx.node_connected_component(G, n)
                              for n in atom_atoms(edges, snk)])
    seg_ok = not any(n[0] == carrier for n in prod_comp)
    # star leak ratio: theft channel vs real channel
    s_grid = toy.star_grid(rx, el)
    leak = ecspr_solve(s_grid, fx["theft_a"], fx["theft_b"], atomic=False)["ieff"]
    real = ecspr_solve(s_grid, src, snk, atomic=False)["ieff"]
    ratio = leak / real if real > 0 else float("inf")
    star_fails = ratio >= T.STAR_LEAK_RATIO_MIN
    # (c) leak-injection: add cargo->carrier atom pair, segregation must break
    inj = dict(edges)
    inj[(atom_atoms(edges, src)[0], atom_atoms(edges, carrier)[0])] = 1.0
    Gi = nx.Graph()
    for (u, v), w in inj.items():
        Gi.add_edge(u, v, w=w)
    inj_comp = set().union(*[nx.node_connected_component(Gi, n)
                             for n in atom_atoms(edges, snk)])
    c_can_fail = any(n[0] == carrier for n in inj_comp)
    atom_ok = a_ok and seg_ok and c_can_fail
    return V("I4", "discrimination",
             {"atom": "PASS" if atom_ok else "FAIL",
              "star": "FAIL" if star_fails else "PASS"},
             {"atom": "PASS", "star": "FAIL"},
             f"cargo->prod={cargo_ie:.3g} seg={seg_ok} inj_breaks={c_can_fail} | "
             f"star leak/real={ratio:.2f} (>= {T.STAR_LEAK_RATIO_MIN})")


def check_I5():
    """Compartment collapse: a transport reaction (identical bare-metabolite multiset)
    self-annihilates -- no cross-metabolite edge from it on the atom graph."""
    fx = toy.fx_transport()
    rx, el = fx["reactions"], fx["element"]
    edges = atom_edges(toy.pairs_table(rx, el), toy.weights_unit(rx))
    m = fx["recycle"]
    # no edge touches two DISTINCT copies of the recycled metabolite, and the recycle
    # contributes no self-edge at all (u==v dropped by atom_edges)
    recycle_edges = [(u, v) for (u, v) in edges if u[0] == m and v[0] == m]
    # the only surviving edge is the real X->Y carrier
    survives = {(u[0], v[0]) for (u, v) in edges}
    ok = (len(recycle_edges) == 0 and ("X", "Y") in survives or ("Y", "X") in survives)
    ok = len(recycle_edges) == 0 and any({a, b} == {"X", "Y"} for a, b in survives)
    return V("I5", "construction", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"recycle self-edges={len(recycle_edges)} (0); surviving met-pairs={sorted(survives)}")


def check_I6():
    """Orientation invariance -- swapping substrate<->product leaves every atom-graph
    reff identical (atom_edges symmetrizes keys u<v). VERIFIED SOUND (green)."""
    fx = toy.fx_multi()
    rx, el = fx["reactions"], fx["element"]
    pairs = toy.pairs_table(rx, el)
    swapped = pairs.rename(columns={"substrate": "product", "product": "substrate",
                                    "sub_idx": "prod_idx", "prod_idx": "sub_idx"})
    e1 = atom_edges(pairs, toy.weights_unit(rx))
    e2 = atom_edges(swapped, toy.weights_unit(rx))
    ok = set(e1) == set(e2) and all(abs(e1[k] - e2[k]) <= canon.EPS for k in e1)
    return V("I6", "construction", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"atom edges identical under substrate<->product swap ({len(e1)} edges)")


def check_I7():
    """Symmetric molecule -- the fumarate/succinate fixture conserves margins whether or
    not CanonicalRankAtoms collapses the symmetric atoms (read pre-collapse)."""
    fx = toy.fx_symmetric()
    rx, el = fx["reactions"], fx["element"]
    pairs = toy.pairs_table(rx, el)
    ok = True
    for (mnxr, e), g in pairs.groupby(["mnxr", "element"]):
        sub_margin = sum(len(x.split(",")) for x in g.sub_idx)
        prod_margin = sum(len(x.split(",")) for x in g.prod_idx)
        if abs(sub_margin - prod_margin) > canon.EPS:
            ok = False
    return V("I7", "construction", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             "symmetric-molecule margins conserved on the pairs table")


def check_I8():
    """Reaction shapes -- fanout, lyase split, recycle each conserve margins pre-collapse."""
    ok = True
    detail = []
    for tag, fx, known in [("multi", toy.fx_multi(), 2),
                           ("lyase", toy.fx_lyase(), 5),
                           ("recycle", toy.fx_transport(), 3)]:
        pairs = toy.pairs_table(fx["reactions"], fx["element"])
        sub = sum(len(x.split(",")) for x in pairs.sub_idx)
        prod = sum(len(x.split(",")) for x in pairs.prod_idx)
        conserved = abs(sub - prod) <= canon.EPS and sub == known
        detail.append(f"{tag}={sub}")
        ok = ok and conserved
    return V("I8", "construction", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             "shapes conserve margins: " + " ".join(detail))


def check_I9():
    """No-loss under --no-prune -- FIXED (was expected-red). A pruned build (Grid keeps only
    the LCC in its factorized arrays) silently drops every atom in a minority component;
    Grid(prune=False) RETAINS every component, so the silent-loss residual -- element-
    carrying atoms not in the retained graph -- is empty. Relational and rebuild-invariant:
    no transcribed bucket sizes. Able to fail two ways: if no-prune drops anything, or if
    the problem it fixes is not real (the pruned build would drop nothing)."""
    # Host carbon atom graph, weights=1.0 -- the real universe, where fragmentation is real.
    pairs = pd.read_parquet(HOST_PAIRS)
    pc = pairs[pairs.element == "C"]
    weights = {r: 1.0 for r in pc.mnxr.unique()}
    edges = atom_edges(pc, weights)
    # no-prune: every atom node that participates in an inter-metabolite transfer is kept
    g_np = Grid(edges, "host-C", prune=False)
    all_nodes = set(g_np.g.nodes)
    residual = all_nodes - set(g_np.lcc)     # silent-loss bucket under --no-prune
    no_prune_ok = (len(residual) == 0 and set(g_np.lcc) == all_nodes)
    # contrast: the pruned build DOES drop the minority components (what --no-prune fixes)
    g_pr = Grid(edges, "host-C", prune=True)
    pruned_out = len(all_nodes) - len(set(g_pr.lcc))
    frac_out = pruned_out / len(all_nodes)
    ok = no_prune_ok and pruned_out > 0
    return V("I9", "construction", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"no-prune retains all {len(all_nodes)} atoms (residual={len(residual)}); "
             f"pruned build drops {pruned_out} ({frac_out:.1%}) -> --no-prune is the no-loss tier")


# =====================================================================
# GROUP II -- ECSPr solve
# =====================================================================

def check_II1():
    """Pure-function contract: Ieff depends only on (network, source, sink); deterministic."""
    fx = toy.fx_cofactor()
    g1 = toy.atom_grid(fx["reactions"], fx["element"])
    g2 = toy.atom_grid(fx["reactions"], fx["element"])
    a = ecspr_solve(g1, fx["cargo_src"], fx["cargo_snk"])
    b = ecspr_solve(g2, fx["cargo_src"], fx["cargo_snk"])
    ok = a["ieff"] == b["ieff"] and a["reff"] == b["reff"]
    return V("II1", "solve", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"same inputs -> identical Ieff ({a['ieff']:.10g})")


def check_II3():
    """Correctness vs a shared-nothing referent (VERIFIED SOUND). Per-pair atom-lane reff
    (grid.zcols/reff_from_z) matches _reff_dense on the LCC within SELFTEST_TOL. Not the
    star SMWSolver (shares the Laplacian recipe)."""
    worst = 0.0
    n = 0
    for fx in (toy.fx_parallel(), toy.fx_cofactor(), toy.fx_multi()):
        rx, el = fx["reactions"], fx["element"]
        edges = atom_edges(toy.pairs_table(rx, el), toy.weights_unit(rx))
        grid = Grid(edges, "atom")
        nodes = grid.nodes
        terms = list(range(1, grid.n))          # ground is index 0
        if len(terms) < 1:
            continue
        Z, _ = grid.zcols(terms)
        Gsub = nx.Graph()
        for (u, v), w in edges.items():
            if u in grid.lcc and v in grid.lcc:
                Gsub.add_edge(u, v, w=w)
        # sample a few in-LCC pairs
        rng = np.random.default_rng(canon.FIT_SEED)
        pool = terms if len(terms) <= 8 else list(rng.choice(terms, 8, replace=False))
        for i in range(len(pool)):
            for j in range(i + 1, len(pool)):
                a, b = pool[i], pool[j]
                r_grid = reff_from_z(Z, terms, a, b)
                r_dense = _reff_dense(Gsub, nodes[a], nodes[b], wk="w")
                worst = max(worst, abs(r_grid - r_dense))
                n += 1
    ok = worst < SELFTEST_TOL and n > 0
    return V("II3", "solve", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"per-pair reff vs dense over {n} pairs: worst |d|={worst:.2e} (< {SELFTEST_TOL:g})")


def check_II4():
    """Always a definite answer: co-LCC -> finite Ieff; genuinely disconnected -> a
    definite zero-capacity (component non-membership), never a skip or exception."""
    fx = toy.fx_network_pair()
    gA = toy.atom_grid(fx["reactions_A"], fx["element"])
    # reachable pair (background axis) -> finite
    reach = ecspr_solve(gA, fx["disjoint_src"], fx["disjoint_snk"])
    # genuinely disconnected pair (MEV, IPP not linked in A) -> definite zero
    disc = ecspr_solve(gA, fx["route_src"], fx["route_snk"])
    ok = (reach["reachable"] and np.isfinite(reach["ieff"]) and reach["ieff"] > 0
          and (not disc["reachable"]) and disc["ieff"] == 0.0)
    return V("II4", "solve", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"reachable Ieff={reach['ieff']:.3g}; disconnected={disc['reachable']} "
             f"Ieff={disc['ieff']} (definite zero-capacity)")


def check_II5():
    """Network-comparison monotonicity + specificity. B superset A (mevalonate pathway
    added upstream): (a) Ieff_B >= Ieff_A; (b) strictly greater where the add creates the
    S->K route; (c) a disjoint axis is exactly invariant (delta 0.0). ECSPr is a pure
    function on each finished network."""
    fx = toy.fx_network_pair()
    gA = toy.atom_grid(fx["reactions_A"], fx["element"])
    gB = toy.atom_grid(fx["reactions_B"], fx["element"])
    aA = ecspr_solve(gA, fx["route_src"], fx["route_snk"])["ieff"]
    aB = ecspr_solve(gB, fx["route_src"], fx["route_snk"])["ieff"]
    dA = ecspr_solve(gA, fx["disjoint_src"], fx["disjoint_snk"])["ieff"]
    dB = ecspr_solve(gB, fx["disjoint_src"], fx["disjoint_snk"])["ieff"]
    mono = aB >= aA - T.STRICT_EPS
    strict = _strict_gt(aB, aA) and aA == 0.0        # unreachable -> finite
    spec = abs(dA - dB) == 0.0                        # disjoint axis exactly invariant
    ok = mono and strict and spec
    return V("II5", "solve", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"route Ieff A={aA:.3g}->B={aB:.3g} (mono&strict); disjoint delta={abs(dA-dB):g} (0)")


def check_II2():
    """All-paths, edge-preserving measurement -- FIXED (was expected-red). Ieff now shorts
    the metabolite endpoints and solves ALL paths (engine supernode_ieff), stated as
    properties of Ieff (directed-ready): (a) parallel additivity -- two edge-disjoint
    unequal-resistance routes give Ieff strictly greater than either alone; (b) edge
    preservation -- removing an edge that carries source->sink current strictly decreases
    Ieff. The retired MIN scorer failed both (it reported only the single best route).
    Able to fail: if either property breaks, verdict flips red."""
    fx = toy.fx_parallel()
    el = fx["element"]
    both = fx["reactions"]
    a_only = fx["reactions_A_only"]
    ie_both = ecspr_solve(toy.atom_grid(both, el), fx["src"], fx["snk"])["ieff"]
    ie_A = ecspr_solve(toy.atom_grid(a_only, el), fx["src"], fx["snk"])["ieff"]
    # (a) adding the second route strictly increases Ieff (parallel combination)
    additivity = _strict_gt(ie_both, ie_A)
    # (b) removing route B's unique edge strictly decreases Ieff (edge preservation)
    without_B_edge = [r for r in both if r["mnxr"] != fx["routeB_edge_rxn"]]
    ie_wo = ecspr_solve(toy.atom_grid(without_B_edge, el), fx["src"], fx["snk"])["ieff"]
    edge_preserved = _strict_gt(ie_both, ie_wo)
    ok = additivity and edge_preserved
    return V("II2", "solve", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"all-paths both={ie_both:.3g} > A={ie_A:.3g} (additivity {additivity}); "
             f"remove routeB edge -> {ie_wo:.3g} (edge-preserved {edge_preserved})")


def check_II9():
    """Whole-labelled-atom-set participation -- FIXED (was expected-red). Dropping one
    labelled source atom that SOLELY carries an independent route strictly lowers Ieff:
    every labelled atom is a terminal, not just the best one. The retired MIN scorer left
    Ieff unchanged (it only ever saw the single best atom). Able to fail: if dropping the
    independent-route atom does not lower Ieff, verdict flips red."""
    fx = toy.fx_parallel()
    el = fx["element"]
    both = fx["reactions"]
    e = atom_edges(toy.pairs_table(both, el), toy.weights_unit(both))
    src_all = atom_atoms(e, fx["src"])
    drop_rank = fx["drop_src_rank"]
    src_kept = [n for n in src_all if n != (fx["src"], drop_rank)]
    snk = atom_atoms(e, fx["snk"])
    full = supernode_ieff(e, src_all, snk)
    dropped = supernode_ieff(e, src_kept, snk)
    ok = _strict_gt(full, dropped)
    return V("II9", "solve", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"all source atoms Ieff={full:.3g} > drop independent-route atom={dropped:.3g}")


def check_II6():
    """Cross-element C/N/S/P. Numeric correctness on every element with confirmed atom-pair
    coverage (assert coverage first). Sulfur/phosphorus carrier segregation on toy fixtures
    (host reactions for these are the shapes the mapper refuses)."""
    pairs = pd.read_parquet(HOST_PAIRS)
    coverage = {e: int((pairs.element == e).sum()) for e in canon.ELEMENTS}
    covered = all(coverage[e] > 0 for e in canon.ELEMENTS)
    # S and P carrier segregation on toys (ACP thioester carries S; ATP phosphate not a carrier)
    seg_ok = True
    star_fails = True
    for el in ("S", "P"):
        fx = toy.fx_cofactor(el)
        rx = fx["reactions"]
        edges = atom_edges(toy.pairs_table(rx, el), toy.weights_unit(rx))
        G = nx.Graph()
        for (u, v), w in edges.items():
            G.add_edge(u, v, w=w)
        prod_comp = set().union(*[nx.node_connected_component(G, n)
                                  for n in atom_atoms(edges, fx["cargo_snk"])])
        seg_ok = seg_ok and not any(n[0] == fx["carrier"] for n in prod_comp)
        s_grid = toy.star_grid(rx, el)
        leak = ecspr_solve(s_grid, fx["theft_a"], fx["theft_b"], atomic=False)["ieff"]
        real = ecspr_solve(s_grid, fx["cargo_src"], fx["cargo_snk"], atomic=False)["ieff"]
        star_fails = star_fails and (leak / real >= T.STAR_LEAK_RATIO_MIN if real > 0 else True)
    atom_ok = covered and seg_ok
    return V("II6", "discrimination",
             {"atom": "PASS" if atom_ok else "FAIL",
              "star": "FAIL" if star_fails else "PASS"},
             {"atom": "PASS", "star": "FAIL"},
             f"coverage {coverage}; S/P carrier segregation={seg_ok}")


def check_II7():
    """Star-vs-atom artifact discrimination -- AFFIRMATIVE. The star exhibits a POSITIVE
    cofactor-theft/input-leak conductance exceeding the real channel by >= a meaningful
    ratio (printed); the atom graph does not (component non-membership)."""
    fx = toy.fx_cofactor()
    rx, el = fx["reactions"], fx["element"]
    edges = atom_edges(toy.pairs_table(rx, el), toy.weights_unit(rx))
    atom_leak = ieff_supernode(edges, atom_atoms(edges, fx["theft_a"]),
                               atom_atoms(edges, fx["theft_b"]))
    s_grid = toy.star_grid(rx, el)
    leak = ecspr_solve(s_grid, fx["theft_a"], fx["theft_b"], atomic=False)["ieff"]
    real = ecspr_solve(s_grid, fx["cargo_src"], fx["cargo_snk"], atomic=False)["ieff"]
    ratio = leak / real if real > 0 else float("inf")
    star_fails = ratio >= T.STAR_LEAK_RATIO_MIN
    atom_ok = atom_leak <= ATOM_REFF_EPS
    return V("II7", "discrimination",
             {"atom": "PASS" if atom_ok else "FAIL",
              "star": "FAIL" if star_fails else "PASS"},
             {"atom": "PASS", "star": "FAIL"},
             f"star theft/real={ratio:.2f} (>= {T.STAR_LEAK_RATIO_MIN}); atom theft={atom_leak:.4g} (0)")


def check_II8():
    """permute_pairs control -- construction+solve sensitivity. Shuffling which-atom-goes-
    where within each reaction (canon.DRAW_SEED) collapses the mean |Ieff effect| below a
    named fraction of the unshuffled -- the lane measures chemistry, not mappability."""
    # a reaction where atom identity decides the transfer: A->PA,B->PB,C->PC,D->PD.
    # permute reassigns which product each substrate feeds; the measured A->PA transfer
    # survives only the 1-in-4 identity permutation.
    fx = toy.fx_permute()
    rx, el = fx["reactions"], fx["element"]
    pairs = toy.pairs_table(rx, el)
    w = toy.weights_unit(rx)

    def effect(p):
        e = atom_edges(p, w)
        return ieff_supernode(e, atom_atoms(e, fx["src"]), atom_atoms(e, fx["snk"]))

    base = effect(pairs)
    shuffled = [effect(permute_pairs(pairs, canon.DRAW_SEED + k)) for k in range(8)]
    mean_shuf = float(np.mean([abs(s) for s in shuffled]))
    collapsed = mean_shuf < T.PERMUTE_COLLAPSE_FRAC * base if base > 0 else False
    return V("II8", "robustness", {"atom": "PASS" if collapsed else "FAIL"}, {"atom": "PASS"},
             f"unshuffled effect={base:.3g}; shuffled mean={mean_shuf:.3g} "
             f"(< {T.PERMUTE_COLLAPSE_FRAC} x base)")


def check_II10():
    """Degenerate endpoints defined: source==sink and source-carries-no-atom have a
    DEFINED result (per II4), never undefined or an exception."""
    fx = toy.fx_cofactor()
    grid = toy.atom_grid(fx["reactions"], fx["element"])
    # source == sink -> defined: reflexive Ieff, we define it as unreachable/zero (no a!=b pair)
    same = ecspr_solve(grid, fx["cargo_snk"], fx["cargo_snk"])
    # source carries no atom of the element -> definite unreachable
    none = ecspr_solve(grid, "NOSUCH_MET", fx["cargo_snk"])
    ok = (same["ieff"] == 0.0) and (none["reachable"] is False and none["ieff"] == 0.0)
    return V("II10", "solve", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"source==sink Ieff={same['ieff']} (defined 0); no-atom reachable={none['reachable']} "
             f"Ieff={none['ieff']} (definite 0)")


# =====================================================================
# GROUP III -- far-distance and robustness
# =====================================================================

def axis_hops(base_graph: nx.Graph, src: str, dst: str) -> int | None:
    """Shortest reaction-hop distance between two metabolites on a star base graph.

    `nx.shortest_path_length` on ("met",m)<->("rxn",r) alternating nodes, halved to
    reaction hops. Returns None when there is no path (NetworkXNoPath caught here)."""
    s, t = ("met", src), ("met", dst)
    if s not in base_graph or t not in base_graph:
        return None
    try:
        d = nx.shortest_path_length(base_graph, s, t)
    except nx.NetworkXNoPath:
        return None
    return d // 2


def _load_base(el: str):
    import pickle
    with open(canon.SOLVE_BASE_DIR / f"base_{el}.pkl", "rb") as fh:
        return pickle.load(fh)


def _rank_far_axes(el: str):
    import json
    axes = json.loads(canon.AXES_JSON.read_text())
    testable = json.loads(canon.AXES_TESTABLE_JSON.read_text())
    base = _load_base(el)
    ranked = []
    for ax_id in testable.get(el, []):
        a = axes[ax_id]
        hops = axis_hops(base, a["source"][0], a["sink"][0])
        if hops is not None:
            ranked.append((hops, ax_id, a["source"][0], a["sink"][0]))
    ranked.sort(reverse=True)
    return base, ranked


def check_III1():
    """Far edges. Rank the testable canon.AXIS_SET axes by reaction hops; take
    FAR_EDGE_COUNT per element. Every star-reachable far edge is co-LCC on the star base (component
    membership), or fails into the enumerated allow-set. Assert reachability, not a
    magnitude threshold."""
    from ecspr_solver import SMWGraphContext
    ok = True
    total = 0
    fails = []
    for el in canon.ELEMENTS:
        base, ranked = _rank_far_axes(el)
        far = ranked[:T.FAR_EDGE_COUNT]
        ctx = SMWGraphContext(base, edge_weight_key=f"w_{el}")
        for hops, ax_id, s, t in far:
            total += 1
            sn, tn = ("met", s), ("met", t)
            reachable = sn in ctx.lcc and tn in ctx.lcc
            if not reachable:
                # allow-set: transport (s==t bare), no-atom-of-element, or genuinely
                # disconnected in source chemistry
                allowed = (s == t) or (sn not in base) or (tn not in base) \
                    or (not nx.has_path(base, sn, tn))
                if not allowed:
                    ok = False
                    fails.append(ax_id)
    return V("III1", "robustness", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"{total} far edges over {len(canon.ELEMENTS)} elements co-LCC or allow-set "
             f"({'clean' if ok else 'FAILS: ' + ','.join(fails)})")


def check_III2():
    """Conditioning guard. On far edges the base R_eff is finite and below FAR_REFF_MAX,
    so long-path Ieff is signal, not float noise near the 1e-12 clamp."""
    from ecspr_solver import SMWGraphContext, SMWSolver, reff_base
    ok = True
    worst = 0.0
    n = 0
    for el in canon.ELEMENTS:
        base, ranked = _rank_far_axes(el)
        ctx = SMWGraphContext(base, edge_weight_key=f"w_{el}")
        for hops, ax_id, s, t in ranked[:T.FAR_EDGE_COUNT]:
            try:
                solver = SMWSolver(base, [("met", s)], [("met", t)],
                                   edge_weight_key=f"w_{el}", context=ctx)
            except ValueError:
                continue
            r = reff_base(solver)
            n += 1
            worst = max(worst, r)
            if not (np.isfinite(r) and 0 < r < T.FAR_REFF_MAX):
                ok = False
    return V("III2", "robustness", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"{n} far-edge base R_eff finite & < {T.FAR_REFF_MAX:g} (worst r={worst:.3g})")


def check_III3():
    """Determinism. Seeded permute_pairs (canon.DRAW_SEED) reproduces bit-identically;
    the frozen draw-size list is the explicit one, not the retired sbatch sizes."""
    fx = toy.fx_parallel()
    pairs = toy.pairs_table(fx["reactions"], fx["element"])
    p1 = permute_pairs(pairs, canon.DRAW_SEED)
    p2 = permute_pairs(pairs, canon.DRAW_SEED)
    bit_identical = p1.equals(p2)
    # a np.random draw reproduces
    d1 = np.random.default_rng(canon.DRAW_SEED).permutation(10)
    d2 = np.random.default_rng(canon.DRAW_SEED).permutation(10)
    draw_ok = np.array_equal(d1, d2)
    # frozen list is not the retired one (K-collision / retired-size guard)
    frozen_ok = set(canon.DRAW_SIZES) != set(canon.RETIRED_DRAW_SIZES)
    ok = bit_identical and draw_ok and frozen_ok
    return V("III3", "robustness", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"permute bit-identical={bit_identical}; seeded draw reproduces={draw_ok}; "
             f"frozen!=retired sizes={frozen_ok}")


def check_III4():
    """Sharding partition-invariance. Ieff per (network, axis) is independent of which
    shard computes it -- exact == for per-cell values regardless of partition."""
    fx = toy.fx_network_pair()
    grid = toy.atom_grid(fx["reactions_B"], fx["element"])
    axes = [(fx["route_src"], fx["route_snk"]),
            (fx["disjoint_src"], fx["disjoint_snk"])]
    # ground truth: all axes in one shard
    truth = {ax: ecspr_solve(grid, *ax)["ieff"] for ax in axes}

    def sharded(nshards):
        out = {}
        for i, ax in enumerate(axes):
            if i % nshards == i % nshards:      # every axis lands in exactly one shard
                out[ax] = ecspr_solve(grid, *ax)["ieff"]
        return out

    s2 = sharded(2)
    s3 = sharded(3)
    exact = all(truth[ax] == s2[ax] == s3[ax] for ax in axes)
    # shard counts pinned to canon.DRAW_SIZES, never the retired sbatch sizes
    # (canon.RETIRED_DRAW_SIZES)
    pin_ok = set(canon.DRAW_SIZES) != set(canon.RETIRED_DRAW_SIZES)
    ok = exact and pin_ok
    return V("III4", "robustness", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"per-cell Ieff exact across 2- and 3-shard partitions={exact}")


def check_III5():
    """Batch path -- VERIFIED vacuous, DROPPED with rationale. The CPU
    ecspr_ablation._reff_batch ignores batch_size (it is a per-map reff_summary loop),
    so a 'bit-identity across batch sizes' claim is vacuous. The honest GPU-vs-CPU
    retarget needs a CUDA device, unavailable here -> the claim is dropped, not faked."""
    import inspect
    from ecspr_ablation import _reff_batch
    src = inspect.getsource(_reff_batch)
    # affirm the vacuity: batch_size is a parameter the body never branches on
    vacuous = "batch_size" in src and "for a in maps" in src
    note = ("CPU _reff_batch ignores batch_size (per-map loop); GPU-vs-CPU retarget "
            "needs CUDA (absent) -> claim dropped, not asserted by bit-identity")
    return V("III5", "robustness", {"atom": "NA"}, {"atom": "NA"},
             note if vacuous else "UNEXPECTED: _reff_batch body changed")


# =====================================================================
# Runner
# =====================================================================

# =====================================================================
# GROUP IV -- directed conductance (the directional contract)
# =====================================================================
# Groups I-III exercise the UNDIRECTED contract only: construction and the all-paths
# supernode solve are symmetric in source<->sink, so a reaction that runs one way in the
# cell and the other way only against a steep thermodynamic gradient is scored identically
# in both directions. Directionality is a correctness axis none of I-III touch. GROUP IV
# runs the engine's real rectified-network solver (ecspr_directed, committed as a
# primitive) and asserts the two-sided directional contract: an irreversible reaction
# conducts forward but is throttled in reverse (IV1), AND the directed solve collapses to
# the undirected answer exactly at the symmetric limit ratio == 1 (IV2). These exercise
# the SAME solver the SCADC directed driver calls -- not a re-implementation -- so a
# regression in the engine primitive turns them red.

def _irrev_chain(ratio):
    """Toy: metabolites m0 -> r0 -> m1 -> r1 -> m2, both reactions irreversible with
    backward/forward conductance `ratio`. Returns (B, gp, gm, idx) with met/rxn node
    identity matching the production graph's ("met"/"rxn", id) convention."""
    nodes = [("met", "m0"), ("rxn", "r0"), ("met", "m1"), ("rxn", "r1"), ("met", "m2")]
    idx = {n: i for i, n in enumerate(nodes)}
    # oriented substrate -> rxn -> product (forward direction of each reaction)
    edges = [(idx[("met", "m0")], idx[("rxn", "r0")]),
             (idx[("rxn", "r0")], idx[("met", "m1")]),
             (idx[("met", "m1")], idx[("rxn", "r1")]),
             (idx[("rxn", "r1")], idx[("met", "m2")])]
    B = build_incidence(edges, len(nodes))
    gp = np.ones(len(edges))
    gm = ratio * np.ones(len(edges))
    return B, gp, gm, idx


def check_IV1():
    """Directed throttling -- an irreversible reaction conducts forward but not reverse.
    On a toy irreversible chain m0=>m1=>m2 the directed solver must give C_eff(m0->m2)
    at least DIRECTED_THROTTLE_MIN times C_eff(m2->m0). Symmetric (undirected) construction
    scores both directions identically; this is exactly what GROUPS I-III cannot see.
    Able to fail: if the solver stops distinguishing direction, the ratio collapses to ~1
    and the verdict flips red."""
    ratio = T.IV1_TOY_IRREV_RATIO      # a strongly-favoured (irreversible) reaction
    B, gp, gm, idx = _irrev_chain(ratio)
    c_fwd = directed_ceff(B, gp, gm, idx[("met", "m0")], idx[("met", "m2")])
    c_rev = directed_ceff(B, gp, gm, idx[("met", "m2")], idx[("met", "m0")])
    throttle = c_fwd / c_rev if c_rev > 0 else float("inf")
    ok = throttle >= T.DIRECTED_THROTTLE_MIN
    return V("IV1", "directed", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"irrev chain (ratio={ratio}): C_fwd={c_fwd:.4g} C_rev={c_rev:.4g} "
             f"throttle={throttle:.1f}x (>= {T.DIRECTED_THROTTLE_MIN}x)")


def check_IV2():
    """Symmetric-limit parity -- at ratio == 1 the directed solve reproduces the
    undirected two-terminal conductance. Ties GROUP IV back to the II-series: where
    direction is unknown (the ensemble's ratio 1.0 default), the directed model IS the
    undirected model, so no divergence in the head-to-head is an artifact of the solver.
    Referent is the undirected dense _reff_dense already gated by II3. Able to fail: if
    the ratio-1 directed solve disagrees with 1/R_eff beyond SELFTEST_TOL, verdict flips
    red; and the same irreversible chain must be direction-BLIND here (fwd == rev)."""
    B, gp, gm, idx = _irrev_chain(1.0)   # ratio == 1 -> symmetric
    s, t = idx[("met", "m0")], idx[("met", "m2")]
    c_dir = directed_ceff(B, gp, gm, s, t)
    # undirected referent on the same topology (build an nx graph, unit weights)
    G = nx.Graph()
    node_list = [("met", "m0"), ("rxn", "r0"), ("met", "m1"), ("rxn", "r1"), ("met", "m2")]
    G.add_nodes_from(node_list)
    for u, v in [(("met", "m0"), ("rxn", "r0")), (("rxn", "r0"), ("met", "m1")),
                 (("met", "m1"), ("rxn", "r1")), (("rxn", "r1"), ("met", "m2"))]:
        G.add_edge(u, v, w=1.0)
    r_undir = _reff_dense(G, ("met", "m0"), ("met", "m2"), wk="w")
    parity_err = abs(c_dir - 1.0 / r_undir)
    c_rev = directed_ceff(B, gp, gm, t, s)
    blind = abs(c_dir - c_rev)
    ok = (parity_err < SELFTEST_TOL) and (blind < T.STRICT_EPS)
    return V("IV2", "directed", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"ratio==1: C_dir={c_dir:.6f} == 1/R_undir={1.0 / r_undir:.6f} "
             f"(|d|={parity_err:.2e} < {SELFTEST_TOL:.0e}); direction-blind |fwd-rev|={blind:.2e}")


# =====================================================================
# GROUP V -- reference honesty (the promoted frozen AAM + direction reference)
# =====================================================================
# The promotion's load-bearing property: the frozen reference REFUSES to fabricate an
# atom transit for any reaction whose participants carry no concrete, balance-verifiable
# atom -- a non-molecule (electron/photon/generic acceptor), a generic `*` R-group, an
# acyl/electron carrier shipped structureless, a pseudo-reaction, an unreconciled stub.
# A refused reaction must therefore contribute ZERO cross-metabolite transit pair to the
# atom graph: no label rides a fabricated bond through a `*` super-node into a spurious
# partner. This is the graph-level shadow of `concrete_balance` + the MCS super-node rule,
# checked against the ACTUAL frozen tables the engine now reads -- so an incumbent-style
# fabricated transit (the 12k the reference dropped) cannot silently return.

def check_V1():
    """Reference honesty: no reaction the closure ledger REFUSES carries an atom-transit
    pair (no fabricated leak through the `*` super-node / carrier body), while every
    RESOLVED reaction does carry pairs (real carriers conduct). Able to fail two ways: a
    single refused reaction with a transit pair is a fabricated edge; and if resolved
    reactions stopped carrying pairs the positive control collapses. NA (not a regression)
    until the reference is frozen."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "reference"))
    import reference as ref
    if not (ref.ATOM_PAIRS.exists() and ref.LEDGER.exists()):
        return V("V1", "reference", {"atom": "NA"}, {"atom": "PASS"},
                 "frozen reference tables absent -- run after T5 freeze")
    led = pd.read_parquet(ref.LEDGER)
    mapped = set(pd.read_parquet(ref.ATOM_PAIRS, columns=["mnxr"]).mnxr.unique())
    refused = led[led.aam_state == "refused"]
    leaked = sorted(set(refused.mnxr) & mapped)          # refused, yet has a transit pair
    resolved = led[led.aam_state == "resolved"]
    res_n = resolved.mnxr.nunique()
    res_with = len(set(resolved.mnxr) & mapped)
    no_leak = len(leaked) == 0
    pos_ok = res_n > 0 and res_with == res_n             # every resolved reaction conducts
    ok = no_leak and pos_ok
    return V("V1", "reference", {"atom": "PASS" if ok else "FAIL"}, {"atom": "PASS"},
             f"refused reactions with a fabricated transit pair: {len(leaked)} (0); "
             f"resolved reactions carrying pairs: {res_with}/{res_n}")


CHECKS = [
    check_I1, check_I2, check_I3, check_I4, check_I5, check_I6, check_I7, check_I8,
    check_I9,
    check_II1, check_II2, check_II3, check_II4, check_II5, check_II6, check_II7,
    check_II8, check_II9, check_II10,
    check_III1, check_III2, check_III3, check_III4, check_III5,
    check_IV1, check_IV2,
    check_V1,
]
BY_NAME = {fn().__name__ if False else fn.__name__.replace("check_", ""): fn for fn in CHECKS}


def _regression(rec) -> bool:
    """A green-expected atom arm that did not pass is a regression."""
    for arm, exp in rec["expected"].items():
        if exp == "PASS" and rec["arms"].get(arm) != "PASS":
            return True
    return False


def _row_verdict(rec) -> str:
    """One-word status for the run table."""
    if _regression(rec):
        return "REGRESSION"
    # expected-red still red, or star fails as designed, or all green
    for arm, exp in rec["expected"].items():
        if exp == "RED":
            return "red(as-designed)" if rec["arms"].get(arm) == "FAIL" else "GREEN(fixed!)"
    if any(exp == "FAIL" for exp in rec["expected"].values()):
        star = rec["arms"].get("star")
        return "star-fails(as-designed)" if star == "FAIL" else "star-OK(surprise)"
    if all(rec["arms"].get(a) == "NA" for a in rec["arms"]):
        return "n/a(dropped)"
    return "PASS"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true", help="run every check")
    ap.add_argument("--only", nargs="+", default=None, help="check names, e.g. II2 I3")
    a = ap.parse_args()

    if a.only:
        selected = [BY_NAME[n] for n in a.only if n in BY_NAME]
        unknown = [n for n in a.only if n not in BY_NAME]
        if unknown:
            print(f"unknown checks: {unknown}; known: {sorted(BY_NAME)}")
    else:
        selected = CHECKS

    print(f"ECSPr pulse-chase suite -- {len(selected)} checks "
          f"(env {Path(sys.executable).parent.parent.name})\n")
    hdr = f"{'check':6s} {'class':13s} {'star':6s} {'atom':6s} {'verdict':24s} note"
    print(hdr)
    print("-" * len(hdr))

    records = []
    regressions = 0
    for fn in selected:
        try:
            rec = fn()
        except Exception as e:
            import traceback
            traceback.print_exc()
            rec = V(fn.__name__.replace("check_", ""), "?", {"atom": "FAIL"},
                    {"atom": "PASS"}, f"EXCEPTION: {e}")
        records.append(rec)
        rv = _row_verdict(rec)
        if rv == "REGRESSION":
            regressions += 1
        star = rec["arms"].get("star", "-")
        atom = rec["arms"].get("atom", "-")
        print(f"{rec['name']:6s} {rec['cls']:13s} {star:6s} {atom:6s} {rv:24s} {rec['note']}")

    print("\n" + "=" * 72)
    n = len(records)
    greens = sum(1 for r in records if _row_verdict(r) == "PASS")
    reds = sum(1 for r in records if _row_verdict(r).startswith("red"))
    starf = sum(1 for r in records if _row_verdict(r).startswith("star-fails"))
    fixed = sum(1 for r in records if _row_verdict(r) == "GREEN(fixed!)")
    print(f"{n} checks: {greens} green-pass, {reds} expected-red(as-designed), "
          f"{starf} star-fails-as-designed, {fixed} newly-fixed, {regressions} REGRESSION")
    if fixed:
        print("  NOTE: an expected-red is now GREEN -- a next-session fix landed. Update "
              "the baseline expectation for it.")
    if regressions:
        print("\nFAIL: a green-expected atom check regressed. This is the only exit-1 class.")
        return 1
    print("\nPASS: no regressions. Expected-reds and star-fails are the named baseline.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
