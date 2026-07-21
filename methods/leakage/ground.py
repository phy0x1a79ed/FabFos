"""Ground constructions for the ECSPr media->biomass probe, as comparable modes.

Three modes, all measured through the identical solve so a difference between them is a
difference of GROUND and nothing else:

``hard``
    The incumbent. Every set4 precursor contracted into one grounded sink super-node; the
    per-precursor readout is ``Solution.delivered``. Outside the ground set net current is
    identically zero by KCL, so only precursors report a draw.

``leak_raw``
    Universal leakage as first prototyped: one virtual node Omega, an edge from EVERY atom
    node into it, with a large ``port`` conductance on precursor atom nodes and a small
    background ``leak`` elsewhere. Attached PER ATOM NODE, so a precursor's effective sink
    conductance is ``port * n_atoms``.

``leak_norm``
    The same, except each precursor's total port conductance is ``port`` and is SPLIT across
    its atom nodes. This exists because ``leak_raw`` is not ground-neutral: set4 carbon
    precursors carry 3-94 atoms, a 31x spread, so ``leak_raw`` silently weights a precursor
    by its carbon count and at low port ranks precursors by size rather than by routing.
    ``leak_norm`` removes that confound; comparing the two IS the test of whether it mattered.

The source metabolite never leaks under any mode -- a leak edge on the source is a direct
source->Omega path bypassing the network entirely.

Env: numpy + pandas + scipy, via the ECSPr engine lib. No canon, no absolute paths.
"""
from __future__ import annotations

import numpy as np

from ecspr_graph import AtomGraph, Terminal, solve

GND = ("__OMEGA__", 0)
MODES = ("hard", "leak_raw", "leak_norm")


def add_reactions(pairs, element, base_weights, mnxrs, weight, ratios, builder):
    """A graph rebuilt with ``mnxrs`` added to the weight dict -- the GOF perturbation.

    The engine takes no perturbation argument by design; a gain-of-function is a LARGER
    reaction set, hence a different network, and the delta is a subtraction the caller does
    across two independent solves.
    """
    w = dict(base_weights)
    for r in mnxrs:
        w[r] = w.get(r, 0.0) + float(weight)
    return builder(pairs, element, w, ratios)


def merge_reactions(g: AtomGraph, pairs, element: str, mnxrs, weight: float, ratios):
    """``(graph, n_edges_added)`` -- ``mnxrs``' atom transfers added as PARALLEL edges.

    Electrically identical to rebuilding the whole graph with ``w[r] += weight``: for an
    already-present reaction the existing edge carries ``w0*pair_w`` and the new one
    ``w_add*pair_w``, and parallel conductances sum, giving ``(w0+w_add)*pair_w`` -- exactly
    what re-weighting produces. The backward branch scales by the same ratio, and the
    ratio>1 edge flip is applied identically because it goes through the same builder.

    Used for BOTH lanes on purpose. netA cannot be rebuilt per condition without cobra
    (`ecspr_build._load_model`), and applying the perturbation by two different mechanisms
    would confound a netA-vs-netB difference with a difference in how it was perturbed.
    It is also ~1000x faster than a full rebuild per condition.
    """
    import ecspr_build as _B
    sub = pairs[pairs.mnxr.isin(list(mnxrs))]
    if sub.empty:
        return g, 0
    gadd = _B.graph_from_pairs(sub, element, {r: float(weight) for r in mnxrs}, ratios)
    if not gadd.m:
        return g, 0
    nodes = list(g.nodes)
    idx = dict(g.idx)

    def _i(k):
        j = idx.get(k)
        if j is None:
            j = idx[k] = len(nodes)
            nodes.append(k)
        return j

    edges = list(g.edges)
    for (a, b) in gadd.edges:
        edges.append((_i(gadd.nodes[a]), _i(gadd.nodes[b])))
    return (AtomGraph(nodes, edges, np.concatenate([g.gp, gadd.gp]),
                      np.concatenate([g.gm, gadd.gm]), dict(g.meta)), gadd.m)


def _attach(g: AtomGraph, source_met: str, precursors, port: float, leak: float,
            normalize: bool):
    """``(graph_with_Omega, {edge_index: metabolite})``.

    Leak edges are symmetric (``gm == gp``): dilution is not rectified, and handing them to
    the diode would push them through ``DIODE_BACKWARD_FLOOR`` for no reason.
    """
    prec = set(precursors)
    nodes = list(g.nodes) + [GND]
    gi = len(nodes) - 1
    edges, gp, gm, owner = list(g.edges), list(g.gp), list(g.gm), {}
    natoms = {m: max(1, len(g.atoms_of(m))) for m in prec}
    for i, nd in enumerate(g.nodes):
        met = nd[0]
        if met == source_met:
            continue
        if met in prec:
            c = port / natoms[met] if normalize else port
        else:
            c = leak
        if c <= 0:
            continue
        owner[len(edges)] = met
        edges.append((i, gi))
        gp.append(c)
        gm.append(c)
    return (AtomGraph(nodes, edges, np.asarray(gp, float), np.asarray(gm, float),
                      dict(g.meta)), owner)


def measure(g: AtomGraph, media: str, sets: dict, mode: str, *, port=None, leak=None,
            tol=None) -> dict:
    """One solve. Returns ``{draw: {met: current}, total, leak_frac, converged, ...}``.

    ``draw`` is the per-metabolite current LEAVING the network at that metabolite. Under
    ``hard`` that is ``delivered`` (net inflow at a contracted precursor); under the leak
    modes it is the current in that metabolite's own port edge, because with a finite port
    the metabolite is an interior node and its net inflow is zero by KCL.
    """
    precursors = list(sets["precursors"])
    if mode == "hard":
        src = Terminal.metabolite(g, media, label=media)
        gnd = Terminal.merge(g, precursors, label="ground")
        sol = solve(g, src, gnd, tol=tol)
        draw = {m: (sol.delivered(m) if g.atoms_of(m) else 0.0) for m in precursors}
        prec_share = sum(draw.values())
        return dict(draw=draw, total=sol.total, prec_share=prec_share,
                    leak_frac=0.0, converged=int(sol.converged), sol=sol, graph=g)

    g2, owner = _attach(g, media, precursors, port, leak, normalize=(mode == "leak_norm"))
    sol = solve(g2, Terminal.metabolite(g2, media, label=media),
                Terminal.of_nodes("omega", [GND]), tol=tol)
    oe, cur = sol.edge_currents()
    allc = {}
    for k in range(len(oe)):
        m = owner.get(int(oe[k]))
        if m is not None:
            allc[m] = allc.get(m, 0.0) + float(cur[k])
    draw = {m: allc.get(m, 0.0) for m in precursors}
    prec_share = sum(draw.values())
    return dict(draw=draw, total=sol.total, prec_share=prec_share,
                leak_frac=sum(allc.values()) - prec_share,
                converged=int(sol.converged), sol=sol, graph=g2)


def source_component(g: AtomGraph, media: str) -> set:
    """Atom nodes undirected-reachable from ``media``'s atoms in ``g``.

    The diagnostic that makes a zero attributable. A gain-of-function whose edges land
    outside this set cannot move any readout under ANY ground -- universal leakage included,
    because Omega gives every node a way OUT but never a path IN from the source. Reporting
    such a condition as "no effect" alongside genuine no-effects is the confusion this
    exists to prevent.
    """
    import scipy.sparse as sp
    if not g.m:
        return set()
    t = np.array([e[0] for e in g.edges], np.int64)
    h = np.array([e[1] for e in g.edges], np.int64)
    A = sp.coo_matrix((np.ones(t.size), (t, h)), shape=(g.n, g.n)).tocsr()
    lab = sp.csgraph.connected_components(A, directed=False)[1]
    seeds = {lab[g.idx[nd]] for nd in g.atoms_of(media) if nd in g.idx}
    return set(np.flatnonzero(np.isin(lab, list(seeds))).tolist()) if seeds else set()


def added_edge_report(g_base: AtomGraph, g_cond: AtomGraph, comp: set, media: str) -> dict:
    """How the perturbation actually landed, counted over the edges APPENDED by
    :func:`merge_reactions` -- which appends, so the new edges are exactly the tail.

    A previous version compared endpoint-pair KEYS against the base edge set and skipped any
    pair already present. That silently reported ``n_edges_added=0, reachable=0`` for a
    reaction whose transfers are all PARALLEL to existing ones -- even though a parallel edge
    raises the conductance and moves the answer. It made 27 real effects look like a
    negative-control floor, i.e. it manufactured a noise floor out of true signal. Counting
    the appended tail cannot make that mistake.

    ``touches_source_component`` asks whether an appended edge has an endpoint in the base
    graph's source component. A perturbation entirely outside it cannot move any readout
    under ANY ground -- Omega gives every node a way OUT, never a path IN from the source.
    """
    n_new = g_cond.m - g_base.m
    if n_new <= 0:
        return dict(n_edges_added=0, n_edges_touching_source_component=0, reachable=0,
                    n_parallel_added=0)
    idx_b, n_touch, n_par = g_base.idx, 0, 0
    base_pairs = {(g_base.nodes[a], g_base.nodes[b]) for a, b in g_base.edges}
    for a, b in g_cond.edges[g_base.m:]:
        ka, kb = g_cond.nodes[a], g_cond.nodes[b]
        if (ka, kb) in base_pairs:
            n_par += 1
        if (idx_b.get(ka) in comp) or (idx_b.get(kb) in comp):
            n_touch += 1
    return dict(n_edges_added=n_new, n_edges_touching_source_component=n_touch,
                reachable=int(n_touch > 0), n_parallel_added=n_par)


def throughput_at(res: dict, mets) -> dict:
    """Throughput (``max(inflow, outflow)``) per metabolite -- the quantity to read at a
    CENTRAL metabolite, where net current is 0 by KCL and therefore says nothing.

    Reported as a current, not a conductance: it is a fraction of the unit injection.
    """
    g, sol = res["graph"], res["sol"]
    return {m: (sol.throughput_metabolite(m) if g.atoms_of(m) else 0.0) for m in mets}
