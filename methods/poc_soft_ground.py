"""POC: where does the measurement horizon sit, and can a soft ground move it?

THE DEFECT, STATED PRECISELY
----------------------------
The incumbent probe contracts every set4 biomass precursor into one super-node. Current
leaves the system ONLY there. So for any metabolite outside the ground set, net current is
identically zero -- not small, zero, by KCL, because every path out eventually returns to
the ground. Such a metabolite can carry throughput but cannot DRAW. "How much does GMP get"
is not a question the instrument can express, and a gain-of-function that lands past a
precursor is unmeasurable rather than merely weak.

WHAT THIS TRIES
---------------
A virtual ground node Omega with an edge from every atom node into it, so every metabolite
is a sink and current can stop anywhere. Two configurations:

  * UNIFORM  -- one leak conductance everywhere. What the "just connect everything" idea
    actually does.
  * TWO-TIER -- a large port at the biomass precursors, a small background leak elsewhere.
    Still not a calibration: one conductance per tier, no per-metabolite fitting.

THE TEST
--------
Two numbers per configuration, and they pull against each other:

  * PRECURSOR SHARE -- fraction of injected current that reaches a biomass precursor. If
    this collapses, the leak has replaced the measurement rather than regularised it.
  * d(IMP draw) under a GOF placed STRICTLY BEYOND the horizon -- a synthetic edge from IMP
    to AMP, i.e. downstream of a precursor and nowhere else. Under the hard ground this is
    exactly zero by construction. Anything else means the horizon moved.

IMP is a set4 carbon precursor (a sink of both nucleotide axes, never a source), which the
script asserts rather than assumes. XMP / GMP / AMP / guanosine sit beyond it.

    PYTHONPATH=src mamba run -n scadc-metabolic-model python methods/poc_soft_ground.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent / "src"))
sys.path.insert(0, str(_HERE.parent / "src/metasmith_libraries/resources/lib"))

from fabfos import canon                                            # noqa: E402
import ecspr_build as B                                             # noqa: E402
import ecspr_probe as P                                             # noqa: E402
from ecspr_graph import AtomGraph, Terminal, solve                  # noqa: E402

ELEMENT = "C"
GND = ("__OMEGA__", 0)

IMP, XMP, GMP, AMP, GUO = ("MNXM1101868", "MNXM1104385", "MNXM1101285",
                           "MNXM728294", "MNXM1103400")
LADDER = {IMP: "IMP", XMP: "XMP", GMP: "GMP", AMP: "AMP", GUO: "guanosine"}


# =====================================================================
# Graph surgery
# =====================================================================

def add_edge(g: AtomGraph, tail, head, gp: float, gm: float = None) -> AtomGraph:
    """``g`` plus one explicit atom-transfer edge -- the synthetic GOF.

    Deliberately synthetic rather than a real unused MNXR: adding a reaction adds ALL of its
    atom pairs, and a reaction chosen for one IMP -> AMP pair drags in pairs elsewhere in
    the network, which would confound "strictly beyond the horizon" with "also feeds the
    precursor". One edge has no such ambiguity.
    """
    return AtomGraph(list(g.nodes), list(g.edges) + [(g.idx[tail], g.idx[head])],
                     np.append(g.gp, gp), np.append(g.gm, gp if gm is None else gm),
                     dict(g.meta))


def attach_ground(g: AtomGraph, source_met: str, g_port: dict, g_default: float):
    """``g`` plus Omega and one leak edge per atom node. Returns ``(graph, edge -> met)``.

    Leak edges are SYMMETRIC (``gm == gp``): dilution is not rectified, and handing them to
    the diode would put them through ``DIODE_BACKWARD_FLOOR`` for no reason.

    The SOURCE metabolite never gets a leak edge. One there is a direct source -> Omega path
    bypassing the whole network -- the single short this construction can create.
    """
    nodes = list(g.nodes) + [GND]
    gi = len(nodes) - 1
    edges, gp, gm, owner = list(g.edges), list(g.gp), list(g.gm), {}
    for i, nd in enumerate(g.nodes):
        if nd[0] == source_met:
            continue
        c = g_port.get(nd[0], g_default)
        if c <= 0:
            continue
        owner[len(edges)] = nd[0]
        edges.append((i, gi))
        gp.append(c)
        gm.append(c)
    return (AtomGraph(nodes, edges, np.asarray(gp, float), np.asarray(gm, float),
                      dict(g.meta)), owner)


def soft_solve(g: AtomGraph, media: str, g_port: dict, g_default: float, *, tol=None):
    """Solve ``media -> Omega``; return ``(solution, {metabolite: port current}, graph)``.

    The draw is read off each metabolite's OWN port edge. ``Solution.delivered`` returns ~0
    here and correctly so -- with a finite port the metabolite is an interior node, so its
    net inflow is zero by KCL. That the readout MOVES from net-inflow to port-current is the
    whole difference between a soft sink and a contracted one.
    """
    g2, owner = attach_ground(g, media, g_port, g_default)
    sol = solve(g2, Terminal.metabolite(g2, media, label=media),
                Terminal.of_nodes("omega", [GND]), tol=tol)
    oe, cur = sol.edge_currents()
    draw = {}
    for k in range(len(oe)):
        m = owner.get(int(oe[k]))
        if m is not None:
            draw[m] = draw.get(m, 0.0) + float(cur[k])
    return sol, draw, g2


def hard_draw(g: AtomGraph, media: str, sets: dict):
    """The incumbent probe as a ``{metabolite: net current}`` map over the ladder."""
    _, sol = P.probe(g, media, sets)
    return sol, {m: (sol.delivered(m) if g.atoms_of(m) else 0.0) for m in LADDER}


# =====================================================================

def main() -> int:
    print("=" * 92)
    print("POC: the measurement horizon under a contracted ground vs a distributed one")
    print("=" * 92)

    axes = json.loads(Path(canon.AXES_JSON).read_text())
    sets = P.axis_node_sets(axes, ELEMENT)
    precursors = set(sets["precursors"])
    media = sets["catabolic_sources"][0] if sets["catabolic_sources"] else sets["media"][0]
    assert IMP in precursors, "IMP is expected to be a set4 carbon precursor"
    assert not (precursors & {XMP, GMP, AMP, GUO}), "the ladder must sit OUTSIDE the ground"

    pairs = B.load_pairs(canon.REFERENCE_ATOM_PAIRS, ELEMENT)
    ratios = B.load_direction_ratios(canon.REFERENCE_DIRECTION)
    weights = B.load_evidence_weights(canon.EVIDENCE_WEIGHTS)
    g = B.graph_from_pairs(pairs, ELEMENT, weights, ratios)
    med = float(np.median(g.gp))
    n_prec_atoms = sum(len(g.atoms_of(m)) for m in precursors)
    print(f"\nmedia {media}   graph {g.n:,} nodes / {g.m:,} edges / "
          f"{g.meta['n_reactions_used']:,} reactions")
    print(f"median edge conductance {med:.6g}   precursors {len(precursors)} "
          f"({n_prec_atoms} atom nodes)   ladder beyond the ground: "
          f"{[LADDER[m] for m in (XMP, GMP, AMP, GUO)]}")

    # The GOF: one edge, strictly downstream of the IMP precursor.
    t_nd, h_nd = g.atoms_of(IMP)[0], g.atoms_of(AMP)[0]
    g_gof = add_edge(g, t_nd, h_nd, med)
    print(f"\nsynthetic GOF  : {t_nd} -> {h_nd} at gp = gm = {med:.6g}")
    print("                 strictly beyond the horizon: its tail IS a precursor, so under "
          "the\n                 contracted ground it hangs off a node pinned at phi = 0.")

    # ---- 1. the incumbent ------------------------------------------------
    print("\n" + "-" * 92)
    print("1. HARD GROUND (incumbent)")
    print("-" * 92)
    t0 = time.time()
    hsol, hbase = hard_draw(g, media, sets)
    hsol2, hgof = hard_draw(g_gof, media, sets)
    print(f"   conductance to ground {hsol.total:.6g} -> {hsol2.total:.6g} with the GOF "
          f"(d = {hsol2.total - hsol.total:+.6g})   [{time.time() - t0:.1f}s]")
    print(f"\n   {'name':10s} {'in ground':10s} {'net current':>14s} {'throughput':>13s} "
          f"{'voltage':>10s} {'d under GOF':>14s}")
    for m, nm in LADDER.items():
        tp = hsol.throughput_metabolite(m) if g.atoms_of(m) else 0.0
        v = hsol.voltage_metabolite(m)["weighted_mean"] if g.atoms_of(m) else float("nan")
        print(f"   {nm:10s} {'YES' if m in precursors else 'no':10s} {hbase[m]:14.6g} "
              f"{tp:13.6g} {v:10.5g} {hgof[m] - hbase[m]:14.6g}")
    print("\n   Every off-ground row is a machine zero in BOTH current columns: outside the\n"
          "   ground set net current is identically zero by KCL, so the GOF is not weakly\n"
          "   measured -- it is unmeasurable.")

    # ---- 2 + 3. the distributed grounds ----------------------------------
    configs = [(f"uniform  leak={b:g}", {}, b * med) for b in (1e-6, 1e-4)]
    for a in (1e1, 1e3, 1e5):
        configs.append((f"two-tier port={a:g} leak=1e-07",
                        {m: a * med for m in precursors}, 1e-7 * med))

    # Two placements, both strictly beyond the horizon (tail = the IMP precursor).
    # XMP is one atom-transfer past IMP and that route ALREADY EXISTS in the base graph, so
    # this GOF strengthens a live channel. AMP is not reachable from IMP in the base, so
    # that GOF opens a new one. They are different questions and the instrument should not
    # answer them identically.
    gofs = {"IMP->XMP": add_edge(g, g.atoms_of(IMP)[0], g.atoms_of(XMP)[0], med),
            "IMP->AMP": g_gof}

    print("\n" + "-" * 92)
    print("2. DISTRIBUTED GROUND -- conductances are multiples of the median edge "
          f"({med:.4g})")
    print("-" * 92)
    print(f"   {'configuration':28s} {'prec.share':>11s} {'V(IMP)':>10s} {'V(XMP)':>10s} "
          f"{'IMP draw':>11s} {'GOF':10s} {'d G(src)':>11s} {'d IMP draw':>12s} "
          f"{'rel':>10s}")
    for name, port, dflt in configs:
        t0 = time.time()
        sol, draw, _ = soft_solve(g, media, port, dflt)
        share = sum(draw.get(m, 0.0) for m in precursors)
        b0 = draw.get(IMP, 0.0)
        vi = sol.voltage_metabolite(IMP)["weighted_mean"]
        vx = sol.voltage_metabolite(XMP)["weighted_mean"]
        head = f"   {name:28s} {share:11.6f} {vi:10.5g} {vx:10.5g} {b0:11.5g} "
        for gname, gg in gofs.items():
            sol2, draw2, _ = soft_solve(gg, media, port, dflt)
            d = draw2.get(IMP, 0.0) - b0
            print(head + f"{gname:10s} {sol2.total - sol.total:11.4g} {d:12.4g} "
                  f"{(d / b0 if b0 else float('nan')):10.2e}")
            head = " " * 76
        print(f"{'':76s}[{time.time() - t0:.0f}s]")

    print("\n   prec.share -> 1 means the biomass measurement survived; -> 0 means the leak\n"
          "   replaced it. d IMP draw away from zero means the horizon actually moved.\n"
          "   V(IMP) vs V(XMP) is the diagnostic for WHY a configuration fails: an edge\n"
          "   between two nodes at the same potential carries no current, so if the ground\n"
          "   has flattened the far field the new edge is invisible for that reason instead.")
    print("\n" + "=" * 92)
    return 0


if __name__ == "__main__":
    sys.exit(main())
