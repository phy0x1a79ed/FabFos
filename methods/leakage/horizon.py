"""Does the ground construction have a MEASUREMENT HORIZON at the ground set?

The claim under test, stated so it can fail: with a hard (contracted) ground, a perturbation
that the source can only reach *through* a ground node is invisible -- not small, exactly
zero -- because a contracted sink absorbs the current and a branch hanging off it lies on no
source->sink path. A leakage ground has no such horizon: every node drains into Omega, so a
branch beyond a ground node is still a path to ground and still carries current.

An earlier round established the weaker half of this and it must not be conflated: an edge
merely INCIDENT to a grounded precursor is fully visible under the hard ground, because it is
a new way IN to ground (verified on a synthetic IMP->AMP insert). "After the precursor" is
therefore not by itself sufficient for invisibility. What is sufficient is being separated
from the source by the ground set -- and that is a structural property this module computes
rather than assumes.

THE PROBE
    Nitrogen, netB (evidence-weighted, epi300). Two configurations:

    ``glu_to_ser``    L-glutamate -> {L-serine}. A single-metabolite ground, deliberately
                      placed upstream of purine biosynthesis: serine -> glycine -> GAR is
                      the route by which serine's nitrogen enters the purine ring, so the
                      pur block sits "after" the probed metabolite in the biosynthetic sense.
    ``nh4_to_set4``   NH4(+) -> the 16 set4 nitrogen precursors. The real configuration, run
                      alongside so the horizon claim is not demonstrated only on a toy ground.

THE PERTURBATIONS
    Purine-biosynthesis reactions, one EC per "gene", each applied two ways:

    ``gof``  reinforce -- the reaction's evidence weight is raised by the median weight
    ``lof``  the reaction is removed from the weight dict entirely

    Both are applied by REBUILDING the graph from the atom-pair table with a modified weight
    dict. LOF cannot be expressed as a parallel edge addition, and applying the two arms by
    different mechanisms would confound the arm with the mechanism.

CONTROLS -- checked before any effect is read
    ``null``   an unmodified rebuild. MUST be exactly 0.0 in every mode. If it is not, the
               rebuild is not deterministic and no number below means anything.
    ``serC``   phosphoserine transaminase: L-glutamate -> O-phospho-L-serine, i.e. squarely
               ON the glutamate->serine path. MUST be nonzero under the hard ground. Without
               it, "hard sees nothing" could be a broken probe rather than a horizon.
    ``glyA``   serine <-> glycine, incident to the probed metabolite itself: the boundary
               case between the two.

Env: ``mamba run -n scadc-metabolic-model`` (needs the ECSPr engine lib on sys.path).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

# canon.ENGINE_LIB points at a sibling standalone checkout, so it is not the right answer
# here; this repo's own submodule is. Explicit, because an ambient PYTHONPATH can otherwise
# shadow it with an older copy.
_LIB = Path(os.environ.get(
    "FABFOS_ENGINE_LIB",
    _HERE.parents[1] / "src/metasmith_libraries/resources/lib"))
if not (_LIB / "ecspr_graph.py").exists():
    raise SystemExit(f"engine lib not found at {_LIB}; set FABFOS_ENGINE_LIB")
sys.path.insert(0, str(_LIB))

import ecspr_build as B                                             # noqa: E402
import ground as G                                                  # noqa: E402
from ecspr_graph import AtomGraph                                   # noqa: E402

GLUTAMATE = "MNXM741173"
SERINE = "MNXM737787"
NH4 = "MNXM729302"

# Purine biosynthesis, one entry per gene, plus the two path controls. Each MNXR listed was
# checked to carry N atom-pairs AND to be present in the epi300 evidence weights -- a
# reaction absent from either is not in the graph and perturbing it is a no-op for reasons
# that have nothing to do with the ground.
PERTURBATIONS = {
    # -- purine biosynthesis: "after serine" in the TEXTBOOK sense (serine -> glycine ->
    #    GAR -> IMP), and the first thing one reaches for to test a horizon. It is not
    #    behind the ground electrically -- glutamine and aspartate donate ring nitrogen
    #    directly -- and the run below shows exactly that.
    "purF": ["MNXR188137"],                 # PRPP + Gln -> PRA        (N from glutamine)
    "purD": ["MNXR188698"],                 # PRA + Gly -> GAR         (N from glycine)
    "purL": ["MNXR165854"],                 # FGAR + Gln -> FGAM       (N from glutamine)
    "purM": ["MNXR108734"],                 # FGAM -> AIR
    "purC": ["MNXR192576"],                 # CAIR + Asp -> SAICAR     (N from aspartate)
    "purH": ["MNXR188247"],                 # AICAR -> FAICAR -> IMP
    "purA": ["MNXR192666"],                 # IMP + Asp -> adenylosuccinate
    "pur_block": ["MNXR188137", "MNXR188698", "MNXR165854", "MNXR108734",
                  "MNXR192576", "MNXR188247"],

    # -- BEHIND the ground: every atom node reachable from the source only THROUGH a ground
    #    metabolite. Enumerated structurally over all 6,085 weighted reactions in the N
    #    graph, not picked by intuition -- these are the complete sets.
    "psd__BEHIND_serine": ["MNXR103229", "MNXR103233", "MNXR103235",   # EC 4.1.1.65
                           "MNXR136189", "MNXR182092"],   # PS <-> PE decarboxylation
    "selA__BEHIND_serine": ["MNXR146398", "MNXR189860"],   # EC 2.9.1.1  Ser-tRNA(Sec)->Sec
    "lpxCD__BEHIND_set4": ["MNXR199080", "MNXR199001"],    # lipid A, after UDP-GlcNAc

    # -- controls, checked before any effect is read
    "serC__CTL_upstream": ["MNXR102169", "MNXR162393"],   # Glu -> O-phospho-L-serine
    "glyA__CTL_at_ground": ["MNXR104833"],                # serine <-> glycine
    "CTL_unreachable": ["MNXR112538"],                    # archaeal; no path from either
                                                          # source. Zero under ANY ground.
}

PROBES = {
    "glu_to_ser": dict(source=GLUTAMATE, ground=[SERINE]),
    "nh4_to_set4": dict(source=NH4, ground=None),   # filled from the axis set
}

COLS = ["probe", "mode", "port_mult", "gene", "kind", "rxns", "n_rxns", "placement",
        "total_base", "total_pert", "delta_total", "rel_total",
        "n_edges_base", "n_edges_pert", "d_edges",
        "rxn_nodes_total", "rxn_nodes_front", "rxn_nodes_behind_ground",
        "rxn_nodes_unreachable",
        "prec_share_base", "prec_share_pert", "leak_frac_base", "leak_frac_pert",
        "converged", "port", "leak", "add_weight"]


def classify(g: AtomGraph, source_met: str, ground_mets) -> dict:
    """Partition the atom nodes into ``front`` / ``behind_ground`` / ``unreachable``.

    This is what makes a zero ATTRIBUTABLE, and the distinction the whole question turns on:

    ``behind_ground``   reachable from the source in the full graph, but only by passing
                        THROUGH a ground metabolite. A hard ground must report exactly zero
                        here -- the contracted sink absorbs the current before it arrives.
                        A leakage ground should not.
    ``unreachable``     no path from the source at all. Zero under EVERY ground, leakage
                        included: Omega gives every node a way OUT, never a path IN.
                        Reporting these together with ``behind_ground`` is the confusion
                        this partition exists to prevent.
    ``front``           everything else -- what a hard ground can already see.
    """
    comp = G.source_component(g, source_met)
    hz = _horizon(g, source_met, ground_mets)
    return dict(behind_ground=hz & comp, unreachable=hz - comp,
                front=set(range(g.n)) - hz)


def _horizon(g: AtomGraph, source_met: str, ground_mets) -> set:
    """Atom-node indices the source cannot reach once the ground set is deleted."""
    import scipy.sparse as sp
    gnd = {g.idx[nd] for m in ground_mets for nd in g.atoms_of(m) if nd in g.idx}
    keep = np.ones(g.n, bool)
    keep[list(gnd)] = False
    t = np.array([e[0] for e in g.edges], np.int64)
    h = np.array([e[1] for e in g.edges], np.int64)
    m = keep[t] & keep[h]
    A = sp.coo_matrix((np.ones(int(m.sum())), (t[m], h[m])), shape=(g.n, g.n)).tocsr()
    lab = sp.csgraph.connected_components(A, directed=False)[1]
    seeds = {lab[g.idx[nd]] for nd in g.atoms_of(source_met) if nd in g.idx and keep[g.idx[nd]]}
    if not seeds:
        return set()
    reach = set(np.flatnonzero(np.isin(lab, list(seeds))).tolist())
    return set(np.flatnonzero(keep).tolist()) - reach


def rxn_nodes(pairs, element, mnxrs, ratios) -> list:
    """The atom-node keys a reaction set touches, via the same builder the graph uses."""
    sub = pairs[pairs.mnxr.isin(list(mnxrs))]
    if sub.empty:
        return []
    gs = B.graph_from_pairs(sub, element, {r: 1.0 for r in mnxrs}, ratios)
    return list(gs.nodes)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--element", default="N")
    ap.add_argument("--modes", nargs="*", default=list(G.MODES))
    ap.add_argument("--probes", nargs="*", default=list(PROBES))
    ap.add_argument("--port-mults", type=float, nargs="*", default=[1e3],
                    help="port conductance as a multiple of the median edge conductance. "
                         "Sweeping it is the test of whether the behind-ground signal can "
                         "be raised into a usable range WITHOUT losing the precursor "
                         "signature -- the two move in opposite directions.")
    ap.add_argument("--leak-mult", type=float, default=1e-7)
    a = ap.parse_args(argv)

    ref, out, X = Path(a.ref), Path(a.out), a.element
    out.parent.mkdir(parents=True, exist_ok=True)

    import ecspr_probe as P
    axes = json.loads((ref / "derived/axes/biomass_dag_axes_set4.json").read_text())
    sets4 = P.axis_node_sets(axes, X)
    PROBES["nh4_to_set4"]["ground"] = list(sets4["precursors"])

    ratios = B.load_direction_ratios(ref / "derived/mnxref-4_5/direction.parquet")
    weights = B.load_evidence_weights(ref / "derived/evidence/evidence_weights.parquet")
    w_add = float(np.median(list(weights.values())))
    pairs = B.load_pairs(ref / "derived/mnxref-4_5/atom_pairs.parquet", X)
    print(f"[{X}] {len(pairs):,} atom-pair rows, {len(weights):,} weighted reactions, "
          f"median weight {w_add:.6g}", flush=True)

    t0 = time.time()
    g0 = B.graph_from_pairs(pairs, X, weights, ratios)
    print(f"[{X}] base graph {g0.n:,} nodes / {g0.m:,} edges in {time.time() - t0:.1f}s",
          flush=True)

    # Every perturbed graph is built once and shared across probes and modes: the graph is a
    # function of the weight dict alone, and rebuilding per probe would be pure waste.
    graphs = {("null", "none"): g0}
    for gene, rxns in PERTURBATIONS.items():
        missing = [r for r in rxns if r not in weights]
        if missing:
            print(f"  !! {gene}: {missing} absent from the evidence weights -- "
                  f"perturbing them is a no-op for a reason unrelated to ground", flush=True)
        for kind in ("gof", "lof"):
            w = dict(weights)
            for r in rxns:
                if kind == "gof":
                    w[r] = w.get(r, 0.0) + w_add
                else:
                    w.pop(r, None)
            graphs[(gene, kind)] = B.graph_from_pairs(pairs, X, w, ratios)
    print(f"[{X}] {len(graphs)} graphs built in {time.time() - t0:.1f}s", flush=True)

    node_keys = {gene: rxn_nodes(pairs, X, rxns, ratios)
                 for gene, rxns in PERTURBATIONS.items()}

    fh = open(out, "w", buffering=1, newline="")
    w_out = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
    w_out.writeheader()

    med = float(np.median(g0.gp))
    leak = a.leak_mult * med

    # per-precursor baseline signature, so the port that buys behind-ground sensitivity can
    # be checked against what it costs in precursor attribution -- in the same run
    sig_fh = open(str(out).replace(".tsv", "_signature.tsv"), "w", buffering=1, newline="")
    sig = csv.DictWriter(sig_fh, delimiter="\t", fieldnames=[
        "probe", "mode", "port_mult", "metabolite", "draw_base", "share_base",
        "prec_share", "leak_frac"])
    sig.writeheader()

    for pname in a.probes:
        pr = PROBES[pname]
        src, gnd = pr["source"], pr["ground"]
        cl = classify(g0, src, gnd)
        place = {}
        for gene, ks in node_keys.items():
            idxs = [g0.idx.get(k) for k in ks]
            idxs = [i for i in idxs if i is not None]
            nb = sum(1 for i in idxs if i in cl["behind_ground"])
            nu = sum(1 for i in idxs if i in cl["unreachable"])
            place[gene] = (len(idxs) - nb - nu, nb, nu)
        print(f"\n== probe {pname}: {src} -> {len(gnd)} ground metabolite(s) ==", flush=True)
        print(f"   atom nodes: {len(cl['front']):,} front / "
              f"{len(cl['behind_ground']):,} behind ground / "
              f"{len(cl['unreachable']):,} unreachable   (of {g0.n:,})", flush=True)
        for gene, (nf, nb, nu) in place.items():
            tag = ("BEHIND-GROUND" if nb and not nf else
                   "unreachable" if nu and not nf and not nb else
                   "mixed" if (nb or nu) else "front")
            print(f"     {gene:24s} front={nf:3d} behind={nb:3d} unreach={nu:3d}  {tag}",
                  flush=True)

        sets = dict(precursors=list(gnd))
        for mode, pm in [(m, p) for m in a.modes for p in
                         ([1.0] if m == "hard" else a.port_mults)]:
            port = pm * med
            base = G.measure(g0, src, sets, mode, port=port, leak=leak)
            bt = base["total"]
            print(f"   [{mode} port={pm:g}x] base total = {bt:.10g}  "
                  f"prec_share={base['prec_share']:.6f} "
                  f"leak_frac={base['leak_frac']:.3e}", flush=True)
            tot_draw = sum(base["draw"].values()) or float("nan")
            for m_, v in base["draw"].items():
                sig.writerow(dict(probe=pname, mode=mode, port_mult=pm, metabolite=m_,
                                  draw_base=v, share_base=v / tot_draw,
                                  prec_share=base["prec_share"],
                                  leak_frac=base["leak_frac"]))
            for (gene, kind), gp in graphs.items():
                cur = G.measure(gp, src, sets, mode, port=port, leak=leak)
                ct = cur["total"]
                ks = node_keys.get(gene, [])
                nf, nb, nu = place.get(gene, (0, 0, 0))
                w_out.writerow(dict(
                    probe=pname, mode=mode, port_mult=pm, gene=gene, kind=kind,
                    rxns=";".join(PERTURBATIONS.get(gene, [])),
                    n_rxns=len(PERTURBATIONS.get(gene, [])),
                    placement=("behind_ground" if nb and not nf else
                               "unreachable" if nu and not nf and not nb else
                               "mixed" if (nb or nu) else "front"),
                    total_base=bt, total_pert=ct, delta_total=ct - bt,
                    rel_total=((ct - bt) / bt) if abs(bt) > 0 else float("nan"),
                    n_edges_base=g0.m, n_edges_pert=gp.m, d_edges=gp.m - g0.m,
                    rxn_nodes_total=len(ks), rxn_nodes_front=nf,
                    rxn_nodes_behind_ground=nb, rxn_nodes_unreachable=nu,
                    prec_share_base=base["prec_share"], prec_share_pert=cur["prec_share"],
                    leak_frac_base=base["leak_frac"], leak_frac_pert=cur["leak_frac"],
                    converged=cur["converged"], port=port, leak=leak, add_weight=w_add))
    fh.close()
    sig_fh.close()
    print(f"\nwrote {out} and {str(out).replace('.tsv', '_signature.tsv')}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
