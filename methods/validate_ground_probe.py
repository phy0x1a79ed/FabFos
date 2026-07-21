"""Tiny validation run for the media -> ground probe: does the instrument actually measure?

Builds the CARBON atom graph from the tier-4 reference (atom pairs + the direction
ensemble) restricted to the epi300 evidence weights, runs the glucose -> aggregate-ground
probe, prints the per-precursor share table, then rebuilds WITHOUT one high-traffic
reaction and prints the difference.

The point of the second half is the API claim: the engine takes no perturbation argument.
A loss-of-function is a smaller reaction set, hence a different network, and the comparison
is a subtraction done HERE, by the caller, on two independent solves.

The knockout is chosen from the RUN'S OWN throughput ranking, never named in advance. A
demo that hard-codes its reaction proves nothing the moment the universe moves.

    PYTHONPATH=src mamba run -n scadc-metabolic-model python methods/validate_ground_probe.py

Every path comes from `canon`; this file transcribes none. Every number it prints it
computed.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent / "src"))

from fabfos import canon                                            # noqa: E402

# The IN-REPO submodule, deliberately -- NOT `canon.ENGINE_LIB`, which points at the
# sibling standalone checkout. This run validates the engine this repo ships.
sys.path.insert(0, str(_HERE.parent / "src/metasmith_libraries/resources/lib"))

import ecspr_build as B                                             # noqa: E402
import ecspr_probe as P                                             # noqa: E402
from ecspr_graph import Terminal, solve                             # noqa: E402

ELEMENT = "C"


def _fmt(rows, key="share", n=12):
    out = []
    for r in sorted(rows, key=lambda r: -abs(r[key]))[:n]:
        out.append(f"    {r['metabolite']:14s} {r['role']:9s} "
                   f"share={r['share']:+.6f} thru={r['throughput']:.6f} "
                   f"dV={r['voltage_drop']:.4f}")
    return "\n".join(out)


def main() -> int:
    print("=" * 78)
    print("ECSPr ground-probe validation run -- carbon, epi300 evidence, tier-4 reference")
    print("=" * 78)

    axes = json.loads(Path(canon.AXES_JSON).read_text())
    canon.assert_canonical_axes(canon.AXES_JSON)
    sets = P.axis_node_sets(axes, ELEMENT)
    print(f"\nnode sets derived from {Path(canon.AXES_JSON).name} [{ELEMENT}]:")
    print(f"  media entries : {len(sets['media'])}  {sets['media']}")
    print(f"  central       : {len(sets['central'])}")
    print(f"  precursors    : {len(sets['precursors'])} (merged into the aggregate ground)")
    print(f"  catabolic src : {sets['catabolic_sources']}")

    t0 = time.time()
    pairs = B.load_pairs(canon.REFERENCE_ATOM_PAIRS, ELEMENT)
    ratios = B.load_direction_ratios(canon.REFERENCE_DIRECTION)
    weights = B.load_evidence_weights(canon.EVIDENCE_WEIGHTS)
    print(f"\ninputs loaded in {time.time() - t0:.1f}s: {len(pairs):,} pair rows, "
          f"{len(ratios):,} direction ratios, {len(weights):,} weighted reactions")

    t0 = time.time()
    g = B.graph_from_pairs(pairs, ELEMENT, weights, ratios, with_provenance=True)
    tb = time.time() - t0
    print(f"graph built in {tb:.1f}s: {g.n:,} nodes / {g.m:,} edges over "
          f"{g.meta['n_metabolites']:,} metabolites")
    print(f"  reactions used {g.meta['n_reactions_used']:,} / "
          f"{g.meta['n_reactions_requested']:,}   AAM gap {g.meta['n_aam_gap']:,}")
    print(f"  pair rows folded in {g.meta['n_pair_rows']:,}; "
          f"{g.meta['n_reversed_rows']:,} reversed by a ratio > 1")

    media = sets["catabolic_sources"][0] if sets["catabolic_sources"] else sets["media"][0]
    print(f"\nprobing {media} -> ground ...")
    t0 = time.time()
    rows, sol = P.probe(g, media, sets)
    ts = time.time() - t0
    drawn = sum(r["share"] for r in rows if r["in_ground"])
    print(f"  solved in {ts:.1f}s ({sol.iters} Newton iterations, "
          f"converged={sol.converged})")
    print(f"  total conductance to ground = {sol.total:.6f}")
    print(f"  conservation: sum of per-precursor shares = {drawn:.12f} "
          f"(|err| = {abs(drawn - 1.0):.2e})")
    missing = [r["metabolite"] for r in rows if r["in_ground"] and not r["in_graph"]]
    print(f"  precursors absent from the graph: {len(missing)} {missing}")

    print("\n  per-precursor draw (top 12 by share):")
    print(_fmt([r for r in rows if r["in_ground"]]))
    print("\n  central metabolites (net current is 0 by KCL -- read throughput):")
    print(_fmt([r for r in rows if r["role"] == "central"], key="throughput"))

    # ---- the caller-side loss of function -------------------------------
    rank = B.reaction_currents(g, sol)
    ko = str(rank.index[0])
    print(f"\nknockout chosen from this run's own throughput ranking: {ko} "
          f"(attributed current {rank.iloc[0]:.4f}; next {rank.index[1]} "
          f"{rank.iloc[1]:.4f})")

    t0 = time.time()
    g2 = B.graph_from_pairs(pairs, ELEMENT, {k: v for k, v in weights.items() if k != ko},
                            ratios)
    rows2, sol2 = P.probe(g2, media, sets)
    print(f"  rebuilt and re-solved in {time.time() - t0:.1f}s: {g2.n:,} nodes / "
          f"{g2.m:,} edges; total = {sol2.total:.6f}")

    a = {r["metabolite"]: r for r in rows}
    b = {r["metabolite"]: r for r in rows2}
    diff = [dict(metabolite=m, role=a[m]["role"],
                 share=b[m]["share"] - a[m]["share"],
                 throughput=b[m]["throughput"] - a[m]["throughput"],
                 voltage_drop=b[m]["voltage_drop"] - a[m]["voltage_drop"])
            for m in a if m in b]
    moved = sum(abs(d["share"]) for d in diff if a[d["metabolite"]]["in_ground"]) / 2.0
    print(f"\n  DELTA (knockout minus base) -- computed by the caller, not the engine")
    print(f"  d(total) = {sol2.total - sol.total:+.6f} "
          f"({100 * (sol2.total - sol.total) / sol.total:+.2f}%)")
    print(f"  total share redistributed among precursors = {moved:.6f} "
          f"of the injected current")
    print("\n  largest per-precursor share shifts:")
    print(_fmt([d for d in diff if a[d["metabolite"]]["in_ground"]]))

    ok = abs(drawn - 1.0) < 1e-6 and sol.converged and sol.total > 0
    print("\n" + "=" * 78)
    print(f"VERDICT: {'PASS' if ok else 'FAIL'} -- conservation "
          f"{'holds' if abs(drawn - 1.0) < 1e-6 else 'BROKEN'}, solve "
          f"{'converged' if sol.converged else 'DID NOT CONVERGE'}, "
          f"redistribution {'observed' if moved > 0 else 'NOT observed'}")
    print("=" * 78)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
