"""Is the behind-ground effect size a MEASUREMENT, or a restatement of the leak parameter?

Behind the ground, current only reaches the perturbed region because the background leak
gives it somewhere to go once it is past a precursor. That makes the leak conductance
load-bearing in a way it is not for a front-of-ground perturbation, where the leak is a
regulariser carrying ~1e-5 of the injected current.

If a behind-ground effect scales linearly with ``leak`` while a front-of-ground effect does
not, then the behind-ground number is parameterised rather than measured: halving an
arbitrary tuning constant halves the finding. That is a different objection from "it is
small" and it is not answered by the noise floor.

Sweeps ``leak`` at fixed ``port``. Env: ``mamba run -n scadc-metabolic-model``.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
_LIB = Path(os.environ.get(
    "FABFOS_ENGINE_LIB", _HERE.parents[1] / "src/metasmith_libraries/resources/lib"))
sys.path.insert(0, str(_LIB))

import ecspr_build as B                                             # noqa: E402
import ecspr_probe as P                                             # noqa: E402
import ground as G                                                  # noqa: E402
from horizon import PERTURBATIONS, PROBES                           # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--element", default="N")
    ap.add_argument("--port-mult", type=float, default=1e3)
    ap.add_argument("--leak-mults", type=float, nargs="*",
                    default=[1e-9, 1e-8, 1e-7, 1e-6, 1e-5])
    ap.add_argument("--genes", nargs="*",
                    default=["psd__BEHIND_serine", "lpxCD__BEHIND_set4",
                             "serC__CTL_upstream", "CTL_unreachable"])
    a = ap.parse_args(argv)
    ref, X = Path(a.ref), a.element

    axes = json.loads((ref / "derived/axes/biomass_dag_axes_set4.json").read_text())
    PROBES["nh4_to_set4"]["ground"] = list(P.axis_node_sets(axes, X)["precursors"])
    ratios = B.load_direction_ratios(ref / "derived/mnxref-4_5/direction.parquet")
    weights = B.load_evidence_weights(ref / "derived/evidence/evidence_weights.parquet")
    pairs = B.load_pairs(ref / "derived/mnxref-4_5/atom_pairs.parquet", X)
    g0 = B.graph_from_pairs(pairs, X, weights, ratios)
    med = float(np.median(g0.gp))
    port = a.port_mult * med

    # base = the network LACKING the reaction, so the delta is an insertion
    lacking = {g: B.graph_from_pairs(
        pairs, X, {k: v for k, v in weights.items() if k not in PERTURBATIONS[g]}, ratios)
        for g in a.genes}

    for pname, pr in PROBES.items():
        src, gnd = pr["source"], pr["ground"]
        sets = dict(precursors=list(gnd))
        print(f"\n== {pname}  port={a.port_mult:g}x  (leak_raw; insertion deltas)")
        print(f"   {'leak':>8s} {'leak_frac':>11s} " +
              " ".join(f"{g.split('__')[0]:>13s}" for g in a.genes) +
              "   scaling vs leak")
        prev = {}
        for lm in a.leak_mults:
            leak = lm * med
            full = G.measure(g0, src, sets, "leak_raw", port=port, leak=leak)
            rels = {}
            for g in a.genes:
                t0 = G.measure(lacking[g], src, sets, "leak_raw",
                               port=port, leak=leak)["total"]
                rels[g] = (full["total"] - t0) / t0 if t0 else float("nan")
            # slope in log-log against the previous leak decade: 1.0 == proportional to leak
            slope = []
            for g in a.genes:
                if g in prev and prev[g] and rels[g] and abs(rels[g]) > 1e-300:
                    slope.append(f"{g.split('__')[0]}={np.log10(abs(rels[g]/prev[g])) / np.log10(lm/prev['_lm']):.2f}")
            print(f"   {lm:8.0e} {full['leak_frac']:11.3e} " +
                  " ".join(f"{rels[g]:13.4e}" for g in a.genes) +
                  ("   " + " ".join(slope) if slope else ""))
            prev = dict(rels)
            prev["_lm"] = lm
    return 0


if __name__ == "__main__":
    sys.exit(main())
