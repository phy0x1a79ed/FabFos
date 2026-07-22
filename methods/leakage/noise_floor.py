"""What is the NUMERICAL floor of a probe delta -- and does the behind-ground signal clear it?

The horizon run reports behind-ground inserts at 1e-11 to 1e-14 relative against a bit-exact
zero control. A bit-exact control proves the *hard* ground is blind; it does NOT establish
that the leak's 1e-11 is a real number, because both grounds' controls are exact for the
same trivial reason (an identical matrix gives an identical factorisation). The question is
what happens when the matrix is mathematically identical but numerically different.

So: rebuild the graph from a ROW-PERMUTED atom-pair table. Same multiset of edges, same
conductances, same answer in exact arithmetic -- but a different node ordering, hence a
different sparse LU pivot order, hence different round-off. The spread of ``total`` across
permutations is the floor a delta has to clear to mean anything.

This is the control the previous run was missing: the null rebuild reuses the same ordering
and so cannot see this. A behind-ground effect below this spread is arithmetic, not biology.

Env: ``mamba run -n scadc-metabolic-model``.
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
from horizon import PERTURBATIONS, PROBES, GLUTAMATE, NH4, SERINE   # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--element", default="N")
    ap.add_argument("--nperm", type=int, default=6)
    ap.add_argument("--port-mults", type=float, nargs="*", default=[1.0, 100.0, 1000.0])
    ap.add_argument("--leak-mult", type=float, default=1e-7)
    ap.add_argument("--genes", nargs="*",
                    default=["psd__BEHIND_serine", "lpxCD__BEHIND_set4",
                             "serC__CTL_upstream"])
    a = ap.parse_args(argv)
    ref, X = Path(a.ref), a.element

    axes = json.loads((ref / "derived/axes/biomass_dag_axes_set4.json").read_text())
    PROBES["nh4_to_set4"]["ground"] = list(P.axis_node_sets(axes, X)["precursors"])
    ratios = B.load_direction_ratios(ref / "derived/mnxref-4_5/direction.parquet")
    weights = B.load_evidence_weights(ref / "derived/evidence/evidence_weights.parquet")
    pairs = B.load_pairs(ref / "derived/mnxref-4_5/atom_pairs.parquet", X)
    g0 = B.graph_from_pairs(pairs, X, weights, ratios)
    med = float(np.median(g0.gp))
    leak = a.leak_mult * med

    # the permuted rebuilds -- mathematically the same network, numerically a different one
    perms = [g0] + [B.graph_from_pairs(pairs.sample(frac=1.0, random_state=k), X,
                                       weights, ratios) for k in range(a.nperm)]
    for g in perms[1:]:
        assert g.m == g0.m and g.n == g0.n, "permutation changed the graph, not its order"
        assert abs(float(g.gp.sum()) - float(g0.gp.sum())) <= 1e-9 * float(g0.gp.sum())
    print(f"[{X}] {len(perms)} orderings of a {g0.n:,}n/{g0.m:,}e graph\n", flush=True)

    # the behind-ground perturbations, as INSERTS: base is the network lacking the reaction
    lof = {}
    for gene in a.genes:
        w = {k: v for k, v in weights.items() if k not in PERTURBATIONS[gene]}
        lof[gene] = [B.graph_from_pairs(p_.sample(frac=1.0, random_state=k) if k >= 0 else p_,
                                        X, w, ratios)
                     for k, p_ in [(-1, pairs)] + [(k, pairs) for k in range(a.nperm)]]

    for pname, pr in PROBES.items():
        src, gnd = pr["source"], pr["ground"]
        sets = dict(precursors=list(gnd))
        for mode in ("hard", "leak_raw"):
            for pm in ([1.0] if mode == "hard" else a.port_mults):
                port = pm * med
                tot = np.array([G.measure(g, src, sets, mode, port=port, leak=leak)["total"]
                                for g in perms], float)
                spread = float(tot.max() - tot.min()) / float(tot[0])
                print(f"== {pname} / {mode} port={pm:g}x")
                print(f"   total = {tot[0]:.12g}   REORDERING FLOOR = "
                      f"{spread:.3e} relative  (n={len(tot)})")
                for gene in a.genes:
                    d = np.array([(tot[i] - G.measure(lof[gene][i], src, sets, mode,
                                                      port=port, leak=leak)["total"])
                                  / G.measure(lof[gene][i], src, sets, mode, port=port,
                                              leak=leak)["total"]
                                  for i in range(len(perms))], float)
                    med_d = float(np.median(np.abs(d)))
                    verdict = ("BELOW FLOOR -- arithmetic" if med_d <= spread else
                               f"{med_d / spread:.3g}x floor" if spread > 0 else "floor=0")
                    print(f"     insert {gene:22s} |rel| = {med_d:.4e}   "
                          f"spread over orderings = {float(d.max() - d.min()):.3e}   "
                          f"{verdict}")
                print(flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
