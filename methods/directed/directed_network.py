"""SCADC directed-network driver -- feeds the real directionality ensemble into the
engine's rectified-network solver (``ecspr_directed``).

This is the *experiment-specific* half the engine primitive deliberately does not carry:
it holds the SCADC paths, reads the base reaction graph and MetaNetX reaction roles, and
joins the directionality ensemble's per-reaction ``g_rev/g_fwd = exp(dG'/RT)`` ratio onto
the edges. The engine primitive stays path-free and canon-free; everything that fixes a
*number* enters here and is recorded with its content hash.

Orientation. Each base-graph edge joins a metabolite node ``("met", id)`` to a reaction
node ``("rxn", mnxr)``. MetaNetX ``reac_prop`` gives, per reaction, the substrate set S
(LHS) and product set P (RHS). We orient a substrate edge tail=met -> head=rxn and a
product edge tail=rxn -> head=met, so a *forward*-flowing reaction carries current
substrate -> rxn -> product. That is the same convention in which the ensemble's ``dG'``
(and hence ``ratio``) is written, so no sign bookkeeping is needed: for a forward-favoured
reaction ``dG' < 0`` gives ``ratio < 1`` and the backward branch is throttled.

Directionality. ``gp = w`` (the element evidence weight, unchanged from the undirected
solve) on every oriented edge; ``gm = ratio * w`` on both the reaction's substrate and
product edges, because reversing flow through the reaction reverses both. A reaction with
no direction evidence has ``ratio == 1.0`` -> ``gm == gp`` -> a symmetric edge -> a
*provable no-op*: the directed solve reproduces the undirected ``R_eff`` exactly there.
Role-unknown edges (the metabolite is in neither S nor P, e.g. a cofactor the roles table
does not resolve) are left symmetric for the same reason -- we cannot orient them.

The engine's backward floor keeps ``gm > 0`` so the grounded potential is unique even for
a reaction the ensemble clamps to ``ratio ~ 3e-18``.
"""
from __future__ import annotations

import hashlib
import pickle
import re
import sys
from pathlib import Path

import numpy as np
import networkx as nx
import pandas as pd

HERE = Path(__file__).resolve().parent
from fabfos import canon                                              # noqa: E402
sys.path.insert(0, str(canon.ENGINE_LIB / "resources" / "lib"))
from ecspr_directed import (build_incidence, directed_ceff, orient_and_weight,  # noqa: E402
                            DIODE_BACKWARD_FLOOR)
from ecspr_solver import _reff_dense                                            # noqa: E402

# Base reaction graphs (met/rxn bipartite, per element X, edge weight key ``w_{X}``).
# These are staged inputs in the DAG; the absolute paths live here in the driver, never
# in the engine.
NET_CACHE = {
    "netA": Path("/home/tony/agentic_workspace/projects/scadc/metabolic-modelling/main/"
                 "metabolic-modelling/04_reaction_network/cache/netA_iECDH10B"),
    "netB": Path("/home/tony/agentic_workspace/projects/scadc/metabolic-modelling/main/"
                 "metabolic-modelling/04_reaction_network/cache"),
}
DIRECTION_ANNOTATION = Path(canon.DATA) / "direction" / "direction_annotation.parquet"
METANETX_REAC_PROP = Path(canon.DATA) / "references" / "metanetx" / "reac_prop.tsv"

_EQ_TERM = re.compile(r"(\d+(?:\.\d+)?)\s+(MNXM\w+)@\w+")


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_roles(reac_prop: Path = METANETX_REAC_PROP) -> dict:
    """MNXR -> (substrate set, product set) parsed from MetaNetX reac_prop equations."""
    roles = {}
    with open(reac_prop) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            p = ln.rstrip("\n").split("\t")
            if len(p) < 2 or p[0] == "EMPTY" or "=" not in p[1]:
                continue
            lhs, rhs = p[1].split("=", 1)
            S = {m for _, m in _EQ_TERM.findall(lhs)}
            P = {m for _, m in _EQ_TERM.findall(rhs)}
            if S and P:
                roles[p[0]] = (S, P)
    return roles


def load_ratios(path: Path = DIRECTION_ANNOTATION) -> dict:
    """MNXR -> g_rev/g_fwd ratio from the directionality ensemble. Missing = no evidence."""
    df = pd.read_parquet(path)
    return dict(zip(df.mnxr.astype(str), df.ratio.astype(float)))


class DirectedNet:
    """An oriented base reaction graph ready for the rectified solve.

    Build once per (network, element); solve any number of (source, sink) pairs. The
    incidence and forward weights are fixed; only the terminals change per axis. Ratios
    are baked into ``gm`` at build time.
    """

    def __init__(self, tag: str, X: str, roles: dict, ratios: dict,
                 floor: float = DIODE_BACKWARD_FLOOR, force_noop: bool = False):
        G = pickle.load(open(NET_CACHE[tag] / f"base_{X}.pkl", "rb"))
        G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
        self.tag = tag
        self._build(G, X, roles, ratios, floor, force_noop)

    @classmethod
    def from_graph(cls, G, X: str, roles: dict, ratios: dict,
                   floor: float = DIODE_BACKWARD_FLOOR, force_noop: bool = False,
                   take_lcc: bool = True):
        """Build directly from an in-memory graph (e.g. a perturbed base from a
        knockout/add-back), bypassing the on-disk cache. Used by the benchmark bridge,
        which perturbs the base graph per condition and re-solves."""
        self = cls.__new__(cls)
        self.tag = None
        if take_lcc and G.number_of_nodes():
            G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
        self._build(G, X, roles, ratios, floor, force_noop)
        return self

    def _build(self, G, X, roles, ratios, floor, force_noop):
        # The orientation logic now lives in the ENGINE (ecspr_directed.orient_and_weight),
        # path-free and canon-free; this driver only supplies the graph object and the
        # roles/ratios dicts. Keeping the DirectedNet field surface (B/gp/gm/idx/nodes/
        # stats) so the ceff/reff_undirected methods and the benchmark bridge are unchanged.
        self.X, self.G = X, G
        onet = orient_and_weight(G, X, roles, ratios, floor=floor, force_noop=force_noop)
        self._onet = onet
        self.nodes, self.idx = onet.nodes, onet.idx
        self.B, self.gp, self.gm, self.floor = onet.B, onet.gp, onet.gm, onet.floor
        self.stats = onet.stats

    def ceff(self, s_node, t_node, tol=None) -> float:
        """Directed two-terminal effective conductance between two metabolite nodes.

        ``tol`` overrides the engine's Newton convergence tolerance. The default (None ->
        engine 1e-10) is exact; a looser tol (e.g. 1e-6) is safe when the result feeds a
        RANK-only metric, and cuts the semismooth-Newton iteration count materially."""
        kw = dict(floor=self.floor)
        if tol is not None:
            kw["tol"] = tol
        return float(directed_ceff(self.B, self.gp, self.gm,
                                   self.idx[s_node], self.idx[t_node], **kw))

    def reff_undirected(self, s_node, t_node) -> float:
        """Undirected R_eff on the same graph (the incumbent's quantity), for parity."""
        return _reff_dense(self.G, s_node, t_node, wk=f"w_{self.X}")


def _prove_noop(tag: str, X: str, roles: dict, ratios: dict, n_pairs: int = 20) -> float:
    """No-evidence no-op: with every ratio forced to 1.0 the directed C_eff must equal
    the undirected 1/R_eff to solver tolerance -- the directed model degrades to the
    undirected one exactly where direction is unknown."""
    net = DirectedNet(tag, X, roles, ratios, force_noop=True)
    mets = [n for n in net.nodes if n[0] == "met"]
    rng = np.random.default_rng(7)
    worst = 0.0
    tested = 0
    for _ in range(n_pairs * 3):
        s, t = rng.choice(len(mets), 2, replace=False)
        s_node, t_node = mets[s], mets[t]
        r_undir = net.reff_undirected(s_node, t_node)
        if not np.isfinite(r_undir) or r_undir <= 0:
            continue
        c_dir = net.ceff(s_node, t_node)
        err = abs(c_dir - 1.0 / r_undir)
        worst = max(worst, err)
        tested += 1
        if tested >= n_pairs:
            break
    print(f"  no-op {tag}/{X}: forced ratio==1.0 -> |C_dir - 1/R_undir| worst={worst:.2e} "
          f"over {tested} pairs")
    return worst


def main():
    import argparse
    ap = argparse.ArgumentParser(description="directed-network driver: coverage + no-op proof")
    ap.add_argument("--elements", default="C")
    ap.add_argument("--nets", default="netA,netB")
    args = ap.parse_args()

    print(f"direction_annotation: {DIRECTION_ANNOTATION}")
    print(f"  sha256 = {sha256_of(DIRECTION_ANNOTATION)}")
    roles = load_roles()
    ratios = load_ratios()
    print(f"roles: {len(roles):,} reactions; ratios: {len(ratios):,} reactions\n")

    worst_noop = 0.0
    for tag in args.nets.split(","):
        for X in args.elements.split(","):
            net = DirectedNet(tag, X, roles, ratios)
            s = net.stats
            frac_dir = s["n_directed_edges"] / s["n_edges"]
            print(f"{tag}/{X}: edges={s['n_edges']} directed={s['n_directed_edges']} "
                  f"({frac_dir:.1%}) role_unknown={s['n_role_unknown']} "
                  f"no_evidence={s['n_no_evidence_edges']}")
            worst_noop = max(worst_noop, _prove_noop(tag, X, roles, ratios))
    print(f"\nworst no-op residual across nets: {worst_noop:.2e}")
    return 0 if worst_noop < 1e-6 else 1


if __name__ == "__main__":
    raise SystemExit(main())
