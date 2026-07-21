"""Toy fixtures for the ECSPr pulse-chase suite -- hand-checked micro-universes.

Each fixture is a small set of reactions with hand-assigned element atom counts and
explicit atom-transfer bijections. From one reaction spec we build BOTH graphs over
the SAME reactions:

  * the ATOM graph -- via the engine's real `ecspr_atom_graph.atom_edges`, nodes
    (metabolite, canonical-atom-rank), edges = atom transfers.
  * the STAR graph -- via the engine's real `ecspr_atom_graph.star_edges` on a
    bipartite ("rxn", r)-("met", m) graph whose edge weight w_X = the metabolite's
    element atom count.

Building both from one spec is the only way a star-vs-atom difference is attributable
to topology rather than to a different reaction set (the engine's own head-to-head
makes the same argument). The suite reads the PAIRS TABLE directly (pre-`atom_edges`)
for the conservation/emission checks, because `atom_edges` drops `u==v` self-mapped
atoms and would break the identity if read post-collapse.

A reaction spec:
    {"mnxr": "R1",
     "sub":  {"SUB": 1, "COFR": 10},     # metabolite -> element-X atom count
     "prod": {"PRD": 1, "COFO": 10},
     "xfer": [("SUB", "PRD", 1), ("COFR", "COFO", 10)]}   # (from, to, n_atoms)

`xfer` names which substrate atoms land in which product; ranks 0..n-1 on each side
(a toy metabolite's atoms are ranks 0..count-1). Everything downstream is the engine's
own code.

No SCADC paths, no canon values. rdkit is used only by the ambiguous-emission fixture.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import networkx as nx

# The engine library (canon.ENGINE_LIB/resources/lib). The suite's runner puts it on
# sys.path; when toy.py is imported standalone we add it here too so the fixtures are
# usable in a REPL. The lib never imports back.
_LIB = Path("/home/tony/agentic_workspace/projects/metasmith-libraries/fabfos"
            "/resources/lib")
if str(_LIB) not in sys.path:
    sys.path.insert(0, str(_LIB))

from ecspr_atom_graph import atom_edges, star_edges, Grid  # noqa: E402
from ecspr_atom_pairs import PAIR_COLS  # noqa: E402


# =====================================================================
# spec -> tables and graphs
# =====================================================================

def _xfer_ranks(entry):
    """A transfer entry is (sm, pm, k) -> ranks 0..k-1 on each side, or
    (sm, pm, sub_ranks, prod_ranks) with explicit (possibly non-zero-based) ranks so a
    multi-atom metabolite can route different atoms to different destinations."""
    if len(entry) == 3:
        sm, pm, k = entry
        return sm, pm, list(range(k)), list(range(k))
    sm, pm, sub_ranks, prod_ranks = entry
    assert len(sub_ranks) == len(prod_ranks), (sm, pm, sub_ranks, prod_ranks)
    return sm, pm, list(sub_ranks), list(prod_ranks)


def pairs_table(reactions: list, element: str) -> pd.DataFrame:
    """The atom-transfer PAIRS table, exactly the schema `atom_edges` consumes and the
    real extractor emits. One row per (reaction, sub metabolite, prod metabolite)."""
    rows = []
    for rx in reactions:
        for entry in rx["xfer"]:
            sm, pm, si, pi = _xfer_ranks(entry)
            rows.append(dict(
                mnxr=rx["mnxr"], element=element, substrate=sm, product=pm,
                n_atoms=len(si),
                sub_idx=",".join(str(i) for i in si),
                prod_idx=",".join(str(j) for j in pi),
                # toy transfers are confident bijections -> unit dilution weight
                pair_w=",".join("1.0" for _ in si),
            ))
    return pd.DataFrame(rows, columns=list(PAIR_COLS))


def weights_unit(reactions: list) -> dict:
    """weights=1.0 for every reaction -- the build the conservation identities need."""
    return {rx["mnxr"]: 1.0 for rx in reactions}


def bipartite(reactions: list, element: str) -> nx.Graph:
    """The star's element-X bipartite universe: ("rxn", r)-("met", m), w_X = atom count."""
    wk = f"w_{element}"
    G = nx.Graph()
    for rx in reactions:
        rn = ("rxn", rx["mnxr"])
        for side in ("sub", "prod"):
            for m, cnt in rx[side].items():
                if cnt <= 0:
                    continue
                mn = ("met", m)
                # a metabolite on both sides keeps the larger atom count (the engine's
                # metabolite_atoms does the same max over incident edges)
                w = max(cnt, G[rn][mn][wk]) if G.has_edge(rn, mn) else cnt
                G.add_edge(rn, mn, **{wk: float(w)})
    return G


def atom_nx(reactions: list, element: str, weights: dict | None = None) -> nx.Graph:
    """The atom graph as a plain nx.Graph (ALL components, not just the LCC) so
    component membership can be queried. Nodes (metabolite, rank); edge attr 'w'."""
    weights = weights_unit(reactions) if weights is None else weights
    pairs = pairs_table(reactions, element)
    G = nx.Graph()
    for (u, v), w in atom_edges(pairs, weights).items():
        G.add_edge(u, v, w=w)
    return G


def atom_grid(reactions: list, element: str, weights: dict | None = None,
              name: str = "atom") -> Grid:
    weights = weights_unit(reactions) if weights is None else weights
    return Grid(atom_edges(pairs_table(reactions, element), weights), name=name)


def star_grid(reactions: list, element: str, weights: dict | None = None,
              name: str = "star") -> Grid:
    weights = weights_unit(reactions) if weights is None else weights
    G_X = bipartite(reactions, element)
    rxns = [rx["mnxr"] for rx in reactions]
    return Grid(star_edges(G_X, element, rxns, weights), name=name)


# =====================================================================
# Fixtures. Each returns a dict: reactions, element, and named endpoints.
# =====================================================================

def fx_cofactor(element: str = "C") -> dict:
    """Cofactor theft (I4, II7, II6-toy) -- MNXR106432 in miniature:

        ACOA(23 X) + CO2(1 X) + NADPH(21 X)  ->  PYR(3 X) + NADP(21 X) + COA(21 X)

    Carbon transfers: the acetyl group ACOA->PYR (2), the CoA scaffold ACOA->COA (21),
    the added carbon CO2->PYR (1), and the cofactor recycle NADPH->NADP (21, a SEPARATE
    component). NADPH transfers ZERO carbon to PYR. On the atom graph PYR and NADPH are
    in different components; the real channel to PYR is the 1-carbon CO2->PYR. The star
    joins all six metabolites to one hub, so eliminating it makes PYR<->NADPH conductance
    set by NADPH's 21 carbons -- exceeding the real CO2->PYR channel many-fold although
    it transfers nothing. `carrier` (NADPH) is the theft partner; `cargo_src`/`cargo_snk`
    the real channel.
    """
    rx = {"mnxr": "RCOF",
          "sub": {"ACOA": 23, "CO2": 1, "NADPH": 21},
          "prod": {"PYR": 3, "NADP": 21, "COA": 21},
          # acetyl carbons -> pyruvate ranks {0,1}; CoA scaffold -> CoA; CO2 -> pyruvate
          # rank {2}; NADPH -> NADP (cofactor recycle, its own component)
          "xfer": [("ACOA", "PYR", [0, 1], [0, 1]),
                   ("ACOA", "COA", list(range(2, 23)), list(range(21))),
                   ("CO2", "PYR", [0], [2]),
                   ("NADPH", "NADP", list(range(21)), list(range(21)))]}
    return {"reactions": [rx], "element": element,
            "cargo_src": "CO2", "cargo_snk": "PYR", "carrier": "NADPH",
            "theft_a": "PYR", "theft_b": "NADPH"}


def fx_condensation(element: str = "C") -> dict:
    """Input->input leak (I2, II7). Two substrates condense into one product:

        A(2 X) + B(3 X)  ->  P(5 X)

    A's atoms land in P ranks {0,1}; B's atoms land in P ranks {2,3,4}. A and B share
    the metabolite P but NOT a single atom node, so on the atom graph they are in
    different components -- there is no A<->B path. The star joins A and B to one hub,
    manufacturing a co-substrate conductance the chemistry forbids.
    """
    rx = {"mnxr": "RCOND",
          "sub": {"A": 2, "B": 3}, "prod": {"P": 5},
          # A's carbons land in P ranks {0,1}; B's in P ranks {2,3,4} -- DISTINCT product
          # atoms, so A and B share no atom node (different components)
          "xfer": [("A", "P", [0, 1], [0, 1]), ("B", "P", [0, 1, 2], [2, 3, 4])]}
    return {"reactions": [rx], "element": element,
            "in1": "A", "in2": "B", "prod": "P"}


def fx_transport(element: str = "C") -> dict:
    """Compartment collapse / recycle (I5, I8). A transport reaction whose bare
    metabolite multiset is identical on both sides:

        M(2 X)  ->  M(2 X)          plus a real X(1)->Y(1) carrier alongside

    atom_edges drops the M->M self-mapped atoms (u==v), so the transport self-
    annihilates: no cross-metabolite edge from it. The pairs table still records it
    (pre-collapse), which is where the conservation read happens.
    """
    rx = {"mnxr": "RTRANS",
          "sub": {"M": 2, "X": 1}, "prod": {"M": 2, "Y": 1},
          "xfer": [("M", "M", 2), ("X", "Y", 1)]}
    return {"reactions": [rx], "element": element, "recycle": "M"}


def fx_lyase(element: str = "C") -> dict:
    """One substrate split across two products (I8):  S(5 X) -> P1(2 X) + P2(3 X)."""
    rx = {"mnxr": "RLY",
          "sub": {"S": 5}, "prod": {"P1": 2, "P2": 3},
          "xfer": [("S", "P1", 2), ("S", "P2", 3)]}
    return {"reactions": [rx], "element": element}


def fx_multi(element: str = "C") -> dict:
    """Multi-substrate / multi-product fanout (I8): A+B -> C+D, A->C, B->D."""
    rx = {"mnxr": "RMULTI",
          "sub": {"A": 1, "B": 1}, "prod": {"C": 1, "D": 1},
          "xfer": [("A", "C", 1), ("B", "D", 1)]}
    return {"reactions": [rx], "element": element}


def fx_symmetric(element: str = "C") -> dict:
    """Symmetric molecule (I7). A fumarate-like 4-carbon symmetric substrate mapping to
    a succinate-like product; the label can spread over either symmetric route. Ranks
    may collapse under CanonicalRankAtoms, so the conservation read is off the pairs
    table (pre-collapse)."""
    rx = {"mnxr": "RSYM",
          "sub": {"FUM": 4}, "prod": {"SUC": 4},
          "xfer": [("FUM", "SUC", 4)]}
    return {"reactions": [rx], "element": element}


def fx_parallel(element: str = "C") -> dict:
    """Two edge-disjoint, UNEQUAL-resistance routes through DIFFERENT atom pairs (II2, II9).

    SRC carries two labelled X atoms (ranks 0, 1); SNK likewise. The two routes end on
    different sink atoms, so no single (source-atom, sink-atom) pair sees both -- which is
    what makes the MIN-over-pairs scorer unable to combine them:

      route A (3 hops, r=3, the LARGER component so it is the engine's LCC):
          (SRC,0) -> PA1 -> PA2 -> (SNK,0)
      route B (2 hops, r=2, the extra route the engine ignores):
          (SRC,1) -> PB1 -> (SNK,1)

    Deliberately NON-balanced (r=3 vs r=2), so neither is a zero-current Wheatstone arm.
    Under MIN, Ieff(both) == Ieff(A alone) -- adding route B changes nothing (II2 parallel
    additivity fails) and removing route B's edge changes nothing (II2 edge preservation
    fails), and dropping source atom 1 -- which solely carries route B -- does not lower
    Ieff (II9). The super-node all-paths reference combines both in parallel and exceeds
    either alone.
    """
    route_A = [
        {"mnxr": "PA_s_a1", "sub": {"SRC": 2}, "prod": {"PA1": 1}, "xfer": [("SRC", "PA1", [0], [0])]},
        {"mnxr": "PA_a1_a2", "sub": {"PA1": 1}, "prod": {"PA2": 1}, "xfer": [("PA1", "PA2", [0], [0])]},
        {"mnxr": "PA_a2_k", "sub": {"PA2": 1}, "prod": {"SNK": 2}, "xfer": [("PA2", "SNK", [0], [0])]},
    ]
    route_B = [
        {"mnxr": "PB_s_b1", "sub": {"SRC": 2}, "prod": {"PB1": 1}, "xfer": [("SRC", "PB1", [1], [0])]},
        {"mnxr": "PB_b1_k", "sub": {"PB1": 1}, "prod": {"SNK": 2}, "xfer": [("PB1", "SNK", [0], [1])]},
    ]
    return {"reactions": route_A + route_B, "element": element,
            "reactions_A_only": route_A, "src": "SRC", "snk": "SNK",
            "drop_src_rank": 1,                 # the atom that solely carries route B
            "routeB_edge_rxn": "PB_b1_k"}       # an edge unique to route B


def fx_network_pair(element: str = "C") -> dict:
    """Two finished networks B superset A for the comparison contract (II5, II4).

    Everything is ONE connected backbone (as the real host is), so adding reactions
    keeps the main component the LCC and cannot flip which axis is measured -- the trap
    a two-separate-components toy falls into. The backbone C - XM - D is the disjoint
    axis; MEV hangs off XM (so it stays in the main LCC); IPP sits in a SMALL separate
    2-node component that A cannot reach. B adds MEV -> IPPI -> IPP, pulling IPP into the
    main component: MEV->IPP goes from unreachable to finite, while C->D is untouched and
    exactly invariant.

    The suite builds each into a finished graph and calls ECSPr as a pure function on
    both -- no build_ar2m under test.
    """
    backbone = [
        {"mnxr": "BG_c_x", "sub": {"C": 1}, "prod": {"XM": 1}, "xfer": [("C", "XM", 1)]},
        {"mnxr": "BG_x_d", "sub": {"XM": 1}, "prod": {"D": 1}, "xfer": [("XM", "D", 1)]},
        # MEV hangs off the backbone so it is in the main LCC in BOTH networks
        {"mnxr": "BG_x_mev", "sub": {"XM": 1}, "prod": {"MEV": 1}, "xfer": [("XM", "MEV", 1)]},
        # IPP seeded in a small separate component (2 nodes) that A cannot reach from MEV
        {"mnxr": "A_seed_ipp", "sub": {"IPP": 1}, "prod": {"ISEED": 1}, "xfer": [("IPP", "ISEED", 1)]},
    ]
    add = [
        {"mnxr": "MEV_step1", "sub": {"MEV": 1}, "prod": {"IPPI": 1}, "xfer": [("MEV", "IPPI", 1)]},
        {"mnxr": "MEV_step2", "sub": {"IPPI": 1}, "prod": {"IPP": 1}, "xfer": [("IPPI", "IPP", 1)]},
    ]
    return {"reactions_A": backbone, "reactions_B": backbone + add, "element": element,
            "route_src": "MEV", "route_snk": "IPP",
            "disjoint_src": "C", "disjoint_snk": "D"}


def fx_permute(element: str = "C") -> dict:
    """A reaction whose atom identity decides the transfer (II8). One reaction maps four
    substrates to four products, atom-for-atom: A->PA, B->PB, C->PC, D->PD. `permute_pairs`
    reassigns which product each substrate feeds within the reaction, so the measured
    transfer A->PA survives only the 1-in-4 identity permutation. The mean effect
    collapses far below the unshuffled value -- proving the lane measures which atom goes
    where, not merely that the reaction was mappable."""
    rx = {"mnxr": "RPERM",
          "sub": {"A": 1, "B": 1, "C": 1, "D": 1},
          "prod": {"PA": 1, "PB": 1, "PC": 1, "PD": 1},
          "xfer": [("A", "PA", 1), ("B", "PB", 1), ("C", "PC", 1), ("D", "PD", 1)]}
    return {"reactions": [rx], "element": element, "src": "A", "snk": "PA"}


# =====================================================================
# I3 -- the ambiguous-emission fixture (exercises the real refusal path)
# =====================================================================

def ambiguous_forced_inputs():
    """Inputs that put `forced_pairs` on its n>1 refusal branch.

    One substrate S (2 carbons, inequivalent) and one product P (2 carbons,
    inequivalent), carbon conserved. `forced_pairs` emits ONLY n==1 correspondences,
    so this conserved 2-carbon reaction is refused -> zero pairs today. The fix would
    emit M=2 candidate atom correspondences at weight er/2. The candidates have DISTINCT
    canonical ranks (both molecules have two inequivalent carbons), so M is readable.

    Returns (sub_mnxms, prod_mnxms, formulas, ranks_of) ready for `forced_pairs`, plus
    the expected candidate count M.
    """
    from rdkit import Chem
    from ecspr_atom_pairs import canonical_ranks
    # inequivalent carbons: ethanol CCO (C0 methyl, C1 CH2-OH) and acetaldehyde CC=O
    smis = {"S_ETOH": "CCO", "P_ACAL": "CC=O"}
    ranks_of = {}
    for m, smi in smis.items():
        mol = Chem.MolFromSmiles(smi)
        rk = canonical_ranks(mol)
        cidx = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C"]
        ranks_of[(m, "C")] = [rk[i] for i in cidx]
    formulas = {"S_ETOH": "C2H6O", "P_ACAL": "C2H4O"}
    return ["S_ETOH"], ["P_ACAL"], formulas, ranks_of, 2


def ambiguous_mapped_inputs():
    """A mapped reaction whose two product templates are DISTINCT metabolites sharing a
    canonical SMILES -- `pairs_from_mapped` refuses it as `ambiguous_duplicate` today.

    Returns (mapped_smiles, sub_mnxms, prod_mnxms, canon_map).
    """
    from ecspr_atom_pairs import canon_smiles
    # substrate C-C-O mapped; two products both C=O (formaldehyde) but labelled as two
    # different MNXMs -> canonical-SMILES collision -> ambiguous.
    mapped = "[CH3:1][CH2:2][OH:3]>>[CH2:1]=[O:3].[CH2:2]=[O]"
    sub = ["S1"]
    prod = ["P_A", "P_B"]
    canon = {"S1": canon_smiles("CCO"),
             "P_A": canon_smiles("C=O"), "P_B": canon_smiles("C=O")}
    return mapped, sub, prod, canon


if __name__ == "__main__":
    # smoke: every fixture builds both graphs without raising
    for fn in (fx_cofactor, fx_condensation, fx_transport, fx_lyase, fx_multi,
               fx_symmetric, fx_parallel):
        fx = fn()
        rx, el = fx["reactions"], fx["element"]
        a = atom_grid(rx, el)
        try:
            s = star_grid(rx, el)
            sn = s.n
        except ValueError:
            sn = 0
        print(f"{fn.__name__:16s} atom_n={a.n:3d} star_n={sn:3d} "
              f"pairs={len(pairs_table(rx, el))}")
    print("toy smoke ok")
