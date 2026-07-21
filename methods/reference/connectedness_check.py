"""Connectedness sanity check on the frozen reference atom graph.

This is NOT the validation suite (specificity AUC / GOF-LOF -- that measures whether
the graph RANKS biology correctly, and another run owns it). This asks the prior
question: is the graph TOPOLOGICALLY SENSIBLE at all? Two independent tests, neither
of which reads any experiment outcome:

  TEST 1 -- E. coli model connectedness (the user's test).
    "For our model E. coli, all 4 elements should be connected to all the metabolites
     that carry them." Induce the reference element-X graph on the curated GEM reactome
     (the same crosswalk the validation netA build uses) and check: of the GEM
     metabolites whose CHEMICAL FORMULA contains element X, what fraction land in the
     one connected component? Do it for the HONEST reference AND the incumbent
     (fabricated-transit) graph, and report the delta -- the decisive question is
     whether refusing ~27k fabricated transits DISCONNECTED any real metabolite, or
     only pruned edges that were never chemistry.

  TEST 2 -- edge/formula consistency (a universe-wide bug check).
    Every w_X>0 edge asserts "an atom of element X transits this reaction through this
    metabolite." That metabolite MUST contain element X in its formula. A carbon-transit
    edge through a carbon-free metabolite is a mapping bug. Count violations over the
    whole reference graph. Sensible => zero (modulo formula-less nodes).

Isolated element-carrying metabolites are then attributed against the closure ledger:
a metabolite is allowed to be off-graph only for an adjudicated reason (its reactions
are pseudo/transport, or refused as structureless), never silently.

Env: p312 (the graphs are numpy-2.x pickled). Reads only frozen inputs + the incumbent
universe graph; writes a report + tables under REF/sanity/. Read-only on all graphs.
"""
from __future__ import annotations

import argparse
import pickle
import re
import sys
from collections import defaultdict
from pathlib import Path

import networkx as nx
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import reference as ref  # noqa: E402

ELEMENTS = list(ref.BASE_GRAPH_ELEMENTS)  # C, N, S, P
INCUMBENT_UNIVERSE = {  # the fabricated-transit graph the reference supersedes
    X: ref.BASE_GRAPH_DIR / f"mnx_bipartite_{X}.pkl" for X in ELEMENTS
}
CROSSWALK = {
    "iML1515": ref.BASE_GRAPH_DIR / "crosswalk" / "gem_rxn_to_mnxr_iML1515.parquet",
    "iECDH10B": ref.BASE_GRAPH_DIR / "crosswalk" / "gem_rxn_to_mnxr_iECDH10B.parquet",
}
OUT = ref.REF_ROOT / "sanity"

_ELEM_TOKEN = re.compile(r"([A-Z][a-z]?)")
_SMI_BRACKET = re.compile(r"\[([A-Za-z][a-z]?)")
_SMI_ORGANIC = re.compile(r"(Cl|Br|[BCNOPSFIcnops])")


def formula_elements(formula: str) -> set[str]:
    """The element symbols present in a MetaNetX formula. Tokenises on real element
    boundaries so 'Cl'/'Na' never masquerade as C/N, and polymer 'R'/'*' are ignored."""
    if not formula:
        return set()
    return set(_ELEM_TOKEN.findall(formula))


def smiles_elements(smiles: str) -> set[str]:
    """Element symbols present in a SMILES: bracket atoms + the organic-subset shorthand
    (aromatic lower-case folded to its element). Used as a fallback ground truth when the
    chem_prop FORMULA column is blank (MetaNetX ships many SMILES with no formula)."""
    if not smiles:
        return set()
    toks = _SMI_BRACKET.findall(smiles) + _SMI_ORGANIC.findall(re.sub(r"\[[^\]]*\]", "", smiles))
    return {t.capitalize() for t in toks}


# =====================================================================
# Inputs
# =====================================================================
def load_reac_participants(reac_prop: Path) -> dict[str, set[str]]:
    """mnxr -> set of participant MNXMs (both sides), parsed from the equation column.
    Coefficients/compartments stripped; used to enumerate a GEM's metabolite set."""
    mnxm = re.compile(r"(MNXM\d+|MNXM[A-Za-z0-9]+|WATER|BIOMASS)")
    out: dict[str, set[str]] = {}
    with open(reac_prop) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 2 or not p[0].startswith("MNXR"):
                continue
            out[p[0]] = set(mnxm.findall(p[1]))
    return out


def load_structure(chem_prop: Path, want: set[str]) -> dict[str, tuple[str, str]]:
    """mnxm -> (formula, smiles), only for the metabolites we need (chem_prop is ~0.8 GB).
    Both are kept because MetaNetX often populates one and not the other; 'carries element
    X' is true if EITHER attests X, and only false when a populated field lacks X."""
    out: dict[str, tuple[str, str]] = {}
    with open(chem_prop) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if not p or p[0] not in want:
                continue
            out[p[0]] = (p[3] if len(p) > 3 else "", p[8] if len(p) > 8 else "")
            if len(out) == len(want):
                break
    return out


def carrier_verdict(struct: tuple[str, str] | None, X: str) -> str:
    """'yes' / 'no' / 'unknown' -- does the metabolite carry element X? 'unknown' when
    chem_prop gives neither a formula nor a SMILES (structureless in this release: the
    ensemble may still have mapped it from an external structure, so a transit edge on it
    is not a contradiction, merely un-checkable here)."""
    if struct is None:
        return "unknown"
    formula, smiles = struct
    if not formula and not smiles:
        return "unknown"
    return "yes" if (X in formula_elements(formula) or X in smiles_elements(smiles)) else "no"


def load_ledger() -> pd.DataFrame:
    return pd.read_parquet(ref.LEDGER)


# =====================================================================
# Graph induction (mirrors validation/60_build_netA_base.build_base_graph)
# =====================================================================
def induce(G: nx.Graph, mnxrs: set[str], X: str) -> nx.Graph:
    """GEM-induced element-X base: keep only w_X>0 rxn<->met edges. E=1.0 (topology only)."""
    wk = f"w_{X}"
    base = nx.Graph()
    for r in mnxrs:
        rn = ("rxn", r)
        if rn not in G:
            continue
        for _, m, d in G.edges(rn, data=True):
            if d.get(wk, 0.0) > 0:
                base.add_edge(rn, m, **{wk: float(d[wk])})
    return base


def met_sets(base: nx.Graph) -> tuple[set[str], set[str]]:
    """(all connected metabolites, metabolites in the largest connected component)."""
    conn = {n[1] for n in base.nodes if n[0] == "met"}
    if not base.number_of_nodes():
        return conn, set()
    lcc = max(nx.connected_components(base), key=len)
    return conn, {n[1] for n in lcc if n[0] == "met"}


# =====================================================================
# Main
# =====================================================================
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--model", choices=["iECDH10B", "iML1515"], default="iECDH10B",
                    help="the curated GEM whose reactome defines 'our model E. coli'")
    ap.add_argument("--out", type=Path, default=OUT)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    # the GEM reactome, in MetaNetX space (same crosswalk the validation netA build uses)
    xw = pd.read_parquet(CROSSWALK[a.model])
    mnxrs = set(xw.mnxr.dropna().astype(str).unique())
    print(f"[check] {a.model}: {len(xw)} GEM rxn rows -> {len(mnxrs)} unique MNXR", flush=True)

    participants = load_reac_participants(ref.REAC_PROP)
    gem_mets: set[str] = set()
    for r in mnxrs:
        gem_mets |= participants.get(r, set())
    # also need structures for every metabolite that carries a transit edge (Test 2)
    edge_mets: set[str] = set()
    for X in ELEMENTS:
        G = pickle.load(open(ref.GRAPH_DIR / f"mnx_bipartite_{X}.pkl", "rb"))
        for u, v, d in G.edges(data=True):
            if d.get(f"w_{X}", 0.0) > 0:
                edge_mets.add(v[1] if v[0] == "met" else u[1])
    struct = load_structure(ref.CHEM_PROP, gem_mets | edge_mets)
    with_struct = sum(1 for m in gem_mets if struct.get(m, ("", ""))[0] or struct.get(m, ("", ""))[1])
    print(f"[check] GEM touches {len(gem_mets)} metabolites; "
          f"{with_struct} have a chem_prop formula or SMILES", flush=True)

    ledger = load_ledger()
    reason_of = dict(zip(ledger.mnxr.astype(str), ledger.aam_reason.fillna("")))
    state_of = dict(zip(ledger.mnxr.astype(str), ledger.aam_state.fillna("")))

    rows = []
    edge_violations = {}
    isolated_detail = defaultdict(list)
    for X in ELEMENTS:
        Gref = pickle.load(open(ref.GRAPH_DIR / f"mnx_bipartite_{X}.pkl", "rb"))
        Ginc = pickle.load(open(INCUMBENT_UNIVERSE[X], "rb"))

        # metabolites the GEM carries element X for, by formula-or-SMILES (ground truth)
        carry = {m for m in gem_mets if carrier_verdict(struct.get(m), X) == "yes"}

        base_ref = induce(Gref, mnxrs, X)
        base_inc = induce(Ginc, mnxrs, X)
        conn_ref, lcc_ref = met_sets(base_ref)
        conn_inc, lcc_inc = met_sets(base_inc)

        carry_conn_ref = carry & conn_ref
        carry_lcc_ref = carry & lcc_ref
        # metabolites the incumbent connected but the honest graph drops
        lost = (carry & conn_inc) - conn_ref
        # metabolites that carry X yet sit off the reference graph entirely
        isolated = carry - conn_ref

        rows.append(dict(
            element=X,
            gem_carry_X=len(carry),
            connected_ref=len(carry_conn_ref),
            in_LCC_ref=len(carry_lcc_ref),
            pct_LCC_ref=round(100 * len(carry_lcc_ref) / max(1, len(carry)), 2),
            connected_incumbent=len(carry & conn_inc),
            lost_vs_incumbent=len(lost),
            isolated_carry_X=len(isolated),
        ))

        # attribute every lost/isolated metabolite to a ledger verdict on its GEM reactions
        for m in sorted(isolated):
            rxn_states = []
            for r in mnxrs:
                if m in participants.get(r, set()):
                    rxn_states.append((r, state_of.get(r, "?"), reason_of.get(r, "")))
            # dominant refusal reason among this metabolite's GEM reactions
            reasons = [rs[2] for rs in rxn_states if rs[1] == "refused" and rs[2]]
            top = max(set(reasons), key=reasons.count) if reasons else "(all resolved/diluted -- currency-only, no cross-transit)"
            isolated_detail[X].append(dict(
                mnxm=m, formula=struct.get(m, ("", ""))[0], n_gem_rxn=len(rxn_states),
                dropped_by_incumbent=m in lost, dominant_reason=top))

        # TEST 2: edge/formula consistency over the WHOLE reference graph for X.
        # A CONTRADICTION = an edge on a metabolite whose populated formula/SMILES lacks X
        # (a fabricated element). 'unknown' = structureless in chem_prop, un-checkable here.
        wk = f"w_{X}"
        contradiction = 0
        checkable = 0
        unknown = 0
        for u, v, d in Gref.edges(data=True):
            if d.get(wk, 0.0) <= 0:
                continue
            met = v[1] if v[0] == "met" else (u[1] if u[0] == "met" else None)
            if met is None:
                continue
            verd = carrier_verdict(struct.get(met), X)
            if verd == "unknown":
                unknown += 1
            else:
                checkable += 1
                if verd == "no":
                    contradiction += 1
        edge_violations[X] = (contradiction, checkable, unknown)
        print(f"[check] {X}: carry={len(carry)} connected_ref={len(carry_conn_ref)} "
              f"LCC={len(carry_lcc_ref)} lost_vs_incumbent={len(lost)} "
              f"edge_contradiction={contradiction}/{checkable} (unknown {unknown})", flush=True)

    summ = pd.DataFrame(rows)
    summ.to_csv(a.out / f"connectedness_{a.model}.tsv", sep="\t", index=False)
    iso_rows = [dict(element=X, **r) for X, rs in isolated_detail.items() for r in rs]
    pd.DataFrame(iso_rows).to_csv(a.out / f"isolated_metabolites_{a.model}.tsv",
                                  sep="\t", index=False)

    # markdown report
    md = [f"# Reference graph connectedness sanity check -- {a.model}", "",
          f"Reference: MetaNetX {ref.MNXREF_VERSION} ({ref.MNXREF_DATE}), "
          f"honest atom graph at `REF/graph/`.  GEM reactome: {len(mnxrs)} MNXR.", "",
          "## TEST 1 -- element connectedness of the model's metabolites", "",
          "'connected' = has >=1 element-X transit edge; 'in LCC' = in the single "
          "largest connected component.  Ground truth for 'carries X' = the metabolite's "
          "chem_prop formula contains element X.", "",
          "| element | GEM mets carrying X | connected (honest) | in LCC (honest) | "
          "% in LCC | connected (incumbent) | lost vs incumbent | isolated |",
          "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['element']} | {r['gem_carry_X']} | {r['connected_ref']} | "
                  f"{r['in_LCC_ref']} | {r['pct_LCC_ref']} | {r['connected_incumbent']} | "
                  f"{r['lost_vs_incumbent']} | {r['isolated_carry_X']} |")
    md += ["", "**lost vs incumbent** = metabolites the fabricated-transit graph connected "
           "that the honest graph drops.  If these are all structureless stubs / "
           "non-molecules (see the isolated table), the refusals pruned non-chemistry, "
           "not real metabolism.", "",
           "## TEST 2 -- edge / structure consistency (whole reference graph)", "",
           "Every w_X>0 edge must touch a metabolite whose formula OR SMILES contains X.  "
           "A contradiction = a fabricated element (transit of X through a metabolite that "
           "provably lacks X).  'unknown' = the metabolite is structureless in chem_prop "
           "(no formula, no SMILES) -- the ensemble mapped it from an external structure, so "
           "the edge is un-checkable here, not wrong.", "",
           "| element | contradictions | edges checkable | unknown (structureless) |",
           "|---|---|---|---|"]
    for X in ELEMENTS:
        c, ck, un = edge_violations[X]
        md.append(f"| {X} | {c} | {ck} | {un} |")
    md += ["", f"Isolated-metabolite attribution: `isolated_metabolites_{a.model}.tsv` "
           "(each off-graph carrier tied to the ledger verdict on its GEM reactions).", ""]
    (a.out / f"REPORT_{a.model}.md").write_text("\n".join(md) + "\n")
    print(f"\n[check] wrote {a.out}/REPORT_{a.model}.md (+ 2 tsv)")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
