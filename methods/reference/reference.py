"""The frozen ECSPr atom-mapping + directionality reference, as data + assertions.

WHY THIS EXISTS
---------------
ECSPr scores conductance on an atom graph whose edges are reactions and whose nodes
are ``(metabolite, atom-rank)``. Two properties feed that graph: atom-atom mapping
(AAM -- which substrate atom becomes which product atom) and directionality (the
thermodynamic forward/backward conductance ratio). Both were treated as things a
FabFos run computes.

They are not dynamic. Every reaction in any FabFos network is a MetaNetX reaction id
(``("rxn", mnxr)``); its substrates/products come from MetaNetX ``reac_prop``; the
annotation lanes only ever *project into* an MNXR via ``reac_xref``, never mint a
reaction outside the index; even a GOF/LOF reinforcement copy (``("rxn_reinf", mnxr)``)
is still an existing MNXR. So ``AAM(mnxr)`` and ``direction(mnxr)`` are STATIC functions
of the MetaNetX reaction id. The only experiment-dynamic things are *which* MNXRs are
present, their evidence weight, and reinforcement copies -- never the chemistry.

Because the reaction universe is closed and finite (bounded by a MetaNetX release),
completeness is achievable and verifiable. This module names the frozen, version-pinned,
MNXR-keyed reference asset that any consumer reads as static input, and holds the
schema, the verdict vocabulary, and the assertion helper that refuses a stale or
un-pinned reference loudly -- the same discipline ``canon.py`` applies to the basis.

RELATION TO canon.py
--------------------
``canon.py`` holds the SCADC experiment basis (fosmid count, axes, draw sizes, scorer).
This module holds the METHOD reference (AAM + direction over the whole MNXR universe),
which is upstream of any experiment and shared across all of them. Numbers here are the
universe's, not the basis's; the no-transcribed-numbers rule that governs ``canon.py``
governs this file too -- it is another data module, not prose.

The reference is DECOUPLED from FabFos: builders may live in the ``fabfos`` /
``fabfos-directionality`` / ``fabfos-metabolite-graph`` scopes, but the frozen output is
one shared asset under the data tree (``REF_ROOT``), read by everything and owned by
nothing.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

# =====================================================================
# MetaNetX release pin
# =====================================================================
# The reaction universe is bounded by exactly one MetaNetX release. Pin it: a
# reference built against one release is not interchangeable with another, and a
# consumer that reads the reference must be able to prove which universe it names.
MNXREF_VERSION = "4.5"
MNXREF_DATE = "2025-08-13"

DATA = Path("/home/tony/agentic_workspace/data/scadc")
REFERENCES = DATA / "references" / "metanetx"
REAC_PROP = REFERENCES / "reac_prop.tsv"     # the reaction universe (~83.8k rows)
CHEM_PROP = REFERENCES / "chem_prop.tsv"     # per-metabolite name/formula/SMILES
REAC_XREF = REFERENCES / "reac_xref.tsv"     # lane->MNXR crosswalk

# =====================================================================
# Reference home (version-pinned, decoupled from any FabFos experiment)
# =====================================================================
REF_ROOT = DATA / "ecspr_reference" / f"mnxref-{MNXREF_VERSION.replace('.', '_')}"

# The ATOM-PAIR UNIVERSE TIER. Not a runtime flag anywhere in this system -- the tier
# is a DECLARATION, and this constant plus the sibling declaration in the fabfos
# project's provenance/data/_declared.yml are the two places it is made. Tier 4 =
# tier 3 (predicted ensemble) + the strictly ADDITIVE MetaCyc curated increment, frozen
# 2026-07-20 (see REF_ROOT/TIER4_FREEZE.md). The tier3 siblings stay on disk beside
# these as the comparison baseline; nothing reads them except a movement comparison.
#
# WHY A SUFFIX RATHER THAN A SWAP: the tier3 tables are the only evidence of what moved
# between tiers, so they are never overwritten. The cost is that path and tier are
# coupled here -- which is exactly why this constant exists rather than four literals.
TIER = 4
_T = "" if TIER == 3 else f"_tier{TIER}"

ATOM_PAIRS = REF_ROOT / f"atom_pairs{_T}.parquet"  # AAM table (MNXR-keyed)
DIRECTION = REF_ROOT / "direction.parquet"        # direction table (MNXR-keyed)
LEDGER = REF_ROOT / "closure_ledger.parquet"      # the closure ledger (per reaction)
MANIFEST = REF_ROOT / "MANIFEST.json"             # pin + content hashes + counts

# =====================================================================
# Seed inputs (the incumbent, single-source tables the reference supersedes)
# =====================================================================
# Until the reference tables are frozen (T2-T4), the ledger reads these so it can be
# built and the worklist measured from day one. Present in atom_pairs_universe = mapped;
# present in direction_annotation = has a direction verdict.
SEED_ATOM_PAIRS = DATA / "fabfos_2026_199" / "ecspr_atom" / "atom_pairs_universe.parquet"
SEED_DIRECTION = DATA / "direction" / "direction_annotation.parquet"

# The reachable sets, for prioritisation (the target is full coverage regardless).
# EVIDENCE_TABLE names the reactions any host/fosmid ORF projects onto (network B's
# pre-graph reach); the base graphs are the strict in-network subset.
EVIDENCE_TABLE = (DATA / "fabfos_2026_199" / "ecspr_clean" / "evidence_network"
                  / "evidence_table_clean.parquet")
BASE_GRAPH_DIR = (Path("/home/tony/agentic_workspace/projects/scadc/metabolic-modelling")
                  / "main" / "metabolic-modelling" / "04_reaction_network" / "cache")
BASE_GRAPH_ELEMENTS = ("C", "N", "S", "P")

# The unweighted met-rxn topology skeleton -- pure reac_prop participation, AAM-independent.
# The reference-fed atom graph re-weights THIS skeleton's edges with w_X derived from the
# frozen ATOM_PAIRS, exactly as the incumbent 08_weighted_bipartite_expanded.py re-weighted
# it from aam_unified_expanded.tsv. Same topology; the only change is honest weights (the
# reference refuses the fabricated transits -- non-molecule/stub/pseudo -- the incumbent kept).
SKELETON = (Path("/home/tony/agentic_workspace/projects/scadc/metabolic-modelling")
            / "main" / "metabolic-modelling" / "_reference_try1" / "betweenness"
            / "cache" / "mnx_bipartite.pkl")
GRAPH_DIR = REF_ROOT / f"graph{_T}"   # reference-fed per-element bipartite graphs
GRAPH_FILES = tuple(f"mnx_bipartite_{e}.pkl" for e in BASE_GRAPH_ELEMENTS)

# =====================================================================
# Schema
# =====================================================================
# atom-pairs = the pair columns + provenance, mirroring the AAM ensemble combiner
# (aam/combine_aam.py). sub_idx/prod_idx are canonical ranks; method records how the
# members agreed; confidence folds into edge weight at promotion.
PAIR_COLS = ("mnxr", "element", "substrate", "product", "sub_idx", "prod_idx", "pair_w")
PAIR_PROVENANCE = ("method", "source", "confidence")
ATOM_PAIRS_COLS = PAIR_COLS + PAIR_PROVENANCE

# direction = the per-MNXR direction schema (the directionality ensemble's DIR_COLUMNS).
DIRECTION_COLS = ("mnxr", "dG_prime", "sigma", "ratio",
                  "dir_tier", "dir_method", "dir_confidence")

# =====================================================================
# Verdict vocabulary (finite, class-based)
# =====================================================================
# Every reaction carries exactly one AAM verdict and one direction verdict. A verdict
# is one of three STATES; a refusal additionally carries one reason CLASS. This is the
# whole acceptance instrument: "100% adjudicated" means every reaction has a verdict
# whose state is one of these, and none is left un-adjudicated.
V_RESOLVED = "resolved"            # a confident map / a ratio with evidence
V_DILUTED = "diluted-ambiguous"    # known but spread (fanout, symmetry, disagreement)
V_REFUSED = "refused"              # an explicit, reasoned refusal (carries a reason class)
VERDICT_STATES = (V_RESOLVED, V_DILUTED, V_REFUSED)
# The one non-verdict the ledger must drive to zero: no verdict assigned yet.
V_PENDING = "pending"

# =====================================================================
# Refusal reason taxonomy (finite; reused from ecspr_aam_rescue's classes)
# =====================================================================
# The classes ecspr_aam_rescue already names -- electron carrier, acyl carrier,
# non-molecule, SEED lump, unreconciled stub -- plus the structural boundary classes a
# universe-scale ledger needs (pseudo-reaction, generic R-group, still-computable,
# no-transfer). Class-based, so one rule disposes of many reactions in T4.
R_PSEUDO = "pseudo_reaction"       # transport / exchange / degenerate -- non-chemistry
R_ELECTRON_CARRIER = "electron_carrier"   # rescuable via a conserved placeholder super-node
R_ACYL_CARRIER = "acyl_carrier"    # ACP/thioester -- atoms flow through; needs MCS super-node
R_NON_MOLECULE = "non_molecule"    # Unknown/Carbon/enzyme-complex/acceptor -- no structure to map
R_SEED_LUMP = "seed_lump"          # a SEED multi-EC lump -- mapping it manufactures the edge
R_UNRECONCILED_STUB = "unreconciled_stub"  # ordinary metabolite, no structure; needs --resolved
R_GENERIC_RGROUP = "generic_rgroup"  # generic reaction CLASS with R-group tokens
R_UNMAPPED_COMPUTABLE = "unmapped_computable"  # all SMILES present, no map yet (a mapper target)
R_NO_TRANSFER = "no_transfer"      # no C/N/S/P atom crosses -- nothing for the atom lane to map
REASON_CLASSES = (
    R_PSEUDO, R_ELECTRON_CARRIER, R_ACYL_CARRIER, R_NON_MOLECULE, R_SEED_LUMP,
    R_UNRECONCILED_STUB, R_GENERIC_RGROUP, R_UNMAPPED_COMPUTABLE, R_NO_TRANSFER,
)
# T4 ADJUDICATION (MNXref 4.5, closed 2026-07-17). The automatic layers ran to exhaustion --
# the 4-member AAM ensemble (RXNMapper + LocalMapper + MetaCyc + MCS/Indigo over the whole
# universe), the redox-carrier placeholder library, and a MASS-BALANCE-VERIFIED --resolved
# crosswalk (a same-name donor accepted only where it balances the stub's own reactions, the
# reaction mass balance standing in for the formula MetaNetX omits). What remains refused is
# refused for a reason no further layer converts WITHOUT FABRICATION, so each residual class is
# now a terminal verdict -- an honest, reasoned refusal, not a deferral:
#   * unreconciled_stub  -- the participant is structureless in MNXref 4.5 (no SMILES, no
#     formula) and no balance-verifiable donor exists; the atom lane declines to invent one.
#     RECLAIM PATH: a future CURATED crosswalk (a stub sharing a name with a structured
#     metabolite but lacking a verifiable link -- the published reclaimable ceiling).
#   * acyl_carrier -- an ACP thioester whose acyl chain AND carrier body are both shipped
#     structureless; atoms flow through, but resolving them means fabricating the specific
#     acyl geometry from a name. RECLAIM PATH: a curated acyl-ACP structure library + the MCS
#     super-node rule (bounded to standard straight-chain acyls).
#   * generic_rgroup -- a participant SMILES carries `*` (an alkane*, I-antigen*, a cytochrome
#     c*): a generic reaction CLASS, not a concrete instance. Terminal by construction (like
#     no_transfer) -- there is no concrete atom to map, ever.
#   * electron_carrier -- a placeholder-able carrier whose reaction stayed unbuildable/unmapped
#     because ANOTHER participant is structureless; the carrier itself transfers no scored
#     element, so the reaction is refused on the blocking participant's account.
# unmapped_computable is NOT terminal-refused: every structure is present but no member produced
# a map (heme/tetrapyrrole/FAS -- mapper-defeating size/symmetry), so every atom ABSTAINS -- the
# AAM analogue of direction ratio 1.0 -- which is a DILUTED verdict, not a refusal. It keeps its
# reason label for provenance. A future release re-opens any of these by re-running T4.
REASON_STILL_RESOLVABLE = ()
# The reason classes that are terminal -- an honest, adjudicated refusal (see T4 note above).
REASON_TERMINAL = (R_PSEUDO, R_NON_MOLECULE, R_SEED_LUMP, R_NO_TRANSFER,
                   R_ELECTRON_CARRIER, R_ACYL_CARRIER, R_UNRECONCILED_STUB, R_GENERIC_RGROUP)


# =====================================================================
# Content pin
# =====================================================================
def sha256(path: Path, *, chunk: int = 1 << 20) -> str:
    """Streaming SHA-256 of a file (chem_prop is ~0.8 GB; do not read it whole)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def reac_prop_hash() -> str:
    """The content hash of the reaction universe -- what the pin is ABOUT.

    reac_prop (~10 MB) defines every reaction id and equation, so it is the smallest
    file whose hash fixes the universe. chem_prop is hashed too at freeze time (it
    supplies the structures) but it is 80x larger, so the ledger keys on this one.
    """
    return sha256(REAC_PROP)


# =====================================================================
# The reference error and the assertion helper
# =====================================================================
class ReferenceError(AssertionError):
    """The reference is stale, un-pinned, or incomplete. Always loud, never a warning."""


def load_manifest() -> dict:
    if not MANIFEST.exists():
        raise ReferenceError(
            f"no reference manifest at {MANIFEST}. The reference has not been frozen. "
            f"Build and freeze it (T1->T5) before a consumer may read it as static input."
        )
    return json.loads(MANIFEST.read_text())


def assert_reference(*, check_hash: bool = True, require_closed: bool = True) -> dict:
    """Refuse a reference that is not the pinned, complete one. Mirrors
    ``canon.assert_canonical_*``: raises ``ReferenceError`` loudly, returns the manifest.

    Checks, most-useful error first:
      * the manifest exists and pins this module's MetaNetX version;
      * the reac_prop content hash still matches the pin (the universe has not moved
        under the reference) -- skippable with ``check_hash=False`` for a fast path;
      * both tables exist and carry the declared schema, keyed on MNXR;
      * (``require_closed``) the ledger asserts 0 un-adjudicated reactions.

    A consumer calls this before trusting the reference, so a drifted universe or an
    incomplete freeze fails at the point of use instead of silently scoring the wrong
    atom graph.
    """
    man = load_manifest()

    if man.get("mnxref_version") != MNXREF_VERSION:
        raise ReferenceError(
            f"manifest pins MetaNetX {man.get('mnxref_version')!r} but reference.py is "
            f"{MNXREF_VERSION!r}. A reference built against one release does not name "
            f"the universe of another."
        )
    if check_hash:
        got = reac_prop_hash()
        if man.get("reac_prop_sha256") != got:
            raise ReferenceError(
                f"reac_prop hash moved under the reference: manifest "
                f"{man.get('reac_prop_sha256')!r} != current {got!r}. The reaction "
                f"universe changed; the reference must be rebuilt, not trusted."
            )

    import pandas as pd
    for path, cols, label in ((ATOM_PAIRS, ATOM_PAIRS_COLS, "atom-pairs"),
                              (DIRECTION, DIRECTION_COLS, "direction")):
        if not path.exists():
            raise ReferenceError(f"{label} table absent at {path}; reference incomplete.")
        head = pd.read_parquet(path, columns=None).head(0)
        missing = [c for c in cols if c not in head.columns]
        if missing:
            raise ReferenceError(
                f"{label} table {path.name} missing column(s) {missing}; "
                f"expected reference.{'ATOM_PAIRS_COLS' if label=='atom-pairs' else 'DIRECTION_COLS'}."
            )
        if "mnxr" not in head.columns:
            raise ReferenceError(f"{label} table is not MNXR-keyed.")

    if require_closed:
        n_pending = man.get("n_pending", None)
        if n_pending is None:
            raise ReferenceError(
                "manifest does not record n_pending; the closure ledger has not asserted "
                "completeness. Re-run the ledger and refreeze."
            )
        if n_pending != 0:
            raise ReferenceError(
                f"closure ledger reports {n_pending} un-adjudicated reaction(s); the "
                f"reference is not closed. Every MNXR must carry a verdict (resolved / "
                f"diluted-ambiguous / refused) before the reference may be trusted."
            )
    return man
