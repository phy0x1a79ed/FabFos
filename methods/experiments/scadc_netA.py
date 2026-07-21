"""SCADC Network A: the curated-GEM host base.

The counterpart to `scadc_main`, which scores against the annotation-derived host
base (Network B). This experiment asks the same question of a host reconstructed a
completely different way: a curated GEM's reactome, crosswalked to current MNXR and
induced on the atom-mapped universe with uniform conductance. Two reconstructions
that share no derivation path; a finding that holds on both rests on neither.

WHY THIS EXPERIMENT EXISTS AS A SEPARATE FILE
---------------------------------------------
Not for tidiness. `ecspr::base_graphs` has two producers -- the annotation lane and
the GEM lane -- and the planner enumerates producers by properties alone. If both are
loadable at once and both are satisfiable, the planner picks by score, which means
the network a run scored is decided by a tiebreak rather than by the caller. The two
lanes therefore live in two DOMAINS, and an experiment loads exactly one. That is
what makes "which host did this number come from?" answerable from the manifest.

WHY THE SOLVE IS NOT STAGED HERE
--------------------------------
`scadc_main` stages the solve because its job is to gate the spine and the scorer
against a frozen referent, and re-deriving the solve would destroy that signal.
This experiment has no referent -- Network A has never been solved end to end -- so
the solve is the point, not the risk.

WHY THE SOLVE IS DIRECTED, AND WHERE ITS DIRECTION COMES FROM
------------------------------------------------------------
Network A runs the DIRECTED (rectified "diode") solve, not the undirected one. Its
per-reaction direction is NATIVE to the curated GEM: iECDH10B already encodes each
reaction's reversibility as flux bounds (lower_bound >= 0 -> forward-irreversible,
upper_bound <= 0 -> reverse-irreversible, lower < 0 < upper -> reversible). The
reference and Network B instead take direction from a THERMODYNAMIC ensemble ratio
g_rev/g_fwd = exp(dG'/RT); using the GEM's own bounds is the honest, native
directionality for this host, and it shares no derivation path with netB's -- the same
point the two reconstructions make about the base.

The bounds are BINARY (a reaction is reversible or it is not), so the direction table is
a FIXED-MAGNITUDE diode -- reversible -> ratio 1.0, forward-irreversible ->
DIODE_BACKWARD_FLOOR, reverse-irreversible -> its reciprocal -- unlike netB's continuous
exp(dG'/RT). It is built once by the canon-free engine helper
transforms/ecsprNetA/_gem_direction.py and staged here as `ecspr::direction_ratios`
(canon.NETA_DIR_TABLE); the substrate/product roles the solver orients edges by are the
shared, network-agnostic reac_prop.tsv (canon.DIR_REAC_PROP), staged as
`ecspr::reaction_roles`. The directed solve transform lives in the `ecsprDirected` domain
and produces the SAME `ecspr::{reff,ieff}_axes_report` the significance lane consumes.

CAVEAT for the lead: the directed solve and the undirected `solve.py` (currently in the
`ecspr` domain) produce the same axes_report products. Loading both is the very planner
tiebreak the module docstring below warns about -- so this spec loading `ecsprDirected`
REQUIRES the undirected `solve.py` NOT to be co-loaded (a domain split of `ecspr` into
its shared lane vs. the undirected solve). netA's binary bounds also produce EXTREME
ratios (~1e-9 / 1e9), the backflow regime the directed solver is being hardened for, so
these results depend on the solver-stability fix landing.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from fabfos import canon
from spec import DirInput, ExperimentSpec, Value


def assert_axes():
    """Assert the axis set by COUNT, not by filename -- the filename is what drifted."""
    canon.assert_canonical_axes(canon.AXES_JSON)


def assert_gem_is_dh10b():
    """Refuse a GEM that is not the DH10B reconstruction.

    The sibling curated model (iML1515) is K-12 and carries a DIFFERENT genotype;
    scoring the SCADC host against it would be a silent substitution of one strain
    for another. The model id is inside the JSON, so this asserts the content rather
    than the filename.
    """
    import json
    with open(canon.NETA_GEM) as fh:
        model = json.load(fh)
    mid = str(model.get("id", ""))
    if "DH10B" not in mid.upper():
        raise SystemExit(
            f"curated GEM id is {mid!r}; expected the DH10B reconstruction. "
            f"The K-12 sibling is a different strain with a different genotype."
        )


def assert_direction_table():
    """Refuse a GEM-derived direction table that is not the load_direction_ratios contract.

    The directed solve reads this parquet with `ecspr_network.load_direction_ratios`, which
    keys `ratio` by `mnxr`. This is the NETA lane's OWN table (a fixed-magnitude diode from
    GEM bounds), NOT netB's 7-column thermodynamic annotator, so it is checked against the
    minimal ratio contract rather than canon.assert_canonical_direction_table: exactly the
    two columns, one row per MNXR, and a strictly-positive conductance ratio (a missing MNXR
    is the solver's own reversible default -- never a null here).
    """
    import pandas as pd
    df = pd.read_parquet(canon.NETA_DIR_TABLE)
    cols = list(df.columns)
    if cols[:2] != ["mnxr", "ratio"]:
        raise SystemExit(
            f"{canon.NETA_DIR_TABLE.name} columns are {cols}; expected (mnxr, ratio) -- the "
            f"ecspr_network.load_direction_ratios contract."
        )
    if not df["mnxr"].is_unique:
        raise SystemExit("direction table has duplicate MNXR keys; the per-MNXR collapse failed.")
    if not (df["ratio"] > 0).all():
        raise SystemExit("direction ratio must be strictly positive (it is a conductance ratio).")


SPEC = ExperimentSpec(
    name="scadc_netA",
    inputs={
        # the curated host reconstruction + the crosswalk to current MNXR
        "ecspr::curated_gem": canon.NETA_GEM,
        "ecspr::metanetx_reac_xref": canon.REAC_XREF,
        # reused reference: per-element atom-mapped bipartite universe (the honest
        # reference graph, via the canonical pin -- not the incumbent star cache)
        "ecspr::mnx_bipartite": DirInput(
            [canon.BIPARTITE_DIR / f for f in canon.BIPARTITE_FILES]),
        "ecspr::biomass_axes": canon.AXES_JSON,
        # directionality for the DIRECTED solve. roles = the shared, network-agnostic
        # MetaNetX reac_prop (substrate/product per reaction, the edge orientation);
        # direction_ratios = Network A's OWN table, the GEM's native flux bounds mapped to
        # a fixed-magnitude diode (canon.NETA_DIR_TABLE, built by gem_direction.py). This is
        # what makes the solve directed WITHOUT netB's thermodynamic ensemble.
        "ecspr::reaction_roles": canon.DIR_REAC_PROP,
        "ecspr::direction_ratios": canon.NETA_DIR_TABLE,
        # the reactions each fosmid injects onto the host base
        "ecspr::addition_weights": canon.ADDITION_WEIGHTS,
        # ORF counts for null size matching
        "sequences::open_reading_frames": canon.ORFS_FAA,
        "ecspr::frozen_null": DirInput(
            [canon.INCUMBENT_K1000 / f for f in canon.FROZEN_NULL_FILES]),
        "ecspr::compute_profile": Value("compute_profile.yml", canon.COMPUTE_CPU),
    },
    per_experiment=frozenset({
        "ecspr::addition_weights",
        "sequences::open_reading_frames",
    }),
    targets=(
        "ecspr::reff_significance",
        "ecspr::ieff_significance",
    ),
    # The shared lane, EXACTLY ONE base builder (`ecsprNetA`), and the DIRECTED solve
    # (`ecsprDirected`). `ecspr` deliberately holds no base builder, so neither network
    # arrives by default. `ecsprDirected` supersedes the undirected `solve.py` as the
    # producer of `ecspr::{reff,ieff}_axes_report`; see the module docstring's caveat --
    # `solve.py` must not be co-loaded, or the planner faces two producers for one product.
    domains=("ecspr", "ecsprNetA", "ecsprDirected"),
    namespaces=("sequences", "fosmids", "ecspr"),
    preflight=(assert_axes, assert_gem_is_dh10b, assert_direction_table),
)
