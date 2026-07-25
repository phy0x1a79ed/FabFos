"""R6b -- per-reaction directionality, as a continuous ratio.

Three members: eQuilibrator and dGbyG are both TECRDB-fitted and therefore
correlated (floored by a shared-error term so they cannot vote as two independent
predictors); MetaCyc REACTION-DIRECTION is independent.

Two correctness facts this step exists to get right:
  * ORIENTATION. REACTION-DIRECTION is stated in MetaCyc's equation orientation,
    but MNXref re-canonicalises orientation on import, so a naive metacyc->MNXR
    join INVERTS the curated direction on ~60% of reactions. Every curated call is
    re-expressed in MNXR orientation by comparing compound sets; undecidable cases
    are recorded, never guessed.
  * WILDCARD ABSTENTION. dGbyG returns a confident number for R-group wildcards,
    because RDKit gives `*` a valid feature vector. The member abstains on any
    reaction carrying one.

No evidence shrinks toward dG'=0, giving ratio 1.0 -- reversible as a LIMIT, not as
an if-branch, so a reaction the ensemble is silent on is a provable no-op.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image      = model.AddRequirement(lib.GetType("env::equilibrator.env"))
reac_prop  = model.AddRequirement(lib.GetType("raw::metanetx_reac_prop"))
chem_prop  = model.AddRequirement(lib.GetType("raw::metanetx_chem_prop"))
chem_xref  = model.AddRequirement(lib.GetType("raw::metanetx_chem_xref"))
reac_xref  = model.AddRequirement(lib.GetType("raw::metanetx_reac_xref"))
curated    = model.AddRequirement(lib.GetType("raw::metacyc_reactions"))
eq_cache   = model.AddRequirement(lib.GetType("raw::equilibrator_cache"))
annotation = model.AddProduct(lib.GetType("interm::direction_annotation"))

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- compile/direction_ensemble.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the build is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=8, memory=Size.GB(32), duration=Duration(hours=12)),
)
