"""MetaCyc flat-files -- LICENSED. Verifies the drop-in and splits it; never fetches.

This transform must not acquire anything. MetaCyc/BioCyc is licensed and not
redistributable, so the only correct behaviour is to check that the operator has
placed the distribution, verify its digests, and refuse with instructions if not.
A download here would be a licensing violation wearing the shape of convenience.

`raw::metacyc_flatfiles` is THE ONE GIVEN of the reference DAG -- the single input
that is staged rather than produced. Everything else closes from it.

It is split into two typed files because two INDEPENDENT ensemble members read
them, and typing only the directory would make each member depend on the other's
input:

  * atom-mappings-smiles.dat -> the AAM ensemble's curated member. It ships mapped
    reaction SMILES in exactly the [C:n]...>>... form the extractor already
    consumes, so there is no index decoder to re-derive.
  * reactions.dat -> the direction ensemble's curated member (REACTION-DIRECTION).

Both are the independent vote in their ensemble; the other members are correlated
with each other (two transformers, two TECRDB-fitted predictors). Losing this
drop-in does not shrink either ensemble evenly -- it removes the only member that
can break a tie.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
flatfiles = model.AddRequirement(lib.GetType("raw::metacyc_flatfiles"))
smiles    = model.AddProduct(lib.GetType("raw::metacyc_atom_mappings_smiles"))
rxns      = model.AddProduct(lib.GetType("raw::metacyc_reactions"))

DROP_IN = "data/raw/metacyc"
LICENSE_NOTE = (
    "MetaCyc flat-files are licensed and not redistributable. Place the "
    f"distribution at {DROP_IN}/ manually; this transform verifies and splits it, "
    "it never fetches."
)


def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- acquire/metacyc_licensed.py declares what it "
        "consumes and produces so the planner can resolve the DAG; the verify "
        "step is not written yet. It must NEVER grow a download path."
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=flatfiles,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(2), duration=Duration(minutes=10)),
)
