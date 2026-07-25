"""Rhea DR mappings -- the UniProt half of the bridge.

Only the two rhea2uniprot files are typed. The rest of the Rhea distribution
(rhea2ec, rhea2kegg_reaction, rhea2xrefs, the reaction SMILES) is not consumed by
any compiled reference here: UniProt reaches MNXR through MetaNetX's reac_xref
`rhea:` rows, not through Rhea's own cross-references.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image   = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
swiss   = model.AddProduct(lib.GetType("raw::rhea2uniprot"))
trembl  = model.AddProduct(lib.GetType("raw::rhea2uniprot_trembl"))

BASE_URL = "https://ftp.expasy.org/databases/rhea/tsv"

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- acquire/rhea.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the fetch is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(hours=2)),
)
