"""MetaNetX 4.5 -- the identity space everything else is expressed in.

Four files, four types. Typing them separately is what lets the DAG show that
`mnxr_lookup` needs reac_prop + reac_xref while the atom mapping needs chem_prop
and neither needs the other's -- one `metanetx` bundle type would make every
consumer depend on all four and the graph would stop being informative.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image      = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
chem_prop  = model.AddProduct(lib.GetType("raw::metanetx_chem_prop"))
chem_xref  = model.AddProduct(lib.GetType("raw::metanetx_chem_xref"))
reac_prop  = model.AddProduct(lib.GetType("raw::metanetx_reac_prop"))
reac_xref  = model.AddProduct(lib.GetType("raw::metanetx_reac_xref"))

BASE_URL = "https://www.metanetx.org/ftp/4.5"

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- acquire/metanetx.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the fetch is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(hours=2)),
)
