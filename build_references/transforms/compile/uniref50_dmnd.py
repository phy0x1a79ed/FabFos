"""R4 -- the DIAMOND database.

Separate from the fetch so a DIAMOND version bump rebuilds the index without
re-downloading 12 GB.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("env::diamond.env"))
fasta = model.AddRequirement(lib.GetType("raw::uniref50_fasta"))
db    = model.AddProduct(lib.GetType("ref::uniref50_diamond_db"))

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- compile/uniref50_dmnd.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the build is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=8, memory=Size.GB(64), duration=Duration(hours=12)),
)
