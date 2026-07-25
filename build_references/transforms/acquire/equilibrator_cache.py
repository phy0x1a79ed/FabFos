"""eQuilibrator compound cache.

Pinned explicitly rather than left to the package's first-use fetch. Otherwise a
build silently depends on a network call whose result is not recorded anywhere,
and two builds a year apart can disagree with nothing in the provenance to show it.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("env::equilibrator.env"))
cache = model.AddProduct(lib.GetType("raw::equilibrator_cache"))

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- acquire/equilibrator_cache.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the fetch is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(hours=2)),
)
