"""B1 -- the GPR a curated genome-scale model asserts, per host.

Per host (group_by the model), and deliberately a SEPARATE file from the de-novo
table rather than merged with it. The two are independent lines of evidence and the
comparison between them is the point; conflating them into one table destroys the
only thing they are jointly good for.

Needs cobra to open the model, which is also why this cannot run inside the pinned
ecspr image -- that image carries no cobra.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::cobra.env"))
gem       = model.AddRequirement(lib.GetType("raw::host_gem"))
bridge    = model.AddRequirement(lib.GetType("ref::mnxr_lookup"))
reac_xref = model.AddRequirement(lib.GetType("raw::metanetx_reac_xref"))
out       = model.AddProduct(lib.GetType("ref::gpr_table_gem"))

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- benchmark/host_gpr_gem.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the build is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=gem,
    labels=["local"],
    resources=Resources(cpus=2, memory=Size.GB(8), duration=Duration(hours=1)),
)
