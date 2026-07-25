"""B4 -- one row per network to build.

Names each condition's host, arm and tier, and which condition_gpr edges it
applies. This is the table the run-time network construction reads, and it is what
makes X derivable rather than stored.

Takes the condition GPR as input so a condition can never name an edge set that
does not exist -- the alternative, building the two independently and joining them
later, produces a conditions table that looks complete and silently drops rows at
scoring time.

PROVISIONAL contract, as B3.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
cond_gpr  = model.AddRequirement(lib.GetType("bench::condition_gpr"))
laser     = model.AddRequirement(lib.GetType("raw::laser_records"))
keio      = model.AddRequirement(lib.GetType("raw::keio_records"))
eydallin  = model.AddRequirement(lib.GetType("raw::eydallin_records"))
het       = model.AddRequirement(lib.GetType("raw::het_screen_records"))
out       = model.AddProduct(lib.GetType("bench::conditions"))

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- benchmark/conditions.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the build is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=2, memory=Size.GB(16), duration=Duration(hours=1)),
)
