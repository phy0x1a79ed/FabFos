"""B3 -- ONE GPR table over every benchmark condition.

Carries the edges each condition adds or deletes, across four cohorts:
gof_native, gof_het, lof, eydallin. Together with the per-host GPR tables it holds
every edge any condition needs, which is what lets each x in X be CONSTRUCTED at
run time from the conditions table instead of being stored as a network.

The per-protein belief weight is the deployed scheme and is not reinvented here:
each protein carries a total belief of 1.0 split equally across the lanes that
annotated it, so a protein resolved by one lane and one resolved by three carry
identical mass -- the contrast downstream measures topology, not how many lanes
happened to fire. Within a lane the share splits across nominations by raw score,
then spreads evenly across each nomination's fanout. Fanout is left uncapped:
dilution is the designed answer to promiscuous EC fanout, and capping is a
scoring-time policy that does not belong in an annotation artifact.

PROVISIONAL contract -- the literature inputs are confirmed while porting
benchmark_v4.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
laser     = model.AddRequirement(lib.GetType("raw::laser_records"))
keio      = model.AddRequirement(lib.GetType("raw::keio_records"))
eydallin  = model.AddRequirement(lib.GetType("raw::eydallin_records"))
het       = model.AddRequirement(lib.GetType("raw::het_screen_records"))
bridge    = model.AddRequirement(lib.GetType("ref::mnxr_lookup"))
out       = model.AddProduct(lib.GetType("bench::condition_gpr"))

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- benchmark/condition_gpr.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the build is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=4)),
)
