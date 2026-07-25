"""B2 -- the GPR a host's own annotation lanes infer.

Thin by design. The GPR table itself is produced upstream by the shipped 4-lane
mapper running on the host's proteome, exactly as it would on a fosmid ORF set --
this step attaches host attribution and re-emits it under the reference type.

That the benchmark's de-novo evidence comes out of the SHIPPED pipeline rather than
from a frozen intermediate is the point of wiring it this way: it makes the
benchmark a test of the method rather than of a file someone once produced.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
gpr   = model.AddRequirement(lib.GetType("annotation::gpr_table_4lane"))
out   = model.AddProduct(lib.GetType("ref::gpr_table_denovo"))

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- benchmark/host_gpr_denovo.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the build is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=gpr,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(minutes=30)),
)
