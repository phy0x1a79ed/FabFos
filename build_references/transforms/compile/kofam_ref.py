"""R3 -- unpack the KOfam profiles and carry the KO list through.

The ko_list passes through unchanged. It is re-emitted under the run-side type
rather than consumed straight from raw:: so that every input a run tool sees comes
from the compiled tier -- a run tool reaching into raw/ is how a pipeline quietly
acquires an undeclared dependency on the acquisition layout.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image      = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
archive    = model.AddRequirement(lib.GetType("raw::kofam_profiles_archive"))
raw_kolist = model.AddRequirement(lib.GetType("raw::kofam_ko_list"))
profiles   = model.AddProduct(lib.GetType("ref::kofamscan_profiles"))
ko_list    = model.AddProduct(lib.GetType("ref::kofamscan_ko_list"))

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- compile/kofam_ref.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the build is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(hours=1)),
)
