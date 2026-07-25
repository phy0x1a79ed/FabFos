"""KOfam profiles + KO list.

The archive is acquired as shipped and untarred by compile/kofam_ref.py, so the
acquisition stays a byte-for-byte copy of upstream and the unpacking is a separate,
re-runnable step.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
profiles  = model.AddProduct(lib.GetType("raw::kofam_profiles_archive"))
ko_list   = model.AddProduct(lib.GetType("raw::kofam_ko_list"))

PROFILES_URL = "ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz"
KO_LIST_URL  = "ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz"

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- acquire/kofam.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the fetch is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(hours=4)),
)
