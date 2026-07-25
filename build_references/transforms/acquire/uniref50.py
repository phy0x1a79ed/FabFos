"""UniRef50 cluster representatives.

Acquisition only -- the DIAMOND database is built by compile/uniref50_dmnd.py.
The shipped library's logistics/downloadUniRef50DB.py fuses the two into one step;
splitting them is what lets the 12 GB fetch be pinned and reused when the database
is rebuilt for a new DIAMOND version.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
fasta = model.AddProduct(lib.GetType("raw::uniref50_fasta"))

UNIREF50_URL = "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- acquire/uniref50.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the fetch is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(hours=12)),
)
