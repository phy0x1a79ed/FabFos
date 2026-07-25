"""The curated genome-scale model for a host, from BiGG.

Per host, keyed on the same accession the genome is fetched with, so the model and
the genome cannot drift apart in the DAG.

EPI300 has no published model of its own and borrows DH10B's. That is not an
omission to paper over: `check_epi300_identity` re-measures the delta, and its
edit list came back EMPTY, which is what makes the two hosts' GEM tables
legitimately identical rather than accidentally so. Any GEM-side comparison of the
two strains is therefore null by construction.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
acc   = model.AddRequirement(lib.GetType("ncbi::assembly_accession"))
gem   = model.AddProduct(lib.GetType("raw::host_gem"))

BIGG_URL = "http://bigg.ucsd.edu/static/models/{model}.json"

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- acquire/host_gem.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the fetch is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=acc,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(2), duration=Duration(minutes=30)),
)
