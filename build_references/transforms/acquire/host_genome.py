"""The host genome and its proteome, from an NCBI assembly accession.

Deliberately duplicated rather than reused from the shipped library's
logistics/getNcbiAssembly. That directory also carries downloadKofamDB and
downloadUniRef50DB, which produce ref::kofamscan_profiles / ref::kofamscan_ko_list
/ ref::uniref50_diamond_db -- exactly the types compile/kofam_ref.py and
compile/uniref50_dmnd.py produce. Loading it would give the planner two producers
for each and let a TIEBREAK decide which one built a reference, which is not a
thing that should be decided by a tiebreak. Keeping the accession fetch here means
the build library never has to load a directory whose other contents compete with
it.

The proteome half (sequences::orfs) is what the de-novo GPR lanes consume, so a
host's evidence comes out of the shipped annotation pipeline rather than from a
frozen annotation intermediate.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("env::ncbi-datasets.env"))
acc   = model.AddRequirement(lib.GetType("ncbi::assembly_accession"))
fna   = model.AddProduct(lib.GetType("sequences::isolate_assembly"))
faa   = model.AddProduct(lib.GetType("sequences::orfs"))
gff   = model.AddProduct(lib.GetType("sequences::gff"))
gbk   = model.AddProduct(lib.GetType("sequences::gbk"))


def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- acquire/host_genome.py declares what it consumes "
        "and produces so the planner can resolve the DAG; the fetch is not written "
        "yet. See logistics/getNcbiAssembly.py in the shipped library for the "
        "`datasets download genome accession` form this will take."
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=acc,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(hours=1)),
)
