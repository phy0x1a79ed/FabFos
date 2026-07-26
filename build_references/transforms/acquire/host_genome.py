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
import glob
import shutil
from pathlib import Path

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
    # Same body as the shipped logistics/getNcbiAssembly.py. The duplication is the
    # point and is argued in the module docstring: loading that directory would give the
    # planner a second producer for ref::kofamscan_* and ref::uniref50_diamond_db.
    dep_path = context.Input(acc)
    with open(dep_path.local) as f:
        accession = f.readline().strip()

    context.ExecWithContainer(
        image=image,
        cmd=f"""\
            datasets download genome accession {accession} \
                --include gff3,protein,genome,gbff
        """,
    )
    context.LocalShell("unzip -o ncbi_dataset.zip")

    output_manifest = {}

    def fix_out(dep, p: Path):
        op = context.Output(dep)
        shutil.move(p, op.local)
        output_manifest[dep] = op.local

    for f in glob.glob("ncbi_dataset/*/*/*"):
        p = Path(f)
        Log.Info(f"scanning file [{p}]")
        match p.name:
            case "genomic.gff":
                fix_out(gff, p)
            case "genomic.gbff":
                fix_out(gbk, p)
            case "protein.faa":
                fix_out(faa, p)
        # `cds_from_genomic.fna` also ends in genomic.fna and is NOT the assembly; taking
        # it would make the background reference a CDS set and every background filter
        # would silently stop matching intergenic sequence.
        if not p.name.startswith("cds") and p.name.endswith("genomic.fna"):
            fix_out(fna, p)

    return ExecutionResult(
        manifest=[output_manifest],
        success=len(output_manifest) == len(model.produces[0]),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=acc,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(hours=1)),
)
