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
    ifasta = context.Output(fasta)
    # -c so an interrupted 12 GB transfer resumes rather than restarting, and the file
    # lands gzipped: `diamond makedb` reads .gz directly and the label pool streams it,
    # so decompressing would add 60 GB to the acquisition tier for nothing.
    _cmd = f"""
        wget -q -c {UNIREF50_URL} -O {ifasta.container}
    """
    context.ExecWithEnv() \
        .ifContainerDo(env=image, cmd=_cmd) \
        .ifVirtualEnvDo(env=image, cmd=_cmd)
    return ExecutionResult(
        manifest=[{fasta: ifasta.local}],
        success=ifasta.local.exists() and ifasta.local.stat().st_size > 0,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(hours=12)),
)
