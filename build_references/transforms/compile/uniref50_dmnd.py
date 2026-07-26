"""R4 -- the DIAMOND database.

Separate from the fetch so a DIAMOND version bump rebuilds the index without
re-downloading 12 GB.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("env::diamond.env"))
fasta = model.AddRequirement(lib.GetType("raw::uniref50_fasta"))
db    = model.AddProduct(lib.GetType("ref::uniref50_diamond_db"))

def protocol(context: ExecutionContext):
    ifasta = context.Input(fasta)
    idb = context.Output(db)
    # `diamond makedb` reads the gzip directly, so the 12 GB acquisition never has to be
    # expanded to ~60 GB on disk. -d takes a PREFIX and appends .dmnd, so the database is
    # built beside the product and moved onto it.
    _cmd = f"""
        diamond makedb --in {ifasta.container} -d uniref50 --threads ${{SLURM_CPUS_PER_TASK:-8}}
        mv uniref50.dmnd {idb.container}
    """
    context.ExecWithEnv() \
        .ifContainerDo(env=image, cmd=_cmd) \
        .ifVirtualEnvDo(env=image, cmd=_cmd)
    size = idb.local.stat().st_size if idb.local.exists() else 0
    Log.Info(f"uniref50.dmnd {size/1e9:.1f} GB")
    # A DIAMOND database of UniRef50 is tens of GB. A few hundred MB means makedb read a
    # truncated download and exited zero, and the uniref lane would then quietly find
    # nothing for most ORFs.
    return ExecutionResult(
        manifest=[{db: idb.local}],
        success=size > 1_000_000_000,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=8, memory=Size.GB(64), duration=Duration(hours=12)),
)
