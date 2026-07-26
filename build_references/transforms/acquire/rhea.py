"""Rhea DR mappings -- the UniProt half of the bridge.

Only the two rhea2uniprot files are typed. The rest of the Rhea distribution
(rhea2ec, rhea2kegg_reaction, rhea2xrefs, the reaction SMILES) is not consumed by
any compiled reference here: UniProt reaches MNXR through MetaNetX's reac_xref
`rhea:` rows, not through Rhea's own cross-references.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image   = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
swiss   = model.AddProduct(lib.GetType("raw::rhea2uniprot"))
trembl  = model.AddProduct(lib.GetType("raw::rhea2uniprot_trembl"))

BASE_URL = "https://ftp.expasy.org/databases/rhea/tsv"


def protocol(context: ExecutionContext):
    isw = context.Output(swiss)
    itr = context.Output(trembl)
    # The TrEMBL half stays GZIPPED. It is the bulk of the 35.7M bridge rows and the
    # consumer reads it compressed; decompressing here would put ~1 GB of TSV in the
    # acquisition tier to save one pandas argument.
    _cmd = f"""
        wget -q {BASE_URL}/rhea2uniprot.tsv -O {isw.container}
        wget -q {BASE_URL}/rhea2uniprot_trembl.tsv.gz -O {itr.container}
    """
    context.ExecWithEnv() \
        .ifContainerDo(env=image, cmd=_cmd) \
        .ifVirtualEnvDo(env=image, cmd=_cmd)
    return ExecutionResult(
        manifest=[{swiss: isw.local, trembl: itr.local}],
        success=isw.local.exists() and itr.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(hours=2)),
)
