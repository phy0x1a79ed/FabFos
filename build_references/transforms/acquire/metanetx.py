"""MetaNetX 4.5 -- the identity space everything else is expressed in.

Four files, four types. Typing them separately is what lets the DAG show that
`mnxr_lookup` needs reac_prop + reac_xref while the atom mapping needs chem_prop
and neither needs the other's -- one `metanetx` bundle type would make every
consumer depend on all four and the graph would stop being informative.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image      = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
chem_prop  = model.AddProduct(lib.GetType("raw::metanetx_chem_prop"))
chem_xref  = model.AddProduct(lib.GetType("raw::metanetx_chem_xref"))
reac_prop  = model.AddProduct(lib.GetType("raw::metanetx_reac_prop"))
reac_xref  = model.AddProduct(lib.GetType("raw::metanetx_reac_xref"))

BASE_URL = "https://www.metanetx.org/ftp/4.5"

# The four files, in the order the products are declared. Named as a table so the fetch
# loop cannot pair a file with the wrong product -- which would be undetectable, because
# all four are MetaNetX TSVs with the same comment header.
FILES = ("chem_prop.tsv", "chem_xref.tsv", "reac_prop.tsv", "reac_xref.tsv")


def protocol(context: ExecutionContext):
    outs = [context.Output(d) for d in (chem_prop, chem_xref, reac_prop, reac_xref)]
    fetch = "\n".join(
        f"wget -q --show-progress {BASE_URL}/{name} -O {o.container}"
        for name, o in zip(FILES, outs))
    context.ExecWithEnv() \
        .ifContainerDo(env=image, cmd=fetch) \
        .ifVirtualEnvDo(env=image, cmd=fetch)
    return ExecutionResult(
        manifest=[{chem_prop: outs[0].local, chem_xref: outs[1].local,
                   reac_prop: outs[2].local, reac_xref: outs[3].local}],
        success=all(o.local.exists() and o.local.stat().st_size > 0 for o in outs),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(hours=2)),
)
