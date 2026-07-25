"""R6a -- the ensemble atom-atom mapping.

Three members, and the split is the whole design: RXNMapper and LocalMapper are
both transformers over reaction SMILES and therefore CORRELATED, so their agreement
is discounted rather than counted as two independent votes; MetaCyc is a curated
database and therefore INDEPENDENT, so a MetaCyc-plus-neural consensus fuses at
full weight and MetaCyc breaks ties.

Where members disagree the correspondence is NOT decided by majority or by
confidence -- it is spread across the disputed products by member weight, so a
contested atom dilutes its transfer instead of committing to one answer. Provenance
(method / source / confidence) is a record of who spoke, never a usability gate:
every emitted pair is read.

Every member's correspondence is expressed in one identity, (metabolite, canonical
atom rank), which is a property of the metabolite rather than of the mapper -- so
two mappers that map the same physical atom name the same node.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image      = model.AddRequirement(lib.GetType("env::rdkit.env"))
reac_prop  = model.AddRequirement(lib.GetType("raw::metanetx_reac_prop"))
chem_prop  = model.AddRequirement(lib.GetType("raw::metanetx_chem_prop"))
reac_xref  = model.AddRequirement(lib.GetType("raw::metanetx_reac_xref"))
curated    = model.AddRequirement(lib.GetType("raw::metacyc_atom_mappings_smiles"))
pairs      = model.AddProduct(lib.GetType("interm::aam_pairs"))

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- compile/aam_ensemble.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the build is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=8, memory=Size.GB(64), duration=Duration(hours=24)),
)
