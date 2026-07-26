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
# The method, vendored build-side. ~1,300 lines between the extractor, the two neural
# members' runner, the curated member's loader and the combiner -- none of which belongs
# in a driver string, and none of which should ship in the wheel (nothing here runs
# during a fosmid pipeline).
extractor  = model.AddRequirement(lib.GetType("buildlib::ecspr_atom_pairs.py"))
neural     = model.AddRequirement(lib.GetType("buildlib::aam_neural_members.py"))
curated_m  = model.AddRequirement(lib.GetType("buildlib::aam_metacyc_member.py"))
combiner   = model.AddRequirement(lib.GetType("buildlib::aam_combine.py"))
pairs      = model.AddProduct(lib.GetType("interm::aam_pairs"))


def protocol(context: ExecutionContext):
    irp   = context.Input(reac_prop)
    icp   = context.Input(chem_prop)
    irx   = context.Input(reac_xref)
    icur  = context.Input(curated)
    ilib  = context.Input(combiner)
    iout  = context.Output(pairs)
    libdir = ilib.container.parent

    # The two neural members, each over the whole reaction universe. This is the single
    # biggest compute in the reference build. Both caches are written into the work
    # directory and are RESUMABLE off their own TSV, so a killed run costs the reaction it
    # was on rather than the hours before it.
    #
    # Both members see the SAME reaction SMILES, built once by aam_neural_members from
    # reac_prop + chem_prop. That matters because they are the CORRELATED pair whose
    # agreement the combiner discounts: measuring disagreement between two mappers is only
    # meaningful if it is not partly disagreement between two SMILES builders.
    for member in ("rxnmapper", "localmapper"):
        _cmd = f"""
            PYTHONPATH={libdir} python3 {libdir}/aam_neural_members.py \
                --member {member} \
                --reac-prop {irp.container} \
                --chem-prop {icp.container} \
                --out _aam_{member}.tsv
        """
        context.ExecWithEnv() \
            .ifContainerDo(env=image, cmd=_cmd) \
            .ifVirtualEnvDo(env=image, cmd=_cmd)

    # Fuse. MetaCyc joins as the third, INDEPENDENT member -- it is a curated database
    # rather than a transformer, so a MetaCyc-inclusive consensus fuses undiscounted and
    # MetaCyc breaks ties the two neural members cannot break between themselves.
    _cmd = f"""
        PYTHONPATH={libdir} python3 {libdir}/aam_combine.py \
            --rxnmapper _aam_rxnmapper.tsv \
            --localmapper _aam_localmapper.tsv \
            --metacyc {icur.container} \
            --reac-xref {irx.container} \
            --reac-prop {irp.container} \
            --chem-prop {icp.container} \
            --out {iout.container}
    """
    context.ExecWithEnv() \
        .ifContainerDo(env=image, cmd=_cmd) \
        .ifVirtualEnvDo(env=image, cmd=_cmd)

    return ExecutionResult(
        manifest=[{pairs: iout.local}],
        success=iout.local.exists() and iout.local.stat().st_size > 0,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=8, memory=Size.GB(64), duration=Duration(hours=24)),
)
