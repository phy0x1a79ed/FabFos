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

BIGG_URL = "http://bigg.ucsd.edu/static/models/{bigg_model}.json"

# accession -> the BiGG model that host's GPR is asserted by. Keyed on the accession
# rather than on a host name because the accession is what arrives here; the host name
# never crosses this boundary.
#
# EPI300 maps to DH10B's model, which is a measured claim rather than a convenience:
# check_epi300_identity re-measures the delta and its edit list came back EMPTY. That is
# what makes the two hosts' GEM tables legitimately identical instead of accidentally so,
# and it is also why any GEM-side comparison of the two strains is null by construction.
GEM_FOR_ACCESSION = {
    "GCF_000005845.2": "iML1515",          # E. coli K-12 MG1655
    "GCF_000019425.1": "iECDH10B_1368",    # E. coli DH10B
    "GCF_051228345.1": "iECDH10B_1368",    # E. coli EPI300 -- borrows DH10B's, see above
}


def protocol(context: ExecutionContext):
    with open(context.Input(acc).local) as f:
        accession = f.readline().strip()
    bigg_model = GEM_FOR_ACCESSION.get(accession)
    if bigg_model is None:
        # Not a fetch failure: it means a host was added to acquire/host_accessions.py
        # without deciding which curated model asserts its GPR. Falling back to any model
        # would attribute one strain's biochemistry to another.
        raise SystemExit(
            f"no curated GEM declared for accession {accession}. Add it to "
            f"GEM_FOR_ACCESSION here -- and if the host has no published model, record "
            f"which host's model it borrows and the measurement that licenses that, the "
            f"way EPI300 borrows DH10B's.")

    igem = context.Output(gem)
    Log.Info(f"accession {accession} -> BiGG model {bigg_model}")
    context.ExecWithContainer(image=image, cmd=f"""
        wget -q {BIGG_URL.format(bigg_model=bigg_model)} -O {igem.container}
    """)
    return ExecutionResult(
        manifest=[{gem: igem.local}],
        success=igem.local.exists() and igem.local.stat().st_size > 0,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=acc,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(2), duration=Duration(minutes=30)),
)
