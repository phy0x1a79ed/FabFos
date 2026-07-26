"""The fosmid host set, as a declaration.

Emits one `ncbi::assembly_accession` per host, so which hosts the references are
built for is a fact IN THE LIBRARY rather than an argument the caller supplies.
That is the point of putting it here: a driver that passes the accessions in can
pass a different three, and then two runs of "the reference build" are not the same
build. This way the host set moves only when this file moves.

Scatter pattern -- one step, N outputs via `context.Output(acc, i=i)`, the same
shape logistics/scatterNcbiAccession uses. Downstream steps (host_genome, host_gem)
then fan out per accession.

  e_coli_k12     GCF_000005845.2   K-12 MG1655; curated model iML1515
  e_coli_dh10b   GCF_000019425.1   DH10B;       curated model iECDH10B_1368
  e_coli_epi300  GCF_051228345.1   EPI300;      no model of its own -- borrows
                                                DH10B's, which a measured, EMPTY
                                                edit list is what licenses. Any
                                                GEM-side comparison of the two is
                                                therefore null by construction.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
acc   = model.AddProduct(lib.GetType("ncbi::assembly_accession"))

HOSTS = {
    "e_coli_k12":    "GCF_000005845.2",
    "e_coli_dh10b":  "GCF_000019425.1",
    "e_coli_epi300": "GCF_051228345.1",
}


def protocol(context: ExecutionContext):
    # One output per host -- the scatter. Everything downstream (host_genome, host_gem,
    # the four annotation lanes, both GPR tables) fans out from here, so the host set is
    # a fact in this file and nowhere else.
    outputs = []
    for i, (host, accession) in enumerate(sorted(HOSTS.items())):
        outf = context.Output(acc, i=i)
        # The accession alone, one line, because that is what the shipped
        # logistics/getNcbiAssembly.py reads. The host NAME is deliberately not written
        # here: it would make this file a two-column format that only these transforms
        # understand, and the host name is recoverable from the accession downstream.
        with open(outf.local, "w") as f:
            f.write(accession + "\n")
        Log.Info(f"host {host} -> {accession}")
        outputs.append({acc: outf.local})
    return ExecutionResult(manifest=outputs, success=len(outputs) == len(HOSTS))


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(1), duration=Duration(minutes=5)),
)
