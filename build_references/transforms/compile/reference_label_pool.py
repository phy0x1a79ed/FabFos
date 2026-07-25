"""R7 -- the labelled embedding pool the kNN transfer lane votes against.

No separate labelled-proteome acquisition. The pool is built from two things
already in the graph: UniRef50 sequences, and the labels the bridge already
carries. The bridge's UniProt rows come from rhea2uniprot, so a protein with a
curated reaction assignment is exactly a protein the bridge can label -- which
makes the pool a SUBSET of the reference set rather than a new source beside it.

Scale falls out of the `reviewed` slice rather than being chosen: the bridge is
365,240 reviewed against 35,397,466 unreviewed rows, and reviewed is the right cut
for a set that other proteins get classified against. That lands in the same order
as the deployed pool (273,764 rows / 270,454 labelled / 8,590 distinct MNXR).

WHAT THIS CHANGES. The deployed pool is not this. It is KEGG-derived -- 54,005
sequences keyed on KEGG gene ids (`dme:Dmel_CG3481`), labelled by KO, projecting
KO -> MNXR. Building from the Rhea route instead means:
  + no new acquisition, and no KEGG-licensed sequences in the tree
  + one label source instead of two, so a protein cannot be labelled one way here
    and a different way in the GPR mapper
  - a different set, so the pbert lane's numbers move. This is not a reproduction
    of the deployed lane and must not be reported as one.

Two traps worth stating because both are silent:
  * SAME MODEL. A pool embedded with a different model from the query is not a
    weaker pool, it is a meaningless one -- cosine distance between two embedding
    spaces is a number with no referent.
  * ONE PASS. The consumer indexes the embedding stack by row, so pairing an index
    from one build with a stack from another misindexes every row and emits a full,
    confident, WRONG table with nothing raised. Index and stack come out together
    or not at all.

GATED -- see REFERENCES.md "Open decisions".
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image  = model.AddRequirement(lib.GetType("env::proteinbert.env"))
fasta  = model.AddRequirement(lib.GetType("raw::uniref50_fasta"))
bridge = model.AddRequirement(lib.GetType("ref::mnxr_lookup"))
pool   = model.AddProduct(lib.GetType("ref::reference_label_pool"))

# The cut that defines the pool: bridge rows whose id_source is uniprot and whose
# evidence_quality is reviewed, intersected with UniRef50's representatives.
POOL_ID_SOURCE = "uniprot"
POOL_EVIDENCE = "reviewed"


def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- compile/reference_label_pool.py declares what it "
        "consumes and produces so the planner can resolve the DAG; the build is "
        "not written yet."
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=8)),
)
