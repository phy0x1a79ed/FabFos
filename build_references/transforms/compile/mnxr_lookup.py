"""R5 -- the consolidated intermediate-id -> MNXR bridge.

One table replacing the former ko/ec/uniprot trio: `id, id_source, mnxr,
evidence_quality`. Measured before consolidating -- the three id spaces share no
ids, so a single `id` column is unambiguous and `id_source` is a label rather than
a disambiguator; deduping to distinct (id, mnxr) went 35,762,706 -> 30,467,712 rows
and is behaviour-preserving because the consumer already ran that dedup after its
join, and the duplicate rows differed only in columns nothing downstream reads.
82.5 MB against 251.5 MB across the three files it replaces.

`evidence_quality` is carried at a cost of 0.2 MB because "reviewed" sorts before
"unreviewed", so an existing keep-first dedup already retains the stronger claim.

Three routes into one table:
  ec       reac_prop classifs
  ko       ko_to_kegg_r -> reac_xref `kegg.reaction:` rows
  uniprot  rhea2uniprot{,_trembl} -> reac_xref `rhea:` rows
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image      = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
reac_prop  = model.AddRequirement(lib.GetType("raw::metanetx_reac_prop"))
reac_xref  = model.AddRequirement(lib.GetType("raw::metanetx_reac_xref"))
ko_kegg_r  = model.AddRequirement(lib.GetType("raw::kegg_ko_to_kegg_r"))
rhea_sp    = model.AddRequirement(lib.GetType("raw::rhea2uniprot"))
rhea_tr    = model.AddRequirement(lib.GetType("raw::rhea2uniprot_trembl"))
bridge     = model.AddProduct(lib.GetType("ref::mnxr_lookup"))

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- compile/mnxr_lookup.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the build is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=2, memory=Size.GB(32), duration=Duration(hours=2)),
)
