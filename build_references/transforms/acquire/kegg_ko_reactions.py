"""KO -> KEGG reaction, from each KO flat-file's REACTION block.

An ACQUISITION, not a compile step: every row comes from KEGG REST, one
`get/<ko>` per KO at ~3 requests/second. Classifying it here is what retires the
dependency on the sibling project's `kegg_requests.db` -- that cache was a
licensed artifact being read at build time, which is a licensing question rather
than a path fix.

Takes the KO list as input because that is the KO universe to fetch; it does not
read its EC tags. Routing KO through its `[EC:...]` tag instead of through KEGG
reactions is easy and wrong -- measured on the three hosts it takes the kofam lane
from 646 to 8,548 reactions and makes 91% of them reactions the EC lane already
reaches, so the two lanes stop being independent evidence.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image    = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
ko_list  = model.AddRequirement(lib.GetType("raw::kofam_ko_list"))
out      = model.AddProduct(lib.GetType("raw::kegg_ko_to_kegg_r"))

KEGG_URL   = "https://rest.kegg.jp/get/{ko}"
KEGG_PAUSE = 0.4

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- acquire/kegg_ko_reactions.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the fetch is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(hours=8)),
)
