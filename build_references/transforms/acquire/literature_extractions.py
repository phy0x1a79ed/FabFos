"""Benchmark condition sources that CANNOT be fetched -- Eydallin and the het screen.

Both are tabular extractions a human made once from a paper supplement. There is no URL
for either: Eydallin 2010 ships its hits as a PDF, and the heterologous screening cohort's
protein list was assembled from the screening papers. Under the tier rule that makes them
hand-authored judgement rather than acquired bytes -- `curated/`, not `raw/` -- and
`build_references/REFERENCES.md` records that as an open decision.

Pending that decision they are treated as raw drop-ins, and this transform exists to make
the difference between "absent" and "empty" a sentence rather than a silence. It has one
outcome when the tables are not on disk: REFUSE, by name, saying which paper's extraction
is missing and what it is.

    THIS TRANSFORM CANNOT SUCCEED. That is the design, not a defect.

It is reachable only when the drop-ins are absent, because when they are present the
driver stages them and the planner has nothing to schedule. So its running at all IS the
error condition, and refusing is the whole of its behaviour. The alternative -- leaving
these two types with no producer -- makes the same situation surface as a planner failure
naming a type, which tells a reader that a graph did not close but not that a table needs
extracting.

Building without them would ship a `conditions` table that looks complete and is missing
an arm: the `eydallin` cohort and the entire `gof_het` cohort.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
eydallin  = model.AddProduct(lib.GetType("raw::eydallin_records"))
het       = model.AddProduct(lib.GetType("raw::het_screen_records"))

# What each product must contain, and where the build expects to find it. The path is
# named because a refusal that says "place it in the product directory" is useless when
# the product directory is inside a run that no longer exists.
EXTRACTIONS = (
    ("eydallin", "eydallin2010_hits.tsv", "data/raw/literature/eydallin/",
     "Eydallin 2010's screen hits. The publisher ships them as a PDF supplement; the "
     "tabular form was extracted once by hand in the sibling project."),
    ("het_screen", "heterologous_uniprot.tsv", "data/raw/literature/het_screen/",
     "The heterologous / metagenomic screening cohort's protein list -- one row per "
     "protein with the UniProt accession the screen reported. Itself an extraction from "
     "the screening papers."),
)


def protocol(context: ExecutionContext):
    ieyd = context.Output(eydallin)
    ihet = context.Output(het)

    lines = ["the two extracted literature tables are not on disk.", ""]
    for name, filename, where, what in EXTRACTIONS:
        lines += [f"  {name}: place {filename} in {where}",
                  f"      {what}", ""]
    lines += [
        "This is a NAMED refusal, not a crash. Both are hand extractions with no URL, so",
        "there is nothing for this step to fetch. Building without them ships a",
        "conditions table that looks complete and is missing the eydallin and gof_het",
        "arms entirely.",
        "",
        "Once they are in place the stage 1 driver stages them and this step is never",
        "scheduled -- reaching it at all means they are absent.",
    ]
    raise SystemExit("\n".join(lines))

    # Unreachable, and kept so the contract reads the same as every other transform's:
    # these are the products a reader should expect at these paths.
    return ExecutionResult(                                      # pragma: no cover
        manifest=[{eydallin: ieyd.local, het: ihet.local}],
        success=False,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(1), duration=Duration(minutes=5)),
)
