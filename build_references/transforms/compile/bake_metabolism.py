"""R6c -- bake the ensemble outputs into the three compiled tables.

ONE artifact in three files. Each carries an identical bake-identity block (vocab
digest, bit widths, element order, source hashes) and a reader refuses a mismatched
trio -- reading atom_pairs against another bake's vocab decodes every node to the
wrong metabolite SILENTLY. That is why all three are produced by one transform in
one pass rather than by three.

This is a re-encoding, not a re-derivation: nothing here recomputes a mapping or a
free energy. The failure mode that matters is a too-narrow rank field merging two
distinct atoms onto one node, which RAISES the network's conductance and so reads
as an improvement rather than as a bug -- hence the row-by-row round-trip check
over every row rather than a sample.

Bit widths are measured from the data, not chosen, and the maximum atom rank is
taken over ALL elements: measuring it on the carbon slice alone gives 492 where the
true maximum is 496, and that off-by-a-slice is exactly what truncates a node key.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image      = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
pairs      = model.AddRequirement(lib.GetType("interm::aam_pairs"))
annotation = model.AddRequirement(lib.GetType("interm::direction_annotation"))
out_pairs  = model.AddProduct(lib.GetType("ref::atom_pairs"))
out_vocab  = model.AddProduct(lib.GetType("ref::metabolism_vocab"))
out_dir    = model.AddProduct(lib.GetType("ref::direction_ratios"))

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- compile/bake_metabolism.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the build is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=2)),
)
