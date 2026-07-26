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
encoding   = model.AddRequirement(lib.GetType("buildlib::refs_encoding.py"))
baker      = model.AddRequirement(lib.GetType("buildlib::bake_metabolism.py"))
out_pairs  = model.AddProduct(lib.GetType("ref::atom_pairs"))
out_vocab  = model.AddProduct(lib.GetType("ref::metabolism_vocab"))
out_dir    = model.AddProduct(lib.GetType("ref::direction_ratios"))

# Bake and then verify, in one invocation, against the SAME in-memory paths. The selftest
# is not optional and is not a separate transform: its whole job is to prove this step
# changed nothing, and a bake that is written now and checked later is a bake that can
# ship unchecked.
DRIVER = r'''
import sys
from pathlib import Path

from bake_metabolism import bake, selftest

SRC_PAIRS = Path("{src_pairs}")
SRC_DIR   = Path("{src_direction}")
OUT_VOCAB = Path("{out_vocab}")
OUT_PAIRS = Path("{out_pairs}")
OUT_DIR   = Path("{out_direction}")

bake(SRC_PAIRS, SRC_DIR, OUT_VOCAB, OUT_PAIRS, OUT_DIR)
print()
rc = selftest(SRC_PAIRS, SRC_DIR, OUT_VOCAB, OUT_PAIRS, OUT_DIR)
if rc != 0:
    # A failed round trip means the encoding lost or merged something. The merge case
    # RAISES the network's conductance, so it reads as an improvement downstream -- which
    # is precisely why this exits non-zero instead of warning.
    raise SystemExit(rc)
'''


def protocol(context: ExecutionContext):
    ip   = context.Input(pairs)
    ia   = context.Input(annotation)
    ilib = context.Input(baker)
    iv   = context.Output(out_vocab)
    ipr  = context.Output(out_pairs)
    idr  = context.Output(out_dir)
    libdir = ilib.container.parent

    driver = DRIVER.format(src_pairs=ip.container, src_direction=ia.container,
                           out_vocab=iv.container, out_pairs=ipr.container,
                           out_direction=idr.container)
    context.LocalShell("cat > _bake_metabolism.py << 'PYEOF'\n" + driver + "\nPYEOF\n")
    context.ExecWithContainer(
        image=image, cmd=f"PYTHONPATH={libdir} python3 _bake_metabolism.py")

    return ExecutionResult(
        manifest=[{out_pairs: ipr.local, out_vocab: iv.local, out_dir: idr.local}],
        success=all(p.exists() and p.stat().st_size > 0
                    for p in (ipr.local, iv.local, idr.local)),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=2)),
)
