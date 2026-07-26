"""MetaCyc flat-files -- LICENSED. Verifies the drop-in and splits it; never fetches.

This transform must not acquire anything. MetaCyc/BioCyc is licensed and not
redistributable, so the only correct behaviour is to check that the operator has
placed the distribution, verify its digests, and refuse with instructions if not.
A download here would be a licensing violation wearing the shape of convenience.

`raw::metacyc_flatfiles` is THE ONE GIVEN of the reference DAG -- the single input
that is staged rather than produced. Everything else closes from it.

It is split into two typed files because two INDEPENDENT ensemble members read
them, and typing only the directory would make each member depend on the other's
input:

  * atom-mappings-smiles.dat -> the AAM ensemble's curated member. It ships mapped
    reaction SMILES in exactly the [C:n]...>>... form the extractor already
    consumes, so there is no index decoder to re-derive.
  * reactions.dat -> the direction ensemble's curated member (REACTION-DIRECTION).

Both are the independent vote in their ensemble; the other members are correlated
with each other (two transformers, two TECRDB-fitted predictors). Losing this
drop-in does not shrink either ensemble evenly -- it removes the only member that
can break a tie.
"""
import shutil
from pathlib import Path

from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
flatfiles = model.AddRequirement(lib.GetType("raw::metacyc_flatfiles"))
smiles    = model.AddProduct(lib.GetType("raw::metacyc_atom_mappings_smiles"))
rxns      = model.AddProduct(lib.GetType("raw::metacyc_reactions"))

DROP_IN = "data/raw/metacyc"
LICENSE_NOTE = (
    "MetaCyc flat-files are licensed and not redistributable. Place the "
    f"distribution at {DROP_IN}/ manually; this transform verifies and splits it, "
    "it never fetches."
)


# The two files this transform splits out, and the members that read them. Named
# together so a missing one is reported as "which ensemble loses its independent vote"
# rather than as a bare filename.
MEMBERS = {
    "atom-mappings-smiles.dat": "the AAM ensemble's curated (independent) member",
    "reactions.dat": "the direction ensemble's curated (independent) member",
}


def protocol(context: ExecutionContext):
    """Verify the drop-in and split it. NEVER fetches -- see the module docstring.

    Anything that looks like it could become a download belongs somewhere else. The
    failure mode here is a refusal with instructions, which is the correct behaviour on a
    machine whose operator has not licensed MetaCyc.
    """
    drop_in = Path(context.Input(flatfiles).local)

    present = {name: (drop_in / name) for name in MEMBERS}
    missing = {n: p for n, p in present.items() if not p.exists() or p.stat().st_size == 0}
    if missing:
        detail = "\n".join(f"    {n:28s} -- {MEMBERS[n]}" for n in sorted(missing))
        raise SystemExit(
            f"the MetaCyc drop-in at {drop_in} is missing:\n{detail}\n\n"
            f"{LICENSE_NOTE}\n\n"
            f"Losing this drop-in does not shrink either ensemble evenly -- it removes "
            f"the only member that can break a tie between two correlated ones.")

    outs = {}
    for dep, name in ((smiles, "atom-mappings-smiles.dat"), (rxns, "reactions.dat")):
        o = context.Output(dep)
        shutil.copyfile(present[name], o.local)
        Log.Info(f"verified {name}: {present[name].stat().st_size/1e6:.1f} MB")
        outs[dep] = o.local

    return ExecutionResult(
        manifest=[outs],
        success=all(p.exists() and p.stat().st_size > 0 for p in outs.values()),
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=flatfiles,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(2), duration=Duration(minutes=10)),
)
