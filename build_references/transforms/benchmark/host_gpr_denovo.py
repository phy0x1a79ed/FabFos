"""B2 -- the GPR a host's own annotation lanes infer.

Thin by design. The GPR table itself is produced upstream by the shipped 4-lane
mapper running on the host's proteome, exactly as it would on a fosmid ORF set --
this step attaches host attribution and re-emits it under the reference type.

That the benchmark's de-novo evidence comes out of the SHIPPED pipeline rather than
from a frozen intermediate is the point of wiring it this way: it makes the
benchmark a test of the method rather than of a file someone once produced.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
# The chosen-4 lane set IS `annotation::gpr_table` -- `gpr_4lane` produces the
# type directly rather than a subtype, so there is exactly one producer and no
# tiebreak to express here. (This was `gpr_table_4lane` before the lane sets
# were made canonical; `gpr_table_7lane` remains a subtype.)
gpr   = model.AddRequirement(lib.GetType("annotation::gpr_table"))
# The proteome this table describes, only so its accession can name the host. The mapper
# upstream is host-agnostic by design -- it runs on a fosmid ORF set exactly the same way
# -- so the attribution has to be attached here rather than inside it.
acc   = model.AddRequirement(lib.GetType("ncbi::assembly_accession"))
out   = model.AddProduct(lib.GetType("ref::gpr_table_denovo"))

CHANNEL_PREFIX = "denovo"

# The same frozen 14-column schema the GEM table uses, so the two lines of evidence line
# up column for column. See benchmark/host_gpr_gem.py.
GPR_COLS = (
    "build_id", "host", "unit_id", "feature_id", "feature_kind", "feature_name",
    "mnxr", "channel", "evidence_id", "evidence_name", "raw_score",
    "projection_via", "in_atom_universe", "gpr_rule",
)

HOST_FOR_ACCESSION = {
    "GCF_000005845.2": "e_coli_k12",
    "GCF_000019425.1": "e_coli_dh10b",
    "GCF_051228345.1": "e_coli_epi300",
}

# The mapper's 8-column long format, which this re-expresses in the host schema.
DRIVER = r'''
import numpy as np
import pandas as pd

HOST = "{host}"
GPR_COLS = {gpr_cols}

g = pd.read_parquet("{gpr}")
print(f"[denovo_gpr] {{len(g):,}} mapper rows, {{g['orf'].nunique():,}} ORFs, "
      f"{{g['mnxr'].nunique():,}} MNXR, channels {{sorted(g['channel'].unique())}}",
      flush=True)

out = pd.DataFrame({{
    "build_id": "denovo_" + HOST,
    "host": HOST,
    # The unit is the proteome, not a model: this table's claim is "this host's own
    # annotation lanes infer these reactions", and naming a GEM here would imply a
    # curated model was consulted, which is the whole thing the de-novo line is not.
    "unit_id": "proteome",
    "feature_id": g["orf"],
    "feature_kind": "orf",
    "feature_name": g["intermediate_name"],
    "mnxr": g["mnxr"],
    # The lane stays in the channel, prefixed, so a row's provenance survives the merge
    # with the GEM table (whose single channel is `gem_gpr`).
    "channel": "{prefix}_" + g["channel"].astype(str),
    "evidence_id": g["intermediate_id"],
    "evidence_name": g["intermediate_name"],
    # Carried through from the lane, unlike the GEM table's uniform 1.0: here the score
    # IS evidence strength, and it is what the condition GPR's belief split reads.
    "raw_score": g["raw_score"].astype(np.float32),
    "projection_via": g["projection_via"],
    # Not computable here without the bake, and NOT defaulted to True: a row wrongly
    # marked in-universe claims an edge can exist for a reaction that has no atom pairs.
    # Null means "not asserted", which a consumer can see.
    "in_atom_universe": pd.Series([None] * len(g), dtype="object"),
    # There is no boolean rule: a de-novo call is per ORF, and inventing "orf" as a
    # one-gene rule would make the two tables look like the same kind of claim.
    "gpr_rule": None,
}})[list(GPR_COLS)]

out = out.sort_values(["feature_kind", "feature_id", "mnxr", "channel"],
                      kind="mergesort", na_position="last").reset_index(drop=True)
out.to_parquet("{out}", index=False, compression="zstd")
print(f"[denovo_gpr] wrote {{len(out):,}} rows for {{HOST}} "
      f"({{out['feature_id'].nunique():,}} ORFs, {{out['mnxr'].nunique():,}} MNXR)",
      flush=True)
'''


def protocol(context: ExecutionContext):
    igpr = context.Input(gpr)
    with open(context.Input(acc).local) as f:
        accession = f.readline().strip()
    host = HOST_FOR_ACCESSION.get(accession)
    if host is None:
        raise SystemExit(
            f"no host name for accession {accession}; add it here and to "
            f"acquire/host_accessions.py, which is where the host set is declared")

    iout = context.Output(out)
    driver = DRIVER.format(gpr=igpr.container, host=host, prefix=CHANNEL_PREFIX,
                           gpr_cols=repr(GPR_COLS), out=iout.container)
    context.LocalShell("cat > _host_gpr_denovo.py << 'PYEOF'\n" + driver + "\nPYEOF\n")
    context.ExecWithEnv() \
        .ifContainerDo(env=image, cmd="python3 _host_gpr_denovo.py") \
        .ifVirtualEnvDo(env=image, cmd="python3 _host_gpr_denovo.py")

    return ExecutionResult(
        manifest=[{out: iout.local}],
        success=iout.local.exists() and iout.local.stat().st_size > 0,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=gpr,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(minutes=30)),
)
