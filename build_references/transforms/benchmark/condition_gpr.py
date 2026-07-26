"""B3 -- ONE GPR table over every benchmark condition.

Carries the edges each condition adds or deletes, across four cohorts:
gof_native, gof_het, lof, eydallin. Together with the per-host GPR tables it holds
every edge any condition needs, which is what lets each x in X be CONSTRUCTED at
run time from the conditions table instead of being stored as a network.

The per-protein belief weight is the deployed scheme and is not reinvented here:
each protein carries a total belief of 1.0 split equally across the lanes that
annotated it, so a protein resolved by one lane and one resolved by three carry
identical mass -- the contrast downstream measures topology, not how many lanes
happened to fire. Within a lane the share splits across nominations by raw score,
then spreads evenly across each nomination's fanout. Fanout is left uncapped:
dilution is the designed answer to promiscuous EC fanout, and capping is a
scoring-time policy that does not belong in an annotation artifact.

PROVISIONAL contract -- the literature inputs are confirmed while porting
benchmark_v4.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
laser     = model.AddRequirement(lib.GetType("raw::laser_records"))
keio      = model.AddRequirement(lib.GetType("raw::keio_records"))
eydallin  = model.AddRequirement(lib.GetType("raw::eydallin_records"))
het       = model.AddRequirement(lib.GetType("raw::het_screen_records"))
bridge    = model.AddRequirement(lib.GetType("ref::mnxr_lookup"))
weights   = model.AddRequirement(lib.GetType("buildlib::bench_evidence_weights.py"))
cohorts_m = model.AddRequirement(lib.GetType("buildlib::bench_cohorts.py"))
out       = model.AddProduct(lib.GetType("bench::condition_gpr"))

# The four arms. gof_het is the heterologous / metagenomic screening arm and is the only
# one whose proteins are not the host's own, which is why it is the only one that
# projects through the bridge here rather than resolving against a host GPR table.
COHORTS = ("gof_native", "gof_het", "lof", "eydallin")

DRIVER = r'''
import os, sys
import pandas as pd

sys.path.insert(0, os.path.dirname("{lib}"))
from bench_cohorts import load_all_cohorts
from bench_evidence_weights import nomination_contributions, assert_conservation

genes = load_all_cohorts(laser="{laser}", keio="{keio}",
                         eydallin="{eydallin}", het="{het}")
print(f"[cond_gpr] {{len(genes):,}} cohort entries over "
      f"{{genes['cohort'].nunique()}} cohorts", flush=True)
for c, g in genes.groupby("cohort"):
    print(f"           {{c:12s}} {{len(g):5,}} named, "
          f"{{int(g['uniprot'].notna().sum()):5,}} with a UniProt accession", flush=True)

# The bridge's uniprot slice. Only the heterologous arm carries accessions -- a native
# gene is named by symbol and its edges come from that host's own GPR table at network
# construction time, which is why this table and the host tables between them carry every
# edge any condition needs.
b = pd.read_parquet("{bridge}", columns=["id", "id_source", "mnxr"])
up = b[b["id_source"] == "uniprot"][["id", "mnxr"]].drop_duplicates()
up = up.rename(columns={{"id": "uniprot"}})

with_acc = genes[genes["uniprot"].notna()].copy()
ev = with_acc.merge(up, on="uniprot", how="inner")
print(f"[cond_gpr] {{ev['uniprot'].nunique():,}} accessions project to "
      f"{{ev['mnxr'].nunique():,}} MNXR ({{len(ev):,}} rows)", flush=True)

# The deployed belief scheme, imported rather than reimplemented. The normalisation unit
# is the PROTEIN, so `orf` is the protein key; `channel` is the lane, and here there is
# one -- the bridge -- so each protein's 1.0 splits across its own nominations and then
# evenly across each nomination's fanout. Fanout is left uncapped by design.
ev = ev.rename(columns={{"uniprot": "orf"}})
ev["channel"] = "bridge_uniprot"
ev["intermediate_id"] = ev["orf"]
# A bridge row is an assertion of the mapping, not a scored hit, so every nomination
# enters at the same score and the split is by fanout alone. Inventing a score here would
# be inventing evidence strength for a curated crosswalk.
ev["raw_score"] = 1.0

contrib = nomination_contributions(ev)
assert_conservation(contrib, "condition_gpr")

gpr = (contrib.groupby(["orf", "mnxr"], sort=False)
       .agg(E_full=("contrib", "sum"),
            lanes=("channel", lambda s: ";".join(sorted(set(s)))),
            n_lanes=("channel", "nunique"))
       .reset_index().rename(columns={{"orf": "protein_uid"}}))

meta_cols = ["uniprot", "cohort", "condition_id", "gene_label", "action",
             "host_hint", "source_organism"]
meta = (with_acc[[c for c in meta_cols if c in with_acc.columns]]
        .rename(columns={{"uniprot": "protein_uid"}}).drop_duplicates())
gpr = gpr.merge(meta, on="protein_uid", how="left")

# Native-arm rows carry no MNXR of their own: the edge set they add or delete is the
# host's, named by gene symbol and resolved against ref::gpr_table_{{gem,denovo}} when the
# network is constructed. They are EMITTED anyway, with a null mnxr, so the conditions
# table can only ever name an edge set this table knows about -- an absent row would make
# a condition look edge-less rather than host-resolved.
native = genes[genes["uniprot"].isna()].copy()
native = native.rename(columns={{"gene_label": "gene_label"}})
native_rows = pd.DataFrame({{
    "protein_uid": native["gene_label"],
    "mnxr": None,
    "E_full": float("nan"),
    "lanes": "host_gpr",
    "n_lanes": 1,
}})
for c in meta_cols[1:]:
    native_rows[c] = native[c] if c in native.columns else None

allrows = pd.concat([gpr, native_rows], ignore_index=True)
allrows = allrows.sort_values(["cohort", "condition_id", "protein_uid", "mnxr"],
                              na_position="last").reset_index(drop=True)
allrows.to_parquet("{out}", index=False, compression="zstd")
print(f"[cond_gpr] wrote {{len(allrows):,}} rows "
      f"({{int(allrows['mnxr'].notna().sum()):,}} carry an MNXR directly, "
      f"{{int(allrows['mnxr'].isna().sum()):,}} resolve against a host GPR table)",
      flush=True)
'''


def protocol(context: ExecutionContext):
    iout = context.Output(out)
    driver = DRIVER.format(
        lib=context.Input(weights).container,
        laser=context.Input(laser).container,
        keio=context.Input(keio).container,
        eydallin=context.Input(eydallin).container,
        het=context.Input(het).container,
        bridge=context.Input(bridge).container,
        out=iout.container,
    )
    context.LocalShell("cat > _condition_gpr.py << 'PYEOF'\n" + driver + "\nPYEOF\n")
    context.ExecWithContainer(image=image, cmd="python3 _condition_gpr.py")
    return ExecutionResult(
        manifest=[{out: iout.local}],
        success=iout.local.exists() and iout.local.stat().st_size > 0,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=4)),
)
