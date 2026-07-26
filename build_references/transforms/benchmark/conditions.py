"""B4 -- one row per network to build.

Names each condition's host, arm and tier, and which condition_gpr edges it
applies. This is the table the run-time network construction reads, and it is what
makes X derivable rather than stored.

Takes the condition GPR as input so a condition can never name an edge set that
does not exist -- the alternative, building the two independently and joining them
later, produces a conditions table that looks complete and silently drops rows at
scoring time.

PROVISIONAL contract, as B3.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
cond_gpr  = model.AddRequirement(lib.GetType("bench::condition_gpr"))
laser     = model.AddRequirement(lib.GetType("raw::laser_records"))
keio      = model.AddRequirement(lib.GetType("raw::keio_records"))
eydallin  = model.AddRequirement(lib.GetType("raw::eydallin_records"))
het       = model.AddRequirement(lib.GetType("raw::het_screen_records"))
cohorts_m = model.AddRequirement(lib.GetType("buildlib::bench_cohorts.py"))
out       = model.AddProduct(lib.GetType("bench::conditions"))

# The deployed 20-column schema, from the incumbent conditions.tsv. Frozen here because
# the run-time network construction reads it positionally by name; a column added or
# renamed silently changes what a condition means.
COLUMNS = (
    "condition_id", "arm", "tier", "ptype", "host", "gem", "n_units", "is_control",
    "citation", "note", "gene", "element", "target_mnxm", "target_name",
    "target_basis", "expected_dir", "essential_on_glucose_minimal", "is_neg",
    "obs_id", "host_gem_is_proxy",
)

# One row per (observation, element): a condition is scored per element, because the
# atom-resolved network is built per element and a condition's expected direction is an
# element-specific claim.
ELEMENTS = ("C", "N", "P", "S")

DRIVER = r'''
import os, sys
import pandas as pd

sys.path.insert(0, os.path.dirname("{lib}"))
from bench_cohorts import load_all_cohorts

COLUMNS = {columns}
ELEMENTS = {elements}

cond_gpr = pd.read_parquet("{cond_gpr}")
genes = load_all_cohorts(laser="{laser}", keio="{keio}",
                         eydallin="{eydallin}", het="{het}")

# The join that makes X derivable rather than stored. A condition may only be emitted if
# the condition GPR knows the edge set it names -- building the two independently and
# joining later produces a conditions table that looks complete and silently drops rows at
# scoring time.
known = set(cond_gpr["condition_id"].dropna().unique())
named = set(genes["condition_id"].dropna().unique())
missing = named - known
if missing:
    print(f"[conditions] {{len(missing):,}} observations name no edge in the condition "
          f"GPR and are dropped (e.g. {{sorted(missing)[:3]}})", flush=True)
genes = genes[genes["condition_id"].isin(known)]

# One row per observation: its gene set, its arm, and the host it was run in.
per_obs = (genes.groupby(["condition_id", "cohort"], sort=False)
           .agg(gene=("gene_label", lambda s: "|".join(sorted(set(s)))),
                n_units=("gene_label", "nunique"),
                host=("host_hint", lambda s: next((x for x in s if x), "")))
           .reset_index())

ARM = {{"gof_native": "gof", "gof_het": "gof", "lof": "lof", "eydallin": "lof"}}

rows = []
for r in per_obs.itertuples(index=False):
    edges = cond_gpr[cond_gpr["condition_id"] == r.condition_id]
    for el in ELEMENTS:
        rows.append({{
            "condition_id": f"{{r.condition_id}}__{{el}}",
            "arm": ARM.get(r.cohort, r.cohort),
            "tier": "scored",
            # `ptype` is the shape of the perturbation, read off the actions the cohort
            # recorded rather than assumed from the arm: a LASER record can add and delete
            # in the same observation.
            "ptype": "+".join(sorted({{a for a in
                                      str(edges["action"].dropna().unique().tolist())
                                      .replace("[", "").replace("]", "").replace("'", "")
                                      .split(",") if a.strip()}})) or "unknown",
            "host": r.host,
            "gem": "",
            "n_units": int(r.n_units),
            "is_control": "",
            "citation": "",
            "note": "",
            "gene": r.gene,
            "element": el,
            # The target metabolites are a curated resolution (curated/benchmark_decisions/
            # target_resolution.tsv in the incumbent tree) and are NOT invented here. Left
            # empty rather than guessed: a wrong target silently scores the wrong axis.
            "target_mnxm": "",
            "target_name": "",
            "target_basis": "",
            "expected_dir": "",
            "essential_on_glucose_minimal": "",
            "is_neg": "",
            "obs_id": r.condition_id,
            "host_gem_is_proxy": 0,
        }})

df = pd.DataFrame(rows, columns=list(COLUMNS))
df = df.sort_values(["condition_id"]).reset_index(drop=True)
df.to_csv("{out}", sep="\t", index=False)
print(f"[conditions] {{len(df):,}} rows over {{df['obs_id'].nunique():,}} observations "
      f"x {{len(ELEMENTS)}} elements; arms {{df['arm'].value_counts().to_dict()}}",
      flush=True)
print("[conditions] NOTE target_mnxm/target_name/target_basis/expected_dir are EMPTY: "
      "they come from curated/benchmark_decisions/target_resolution.tsv, which is a "
      "hand-authored input this transform does not declare. See REFERENCES.md C2.",
      flush=True)
'''


def protocol(context: ExecutionContext):
    iout = context.Output(out)
    driver = DRIVER.format(
        lib=context.Input(cohorts_m).container,
        cond_gpr=context.Input(cond_gpr).container,
        laser=context.Input(laser).container,
        keio=context.Input(keio).container,
        eydallin=context.Input(eydallin).container,
        het=context.Input(het).container,
        columns=repr(COLUMNS), elements=repr(ELEMENTS),
        out=iout.container,
    )
    context.LocalShell("cat > _conditions.py << 'PYEOF'\n" + driver + "\nPYEOF\n")
    context.ExecWithContainer(image=image, cmd="python3 _conditions.py")
    return ExecutionResult(
        manifest=[{out: iout.local}],
        success=iout.local.exists() and iout.local.stat().st_size > 0,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=2, memory=Size.GB(16), duration=Duration(hours=1)),
)
