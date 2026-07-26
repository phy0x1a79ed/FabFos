"""R6b -- per-reaction directionality, as a continuous ratio.

Three members: eQuilibrator and dGbyG are both TECRDB-fitted and therefore
correlated (floored by a shared-error term so they cannot vote as two independent
predictors); MetaCyc REACTION-DIRECTION is independent.

Two correctness facts this step exists to get right:
  * ORIENTATION. REACTION-DIRECTION is stated in MetaCyc's equation orientation,
    but MNXref re-canonicalises orientation on import, so a naive metacyc->MNXR
    join INVERTS the curated direction on ~60% of reactions. Every curated call is
    re-expressed in MNXR orientation by comparing compound sets; undecidable cases
    are recorded, never guessed.
  * WILDCARD ABSTENTION. dGbyG returns a confident number for R-group wildcards,
    because RDKit gives `*` a valid feature vector. The member abstains on any
    reaction carrying one.

No evidence shrinks toward dG'=0, giving ratio 1.0 -- reversible as a LIMIT, not as
an if-branch, so a reaction the ensemble is silent on is a provable no-op.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image      = model.AddRequirement(lib.GetType("env::equilibrator.env"))
reac_prop  = model.AddRequirement(lib.GetType("raw::metanetx_reac_prop"))
chem_prop  = model.AddRequirement(lib.GetType("raw::metanetx_chem_prop"))
chem_xref  = model.AddRequirement(lib.GetType("raw::metanetx_chem_xref"))
reac_xref  = model.AddRequirement(lib.GetType("raw::metanetx_reac_xref"))
curated    = model.AddRequirement(lib.GetType("raw::metacyc_reactions"))
eq_cache   = model.AddRequirement(lib.GetType("raw::equilibrator_cache"))
# The method, vendored build-side: the two thermo members, the curated aligner, the
# calibration and the combiner, plus the constants they share.
canon_m    = model.AddRequirement(lib.GetType("buildlib::dir_canon.py"))
refdata    = model.AddRequirement(lib.GetType("buildlib::dir_refdata.py"))
flatfile   = model.AddRequirement(lib.GetType("buildlib::dir_metacyc_flatfile.py"))
member_eq  = model.AddRequirement(lib.GetType("buildlib::dir_thermo_eq.py"))
member_db  = model.AddRequirement(lib.GetType("buildlib::dir_thermo_dgbyg.py"))
curated_m  = model.AddRequirement(lib.GetType("buildlib::dir_curated.py"))
calibrate  = model.AddRequirement(lib.GetType("buildlib::dir_calibrate.py"))
combiner   = model.AddRequirement(lib.GetType("buildlib::dir_combine.py"))
annotation = model.AddProduct(lib.GetType("interm::direction_annotation"))

# eQuilibrator's compound cache is passed in rather than left to the package's first-use
# fetch, so the member reads exactly the cache acquire/equilibrator_cache.py pinned.
CACHE_ENV = "EQUILIBRATOR_CACHE_DIR"

# The reaction universe the ensemble scores: every MNXR in reac_prop. Wider than the
# deployed run (which scored the base graph) and deliberately so -- the bake refuses if
# a reaction has atom pairs but no direction row, and scoring only a subset is how that
# gap appears. A reaction the ensemble is silent on lands at ratio 1.0, which is a real
# physical statement rather than an absence.
UNIVERSE = r'''
import json
mnxrs = []
with open("{reac_prop}") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p and p[0].startswith("MNXR"):
            mnxrs.append(p[0])
json.dump(sorted(set(mnxrs)), open("_universe.json", "w"))
print(f"[direction] universe: {{len(set(mnxrs)):,}} reactions", flush=True)
'''

# Evaluate one thermo member over the universe. Kept as its own step per member, mirroring
# the deployed eval_members.py: each member is staged independently and the combiner
# fuses, so a member that fails is a missing vote rather than a failed build.
EVAL = r'''
import json, sys
import pandas as pd

from dir_refdata import load_mnxr_stoich, load_mnxm_props

member_name = sys.argv[1]
if member_name == "eq":
    from dir_thermo_eq import EquilibratorMember as M
else:
    from dir_thermo_dgbyg import DgbygMember as M

mnxrs = json.load(open("_universe.json"))
stoich = load_mnxr_stoich("{reac_prop}")
props = load_mnxm_props("{chem_prop}")
member = M()

rows = []
for i, mnxr in enumerate(mnxrs, 1):
    s = stoich.get(mnxr)
    if s is None:
        rows.append(dict(mnxr=mnxr, dg=None, sigma=None, flag=None, reason="no_stoich"))
        continue
    st, is_bal, is_tr = s
    dg, sig, flag, reason = member.dgr(st, props)
    rows.append(dict(mnxr=mnxr, dg=dg, sigma=sig, flag=flag, reason=reason))
    if i % 2000 == 0:
        print(f"[eval:{{member_name}}] {{i:,}}/{{len(mnxrs):,}}", flush=True)

df = pd.DataFrame(rows)
df.to_parquet(f"_member_{{member_name}}.parquet", index=False)
ok = df["dg"].notna().sum()
print(f"[eval:{{member_name}}] answered {{ok:,}}/{{len(df):,}} ({{ok/len(df):.1%}})", flush=True)
print(df["reason"].value_counts().to_string(), flush=True)
'''


def protocol(context: ExecutionContext):
    irp  = context.Input(reac_prop)
    icp  = context.Input(chem_prop)
    icx  = context.Input(chem_xref)
    irx  = context.Input(reac_xref)
    icur = context.Input(curated)
    ieq  = context.Input(eq_cache)
    ilib = context.Input(combiner)
    iout = context.Output(annotation)
    libdir = ilib.container.parent

    env = f"export {CACHE_ENV}={ieq.container}\nexport XDG_CACHE_HOME={ieq.container}\n"

    universe = UNIVERSE.format(reac_prop=irp.container)
    context.LocalShell("cat > _dir_universe.py << 'PYEOF'\n" + universe + "\nPYEOF\n")
    ev = EVAL.format(reac_prop=irp.container, chem_prop=icp.container)
    context.LocalShell("cat > _dir_eval.py << 'PYEOF'\n" + ev + "\nPYEOF\n")

    context.ExecWithContainer(image=image, cmd=f"python3 _dir_universe.py")

    # T2 -- the curated member, orientation-aligned. This is the load-bearing step: a
    # naive metacyc->MNXR join INVERTS the curated call on ~60% of reactions, because
    # REACTION-DIRECTION is stated in MetaCyc's equation orientation and MNXref
    # re-canonicalises orientation on import. dir_curated re-expresses every call by
    # comparing compound sets, and records undecidable cases rather than guessing.
    context.ExecWithContainer(image=image, cmd=f"""
        {env}
        PYTHONPATH={libdir} python3 {libdir}/dir_curated.py \
            --metacyc-reactions {icur.container} \
            --reac-xref {irx.container} \
            --reac-prop {irp.container} \
            --chem-xref {icx.container} \
            --out _curated_per_mnxr.parquet \
            --out-per-reaction _curated_per_reaction.parquet
    """)

    # The two thermo members. Correlated (both TECRDB-fitted), which is why the combiner
    # floors their fused uncertainty rather than treating them as two independent votes.
    for member in ("eq", "dgbyg"):
        context.ExecWithContainer(image=image, cmd=f"""
            {env}
            PYTHONPATH={libdir} python3 _dir_eval.py {member}
        """)

    # T3 -- calibrate category -> dG' on the eQuilibrator MEASURED arm only. The
    # group-contribution arm returns identically zero for group-conserving chemistry,
    # which is exactly what dominates the REVERSIBLE bin, so including it manufactures a
    # fictitiously tight zero-centred bin.
    context.ExecWithContainer(image=image, cmd=f"""
        {env}
        PYTHONPATH={libdir} python3 {libdir}/dir_calibrate.py \
            --curated _curated_per_mnxr.parquet \
            --reac-prop {irp.container} \
            --chem-prop {icp.container} \
            --out-calibration _calibration.parquet \
            --out-points _calibration_points.parquet
    """)

    # T4 -- fuse. No evidence shrinks toward dG'=0 giving ratio 1.0, so a reaction the
    # ensemble is silent on is a provable no-op rather than an if-branch.
    context.ExecWithContainer(image=image, cmd=f"""
        {env}
        PYTHONPATH={libdir} python3 {libdir}/dir_combine.py \
            --base-mnxrs _universe.json \
            --eq _member_eq.parquet \
            --dgbyg _member_dgbyg.parquet \
            --curated _curated_per_mnxr.parquet \
            --calibration _calibration.parquet \
            --out {iout.container}
    """)

    return ExecutionResult(
        manifest=[{annotation: iout.local}],
        success=iout.local.exists() and iout.local.stat().st_size > 0,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=8, memory=Size.GB(32), duration=Duration(hours=12)),
)
