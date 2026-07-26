"""B1 -- the GPR a curated genome-scale model asserts, per host.

Per host (group_by the model), and deliberately a SEPARATE file from the de-novo
table rather than merged with it. The two are independent lines of evidence and the
comparison between them is the point; conflating them into one table destroys the
only thing they are jointly good for.

Needs cobra to open the model, which is also why this cannot run inside the pinned
ecspr image -- that image carries no cobra.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::cobra.env"))
gem       = model.AddRequirement(lib.GetType("raw::host_gem"))
# The accession, so the row's `host` can be named. The model alone cannot say: EPI300 and
# DH10B share iECDH10B_1368, so keying the host off the model id would collapse three
# hosts into two tables and silently drop one.
acc       = model.AddRequirement(lib.GetType("ncbi::assembly_accession"))
bridge    = model.AddRequirement(lib.GetType("ref::mnxr_lookup"))
reac_xref = model.AddRequirement(lib.GetType("raw::metanetx_reac_xref"))
# `in_atom_universe` means "this reaction has atom-pair coverage, so an edge can exist for
# it". That is a fact about the BAKE, so the bake is an input. The contract sketch left
# these out; without them the column could only be guessed, and a guessed
# `in_atom_universe` is worse than an absent one because the consumer trusts it.
pairs     = model.AddRequirement(lib.GetType("ref::atom_pairs"))
vocab     = model.AddRequirement(lib.GetType("ref::metabolism_vocab"))
direction = model.AddRequirement(lib.GetType("ref::direction_ratios"))
build_lib = model.AddRequirement(lib.GetType("lib::ecspr_build.py"))
encoding  = model.AddRequirement(lib.GetType("buildlib::refs_encoding.py"))
out       = model.AddProduct(lib.GetType("ref::gpr_table_gem"))

CHANNEL = "gem_gpr"

# The frozen 14-column host GPR schema, from the pre-library hosts.py. Shared with the
# de-novo table so the two lines of evidence can be compared column for column -- which is
# the only thing they are jointly good for, and the reason they stay separate files.
GPR_COLS = (
    "build_id", "host", "unit_id", "feature_id", "feature_kind", "feature_name",
    "mnxr", "channel", "evidence_id", "evidence_name", "raw_score",
    "projection_via", "in_atom_universe", "gpr_rule",
)

# accession -> the host name the tables are keyed by. The host SET is declared in
# acquire/host_accessions.py; this is only the naming, kept here because the accession is
# all that crosses into this step.
HOST_FOR_ACCESSION = {
    "GCF_000005845.2": "e_coli_k12",
    "GCF_000019425.1": "e_coli_dh10b",
    "GCF_051228345.1": "e_coli_epi300",
}

DRIVER = r'''
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname("{build_lib}"))
sys.path.insert(0, os.path.dirname("{encoding}"))
from ecspr_build import crosswalk_gem, load_model
import refs_encoding as refs

HOST = "{host}"
CHANNEL = "{channel}"
GPR_COLS = {gpr_cols}

# Every MNXR carrying at least one atom transfer, read off the COMPILED table rather than
# a source one. The identity assertion is what stops a stale vocab decoding those codes to
# the wrong reactions -- which would mislabel `in_atom_universe` for every row without
# raising.
ident = refs.assert_same_bake("{vocab}", "{pairs}", "{direction}")
V = refs.load_vocab("{vocab}")
codes = pd.read_parquet("{pairs}", columns=["rxn"])["rxn"].to_numpy()
universe = set(V.symbols("rxn")[np.unique(codes)].tolist())
print(f"[gem_gpr] atom universe: {{len(universe):,}} reactions (bake "
      f"{{ident['vocab_sha256'][:16]}})", flush=True)

model = load_model("{gem}")
xw, stats = crosswalk_gem(model, "{reac_xref}", universe)
gem_id = getattr(model, "id", None) or "unknown"
print(f"[gem_gpr] {{HOST}} / {{gem_id}}: {{stats['n_reactions']:,}} reactions -> "
      f"{{stats['n_resolved']:,}} resolved ({{xw['mnxr'].nunique():,}} distinct MNXR), "
      f"{{stats['n_unresolved']:,}} unresolved, {{stats['n_aam_gap']:,}} outside the "
      f"atom universe", flush=True)

by_rxn = dict(zip(xw["rxn_id"], zip(xw["mnxr"], xw["source"], xw["in_universe"])))
build_id = f"gem_{{gem_id}}"
rows = []
n_ruleless = 0
for r in model.reactions:
    hit = by_rxn.get(r.id)
    if hit is None:
        # No candidate MNXR at all: it can never carry an edge, and a row with a null
        # reaction would put an unjoinable key in the table.
        continue
    mnxr, via, in_universe = hit
    rule = (r.gene_reaction_rule or "").strip()
    genes = list(r.genes)
    common = dict(
        build_id=build_id, host=HOST, unit_id=gem_id, mnxr=mnxr, channel=CHANNEL,
        evidence_id=r.id, evidence_name=r.name or None,
        # A curated model asserts that a reaction is PRESENT, not how much evidence there
        # is for it, so weighting it by anything would be inventing a quantity. The
        # evidence-weighted line is gpr_denovo.
        raw_score=1.0, projection_via=via,
        in_atom_universe=bool(in_universe),
        # CARRIED, not evaluated. A reference table has no perturbation, so evaluating the
        # boolean rule here would bake in one condition -- and the consumer must remember
        # that cobra.GPR.eval takes the KNOCKED-OUT set, not the active one. Handing it the
        # active set inverts the question, and a previous run of that logic reported more
        # reactions live after a knockout than before.
        gpr_rule=rule or None,
    )
    if not genes:
        # Exchanges, diffusion and spontaneous chemistry have no gene to attribute, but
        # they are live in every condition. Dropping them would make every gene set look
        # like starvation.
        n_ruleless += 1
        rows.append(dict(common, feature_id=None, feature_kind="ruleless",
                         feature_name=None))
        continue
    for g in genes:
        rows.append(dict(common, feature_id=g.id, feature_kind="gem_gene",
                         feature_name=g.name or None))

df = pd.DataFrame(rows, columns=list(GPR_COLS))
df["raw_score"] = df["raw_score"].astype(np.float32)
df["in_atom_universe"] = df["in_atom_universe"].astype(bool)
df = df.sort_values(["feature_kind", "feature_id", "mnxr", "evidence_id"],
                    kind="mergesort", na_position="last").reset_index(drop=True)
df.to_parquet("{out}", index=False, compression="zstd")

gene_rows = df[df["feature_kind"] == "gem_gene"]
print(f"[gem_gpr] {{len(df):,}} rows  {{gene_rows['feature_id'].nunique():,}} genes  "
      f"{{df['mnxr'].nunique():,}} MNXR  {{n_ruleless:,}} ruleless  "
      f"{{int(df['in_atom_universe'].sum()):,}} rows in the atom universe", flush=True)
if n_ruleless == 0:
    raise SystemExit("no ruleless reactions emitted -- a genome-scale model always has "
                     "exchanges and spontaneous chemistry, so this means they were "
                     "dropped rather than absent")

# How much of this model's reactome the de-novo lanes could even nominate. Reported, not
# enforced: the two tables are independent lines of evidence and the comparison between
# them is the point, so a gap here is a finding about lane reach rather than a build error.
reachable = set(pd.read_parquet("{bridge}", columns=["mnxr"])["mnxr"].unique())
gem_mnxr = set(df["mnxr"].unique())
print(f"[gem_gpr] {{len(gem_mnxr & reachable):,}} of {{len(gem_mnxr):,}} GEM reactions are "
      f"also reachable through the de-novo bridge "
      f"({{len(gem_mnxr & reachable)/max(1,len(gem_mnxr)):.1%}})", flush=True)
'''


def protocol(context: ExecutionContext):
    igem = context.Input(gem)
    with open(context.Input(acc).local) as f:
        accession = f.readline().strip()
    host = HOST_FOR_ACCESSION.get(accession)
    if host is None:
        raise SystemExit(
            f"no host name for accession {accession}; add it here and to "
            f"acquire/host_accessions.py, which is where the host set is declared")

    iout = context.Output(out)
    driver = DRIVER.format(
        build_lib=context.Input(build_lib).container,
        encoding=context.Input(encoding).container,
        gem=igem.container, reac_xref=context.Input(reac_xref).container,
        vocab=context.Input(vocab).container, pairs=context.Input(pairs).container,
        direction=context.Input(direction).container,
        bridge=context.Input(bridge).container,
        host=host, channel=CHANNEL, gpr_cols=repr(GPR_COLS), out=iout.container,
    )
    context.LocalShell("cat > _host_gpr_gem.py << 'PYEOF'\n" + driver + "\nPYEOF\n")
    context.ExecWithEnv() \
        .ifContainerDo(env=image, cmd="python3 _host_gpr_gem.py") \
        .ifVirtualEnvDo(env=image, cmd="python3 _host_gpr_gem.py")

    return ExecutionResult(
        manifest=[{out: iout.local}],
        success=iout.local.exists() and iout.local.stat().st_size > 0,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=gem,
    labels=["local"],
    resources=Resources(cpus=2, memory=Size.GB(8), duration=Duration(hours=1)),
)
