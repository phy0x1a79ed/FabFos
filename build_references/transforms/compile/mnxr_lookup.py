"""R5 -- the consolidated intermediate-id -> MNXR bridge.

One table replacing the former ko/ec/uniprot trio: `id, id_source, mnxr,
evidence_quality`. Measured before consolidating -- the three id spaces share no
ids, so a single `id` column is unambiguous and `id_source` is a label rather than
a disambiguator; deduping to distinct (id, mnxr) went 35,762,706 -> 30,467,712 rows
and is behaviour-preserving because the consumer already ran that dedup after its
join, and the duplicate rows differed only in columns nothing downstream reads.
82.5 MB against 251.5 MB across the three files it replaces.

`evidence_quality` is carried at a cost of 0.2 MB because "reviewed" sorts before
"unreviewed", so an existing keep-first dedup already retains the stronger claim.

Three routes into one table:
  ec       reac_prop classifs
  ko       ko_to_kegg_r -> reac_xref `kegg.reaction:` rows
  uniprot  rhea2uniprot{,_trembl} -> reac_xref `rhea:` rows
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image      = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
reac_prop  = model.AddRequirement(lib.GetType("raw::metanetx_reac_prop"))
reac_xref  = model.AddRequirement(lib.GetType("raw::metanetx_reac_xref"))
ko_kegg_r  = model.AddRequirement(lib.GetType("raw::kegg_ko_to_kegg_r"))
rhea_sp    = model.AddRequirement(lib.GetType("raw::rhea2uniprot"))
rhea_tr    = model.AddRequirement(lib.GetType("raw::rhea2uniprot_trembl"))
# The EC route and the UniProt route already have builders in the shipped evidence
# library, and they are CALLED here rather than reimplemented: the run-side mapper reads
# this table through the same module, so a second copy of the parser is a second thing to
# keep in step with reac_prop's column layout.
ev_lib     = model.AddRequirement(lib.GetType("lib::fabfos_evidence.py"))
bridge     = model.AddProduct(lib.GetType("ref::mnxr_lookup"))


DRIVER = r'''
import os, sys
import pandas as pd

sys.path.insert(0, os.path.dirname("{ev_lib}"))
import fabfos_evidence as fe

REAC_PROP = "{reac_prop}"
REAC_XREF = "{reac_xref}"
KO_KEGG_R = "{ko_kegg_r}"
RHEA_SP   = "{rhea_sp}"
RHEA_TR   = "{rhea_tr}"
OUT       = "{out}"

COLUMNS = ["id", "id_source", "mnxr", "evidence_quality"]


def kegg_r_to_mnxr():
    """`kegg.reaction:R#####` rows of reac_xref. Both prefixes, as the deployed builder
    accepts both."""
    rows = []
    with open(REAC_XREF) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            src, mnxr = parts[0], parts[1]
            if not mnxr.startswith("MNXR") or mnxr == "EMPTY":
                continue
            if src.startswith("kegg.reaction:") or src.startswith("keggR:"):
                rows.append((src.split(":", 1)[1], mnxr))
    return pd.DataFrame(rows, columns=["kegg_r", "mnxr"]).drop_duplicates()


def route_ec():
    """reac_prop classifs. Partial ECs ("3.6.3.-") are kept verbatim because that is how
    reac_prop writes them and the CLEAN lane can emit them too."""
    df = fe.load_ec_to_mnxr(REAC_PROP)
    df = df.rename(columns={{"ec": "id"}})
    df["id_source"] = "ec"
    # MetaNetX's classifs column is a curated assignment with no reviewed/unreviewed
    # split, so the whole route is `reviewed`. Marking it `unreviewed` would let a TrEMBL
    # row outrank a curated EC in the keep-first dedup.
    df["evidence_quality"] = "reviewed"
    return df[COLUMNS]


def route_ko():
    """KO -> KEGG reaction -> MNXR.

    NOT KO -> EC -> MNXR. Measured on the three hosts, routing through the ko_list's
    `[EC:...]` tag takes the kofam lane from 646 to 8,548 reactions and makes 91% of them
    reactions the EC lane already reaches -- the two lanes stop being independent
    evidence, which is the entire reason for having four of them.
    """
    ko_r = pd.read_csv(KO_KEGG_R, sep="\t", dtype=str).dropna()
    r_mnxr = kegg_r_to_mnxr()
    df = (ko_r.merge(r_mnxr, on="kegg_r", how="inner")[["ko", "mnxr"]]
              .drop_duplicates().rename(columns={{"ko": "id"}}))
    n_unmapped = int((~ko_r["kegg_r"].isin(set(r_mnxr["kegg_r"]))).sum())
    print(f"[bridge] ko: {{len(df):,}} rows; reac_xref carries "
          f"{{r_mnxr['kegg_r'].nunique():,}} KEGG reactions, {{n_unmapped:,}} KO links "
          f"point at one it does not have", flush=True)
    df["id_source"] = "ko"
    df["evidence_quality"] = "reviewed"      # KEGG's REACTION block is curated
    return df[COLUMNS]


def route_uniprot(tmp):
    """rhea2uniprot{{,_trembl}} -> reac_xref `rhea:` rows, via the shipped builder.

    TrEMBL is ~10x SwissProt and is what gives the lane its reach on non-model ORFs, so
    it carries `unreviewed` and SwissProt carries `reviewed` -- which is the column the
    label pool later cuts on, and the reason the keep-first dedup below retains the
    stronger claim rather than an arbitrary one.
    """
    fe.build_uniprot_bridge(REAC_XREF, RHEA_SP, RHEA_TR, tmp)
    d = pd.read_parquet(tmp, columns=["uniprot_accession", "mnxr", "evidence_quality"])
    d = d.rename(columns={{"uniprot_accession": "id"}})
    d["id_source"] = "uniprot"
    return d[COLUMNS]


def main():
    frames = [route_ec(), route_ko(), route_uniprot("_uniprot_bridge.parquet")]
    for f, name in zip(frames, ("ec", "ko", "uniprot")):
        print(f"[bridge] {{name:8s}} {{len(f):>12,}} rows  {{f['id'].nunique():>9,}} ids  "
              f"{{f['mnxr'].nunique():>7,}} MNXR", flush=True)
    df = pd.concat(frames, ignore_index=True)
    n0 = len(df)

    # The three id spaces share NO ids (measured before consolidating), so `id` alone is
    # unambiguous and `id_source` is a label rather than a disambiguator. Asserted rather
    # than assumed: if a future source collides, a bare `id` join would silently mix two
    # namespaces' claims into one lane.
    per_id = df.groupby("id")["id_source"].nunique()
    clashes = per_id[per_id > 1]
    if len(clashes):
        raise SystemExit(
            f"{{len(clashes):,}} ids appear in more than one id_source "
            f"(e.g. {{list(clashes.index[:5])}}). `id` is no longer unambiguous, so the "
            f"consumer's per-lane slice would mix namespaces -- this needs a decision, "
            f"not a dedup.")

    # Sort before the dedup so "reviewed" < "unreviewed" alphabetically and keep-first
    # retains the stronger claim. This is the whole reason evidence_quality is carried at
    # all (0.2 MB); dropping it would make the surviving row arbitrary.
    df = (df.sort_values(["id", "mnxr", "evidence_quality"], kind="mergesort")
            .drop_duplicates(subset=["id", "mnxr"], keep="first")
            .reset_index(drop=True))
    df.to_parquet(OUT, index=False)
    print(f"[bridge] {{n0:,}} -> {{len(df):,}} rows after dedup to distinct (id, mnxr)",
          flush=True)
    print(f"[bridge] evidence_quality {{df['evidence_quality'].value_counts().to_dict()}}",
          flush=True)
    print(f"[bridge] wrote {{OUT}} ({{os.path.getsize(OUT)/1e6:.1f}} MB)", flush=True)


main()
'''


def protocol(context: ExecutionContext):
    irp  = context.Input(reac_prop)
    irx  = context.Input(reac_xref)
    ikr  = context.Input(ko_kegg_r)
    isp  = context.Input(rhea_sp)
    itr  = context.Input(rhea_tr)
    iev  = context.Input(ev_lib)
    iout = context.Output(bridge)

    driver = DRIVER.format(
        ev_lib=iev.container, reac_prop=irp.container, reac_xref=irx.container,
        ko_kegg_r=ikr.container, rhea_sp=isp.container, rhea_tr=itr.container,
        out=iout.container,
    )
    context.LocalShell("cat > _mnxr_lookup.py << 'PYEOF'\n" + driver + "\nPYEOF\n")
    context.ExecWithEnv() \
        .ifContainerDo(env=image, cmd="python3 _mnxr_lookup.py") \
        .ifVirtualEnvDo(env=image, cmd="python3 _mnxr_lookup.py")

    return ExecutionResult(
        manifest=[{bridge: iout.local}],
        success=iout.local.exists() and iout.local.stat().st_size > 0,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=2, memory=Size.GB(32), duration=Duration(hours=2)),
)
