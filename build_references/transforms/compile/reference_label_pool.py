"""R7 -- the labelled embedding pool the kNN transfer lane votes against.

No separate labelled-proteome acquisition. The pool is built from two things
already in the graph: UniRef50 sequences, and the labels the bridge already
carries. The bridge's UniProt rows come from rhea2uniprot, so a protein with a
curated reaction assignment is exactly a protein the bridge can label -- which
makes the pool a SUBSET of the reference set rather than a new source beside it.

Scale falls out of the `reviewed` slice rather than being chosen: the bridge is
365,240 reviewed against 35,397,466 unreviewed rows, and reviewed is the right cut
for a set that other proteins get classified against. That lands in the same order
as the deployed pool (273,764 rows / 270,454 labelled / 8,590 distinct MNXR).

WHAT THIS CHANGES. The deployed pool is not this. It is KEGG-derived -- 54,005
sequences keyed on KEGG gene ids (`dme:Dmel_CG3481`), labelled by KO, projecting
KO -> MNXR. Building from the Rhea route instead means:
  + no new acquisition, and no KEGG-licensed sequences in the tree
  + one label source instead of two, so a protein cannot be labelled one way here
    and a different way in the GPR mapper
  - a different set, so the pbert lane's numbers move. This is not a reproduction
    of the deployed lane and must not be reported as one.

Two traps worth stating because both are silent:
  * SAME MODEL. A pool embedded with a different model from the query is not a
    weaker pool, it is a meaningless one -- cosine distance between two embedding
    spaces is a number with no referent.
  * ONE PASS. The consumer indexes the embedding stack by row, so pairing an index
    from one build with a stack from another misindexes every row and emits a full,
    confident, WRONG table with nothing raised. Index and stack come out together
    or not at all.

GATED -- see REFERENCES.md "Open decisions".
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image  = model.AddRequirement(lib.GetType("env::proteinbert.env"))
fasta  = model.AddRequirement(lib.GetType("raw::uniref50_fasta"))
bridge = model.AddRequirement(lib.GetType("ref::mnxr_lookup"))
pool   = model.AddProduct(lib.GetType("ref::reference_label_pool"))

# The cut that defines the pool: bridge rows whose id_source is uniprot and whose
# evidence_quality is reviewed, intersected with UniRef50's representatives.
POOL_ID_SOURCE = "uniprot"
POOL_EVIDENCE = "reviewed"

# The layout the consumer reads. functionalAnnotation/gpr_4lane.py's embedding lane opens
# exactly these two names inside the pool directory and indexes the stack by the index's
# `row` column, which is why they are produced together and never separately.
INDEX_NAME = "orf_index.parquet"
STACK_NAME = "emb_pbert.npy"

# The slice step: pick the pool members and write their sequences out as a FASTA for the
# embedder. Split from the embedding step because they run in different environments --
# this half needs pandas and a 12 GB gzip stream, the other half needs ProteinBERT.
SELECT = r'''
import gzip, os, re
import pandas as pd

BRIDGE = "{bridge}"
FASTA  = "{fasta}"
FASTA_OUT = "_pool.faa"
LABELS_OUT = "_pool_labels.parquet"

b = pd.read_parquet(BRIDGE, columns=["id", "id_source", "mnxr", "evidence_quality"])
sel = b[(b["id_source"] == "{id_source}") & (b["evidence_quality"] == "{evidence}")]
labels = (sel.groupby("id")["mnxr"].apply(lambda s: ";".join(sorted(set(s))))
             .rename("mnxr_list").reset_index().rename(columns={{"id": "accession"}}))
print(f"[pool] bridge slice: {{len(sel):,}} rows -> {{len(labels):,}} labelled accessions",
      flush=True)

wanted = dict(zip(labels["accession"], labels["mnxr_list"]))

# UniRef50 names each cluster by its representative, `>UniRef50_P12345 ...`, so the
# bridge's UniProt accession is the cluster id minus the prefix. A protein that is a
# cluster MEMBER but not its representative has no sequence in this file; that is a
# coverage fact and is counted below rather than substituted for.
written = 0
keep = False
with gzip.open(FASTA, "rt") as fh, open(FASTA_OUT, "w") as out:
    for line in fh:
        if line.startswith(">"):
            acc = line[1:].split(None, 1)[0]
            acc = acc[len("UniRef50_"):] if acc.startswith("UniRef50_") else acc
            keep = acc in wanted
            if keep:
                # Header is the bare accession: the embedder echoes it into its index,
                # and that is the key the labels are re-joined on.
                out.write(f">{{acc}}\n")
                written += 1
        elif keep:
            out.write(line)

print(f"[pool] {{written:,}} of {{len(wanted):,}} labelled accessions are UniRef50 "
      f"representatives and have a sequence", flush=True)
if written == 0:
    raise SystemExit("no pool sequences found -- the accession/representative join broke")
labels[labels["accession"].isin(set(wanted))].to_parquet(LABELS_OUT, index=False)
'''

# The assemble step: stitch the embedder's shards into ONE stack and write the index that
# addresses it, in one pass. The consumer indexes the stack BY ROW, so an index from one
# build against a stack from another misindexes every row and emits a full, confident,
# wrong table with nothing raised -- which is why these two files are written together,
# from the same in-memory arrays, or not at all.
ASSEMBLE = r'''
import numpy as np
import pandas as pd
from pathlib import Path

SHARDS = Path("pbert_output")
POOL = Path("{pool}")
POOL.mkdir(parents=True, exist_ok=True)

# The embedder writes one .npy per shard plus a csv index naming the sequences in order.
# Both are read in the SAME sorted shard order, so row i of the stack is sequence i of
# the index by construction rather than by coincidence.
npys = sorted(SHARDS.glob("*.npy"))
csvs = sorted(SHARDS.glob("*.csv"))
if not npys or not csvs:
    raise SystemExit(f"embedder produced no output under {{SHARDS}}")
stack = np.vstack([np.load(f) for f in npys])[:, -512:].astype(np.float32)
idx = pd.concat([pd.read_csv(f) for f in csvs], ignore_index=True)
if len(idx) != len(stack):
    raise SystemExit(f"index has {{len(idx)}} rows but the stack has {{len(stack)}} -- "
                     f"pairing them would misindex every row silently")

labels = pd.read_parquet("_pool_labels.parquet")
idx = idx.rename(columns={{"sequence_id": "orf"}})
idx["row"] = np.arange(len(idx), dtype=np.int64)
idx["role"] = "reference"
merged = idx.merge(labels.rename(columns={{"accession": "orf"}}), on="orf", how="left")
n_unlabelled = int(merged["mnxr_list"].isna().sum())
merged["mnxr_list"] = merged["mnxr_list"].fillna("")

np.save(POOL / "{stack_name}", stack)
merged[["role", "row", "orf", "mnxr_list"]].to_parquet(POOL / "{index_name}", index=False)
print(f"[pool] {{len(merged):,}} reference embeddings, {{n_unlabelled:,}} unlabelled, "
      f"{{merged['mnxr_list'].str.split(';').explode().replace('', None).nunique():,}} "
      f"distinct MNXR", flush=True)
'''


def protocol(context: ExecutionContext):
    ibridge = context.Input(bridge)
    ifasta  = context.Input(fasta)
    ipool   = context.Output(pool)

    select = SELECT.format(bridge=ibridge.container, fasta=ifasta.container,
                           id_source=POOL_ID_SOURCE, evidence=POOL_EVIDENCE)
    context.LocalShell("cat > _pool_select.py << 'PYEOF'\n" + select + "\nPYEOF\n")

    # Both halves run in the ProteinBERT image. It carries numpy/pandas, and running the
    # slice somewhere else would mean staging a 12 GB gzip across two environments.
    context.ExecWithContainer(image=image, cmd="python3 _pool_select.py")

    threads = context.params.get("cpus", 4)
    # SAME MODEL as the query. A pool embedded with a different model from the one the
    # run-side lane uses is not a weaker pool, it is a meaningless one -- cosine distance
    # between two embedding spaces is a number with no referent. That is enforced by
    # sharing `env::proteinbert.env` with functionalAnnotation/proteinbert.py, and the
    # flags below are that transform's, verbatim.
    context.ExecWithContainer(image=image, cmd=f"""
        pbert run -i _pool.faa -o pbert_output \
            --threads {threads} --protein_size 512 --model_batch 1024 -x 1
    """)

    assemble = ASSEMBLE.format(pool=ipool.container, index_name=INDEX_NAME,
                               stack_name=STACK_NAME)
    context.LocalShell("cat > _pool_assemble.py << 'PYEOF'\n" + assemble + "\nPYEOF\n")
    context.ExecWithContainer(image=image, cmd="python3 _pool_assemble.py")

    ok = (ipool.local / INDEX_NAME).exists() and (ipool.local / STACK_NAME).exists()
    return ExecutionResult(
        manifest=[{pool: ipool.local}],
        success=ok,
    )


TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=4, memory=Size.GB(32), duration=Duration(hours=8)),
)
