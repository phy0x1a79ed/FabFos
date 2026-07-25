# proteinBERT — the labelled reference pool for lane-4 embed transfer

What `reference_label_pool` refers to: the thing `pbert_transfer` transfers labels
*from*. ProteinBERT's own weights are baked into the image; what was missing here was
never the model, it was the labelled reference — reference ORF embeddings plus the MNXR
labels attached to them, so a dark query ORF can be given labels by a kNN vote over its
nearest labelled neighbours.

## Origin

Copied verbatim (2026-07-25) from the deployed scadc build:

    projects/scadc/metabolic-modelling/main/metabolic-modelling/05_embed_transfer/cache/

built 2026-06-14 by the pipeline that `src/metasmith_libraries/resources/lib/
fabfos_embed_transfer.py` is the port of — `10_build_reference_pool.py`, whose label
source is `03_layer2_evidence/cache/evidence_table_dlec.parquet`. Not rebuilt here and
not derivable here: the pool needs a reference proteome embedded with the same model as
the query, and this repo holds only query embeddings.

| file | sha256 | size |
|---|---|---|
| `orf_index.parquet` | `532372ed7d4abf1b841f016eef671284d25aaead1c2aadd68ba84c61ba96268d` | 4.6 MB |
| `emb_pbert.npy` | `798ff2361c1e423f9e33fad72c1de21f167f8ba16ce0a2060378bf0dcf38c7e3` | 560.7 MB |
| `proj_pbert.pt` | `a3950c488a358e53e6c622eaeff5dafe5a80994631fd5c1ac26ccf08dd9847f1` | 1.6 MB |

## Contents

`orf_index.parquet` — 273,764 rows; `orf, source, role, mnxr_list, channels, ec_list,
n_channels_max, is_gold, row`. 270,454 `reference` (all labelled, 8,590 distinct MNXR)
and 3,310 `query`. By source: metag 267,500, fosmid 4,635, epi300 1,629.

`emb_pbert.npy` — (273,764, 512) float32. RAW ProteinBERT, which is the variant the
deployed `pbert_transfer` channel applies, at vote floor 0.20. Not the projected stack:
projection is the ESM-C channel's variant.

`proj_pbert.pt` — the SupCon projection head. **Off the `pbert_transfer` path**, kept
only because the type contract says this directory holds projector weights.

## The invariant that makes it one artifact

`emb_pbert.npy` row *i* is `orf_index` row *i*, and `orf_index.row` is exactly
`arange(273764)` — the consumer indexes the embedding stack by that column. **The two
files are one artifact and must be replaced together.** An index from a different build
against this stack raises nothing: every neighbour is simply the wrong protein, and the
lane emits a full, confident, wrong table.

This is not hypothetical. `data/scadc/fabfos_2026/intermediates/embed_transfer/pool/`
holds a *newer* `orf_index.parquet` — 273,849 rows, 4,720 fosmid ORFs against this
build's 4,635 — from a later run. It was **not** used, because that run kept no pooled
`emb_pbert.npy` (only its own fosmid query embeddings), and pairing its index with this
stack would misindex all 273k rows. If that build's embeddings ever surface, take both
files from it or neither.

The `query` rows are scadc's own dark fosmid ORFs. They are kept verbatim rather than
dropped — `Context` filters to `role == "reference"` itself, and renumbering to strip
3,310 rows would break the `row` invariant above for no gain.

## Status: staged, not yet consumed

The artifact is here; nothing reads it yet. `build_references/lanes.yml` still declares
`reference_label_pool` unavailable, and that is now correct for a **different reason**
than when it was written — the file exists, the code does not:

- `load_bridges()` has no loader branch for `reference_label_pool` (it would hit
  `raise SystemExit(f"no loader for bridge ...")`).
- the `embed_transfer` reader in `build_host_denovo_gpr.py` is `raise NotImplementedError`.

So `pbert_transfer` stays absent from all three hosts' `gpr_denovo.parquet`, reported
unavailable by name. Flipping `available: true` before writing both would fail the build,
which is the intended behaviour, not a bug to route around.
