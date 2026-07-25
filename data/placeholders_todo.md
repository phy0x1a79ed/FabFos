# Placeholders — what is declared but not yet built

**Temporary.** Every path below is a **zero-byte marker** sitting at the exact location its
transform will write to. Zero bytes is deliberate: a zero-byte parquet or tsv fails to parse
loudly, so a placeholder that leaks into a run stops it rather than producing a smaller
answer. Nothing here is data.

Artifact ids (R1…R7, B1…B5, C1) refer to `build_references/REFERENCES.md`, which is the
contract these transforms implement. This file goes away when the list below is empty.

Because each of these paths is (or will be) a DVC chunk, the markers are matched by
`data/.gitignore` and are therefore **invisible to git**. This file is the tracked record
that they exist.

## `reference/`

| path | id | builds from | note |
|---|---|---|---|
| `reference/hosts.fna` | R1 | `raw/hosts/*/genome/*.fna` | concat; cheap, no producer written yet |
| `reference/vector.fna` | R2 | `raw/vector/pcc1.fna` | see *Findings* below — the previous file was not this artifact |
| `reference/kofam/ko_list` (+ `profiles/`) | R3 | `raw/kofam/{ko_list,profiles.tar.gz}` | untar |
| `reference/mnxr_lookup.parquet` | R5 | `raw/metanetx/4.5/{reac_prop,reac_xref}.tsv`, `raw/kegg/ko_to_kegg_r.tsv`, `raw/rhea/rhea2uniprot{,_trembl}.tsv*` | the consolidated bridge; ~30.5M rows, ~82 MB |

Already real, not placeholders: `reference/metabolism/` (11 MB), `reference/proteinBERT/`
(541 MB), `reference/uniref50.dmnd` (24 GB).

## `benchmark/`

| path | id | builds from | note |
|---|---|---|---|
| `benchmark/gpr.parquet` | B3 | `raw/literature/*` + UniProt fetch → lanes → R5 | port of `scadc/fig-screenbench` `benchmark_v4/60_build_gpr.py`; cohorts `gof_native`, `gof_het`, `lof`, `eydallin` |
| `benchmark/hosts/<host>/gpr_gem.parquet` ×3 | B1 | `raw/hosts/*/GEM/*.json` + R5 | epi300 borrows dh10b's model; `check_epi300_identity.py` is what licenses that |
| `benchmark/hosts/<host>/gpr_denovo.parquet` ×3 | B2 | `raw/hosts/*/genome/*.faa` → R3/R4/R5/R7 lanes | heavy compute |

Already real: `benchmark/conditions.tsv` (975 rows, B4), `benchmark/Y/` (3.3 MB, B5).

## `curated/`

| path | id | note |
|---|---|---|
| `curated/ECSPr_axes/biomass_dag_axes_set4.json` | C1 | **hand-authored, 28 KB, no producer and not downloadable.** The 79 biomass source→sink axis definitions. It has to be recovered from the incumbent scadc tree — the pin that used to be here had no content behind it. |

Already real: `curated/benchmark_decisions/*.tsv` (3 files, 92 KB).

## `raw/` — acquisitions not yet made

| path | needed by | source |
|---|---|---|
| `raw/kegg/ko_to_kegg_r.tsv` | R5 | KEGG REST, one `get/<ko>` per KO, ~3 req/s. Classifying this as an acquisition is what retires the dependency on scadc's `kegg_requests.db`. |
| `raw/equilibrator/` | R6 direction | eQuilibrator compound cache. The package fetches it on first use; pin it so a build is not silently dependent on a network call. |
| `raw/literature/laser/` | B3, B4 | LASER (`bitbucket.org/jdwinkler/laser_release`) — an upstream **git checkout**, pinned by url+sha256, never hashed as data bytes. It was referenced by the old `.gitignore` but does not exist in this tree. |
| `raw/literature/eydallin/` | B3 | Eydallin 2010 hits table (`eydallin2010_hits.tsv` in the scadc tree) |
| `raw/literature/het_screen/` | B3 | the heterologous / metagenomic screening cohort's source tables |

Already real: `raw/hosts/` (3 chunks), `raw/vector/`, `raw/kofam/`, `raw/metacyc/` *(licensed)*,
`raw/metanetx/4.5/`, `raw/rhea/`, `raw/uniref/`, `raw/literature/keio/`.

## Findings recorded while reorganizing

**`reference/vector.fna` was not R2.** The file previously at `reference/vector/pcc1.fna`
was 4.5 MB and held **two** records — pCC1fos *and* a 4.7 MB EPI300 assembly contig
(`>C1 from:Epi300_RC_100x_hifiasm.bp.hap1.p.fasta`). That is a combined vector+background
file for the pool-size estimate, not the vector reference. R2 is pCC1fos alone, which is
the 8 KB `raw/vector/pcc1.fna`. Demoted to a placeholder rather than silently shipped as
the vector. The combined file is in the backup if some step turns out to want it.

**`raw/LASER` does not exist** despite the old `.gitignore` carrying a line for it. It is
needed by B3/B4, so it is listed above as an acquisition rather than assumed present.
