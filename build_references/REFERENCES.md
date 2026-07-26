# `data/` — every compiled reference and the raw data it requires

The contract `build_references/` implements. **`build_references/` is now a metasmith
transform instance library**, not a directory of scripts:

```
build_references/
  data_types/     raw:: interm:: bench:: buildlib::   (build-only namespaces)
  resources/
    buildlib/     the ported method -- the two ensembles, the bake encoding, the
                  benchmark cohort readers. Build side only; never in the wheel.
  transforms/
    acquire/      11 -- one per upstream source; nothing here is derived
    compile/       7 -- raw -> the direct refs
    benchmark/     4 -- the evaluation set
  build.sh        compiles _metadata/ for this library AND the shipped one
  check_references.py   the checks over a finished run's results
```

New types a **run-side** transform consumes (`ref::mnxr_lookup`) live in
`src/metasmith_libraries/data_types/ref.yml`, not here — `LoadTypeLibraries` keys a
namespace off the YAML filename stem and raises on a duplicate, so a second `ref.yml`
is not possible. `build.sh` passes both type directories.

**All 22 protocols are implemented.** `examples/build_references_dag.py` remains the
gate — it plans every compiled reference and renders
`tests/artifacts/build_references_dag.svg` (27 steps, every expected transform present),
and what it proves is that the *contracts* compose. Executing the same graph is split in
two at the tier boundary: `tests/build_references_stage1_acquire.py` fills `raw/` (6 steps
on a machine that already holds the DVC-pinned bulk, 11 with `--refetch`) and
`tests/build_references_stage2_curate.py` compiles `reference/` and `benchmark/` from it
(17 steps). Stage 2 loads `acquire/` as well and asserts none of it is scheduled, so
"stage 1 did not deliver" is a named failure rather than an unresolvable plan.

The method behind each table is **ported**, not re-derived — the ensembles, the bake and
the benchmark belief scheme live in `resources/buildlib/` (`buildlib::`, build side only).
Two artifacts are deliberately short of their deployed form and say so where they are
built: B3/B4 read the *extracted* cohort tables rather than re-extracting from the
papers, and `conditions.tsv` leaves the curated target columns empty rather than
inferring them. Both are the open decisions recorded at the end of this file.

**The graph has exactly ONE given**, and it is the MetaCyc drop-in — given only
because it is licensed and cannot be fetched. Everything else is produced, including
which three hosts the references are built for: `acquire/host_accessions.py` declares
the host set so that "the reference build" means one thing rather than whatever
accessions a caller happened to pass. A second given appearing in the render is the
signal that something fetchable is being handed in instead.

Each numbered item below is one artifact produced by one transform. **If a folder
under `data/` is not named in this file, it is deleted** — the backup is at
`projects/fabfos/archive/data/` (hardlinked 2026-07-25, 415 files, zero bytes), and
every DVC-pinned chunk is additionally still pinned on `dev`, `dev1` and `dev3`.

## The tier rule

| tier | rule |
|---|---|
| `raw/` | **acquired only** — a download pinned by url+sha256, or a licensed drop-in. Nothing here is produced by a transform in this repo. |
| `reference/` | **compiled here**, from `raw/` only. These are the direct refs the run pipeline consumes. |
| `benchmark/` | compiled here. The evaluation set. |
| `curated/` | **hand-authored judgement.** No producer, not downloadable. Small text, git-tracked, never DVC bulk. |

Anything derived that is not a named artifact below is **transient** — it lives in a
metasmith work directory for the length of a build and is not a folder in `data/`.
That is what retires `raw/direction/` and `raw/mnxref-4_5/`.

---

## `reference/` — the direct refs

### R1 · host background genome
**Requires:** nothing beyond the host assemblies themselves.
No compile step: the assemblies `acquire/host_genome.py` fetches *are* the
background reference, so R1 is `sequences::isolate_assembly` as a target rather
than a concatenation transform standing between it and the caller.

### R2 · `vector.fna`
**Requires:** nothing. `raw/vector/pcc1.fna` needs no processing at all, so it is not
in the build DAG — a transform that only restates a file is a step that can go wrong
in exchange for nothing. It is a pinned raw file the run pipeline reads directly.

Still worth recording why it is 8 KB: the artifact previously at
`data/reference/vector/` was 4.5 MB and held **two** records, pCC1fos plus a 4.7 MB
EPI300 assembly contig. That is a combined vector+background file for the pool-size
estimate, and shipping it as the vector would have made every vector-vs-insert
comparison quietly also a host comparison.

### R3 · `kofam/` → `ko_list` + `profiles/`
Untarred KOfam HMM profiles + the KO list.
**Requires:** `raw/kofam/{ko_list, profiles.tar.gz}`

### R4 · `uniref50.dmnd`
DIAMOND database.
**Requires:** `raw/uniref/uniref50.fasta.gz`

### R5 · `mnxr_lookup.parquet`
The consolidated bridge — `id, id_source, mnxr, evidence_quality` — replacing the three
separate bridge files. Measured 2026-07-25: 30,467,712 rows after dedup, 82.5 MB against
251.5 MB for the three files it replaces; no id collisions across the three id spaces, so
`id_source` is a label rather than a disambiguator. `evidence_quality` costs 0.2 MB and is
kept because `"reviewed"` sorts before `"unreviewed"`, so the existing `keep="first"` dedup
already retains the stronger claim.

**Requires:**
- `raw/metanetx/4.5/reac_prop.tsv` — EC → MNXR, from the classifs column
- `raw/metanetx/4.5/reac_xref.tsv` — `kegg.reaction:` → MNXR, and Rhea → MNXR
- `raw/kegg/ko_to_kegg_r.tsv` — KO → KEGG reaction, from **one call** to KEGG REST's bulk
  `link/reaction/ko` endpoint (~2 s), which is why it is a `raw/` acquisition and not a
  compiled intermediate. This retires the dependency on scadc's `kegg_requests.db`.
  It is deliberately **not** the deployed method, which crawled `get/<ko>` per KO and was
  only ever run over a cached, host-scoped subset: 12,238 rows over 6,220 KOs against the
  deployed 2,738 over 1,356, verified 99.6% identical per KO on the 3,854 KOs where both
  exist. **The kofam lane therefore reaches ~6.5× the reactions it did in the deployed
  build, so its numbers are not comparable to the archived ones.**
- `raw/rhea/rhea2uniprot.tsv` + `raw/rhea/rhea2uniprot_trembl.tsv.gz` — UniProt → Rhea

*Not* routed KO → EC → MNXR: measured on these three hosts that takes the kofam lane from
646 to 8,548 reactions and makes 91% of them reactions the EC lane already reaches, which
destroys lane independence.

### R6 · `metabolism/{atom_pairs,vocab,direction}.parquet`
One artifact in three files. Each carries the same bake-identity block and
`refs.assert_same_bake` refuses a mismatched trio — reading `atom_pairs` against another
bake's `vocab` decodes every node to the wrong metabolite silently.

**`atom_pairs.parquet`** — the ensemble atom-atom mapping. Three members: RXNMapper and
LocalMapper (neural, correlated, discounted against each other) and MetaCyc curated
(independent, full weight).
**Requires:**
- `raw/metanetx/4.5/reac_prop.tsv` + `chem_prop.tsv` — reaction equations and metabolite
  structures; the neural members' input SMILES are built from these
- `raw/metacyc/atom-mappings-smiles.dat` — the curated member **[LICENSED]**
- `raw/metanetx/4.5/reac_xref.tsv` — the `metacyc.reaction:` → MNXR crosswalk

**`direction.parquet`** — per-reaction `g_rev/g_fwd = exp(ΔG'/RT)`. Three members:
eQuilibrator and dGbyG (both TECRDB-fitted, correlated) and MetaCyc `REACTION-DIRECTION`
(independent).
**Requires:**
- `raw/metanetx/4.5/{reac_prop,chem_prop,chem_xref,reac_xref}.tsv`
- `raw/metacyc/reactions.dat` — the curated member **[LICENSED]**
- `raw/equilibrator/` — eQuilibrator's compound cache. The package fetches it on first
  use; pin it explicitly so a build is not silently dependent on a network call.
- dGbyG weights are **vendored code**, not data — they belong in the library, not `data/`.

**`vocab.parquet`** — int-code vocabulary. Derived from the two above, no additional raw.

### R7 · `proteinBERT/` — the reference label pool  ⚠ CHANGED, see *Open decisions*
Reference embeddings plus their MNXR labels, for the kNN transfer lane.
**Requires:** `raw/uniref/uniref50.fasta.gz` + R5.
No separate labelled-proteome acquisition. The bridge's UniProt rows come from
rhea2uniprot, so a protein with a curated reaction assignment is exactly a protein
the bridge can label — the pool is a *subset* of the reference set, not a new
source beside it, and the `reviewed` slice (365,240 rows against 35,397,466
unreviewed) lands in the same order as the deployed pool's 273,764.

---

## `benchmark/`

Networks are **not stored**. Each x in X is constructed at run time from the conditions
table plus the two GPR tables — which between them carry every edge any condition needs.

### B1 · `hosts/<host>/gpr_gem.parquet`  (3 hosts)
The GPR a curated genome-scale model asserts.
**Requires:** `raw/hosts/e_coli_k12/GEM/iML1515.json`,
`raw/hosts/e_coli_dh10b/GEM/iECDH10B_1368.json`, R5, and
`raw/metanetx/4.5/reac_xref.tsv` for the BiGG → MNXR crosswalk.
EPI300 has no model of its own and borrows DH10B's; `check_epi300_identity.py` is the
measurement that licenses that, and its empty edit list is why the two tables are
legitimately identical.

### B2 · `hosts/<host>/gpr_denovo.parquet`  (3 hosts)  🔲 PLACEHOLDER
The GPR that host's own annotation lanes infer.
**Requires:** `raw/hosts/*/genome/*.faa` → the annotation lanes (R3 kofam, R4 uniref50,
the CLEAN image, R7 pool) → R5.
Heavy compute, deferred. The lanes are the run-side transforms `dev2` already authored
(7 run tools → `gpr_4lane`/`gpr_7lane`), so this is a wiring job, not new method.

### B3 · `gpr.parquet`  🔲 PLACEHOLDER
**One** GPR table over every benchmark condition, carrying the edges each condition adds
or deletes. Ported from `scadc/fig-screenbench` `main/ecspr/benchmark_v4/60_build_gpr.py`,
whose four cohorts are `gof_native`, `gof_het`, `lof`, `eydallin` — `gof_het` is the
heterologous/metagenomic screening arm.
**Requires:** `raw/literature/{laser,keio,eydallin,het_screen}/`, a UniProt fetch for the
heterologous proteins, then the same lanes as B2 → R5.
*Exact source inventory to be confirmed while porting v4 — its `_v4common.py` reaches into
absolute scadc paths that have to be resolved to raw chunks one by one.*

### B4 · `conditions.tsv`
One row per network to build. Names its host, its arm, and which B3 edges it applies.
**Requires:** `raw/literature/*` + `curated/benchmark_decisions/*`

### B5 · `Y/` — `expectations.tsv`, `panel_anchors.tsv`, `panel_edges.tsv`
**Already exists** (`data/benchmark/v3/Y/`). Kept as-is; not rebuilt this pass.

---

## `curated/` — hand-authored, no producer

### C1 · `ECSPr_axes/biomass_dag_axes_set4.json`
The 79 biomass source→sink axis **definitions**. Graph-independent: it is set4 itself, so
it is the same table regardless of which universe the solve runs on. Declared
`producer: hand-curated`, `determinism: frozen` in the incumbent provenance. 28 KB — git
text, not a DVC chunk. Currently a pin (`data/reference/ECSPr_axes.dvc`) with no content;
the file has to be recovered from the incumbent tree.

### C2 · `benchmark_decisions/`
`gene_reaction_overrides.tsv`, `host_lineage.tsv`, `target_resolution.tsv` — curated
judgement calls feeding B4. Currently DVC bulk inside `benchmark/v3/decisions/`; small
enough to be git text.

---

## `raw/` — the acquisition set

Every chunk named by an artifact above, and nothing else.

| chunk | source | note |
|---|---|---|
| `metanetx/4.5/` | MetaNetX 4.5 | 4 tsv; feeds R5, R6, B1 |
| `rhea/` | Rhea FTP | R5 |
| `kofam/` | `ftp.genome.jp/pub/db/kofam` | R3 |
| `uniref/` | UniProt FTP | R4 |
| `kegg/` | KEGG REST | R5; `ko_to_kegg_r.tsv` |
| `metacyc/` | **LICENSED — manual drop-in** | R6 both halves. Verify-and-refuse transform, never a download |
| `equilibrator/` | eQuilibrator compound cache | R6 direction |
| `hosts/` | NCBI assemblies + BiGG models | R1, B1, B2 |
| `vector/` | pCC1 | R2 |
| `literature/` | `laser/`, `keio/`, `eydallin/`, `het_screen/` | B3, B4 |

---

## Deletions

Everything under `data/` not named above.

| deleted | why |
|---|---|
| `raw/AAM/{chebi,bigg,modelseed}` | Declared as external references in the incumbent provenance, but **no code in the tree loads any of their filenames** — verified by grep across `methods/`, `src/` and `transforms/`. `atom_pairs`' declared inputs are `metanetx.reac_prop`, `metanetx.chem_prop`, `metacyc26.flatfiles` — these three are not among them. Carried over wholesale from scadc's `references/`. |
| `raw/AAM/` as a grouping | AAM is a *method*, not a source. It also hid the licensing boundary: `metacyc26_flatfiles` is licensed and the other three are not. `metacyc/` becomes its own chunk. |
| `raw/direction/*.parquet` | Six intermediates of the R6 direction chain — `calibration`, `curated_per_mnxr`, `direction_annotation`, `member_dgbyg_base`, `member_eq_base`, `netA_gem_direction`. Transient by the tier rule. |
| `raw/mnxref-4_5/` | Derived (it *is* a pre-baked R6), and it additionally carries ECSPr **run outputs** — `graph/`, `solve/`, `solve_directed/`, `null/`, `null_directed/`, `ground_probe/`, `worklist*.tsv` — which belong in neither tier. Deleting it means R6 can only be rebuilt through the full ported chain, which is the point. |
| `data/_old/` | Superseded by the archive backup. Exception: `external/pathways_2026/reference_proteins.{faa,tsv}` survives *only* under R7 option A. |
| `reference/functional_annotation/bridges/` | Superseded by R5. |
| `benchmark/v3/{base_graphs,universe}` | Networks are built at run time from B1+B3+B4. |
| `benchmark/v3/{baseline,observations}` | Run outputs, not benchmark definition. |
| `benchmark/v3/{contract,provenance,v1}` | Superseded by B4/B5 and by the library's own index. |
| `reference/_metadata/`, `benchmark/_metadata/` | Byte-identical copies (`cb2abe78…`) of the **old** `.awm/data/ref` library index — 911 items over `external/ derived/ validation/` tiers that describe a layout this tree no longer has. Dropped into two places by the clean-slate commit. |
| `reference/functional_annotation/{kofam,rhea,uniref50}.dvc` | Stale pins naming directories that no longer exist; their content moved to `raw/`. |
| `reference/hosts/` | Emptied already; B1/B2 now land in `benchmark/hosts/`. |

---

## Open decisions

**R7 — the reference label pool. RESOLVED in the contract, open in its consequence.**
The pool is now built from UniRef50 + the bridge rather than from a separate labelled
proteome, which removes an acquisition and the KEGG licensing question with it. What
stays open is that this is *not* the deployed pool: that one is KEGG-derived (54,005
sequences keyed on KEGG gene ids, labelled by KO), so the `pbert_transfer` lane's
numbers will move and must not be reported as a reproduction of the deployed lane.
Porting `benchmark_v4`, which carries its own pool calibration and verification steps,
is where that gets settled.

*The superseded framing, kept because the measurement in it is still the evidence:*
Option A: promote
`_old/external/pathways_2026/reference_proteins.{faa,tsv}` (54,005 KO-labelled KEGG gene
sequences) to `raw/kegg_proteins/` and build the pool here. This genuinely closes the
"cannot be derived here" gap, but produces a *different* pool from the staged scadc one
(273,764 rows), so the `pbert_transfer` lane's numbers move — and the sequences are
KEGG-sourced, the same licensing box as `kegg.requests_db`. Option B: keep the staged
scadc artifact as an unreproducible external import with provenance recorded; numbers stay
comparable to the deployed lane, gap stays open but stops being mysterious. Note that
`benchmark_v4` has its own pool work (`48_pbert_calibrate.py`, `49_pbert_verify_pool.py`,
`49b_pool_chain_check.py`), so porting B3 may settle this on its own.

**`functional_annotation_alt/esmc` (`esmc_600m.tgz`).** Needed only if an ESM-C or EZpred
lane is in the benchmark GPR. `lanes.yml`'s four lanes do not include one; `gpr_7lane` does.
Delete unless the 7-lane mapper is on the benchmark path.

**`benchmark/v3/ground_truth/gof_*.tsv`.** These are the paper tables already extracted into
tabular form — the input to B3/B4 rather than an output. Under the tier rule they are
either a `curated/` artifact (extraction is human judgement) or they get re-extracted from
`raw/literature/` by a transform. v4's cohort structure may supersede them entirely.
