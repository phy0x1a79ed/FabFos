# FabFos

An automated pipeline for resolving inserts from pooled fosmid DNA, rebuilt on
the [metasmith](https://github.com/hallamlab/Metasmith) workflow framework.

FabFos is now a thin front end: it describes the fosmid pipeline as a metasmith
**library** (typed transforms over containers/conda envs) and lets metasmith
plan and execute the workflow.

```bash
fabfos \
    --reads /.../interleaved.fastq.gz --interleaved \
    --background /.../host_background_genome.fasta \
    --vector /.../plasmid_backbone.fasta \
    --endf /.../forward_ends.fasta --endr /.../reverse_ends.fasta \
    --output ./example_out
```

The pipeline resolves: reads → QC/trim (bbduk) → optional host filtering
(minimap2/samtools, forced when `--background` is given) → assembly (megahit) →
non-redundant contigs (length filter + dedup, blast). When `--endf/--endr` are
given, an orthogonal blast step tags those contigs with the fosmid end
mappings; with `--vector`, a pool-size estimate (minimap2/samtools/vsearch) is
also produced. ORF annotation, when requested, runs on the non-redundant
contigs rather than the raw assembly.

## Repository layout

```
setup.py  dev.sh  envs/  conda_recipe/   # project files (src-layout: the package is in src/)
src/
  fabfos/               # THE package: CLI front end + canon.py
  metasmith/            # submodule → Metasmith @ release (the framework)
  metasmith_libraries/  # submodule → MetasmithLibraries @ feat/fabfos (the fosmid transforms)
data_types/             # ecspr:: and ref:: type libraries for the data-dependency library
transforms/{build,run}/ # build = compile data dependencies; run = consume them
containers/<name>/      # only images WE build; public ones are pinned, not vendored
provenance/             # declarations: where every dependency came from + its sha256
_old/                   # the previous Snakemake implementation, kept for reference
```

`src/metasmith_libraries` is pinned to `feat/fabfos`, which is the branch
carrying the `fosmids/` transform domain and the ECSPr prerequisite chain —
`release` has neither. `src/metasmith` is pinned to `release`, the only branch
its remote publishes. Both also carry a secondary `awm` git remote (the sibling
local bare repos) for push-free local sync during co-development; the commits
they are pinned to are on the GitHub remotes, so a plain
`git clone --recurse-submodules` works.

## Data dependencies

Everything the pipeline reads — MetaNetX, UniRef50, KOfam, ESM-C weights, the
frozen MetaNetX-4.5 pre-bake, the direction ensemble, the frozen nulls, the
validation set — lives in one metasmith `DataInstanceLibrary` at
`.awm/data/ref/` (911 items, ~46 GB), tiered `external/`, `external/licensed/`,
`derived/`, `validation/`, `sif/`.

`provenance/data/_declared.yml` is the hand-edited source of truth; the records
beside it are **generated** and carry each item's size, sha256 and per-file
digests. Build and check the library with:

```bash
python transforms/build/_stage/build_ref_library.py --stage --place --index --verify
python transforms/build/_stage/check_provenance.py    # AC3/AC4 + hash cross-checks
python transforms/build/_stage/check_canon.py         # canon resolves through the manifest
```

On a machine where the sources already exist, items are **hardlinked** into the
library, so it costs zero bytes and the originals keep working unchanged.

`src/fabfos/canon.py` is the one code↔data interface: scalars are plain
constants, and every data path resolves through the library manifest lazily via
a module-level `__getattr__`. Importing it for a scalar needs neither the
library nor metasmith. A missing key raises `CanonError` naming the symbol —
there is deliberately no fallback to an absolute path.

## Develop

```bash
git clone --recurse-submodules <repo>
cd <repo>
./dev.sh --ibase                 # create the `fabfos` conda env (python + metasmith)
./dev.sh -r --plan-only \        # resolve the DAG without executing
    -r reads.fq.gz -i -o out --vector backbone.fna
```

`./dev.sh -r` runs the CLI from source against the sibling metasmith +
library submodules. `./dev.sh -b` bundles the library into the package
(`src/fabfos/_library`) for shipping; `-bp`/`-bc` build the wheel / conda
package.

## Status

The core pipeline runs end-to-end on the `fosmids_test` fixtures: reads → QC →
host filter → assembly → non-redundant contigs → pool-size estimate + end tags,
with the metasmith **task cache** promoting each step then re-hitting all of
them on a second run.

The `--ecspr` chain currently stops at `metabolic::fosmid_bipartite_graph`; it
does not reach solve / nulls / significance. Wiring the rest of that chain, and
populating `transforms/{build,run}/`, is the outstanding work.

### Known blockers to running off this machine

- **`external_ecspr` and `external_clean` are not on any registry.** Both exist
  only in a local docker daemon, so the ECSPr numeric core and the CLEAN
  annotation lane have no obtainable runtime elsewhere. See
  `provenance/containers/{ecspr,clean}.yml` for remediation. Five further
  images are unresolvable and recorded with blockers; the other 37 are pinned
  by digest.
- **The licensed BioCyc PGDBs are absent** from the data tree — only the
  MetaCyc flatfiles survive, so the direction ensemble's curated member cannot
  be rebuilt as-is. `provenance/data/biocyc.pgdbs.yml` states what degrades.
- **`references/kegg/` is empty**, so KEGG-keyed rollups are not reproducible
  from what is on disk.
- **Six transforms reference `envs::*.condaenv`**, a format only the unmerged
  mamba-executor branch of metasmith understands; the files carry a bare env
  name with no channels or packages, so nothing is installable from them under
  any runtime. Converting them is a separate project.

See the legacy pipeline in `_old/` for the original Snakemake implementation
and full argument/output reference.
