# FabFos

An automated pipeline for resolving inserts from pooled fosmid DNA, rebuilt on
the [metasmith](https://github.com/hallamlab/Metasmith) workflow framework.

FabFos is now a thin front end: it describes the fosmid pipeline as a metasmith
**library** (typed transforms over containers/conda envs) and lets metasmith
plan and execute the workflow. The default runtime is **mamba**, so a conda
install is self-contained — every tool runs via `mamba run -n <env>` with no
container relay.

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
src/
  metasmith/            # submodule → Metasmith @ feat/fabfos (the framework + mamba executor)
  metasmith_libraries/  # submodule → MetasmithLibraries @ feat/fabfos (the fosmid transforms)
  fabfos/               # this package: CLI front end + conda recipe
_old/                   # the previous Snakemake implementation, kept for reference
```

The two submodules are pinned to their `feat/fabfos` branches and also carry a
secondary `awm` git remote (the sibling local bare repos) for push-free local
sync during co-development.

## Develop

```bash
git clone --recurse-submodules <repo>
cd src/fabfos
./dev.sh --ibase                 # create the `fabfos` conda env (python + metasmith)
./dev.sh -r --plan-only \        # resolve the DAG without executing
    -r reads.fq.gz -i -o out --vector backbone.fna
```

`./dev.sh -r` runs the CLI from source against the sibling metasmith +
library submodules. `./dev.sh -b` bundles the library into the package
(`fabfos/_library`) for shipping; `-bp`/`-bc` build the wheel / conda package.

## Status

The core pipeline runs **end-to-end under the mamba runtime** on the
`fosmids_test` fixtures: reads → QC → host filter → assembly → non-redundant
contigs → pool-size estimate + end tags, all seven steps executing via
`mamba run -n <env>` with no containers, and the metasmith **task cache**
promoting each step then re-hitting all of them on a second run. The tool
environments are provisioned from the unified `*.env.yml` resources (conda
spec + container pointer); FabFos auto-creates any missing ones before a run
(`--provision-only` to just build them, `--no-provision` to skip). The core
tools share the `fabfos-bio` env (diamond in `fabfos-annot`).

The `--ecspr` chain (through the per-fosmid bipartite graph) resolves under
mamba; running it end-to-end additionally needs real ECSPr reference data
(universe bipartite graphs + reaction DB). See the legacy pipeline in `_old/`
for the original Snakemake implementation and full argument/output reference.
