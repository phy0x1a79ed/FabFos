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

The pipeline resolves: reads → QC/trim (bbduk) → assembly (megahit) →
pool-size estimate (minimap2/samtools/vsearch) → end-aware scaffolding (blast).

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

Under active refactor. The package plans the full core pipeline under the
mamba runtime today; running it self-contained additionally requires the
unified `*.env.yml` tool-environment resources (conda spec + container
pointer) and the per-tool conda envs to be provisioned. See the legacy
pipeline in `_old/` for the original Snakemake implementation and full
argument/output reference.
