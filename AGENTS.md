# FabFos dev — agent brief

FabFos is a thin metasmith front end. The pipeline definition (transforms, data
types, resources) lives in the **`src/metasmith_libraries`** submodule; the wheel
bundles it as `fabfos/_library`.

## Worktree topology

Parallel dev is split across four `git worktree` checkouts of this superproject,
each pinning `src/metasmith_libraries` to its own scope branch:

| worktree | role     | submodule branch |
|----------|----------|------------------|
| `dev`    | merger   | `feat/fabfos`    |
| `dev1`   | worker   | `feat/fabfos1`   |
| `dev2`   | worker   | `feat/fabfos2`   |
| `dev3`   | worker   | `feat/fabfos3`   |

The submodule has two remotes:
- **`origin`** → GitHub (`…/MetasmithLibraries.git`) — the clone source; `.gitmodules`
  uses this so `git clone --recursive` works for a fresh checkout.
- **`awm`** → local `metasmith-libraries/.bare` — the canonical scope, whose
  worktrees hold `feat/fabfos{,1,2,3}`. Used for internet-free syncing.

## Syncing the library submodule (local, no internet)

Changes flow between the merger (`feat/fabfos`, "dev") and a worker over local
disk via the `awm` remote. Run from inside the relevant worktree:

- `./dev.sh -sf` (`--sync-from-dev`) — merge the merger branch into **this**
  worktree's submodule and advance the superproject gitlink. Run in a worker to
  pull merger updates; run in `dev` to refresh it after a `-st`.
- `./dev.sh -st` (`--sync-to-dev`) — merge **this** worker's committed submodule
  changes into the merger scope. Refuses to run from `dev`; both sides must be clean.

Commit submodule changes *before* `-st`. Publishing to GitHub is a separate
`git push origin <branch>` (submodule) / superproject push.

## Build

- `./dev.sh -b` — bundle the library into `src/fabfos/_library`.
- `./dev.sh -bp` / `-bc` — build the pip wheel / conda package (bundles first).
- `./dev.sh -r -- …` — run the CLI from source.

`./dev.sh -b` only **bundles**. It does not regenerate the tracked per-library
`_metadata/` snapshots the planner reads, so editing a `data_types/*.yml` or a
transform needs a `metasmith build` first — otherwise the planner throws
`datatype [X] not found in [ns]` against a stale index.

## The two libraries

There are two transform libraries, and which one a transform belongs in is decided
by *when it runs*, not by what it does:

| library | what it holds | ships in the wheel |
|---|---|---|
| `src/metasmith_libraries/` (submodule) | the **run** side — what a fosmid pipeline executes | yes |
| `build_references/` | the **build** side — acquiring raw data and compiling the references a run consumes | no |

`build_references/build.sh` compiles `_metadata/` for **both**, because they share a
type graph: a type a run-side transform consumes (`ref::mnxr_lookup`) must live in the
submodule's `data_types/ref.yml`, and `LoadTypeLibraries` raises on a duplicate
namespace, so a second `ref.yml` in `build_references/` is not possible. Build-only
namespaces (`raw::`, `interm::`, `bench::`) live in `build_references/data_types/`.

`examples/build_references_dag.py` resolves the whole reference build and renders it
to `tests/artifacts/`. It is the gate on the build library: it fails if any expected
transform drops out of the plan. **The graph must have exactly one given** — the
licensed MetaCyc drop-in. A second given means something fetchable is being handed in
instead of produced.

## Data tiers

`build_references/REFERENCES.md` is the contract: every compiled reference, the raw
data it requires, and why everything else under `data/` was deleted. The rule that
keeps it honest is that `raw/` is **acquired only** — a download or a licensed
drop-in, never a transform's output. Anything derived that is not a named artifact in
that document is transient and lives in a metasmith work directory, not in `data/`.

`data/.gitignore` is the **only** `.gitignore` under `data/`; `find data -name
.gitignore | wc -l` must return 1. Its header documents the two traps (no trailing
slashes, and never ignore a container that holds pins).
