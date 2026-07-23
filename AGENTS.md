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
