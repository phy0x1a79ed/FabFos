#!/bin/bash
# FabFos dev/build automation.
#
# FabFos is a thin metasmith front end: the planner/executor comes from the
# `metasmith` conda package, and the fosmid pipeline definition (transforms,
# data types, container/conda env resources) is the metasmith library bundled
# into the wheel as `fabfos/_library`.
set -e
# this script sits at the repo root; the package is in ./src/fabfos (src-layout)
REPO=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
HERE="$REPO"
LIB_SRC="$REPO/src/metasmith_libraries"
LIB_DST="$REPO/src/fabfos/_library"

case $1 in
    --ibase) # create the dev conda env
        mamba env create --no-default-packages -f "$HERE/envs/base.yml"
    ;;
    -b|--bundle-library) # copy the metasmith library into the package for shipping
        echo "bundling metasmith library: $LIB_SRC -> $LIB_DST"
        rm -rf "$LIB_DST"
        mkdir -p "$LIB_DST"
        # only the pieces the planner loads at runtime
        for sub in data_types resources transforms; do
            cp -r "$LIB_SRC/$sub" "$LIB_DST/$sub"
        done
    ;;
    -bp|--build-pip) # build the wheel/sdist (bundle first)
        "$HERE/dev.sh" --bundle-library
        cd "$HERE"
        rm -rf build dist *.egg-info
        python -m build
    ;;
    -bc|--build-conda) # compile + build the conda package
        "$HERE/dev.sh" --bundle-library
        python "$HERE/conda_recipe/compile_recipe.py"
        "$HERE/conda_recipe/call_build.sh"
    ;;
    -r|--run) # run the CLI from source (dev): ./dev.sh -r --plan-only ...
        shift
        PYTHONPATH="$REPO/src:$REPO/src/metasmith/src" \
        FABFOS_LIBRARY="${FABFOS_LIBRARY:-$LIB_SRC}" \
        python -m fabfos "$@"
    ;;

    ###################################################
    # metasmith_libraries submodule sync (local, no internet).
    #
    # Topology: the metasmith_libraries submodule is a separate repo whose
    # `awm` remote points at the canonical metasmith-libraries/.bare. That bare
    # holds the merger branch `feat/fabfos` ("dev") and the worker branches
    # `feat/fabfos{1,2,3}`, each checked out as a worktree. Each fabfos worktree
    # (dev, dev1..3) pins its submodule to the matching branch. These helpers
    # move library changes between the merger and this worktree over local disk;
    # publishing to GitHub is a separate `git push origin`.

    -sf|--sync-from-dev) # merge merger (feat/fabfos) into THIS worktree's submodule, then advance the gitlink
        DEV_BRANCH="feat/fabfos"
        cur=$(git -C "$LIB_SRC" rev-parse --abbrev-ref HEAD)
        echo "sync-from-dev: merging $DEV_BRANCH -> $cur (submodule)"
        git -C "$LIB_SRC" fetch awm
        if git -C "$LIB_SRC" merge --no-edit "awm/$DEV_BRANCH"; then
            if ! git -C "$REPO" diff --quiet -- src/metasmith_libraries; then
                git -C "$REPO" add src/metasmith_libraries
                git -C "$REPO" commit -m "sync: pull metasmith_libraries from dev ($DEV_BRANCH)"
                echo "  gitlink advanced to $(git -C "$LIB_SRC" rev-parse --short HEAD)"
            else
                echo "  already up to date; gitlink unchanged"
            fi
        else
            echo "  !! merge conflict in submodule -- resolve in $LIB_SRC, commit, then: git -C $REPO add src/metasmith_libraries && git commit" >&2
            exit 1
        fi
    ;;

    -st|--sync-to-dev) # merge THIS worktree's committed submodule changes into the merger (feat/fabfos)
        DEV_BRANCH="feat/fabfos"
        cur=$(git -C "$LIB_SRC" rev-parse --abbrev-ref HEAD)
        if [ "$cur" = "$DEV_BRANCH" ]; then
            echo "sync-to-dev: this worktree is already on $DEV_BRANCH (the merger); nothing to push" >&2
            exit 1
        fi
        if [ -n "$(git -C "$LIB_SRC" status --porcelain)" ]; then
            echo "sync-to-dev: submodule has uncommitted changes -- commit them first" >&2
            exit 1
        fi
        bare=$(git -C "$LIB_SRC" remote get-url awm)
        agg=$(git -C "$bare" worktree list --porcelain \
              | awk -v b="refs/heads/$DEV_BRANCH" '$1=="worktree"{w=$2} $1=="branch"&&$2==b{print w; exit}')
        if [ -z "$agg" ]; then
            echo "sync-to-dev: could not find a worktree of $bare on $DEV_BRANCH" >&2
            exit 1
        fi
        if [ -n "$(git -C "$agg" status --porcelain)" ]; then
            echo "sync-to-dev: merger scope ($agg) is dirty -- clean it first" >&2
            exit 1
        fi
        echo "sync-to-dev: merging $cur -> $DEV_BRANCH in $agg"
        git -C "$agg" fetch "$LIB_SRC" "$cur"
        if git -C "$agg" merge --no-edit FETCH_HEAD; then
            echo "  $DEV_BRANCH now at $(git -C "$agg" rev-parse --short HEAD)"
            echo "  (run './dev.sh -sf' in the dev worktree to update its checkout)"
        else
            echo "  !! merge conflict in $agg -- resolve and commit there" >&2
            exit 1
        fi
    ;;

    *)
        echo "usage: dev.sh [--ibase|-b|-bp|-bc|-r ...|-sf|-st]"
        echo "  -sf|--sync-from-dev  merge merger (feat/fabfos) into this worktree's submodule"
        echo "  -st|--sync-to-dev    merge this worktree's submodule changes into the merger"
    ;;
esac
