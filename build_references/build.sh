#!/bin/bash
# Compile the build_references transform library's _metadata/ indexes.
#
# Two type directories, not one. LoadTypeLibraries keys a namespace off the YAML
# filename stem and RAISES on a duplicate, so a `ref.yml` here would collide with
# the shipped library's. The split falls out of that:
#
#   src/metasmith_libraries/data_types/   env:: lib:: ncbi:: sequences:: annotation::
#                                         ref:: fosmids::   -- shared, run-side too
#   build_references/data_types/          raw:: interm:: bench::
#                                                            -- build-only
#
# A new type that a RUN-side transform will consume (ref::mnxr_lookup, say) belongs
# in the shipped library, not here.
#
# The shipped library's own transform dirs are rebuilt too: they carry their own
# pruned copy of the type graph under transforms/*/_metadata/types/, so a change to
# ref.yml that they reference does not reach them until they are recompiled.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
REPO="$(dirname "$HERE")"
MLIB="$REPO/src/metasmith_libraries"

if command -v msm >/dev/null 2>&1; then
    msm=(msm)
else
    msm=(python -m metasmith)
fi

args=(build all
      --types "$MLIB/data_types"
      --types "$HERE/data_types")

for d in "$MLIB"/resources/*/;   do args+=(--uniques    "${d%/}"); done
for d in "$MLIB"/transforms/*/;  do args+=(--transforms "${d%/}"); done
for d in "$HERE"/transforms/*/;  do args+=(--transforms "${d%/}"); done

echo "== ${msm[*]} ${args[*]}"
PYTHONPATH="${PYTHONPATH:-}" "${msm[@]}" "${args[@]}"
