"""Locate the metasmith library that defines the FabFos transforms.

Resolution order:

1. ``$FABFOS_LIBRARY`` — explicit override (a metasmith library root).
2. A copy bundled inside this package at ``fabfos/_library`` — what a conda
   install ships (populated by ``dev.sh`` at build time).
3. The dev sibling submodule ``src/metasmith_libraries`` — when running from
   a source checkout of the FabFos repo.

The library root is the directory that contains ``data_types/``,
``resources/`` and ``transforms/``.
"""
import os
from pathlib import Path

# transform domains the FabFos pipeline draws from; mirrors the library layout.
DOMAINS = [
    "assembly",
    "fosmids",
    "functionalAnnotation",
    "logistics",
    "metabolicModelling",
    "metagenomics",
    "pangenome",
    "responseSurface",
    "transcriptomics",
]

_MODULE = Path(__file__).resolve().parent


def _looks_like_library(root: Path) -> bool:
    return all((root / sub).is_dir() for sub in ("data_types", "resources", "transforms"))


def resolve_library_root() -> Path:
    override = os.environ.get("FABFOS_LIBRARY")
    if override:
        root = Path(override).expanduser().resolve()
        if not _looks_like_library(root):
            raise FileNotFoundError(
                f"FABFOS_LIBRARY=[{root}] is not a metasmith library "
                f"(missing data_types/ resources/ transforms/)"
            )
        return root

    bundled = _MODULE / "_library"
    if _looks_like_library(bundled):
        return bundled

    # src/fabfos/fabfos/library.py -> repo: parents[2] == src/, sibling submodule
    dev_sibling = _MODULE.parents[1] / "metasmith_libraries"
    if _looks_like_library(dev_sibling):
        return dev_sibling

    raise FileNotFoundError(
        "could not locate the FabFos metasmith library. Set FABFOS_LIBRARY, "
        "install the package with a bundled library, or run from a source "
        "checkout with the metasmith_libraries submodule present."
    )
