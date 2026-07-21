"""Build the per-element STAR bipartite graph from the TIER 4 atom-pair universe.

    atom_pairs_tier4.parquet  ->  graph_tier4/mnx_bipartite_{C,N,S,P}.pkl

WHY THIS SCRIPT EXISTS AT ALL
-----------------------------
`notes/TIER4_FREEZE.md` closed by saying the graph rebuild "requires that
harness", because `ecspr_atom_graph.py` is a library with no CLI. That is wrong,
and it is wrong for a traceable reason: the provenance record for the graph
(`provenance/data/_declared.yml`, `mnxref.graph.yml`) named `ecspr_atom_graph.py`
as the producer. It is not. `ecspr_atom_graph.py` has no `main()` and writes no
pickle -- it is the atom-vs-star measurement harness, a different lane entirely.

The actual producer is `fabfos/reference/build_atom_graph.py`, which takes no
arguments and is retargeted by rebinding two module constants on `reference`.
`recovery/build_working_tail.py` already did exactly this for tier 3, in twelve
lines. This is that shim for tier 4, with the producer field corrected in the
declaration alongside it.

WHY IT LIVES IN THE REPO AND NOT BESIDE THE PARQUET
---------------------------------------------------
The tier-1/2/3 shims live in the data tree, untracked, which is how the tier-3
build became unreproducible enough that its own freeze note could misidentify the
producer. The charter's G3 wants these tracked. Output still goes to the data
tree; only the recipe is version-controlled.

WHAT MOVES BETWEEN TIERS, AND WHAT DOES NOT
-------------------------------------------
TOPOLOGY DOES NOT MOVE. The skeleton is fixed at 126,636 nodes / 362,482 edges
and comes from reac_prop participation, AAM-independent. Atom pairs can never
ADD a reaction -- pairs that fall off the skeleton are counted and dropped.

WEIGHTS MOVE, and because every consumer drops `w_X <= 0`, the EFFECTIVE universe
moves with them. Measured live `w_C>0` carbon edges: star 196,018, tier 1 168,939,
tier 3 189,340. Star vs tier 3 share only 172,140 (Jaccard 0.807) -- about 19% of
the carbon edge set turns over. Tier 3 vs tier 4 is far smaller: +223 carbon
reactions on 61,766 (0.36%). So a tier-4 graph is a refinement over tier 3, while
the jump away from the star is the substantive change. Stated here because the
number this build feeds is easy to misattribute to tier 4 when it is mostly
attributable to abandoning the star.

Env: p312 (the skeleton is numpy-2.x pickled). Reads only frozen inputs.
Writes only GRAPH_DIR, which must not already exist -- the library places items
by HARDLINK, so overwriting a staged directory would silently mutate whatever
already points at it.

Usage:  mamba run -n p312 python transforms/build/reference/build_tier4_graph.py
        [--pairs PATH] [--out DIR] [--force]
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

REF_SRC = Path("/home/tony/agentic_workspace/projects/scadc/fabfos/main/fabfos/reference")
DATA_REF = Path("/home/tony/agentic_workspace/data/scadc/ecspr_reference/mnxref-4_5")

# The freeze's own census. A build that does not reproduce these is not tier 4.
TIER4_SHA256 = "7b3b217f91373f3141404767f9be3ad20a87be10d756ccef6a7ecdb1311abe93"
TIER4_PAIRS = 2_455_235
TIER4_RXNS = 63_621


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", default=str(DATA_REF / "atom_pairs_tier4.parquet"))
    ap.add_argument("--out", default=str(DATA_REF / "graph_tier4"))
    ap.add_argument("--force", action="store_true",
                    help="permit writing into an existing output directory")
    a = ap.parse_args()

    pairs, out = Path(a.pairs), Path(a.out)

    if not pairs.exists():
        raise SystemExit(f"tier-4 atom pairs not found: {pairs}")

    # Verify the INPUT is the frozen artifact before spending the build on it.
    # The freeze doc pins this hash; a mismatch means the parquet moved under us
    # and the resulting graph would carry tier 4's name but not its content.
    got = sha256_file(pairs)
    if got != TIER4_SHA256:
        raise SystemExit(
            f"{pairs.name} is not the frozen tier-4 artifact.\n"
            f"  expected sha256 {TIER4_SHA256}\n"
            f"  got            {got}\n"
            "Refusing: a graph built from an unpinned table cannot be called tier 4."
        )
    print(f"[tier4] input verified: {pairs.name} sha256 {got[:12]}...")

    if out.exists() and any(out.iterdir()) and not a.force:
        raise SystemExit(
            f"{out} already exists and is non-empty. The data library places items by "
            "HARDLINK, so overwriting a staged directory mutates every tree that shares "
            "those inodes. Pass --force only if you are certain nothing points here."
        )

    sys.path.insert(0, str(REF_SRC))
    import reference as ref  # noqa: E402

    ref.ATOM_PAIRS = pairs
    ref.GRAPH_DIR = out
    out.mkdir(parents=True, exist_ok=True)

    import build_atom_graph  # noqa: E402

    # Census gate: the builder prints its own row/mnxr counts, but printing is not
    # checking. Assert against the freeze before the graphs are written.
    import pandas as pd  # noqa: E402

    df = pd.read_parquet(pairs, columns=["mnxr"])
    n_pairs, n_rxns = len(df), df["mnxr"].nunique()
    if (n_pairs, n_rxns) != (TIER4_PAIRS, TIER4_RXNS):
        raise SystemExit(
            f"tier-4 census mismatch: got {n_pairs:,} pairs / {n_rxns:,} reactions, "
            f"expected {TIER4_PAIRS:,} / {TIER4_RXNS:,} per notes/TIER4_FREEZE.md"
        )
    print(f"[tier4] census verified: {n_pairs:,} pairs / {n_rxns:,} reactions")

    build_atom_graph.main()

    written = sorted(out.glob("mnx_bipartite_*.pkl"))
    print(f"[tier4] wrote {len(written)} graphs -> {out}")
    for p in written:
        print(f"          {p.name}  {p.stat().st_size:,} bytes  sha256 {sha256_file(p)[:12]}...")
    if len(written) != 4:
        raise SystemExit(f"expected 4 element graphs, got {len(written)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
