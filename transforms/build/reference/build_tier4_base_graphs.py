"""Rebuild the benchmark's four base-graph facets on the TIER 4 universe.

    netA_{iML1515,iECDH10B}/base_{C,N,S,P}.pkl   (curated GEM, uniform E=1.0)
    netB_{iML1515,iECDH10B}/base_{C,N,S,P}.pkl   (4-lane annotation, E_full weighted)

WHY A SHIM AND NOT A FLAG
--------------------------
`60_build_netA_base.py` and `62_build_netB_base.py` live in the metabolic-modelling
tree and read the universe through `_graphs.load_bipartite(X)`, which resolves
`CACHE / mnx_bipartite_{X}.pkl`. `CACHE` is a module constant with no override, and
it is ALSO where they write. So the universe cannot be swapped by argument, and an
in-place rebuild would overwrite the star-derived base graphs.

That overwrite would not be merely untidy. Those files are hardlinked into the data
library as `benchmark.v3.base_graphs`, so writing through them mutates the staged
artifact and every other tree sharing the inode -- silently, since the library
records a sha256 that would no longer describe the bytes.

So this rebinds `_graphs.CACHE` to a private tier-4 directory BEFORE the builders
are imported (their module-level CROSSWALK dict is evaluated at import time against
whatever CACHE is then bound), and mirrors every other cache entry in by symlink so
any input read at runtime that this file did not anticipate still resolves. Only the
universe pickles and the outputs are genuinely new.

WHAT THE OUTPUT MEANS
---------------------
The same four facets over the same reaction sets, re-induced on a universe whose
carbon edge set differs from the star's by 19%. Net A's weight is the raw atom
transit count and Net B's is that count times E_full, so BOTH move with the
universe -- there is no facet that is insulated from this change.

`axes_testable.json` moves with the base graphs, because testability is gated on
which metabolites land in the LCC. That is why the frozen null cannot be reused
across this change: the significance scorer requires observed/null lockstep on
identical testable-axis sets per element.

Env: ml (the builders' own env). Reads frozen inputs; writes only WORK.

Usage:  mamba run -n ml python transforms/build/reference/build_tier4_base_graphs.py
        [--work DIR] [--only netA|netB] [--force]
"""
from __future__ import annotations

import argparse
import os
import pickle
import subprocess
import sys
from pathlib import Path

RN = Path("/home/tony/agentic_workspace/projects/scadc/metabolic-modelling/main"
          "/metabolic-modelling/04_reaction_network")
VALID = RN / "validation"
STAR_CACHE = RN / "cache"
TIER4_GRAPH = Path("/home/tony/agentic_workspace/data/scadc/ecspr_reference"
                   "/mnxref-4_5/graph_tier4")
DEFAULT_WORK = Path("/home/tony/agentic_workspace/data/scadc/ecspr_reference"
                    "/mnxref-4_5/benchmark_tier4")
# Where biomass_edges_set4.tsv actually lives. The builders point at a
# publish/03_model/ecspr_scadc/ path that the figure renumbering retired.
SET4_REAL = Path("/home/tony/agentic_workspace/data/scadc/figures/publish"
                 "/05_scadc/biomass_edges_set4.tsv")

ELEMENTS = ("C", "N", "S", "P")
# Live w_X>0 edge counts of the tier-4 universe, measured at build time. A universe
# that does not match these is not the one this script was written against.
TIER4_LIVE = {"C": 189_846, "N": 121_027, "S": 32_738, "P": 86_833}


def prepare_work(work: Path, force: bool) -> None:
    """A CACHE-shaped directory: tier-4 universe, everything else mirrored in."""
    if work.exists() and any(work.iterdir()) and not force:
        raise SystemExit(
            f"{work} exists and is non-empty. Pass --force to reuse it, but be sure "
            "nothing has hardlinked its outputs into the data library first."
        )
    work.mkdir(parents=True, exist_ok=True)

    # The universe: tier 4, copied by symlink so the source stays authoritative.
    for X in ELEMENTS:
        src = TIER4_GRAPH / f"mnx_bipartite_{X}.pkl"
        if not src.exists():
            raise SystemExit(f"tier-4 universe missing: {src}. Run build_tier4_graph.py first.")
        dst = work / f"mnx_bipartite_{X}.pkl"
        if dst.is_symlink() or dst.exists():
            dst.unlink()
        dst.symlink_to(src)

    # Everything else the builders might reach for, mirrored from the star cache.
    # Deliberately excludes the universe pickles (overridden above) and any prior
    # netA_/netB_ outputs (which are what we are replacing).
    for entry in sorted(STAR_CACHE.iterdir()):
        if entry.name.startswith("mnx_bipartite_"):
            continue
        if entry.name.startswith(("netA_", "netB_")):
            continue
        dst = work / entry.name
        if dst.is_symlink() or dst.exists():
            continue
        dst.symlink_to(entry)

    # netB is stage 62 of a chain: stage 61 writes the annotation evidence into the
    # SAME per-host directory, and 62 reads it back from there. That evidence is
    # universe-INDEPENDENT -- it is lane calls over ORFs, with no graph in its
    # derivation -- so it is seeded in rather than regenerated. Only 62's own
    # outputs are left absent so they are rebuilt on tier 4.
    #
    # The split is not guesswork: in the star cache the base_*.pkl carry nlink=2
    # (hardlinked into the data library as benchmark.v3.base_graphs) while the
    # stage-61 files carry nlink=1. Regenerating what is staged, seeding what is not.
    stage62_outputs = {"evidence_table.parquet", "evidence_weights.parquet",
                       "netB_base_report.md"}
    # Exact host dirs only. A bare netB_* glob also catches sibling variants like
    # netB_iML1515_ref and seeds empty directories for facets nothing will build.
    for host_dir in sorted(STAR_CACHE / f"netB_{m}" for m in ("iML1515", "iECDH10B")):
        if not host_dir.is_dir():
            continue
        dst_dir = work / host_dir.name
        dst_dir.mkdir(exist_ok=True)
        for entry in sorted(host_dir.iterdir()):
            if entry.name.startswith("base_") or entry.name in stage62_outputs:
                continue
            dst = dst_dir / entry.name
            if dst.is_symlink() or dst.exists():
                continue
            dst.symlink_to(entry)
        seeded = sorted(p.name for p in dst_dir.iterdir())
        print(f"[tier4-base] seeded {host_dir.name} stage-61 inputs: {', '.join(seeded)}")

    print(f"[tier4-base] work dir ready: {work}")


def verify_universe(work: Path) -> None:
    for X in ELEMENTS:
        with open(work / f"mnx_bipartite_{X}.pkl", "rb") as fh:
            G = pickle.load(fh)
        wk = f"w_{X}"
        n = sum(1 for _, _, d in G.edges(data=True) if d.get(wk, 0) > 0)
        if n != TIER4_LIVE[X]:
            raise SystemExit(
                f"universe {X} has {n:,} live edges, expected tier-4's {TIER4_LIVE[X]:,}. "
                "The work dir is not pointing at tier 4."
            )
    print(f"[tier4-base] universe verified as tier 4: "
          + "  ".join(f"{X}={TIER4_LIVE[X]:,}" for X in ELEMENTS))


def run_builder(script: str, args: list[str], work: Path) -> None:
    """Run a builder with _graphs.CACHE and its stale SET4 path rebound, in its
    own process.

    A subprocess rather than an in-process monkeypatch because both builders bind
    CACHE-derived constants at import time and 62 additionally exec-loads sibling
    modules that do the same -- one process per (script, model) keeps each import
    graph clean instead of relying on module-cache eviction.

    The module is IMPORTED and then main() is called, rather than runpy'd, so that
    module globals can be corrected between the two. Both builders hardcode

        SET4 = .../figures-main/.awm/data/figures/publish/03_model/ecspr_scadc/
               biomass_edges_set4.tsv

    which no longer exists: the publish tree was renumbered and the file now lives
    at publish/05_scadc/. That is a pre-existing break in the metabolic-modelling
    tree, unrelated to the universe swap, and it is repaired HERE rather than by
    editing another project's checkout or by fabricating a compatibility symlink
    into its data tree.
    """
    shim = (
        "import sys, pathlib, importlib.util\n"
        f"RN = pathlib.Path({str(RN)!r})\n"
        f"VALID = pathlib.Path({str(VALID)!r})\n"
        "sys.path.insert(0, str(RN))\n"
        "sys.path.insert(0, str(VALID))\n"
        "import _graphs\n"
        f"_graphs.CACHE = pathlib.Path({str(work)!r})\n"
        f"SET4 = pathlib.Path({str(SET4_REAL)!r})\n"
        "assert SET4.exists(), f'set4 still missing: {SET4}'\n"  # canon-ok: error text inside an emitted script, not a basis value; SET4_REAL above comes from canon
        f"spec = importlib.util.spec_from_file_location('_builder', str(VALID / {script!r}))\n"
        "mod = importlib.util.module_from_spec(spec)\n"
        "sys.modules['_builder'] = mod\n"
        f"sys.argv = [{script!r}] + {args!r}\n"
        "spec.loader.exec_module(mod)\n"
        "mod.SET4 = SET4\n"
        "mod.main()\n"
    )
    print(f"[tier4-base] === {script} {' '.join(args)} ===", flush=True)
    r = subprocess.run([sys.executable, "-u", "-c", shim], cwd=str(VALID))
    if r.returncode != 0:
        raise SystemExit(f"{script} {' '.join(args)} failed with {r.returncode}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", default=str(DEFAULT_WORK))
    ap.add_argument("--only", choices=["netA", "netB"])
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    work = Path(a.work)
    prepare_work(work, a.force)
    verify_universe(work)

    if a.only != "netB":
        for model in ("iML1515", "iECDH10B"):
            run_builder("60_build_netA_base.py", ["--model", model], work)
    if a.only != "netA":
        for host in ("k12", "dh10b"):
            run_builder("62_build_netB_base.py", ["--host", host], work)

    print(f"\n[tier4-base] outputs under {work}:")
    for d in sorted(work.glob("net[AB]_*")):
        if d.is_symlink():
            continue
        bases = sorted(d.glob("base_*.pkl"))
        print(f"  {d.name}: {len(bases)} base graphs "
              + " ".join(f"{p.name}={p.stat().st_size:,}B" for p in bases))
    return 0


if __name__ == "__main__":
    sys.exit(main())
