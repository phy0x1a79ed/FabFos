"""Stage 2 of the reference build: everything derived here from what stage 1 landed.

    PATH="/home/tony/lib/miniforge3/envs/msm/bin:$PATH" \\
        python tests/build_references_stage2_curate.py                  # plan + DAG
    ... python tests/build_references_stage2_curate.py --run
    ... python tests/build_references_stage2_curate.py --publish

The compute half. Every byte it reads is already on disk under `data/raw/`; every byte it
writes lands in `data/reference/` or `data/benchmark/`. Nothing here touches the network,
and the check below is what proves it rather than what asserts it.

WHAT IT PRODUCES, by artifact id in build_references/REFERENCES.md:

  R1  the host assemblies ARE the background reference -- no transform stands between
      them and the caller, so this target is satisfied by stage 1's output directly
  R3  kofam profiles + ko_list, untarred                          (compile/kofam_ref)
  R4  uniref50.dmnd, ~25 GB                                       (compile/uniref50_dmnd)
  R5  mnxr_lookup, ~30.5M rows over three id spaces               (compile/mnxr_lookup)
  R6  the metabolism trio -- atom_pairs, vocab, direction -- out of ONE transform, so
      they carry one bake identity. The two ensembles feed it:
        aam_ensemble        RXNMapper + LocalMapper (correlated) + MetaCyc (independent)
        direction_ensemble  eQuilibrator + dGbyG (correlated)    + MetaCyc (independent)
      This is the cost centre: both neural mappers over ~57.5k reactions, 24 h declared.
  R7  the ProteinBERT reference label pool                (compile/reference_label_pool)
  B1  per-host GPR from the curated GEMs                  (benchmark/host_gpr_gem)
  B2  per-host GPR from de novo annotation -- the four shipped lanes (kofamscan, CLEAN,
      diamond_uniref50, proteinbert) into gpr_4lane        (benchmark/host_gpr_denovo)
  B3  the condition GPR                                    (benchmark/condition_gpr)
  B4  the conditions table                                 (benchmark/conditions)

THE CHECK THAT MAKES THE SPLIT REAL. This stage loads acquire/ as well, and then asserts
no acquire transform appears in the plan. That is backwards from how a split is usually
enforced and it is deliberate: hiding acquire/ would make an unmet input an unresolvable
plan whose error names a TYPE, and you would be left guessing which fetch was supposed to
have covered it. Loading it means an unmet input shows up as the exact acquisition that
did not happen, in the plan and in the DAG, by name.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_references_common import (                                    # noqa: E402
    DATA, SCRATCH, add_common_args, execute, generate, latest_results, make_standins,
    new_inputs, print_plan, publish_by_type, render, stage_raw,
)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "examples"))
# One target list, not two: the plan-only gate owns it, and this executes the same graph.
from build_references_dag import TARGETS as _TARGETS                     # noqa: E402

TARGETS = [dtype for _id, dtype in _TARGETS]

# The acquisitions. Any of these in stage 2's plan means stage 1 did not deliver.
#
# host_accessions is deliberately NOT here, even though it lives in acquire/. It fetches
# nothing: it writes the three accession strings that are the host set, as a fact in the
# library rather than an argument a caller supplies. Both GPR tables need the accession to
# tag their host, and an accession is not an artifact anyone keeps in data/raw/, so it is
# re-emitted here every time -- free, offline and deterministic. The check is about
# network and licence, not about which directory a transform happens to sit in.
ACQUIRE = {
    "metanetx", "rhea", "kofam", "uniref50", "kegg_ko_reactions",
    "equilibrator_cache", "literature", "literature_extractions", "metacyc_licensed",
    "host_genome", "host_gem",
}

EXPECTED = {
    "kofam_ref", "uniref50_dmnd", "mnxr_lookup",
    "aam_ensemble", "direction_ensemble", "bake_metabolism",
    "reference_label_pool",
    "kofamscan", "clean", "diamond_uniref50", "proteinbert", "gpr_4lane",
    "host_gpr_gem", "host_gpr_denovo", "condition_gpr", "conditions",
}

# Where each compiled artifact lands. R1 is absent on purpose -- REFERENCES.md is explicit
# that the host assemblies are the background reference themselves, so publication
# concatenates the three staged assemblies rather than a transform restating them.
PUBLISH = {
    "ref::kofamscan_profiles":   "reference/kofam/profiles",
    "ref::kofamscan_ko_list":    "reference/kofam/ko_list",
    "ref::uniref50_diamond_db":  "reference/uniref50.dmnd",
    "ref::mnxr_lookup":          "reference/mnxr_lookup.parquet",
    "ref::atom_pairs":           "reference/metabolism/atom_pairs.parquet",
    "ref::metabolism_vocab":     "reference/metabolism/vocab.parquet",
    "ref::direction_ratios":     "reference/metabolism/direction.parquet",
    "ref::reference_label_pool": "reference/proteinBERT",
    "bench::condition_gpr":      "benchmark/gpr.parquet",
    "bench::conditions":         "benchmark/conditions.tsv",
}

DVC_CHUNKS = ("reference/metabolism", "reference/proteinBERT", "reference/uniref50.dmnd",
              "reference/kofam", "reference/mnxr_lookup.parquet", "benchmark/hosts")

WORK = SCRATCH / "stage2_curate"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", type=Path, default=WORK)
    ap.add_argument("--allow-acquisition", action="store_true",
                    help="plan even if stage 1 left something unacquired. The graph then "
                         "fetches inside stage 2, which is the thing the split exists to "
                         "prevent -- so it is a flag, not a fallback")
    ap.add_argument("--standins", action="store_true",
                    help="fill anything stage 1 has not produced yet with empty stand-ins, "
                         "so the pure curation DAG can be rendered and inspected before "
                         "stage 1 has run. Planning never reads bytes; RUNNING would, so "
                         "this refuses to combine with --run")
    add_common_args(ap)
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)

    if a.publish or a.publish_dry_run:
        print("== stage 2 publish -> data/reference/, data/benchmark/ ==")
        n = publish_by_type(latest_results(a.work), PUBLISH, DATA, dry_run=a.publish_dry_run)
        print("\nDVC chunks whose directory hash this changes -- run these yourself:")
        for chunk in DVC_CHUNKS:
            print(f"    dvc add data/{chunk}")
        print("\ndata/.gitignore is the only .gitignore under data/; check that still holds:")
        print("    [ $(find data -name .gitignore | wc -l) -eq 1 ]")
        return 0 if n else 1

    if a.standins and a.run:
        raise SystemExit(
            "--standins with --run would hand a transform an empty file, and most of "
            "them would produce an empty output and SUCCEED. Plan-only.")

    inputs = new_inputs(a.work)
    n, missing, staged_types = stage_raw(inputs)
    print(f"inputs: {n + 1} staged item(s) from data/raw/")
    if a.standins and missing:
        for p, dtype in make_standins(missing, a.work / "standins"):
            inputs.AddItem(p, dtype)
        print(f"--standins: {len(missing)} product type(s) filled with empty stand-ins, "
              f"so what renders below is the curation graph stage 1 makes possible")
        staged_types |= {d for d, _ in missing}
        missing = []

    # R1 is the reason this filter is not just defensive. REFERENCES.md is explicit that
    # the host assemblies ARE the background reference -- no transform stands between them
    # and the caller -- so `sequences::isolate_assembly` is a stage 1 product that stage 2
    # only passes through. Demanding it as a target here asks the planner to build
    # something that is already sitting in the input library.
    wanted = [t for t in TARGETS if t not in staged_types]
    passthrough = [t for t in TARGETS if t in staged_types]
    if passthrough:
        print(f"\nsatisfied by stage 1 directly, so not built here: {passthrough}")

    agent, task = generate(a.work, inputs, wanted, with_curation=True)
    if not task.ok:
        print(f"FAILED to plan:\n{task.plan}")
        return 1

    used = print_plan(task)
    render(task, "build_references_stage2_curate")

    leaked = used & ACQUIRE
    if leaked:
        print(f"\nSTAGE 1 DID NOT DELIVER: {sorted(leaked)} are scheduled here.")
        print("  Each of those is a fetch that would run inside the curation stage.")
        for dtype, where in missing:
            print(f"    missing  {dtype:36s} {where}")
        if not a.allow_acquisition:
            print("\n  Run stage 1 first, or pass --allow-acquisition to fetch here anyway.")
            return 1

    absent = EXPECTED - used
    if absent:
        print(f"\nMISSING expected transforms: {sorted(absent)}")
        return 1
    print("\nevery expected curation transform is in the plan"
          + ("" if leaked else ", and no acquisition is"))

    if not a.run:
        print("\n(plan only -- pass --run to execute)")
        return 0
    execute(agent, task, a.work, a.timeout_hours)
    print("publish with --publish (never automatic: it rewrites DVC directory hashes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
