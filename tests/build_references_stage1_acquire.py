"""Stage 1 of the reference build: acquire every byte that comes from outside.

    PATH="/home/tony/lib/miniforge3/envs/msm/bin:$PATH" \\
        python tests/build_references_stage1_acquire.py                 # plan + DAG
    ... python tests/build_references_stage1_acquire.py --refetch       # ignore data/raw/
    ... python tests/build_references_stage1_acquire.py --run
    ... python tests/build_references_stage1_acquire.py --publish

This is the network half. It fills `data/raw/` and touches nothing else, which is the
tier rule stated as a program: raw/ is acquired only, and this is the only thing allowed
to write it. Stage 2 reads what lands here and never fetches.

WHAT IT PRODUCES -- 11 transforms, and they are three different kinds of thing:

  DOWNLOADS (a URL is in the transform, and it is the provenance)
    metanetx            4 tables from metanetx.org/ftp/4.5      ~1.5 GB
    rhea                2 tables from ftp.expasy.org            ~215 MB
    kofam               profiles.tar.gz + ko_list from GenomeNet ~1.5 GB
    uniref50            uniref50.fasta.gz from ftp.uniprot.org   ~12 GB   <- the long pole
    equilibrator_cache  the compound cache, pinned rather than left to first use
    host_genome         3 assemblies via ncbi datasets (fna/faa/gff/gbk)  ~61 MB
    host_gem            2 curated models from BiGG; EPI300 borrows DH10B's
    literature          the LASER checkout (git, pinned by rev) + the Keio supplement zip

    kegg_ko_reactions   the whole KO -> reaction map in ONE call to KEGG REST's bulk
                        link endpoint -> ko_to_kegg_r.tsv, ~2 s. This replaced a
                        per-KO crawl; see the transform for why, and for the
                        measurement that says the two agree on 99.6% of KOs

  DECLARATIONS AND REFUSALS (no network at all)
    host_accessions     the three hosts, as a fact in the library rather than an argument
    metacyc_licensed    verifies the licensed drop-in and splits out its two members;
                        NEVER fetches -- that is the point of it being the one given
    literature_extractions
                        eydallin and het_screen are HAND EXTRACTIONS from paper
                        supplements, so there is nothing to download. This transform
                        CANNOT succeed: it is reachable only when the tables are absent
                        from data/raw/literature/, and refusing by name is the whole of
                        its behaviour. Building without them ships a conditions table
                        that looks complete and is missing two arms.

On a machine that already has the bulk pinned by DVC, most of the above is skipped: every
chunk already under `data/raw/` is staged as an input, so the planner has no reason to
schedule its fetch. `--refetch` turns that off and proves the fetch path still works,
which is the only way that code is ever exercised here.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_references_common import (                                    # noqa: E402
    DATA, RAW_AT, SCRATCH, add_common_args, generate, latest_results, new_inputs,
    print_plan, publish_by_type, render, stage_raw, execute,
)

# Everything acquire/ produces. Reading the tier rule off the transforms rather than off
# a hand-written list: if a new acquisition lands, its product belongs here, and a target
# that no transform produces fails loudly at plan time.
#
# ncbi::assembly_accession is deliberately absent -- it is host_accessions' product, but
# demanding it as a target as well as a requirement of host_genome/host_gem adds nothing
# and would suggest an accession file is an artifact worth keeping.
TARGETS = [
    "raw::metanetx_chem_prop",
    "raw::metanetx_chem_xref",
    "raw::metanetx_reac_prop",
    "raw::metanetx_reac_xref",
    "raw::rhea2uniprot",
    "raw::rhea2uniprot_trembl",
    "raw::kofam_profiles_archive",
    "raw::kofam_ko_list",
    "raw::kegg_ko_to_kegg_r",
    "raw::uniref50_fasta",
    "raw::equilibrator_cache",
    "raw::laser_records",
    "raw::keio_records",
    "raw::eydallin_records",
    "raw::het_screen_records",
    "raw::metacyc_atom_mappings_smiles",
    "raw::metacyc_reactions",
    "raw::host_gem",
    "sequences::isolate_assembly",
    "sequences::orfs",
    "sequences::gff",
    "sequences::gbk",
]

EXPECTED = {
    "metanetx", "rhea", "kofam", "uniref50", "kegg_ko_reactions",
    "equilibrator_cache", "literature", "literature_extractions", "metacyc_licensed",
    "host_accessions", "host_genome", "host_gem",
}

WORK = SCRATCH / "stage1_acquire"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refetch", action="store_true",
                    help="do NOT stage data/raw/; re-acquire everything from source. The "
                         "only way the fetch path is ever exercised on a machine that "
                         "already holds the data")
    ap.add_argument("--work", type=Path, default=WORK)
    add_common_args(ap)
    a = ap.parse_args()
    a.work.mkdir(parents=True, exist_ok=True)

    if a.publish or a.publish_dry_run:
        print(f"== stage 1 publish -> data/raw/ ==")
        n = publish_by_type(latest_results(a.work), RAW_AT, DATA / "raw",
                            dry_run=a.publish_dry_run)
        print("\nraw/ chunks whose DVC directory hash this changes -- run these yourself:")
        for chunk in ("raw/metanetx", "raw/rhea", "raw/kofam", "raw/uniref", "raw/kegg",
                      "raw/equilibrator", "raw/literature/keio", "raw/literature/eydallin",
                      "raw/literature/het_screen", "raw/hosts/e_coli_k12",
                      "raw/hosts/e_coli_dh10b", "raw/hosts/e_coli_epi300"):
            print(f"    dvc add data/{chunk}")
        return 0 if n else 1

    # plan-only never opens the given; only --run needs the licensed bytes
    inputs = new_inputs(a.work, require_given=a.run)
    if a.refetch:
        print("--refetch: nothing staged; every acquisition runs")
        wanted = list(TARGETS)
    else:
        n, gaps, staged_types = stage_raw(inputs)
        print(f"inputs: {n + 1} staged item(s) already under data/raw/")
        # Ask only for what is NOT on disk -- see stage_raw's docstring for why a target
        # that is already satisfied is a plan failure rather than a no-op.
        wanted = [t for t in TARGETS if t not in staged_types]

    if not wanted:
        print("\ndata/raw/ is complete -- nothing to acquire. Stage 2 can run now.")
        return 0

    print(f"\nacquiring {len(wanted)} product type(s):")
    for t in wanted:
        print(f"    {t}")

    agent, task = generate(a.work, inputs, wanted, with_curation=False)
    if not task.ok:
        print(f"FAILED to plan:\n{task.plan}")
        return 1

    used = print_plan(task)
    idle = EXPECTED - used
    if idle:
        # Not a failure. An acquisition absent from the plan means its product is already
        # on disk, which is exactly what staging is for -- but it is worth naming, because
        # "nothing to do" and "the branch was dropped" render identically in a DAG.
        print(f"\nalready on disk, so not scheduled: {sorted(idle)}")
    render(task, "build_references_stage1_acquire")

    if not a.run:
        print("\n(plan only -- pass --run to execute)")
        return 0
    execute(agent, task, a.work, a.timeout_hours)
    print("publish into data/raw/ with --publish (never automatic: it rewrites DVC "
          "directory hashes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
