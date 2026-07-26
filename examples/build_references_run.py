"""Execute the reference build -- the same DAG ``build_references_dag.py`` renders.

    PATH="/home/tony/lib/miniforge3/envs/msm/bin:$PATH" \\
        python examples/build_references_run.py --local-raw --run

The gate beside this file stays plan-only, and deliberately: it asserts that the
CONTRACTS compose, from exactly one given, on a machine holding none of the 41 GB. This
driver is the other half -- it resolves the same graph against real inputs and runs it.

Three modes, and they are separate on purpose:

  --local-raw   add every chunk under data/raw/ to the input library, so the planner
                treats them as already staged and never schedules a download. On this
                machine that is ~14 GB already pinned by DVC; re-fetching it to prove the
                fetch works is a day of network for no new information. The fetch path
                still exists and is what runs on a machine that has nothing.

  --run         Deploy / StageWorkflow / RunWorkflow / CheckWorkflow, the shape
                src/fabfos/pipeline.py uses. Without it this prints the plan and stops.

  --publish     copy what the run produced to the paths data/ declares. SEPARATE from
                --run, and never a side effect of it: publication rewrites DVC directory
                hashes, which is why the pre-library run.sh stopped short and printed the
                commands instead. See --help for what it will and will not overwrite.

THE ONE-GIVEN ASSERTION BELONGS TO THE GATE, NOT HERE. `--local-raw` legitimately adds
many givens -- that is what it is for -- so running that check in this driver would fail
a correct build. If you want to know whether the graph still closes from one given, run
``build_references_dag.py``; that is the question it exists to answer.
"""
from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

from metasmith.python_api import (
    Agent,
    DataInstanceLibrary,
    Duration,
    Resources,
    Size,
    Source,
    TargetBuilder,
    TransformInstanceLibrary,
)

# The runtime enum was renamed ContainerRuntime -> Runtime when the mamba executor
# landed (a container runtime is no longer the only kind). Accept either, so this runs
# against the pinned submodule engine and against an older installed one.
try:
    from metasmith.python_api import Runtime
except ImportError:                                            # pragma: no cover
    from metasmith.python_api import ContainerRuntime as Runtime

REPO = Path(__file__).resolve().parent.parent
MLIB = REPO / "src" / "metasmith_libraries"
BREF = REPO / "build_references"
DATA = REPO / "data"

sys.path.insert(0, str(REPO / "examples"))
from build_references_dag import TARGETS  # noqa: E402  -- one list of targets, not two

# raw:: type -> where that chunk sits under data/raw/, for --local-raw. One entry per
# acquisition product, because the point of pre-staging is to skip a specific fetch: a
# glob that happened to match nothing would silently re-download instead of failing.
LOCAL_RAW = {
    "raw::metanetx_chem_prop":            "metanetx/4.5/chem_prop.tsv",
    "raw::metanetx_chem_xref":            "metanetx/4.5/chem_xref.tsv",
    "raw::metanetx_reac_prop":            "metanetx/4.5/reac_prop.tsv",
    "raw::metanetx_reac_xref":            "metanetx/4.5/reac_xref.tsv",
    "raw::rhea2uniprot":                  "rhea/rhea2uniprot.tsv",
    "raw::rhea2uniprot_trembl":           "rhea/rhea2uniprot_trembl.tsv.gz",
    "raw::kofam_profiles_archive":        "kofam/profiles.tar.gz",
    "raw::kofam_ko_list":                 "kofam/ko_list",
    "raw::kegg_ko_to_kegg_r":             "kegg/ko_to_kegg_r.tsv",
    "raw::uniref50_fasta":                "uniref/uniref50.fasta.gz",
    "raw::equilibrator_cache":            "equilibrator",
    "raw::laser_records":                 "literature/laser",
    "raw::keio_records":                  "literature/keio",
    "raw::eydallin_records":              "literature/eydallin",
    "raw::het_screen_records":            "literature/het_screen",
}

# The MetaCyc drop-in. Always given -- it is licensed and cannot be fetched, which is the
# whole reason the graph has a given at all.
METACYC_DROP_IN = "metacyc"

# Where each produced artifact lands, for --publish. Keyed by the type the transform
# produces, valued by the path data/ declares in build_references/REFERENCES.md.
PUBLISH = {
    "ref::kofamscan_profiles":     "reference/kofam/profiles",
    "ref::kofamscan_ko_list":      "reference/kofam/ko_list",
    "ref::uniref50_diamond_db":    "reference/uniref50.dmnd",
    "ref::mnxr_lookup":            "reference/mnxr_lookup.parquet",
    "ref::atom_pairs":             "reference/metabolism/atom_pairs.parquet",
    "ref::metabolism_vocab":       "reference/metabolism/vocab.parquet",
    "ref::direction_ratios":       "reference/metabolism/direction.parquet",
    "ref::reference_label_pool":   "reference/proteinBERT",
    "bench::condition_gpr":        "benchmark/gpr.parquet",
    "bench::conditions":           "benchmark/conditions.tsv",
}

# The per-host tables land under a host directory, so they are published by host rather
# than by a single path. R1 is the same shape in reverse: REFERENCES.md is explicit that
# the host assemblies ARE the background reference and no concatenation transform stands
# between them and the caller, so `reference/hosts.fna` is produced HERE, at publication,
# rather than by a transform that would only restate three files.
PUBLISH_PER_HOST = {
    "ref::gpr_table_gem":    "benchmark/hosts/{host}/gpr_gem.parquet",
    "ref::gpr_table_denovo": "benchmark/hosts/{host}/gpr_denovo.parquet",
}
PUBLISH_HOSTS_FNA = "reference/hosts.fna"

# R2. Not in the DAG -- a transform that only restates a file is a step that can go wrong
# in exchange for nothing -- so publication copies the pinned raw file directly.
PUBLISH_VECTOR = ("raw/vector/pcc1.fna", "reference/vector.fna")

# Chunks that are DVC bulk rather than git text. Publishing rewrites their directory
# hashes, which is why this is printed as a command rather than run.
DVC_CHUNKS = ("reference/metabolism", "reference/proteinBERT", "reference/uniref50.dmnd",
              "reference/kofam", "reference/mnxr_lookup.parquet", "benchmark/hosts")


def build_inputs(work: Path, local_raw: bool) -> DataInstanceLibrary:
    inputs = DataInstanceLibrary(work / "inputs.xgdb")
    for tl in ("ncbi.yml", "sequences.yml", "annotation.yml", "ref.yml", "lib.yml"):
        inputs.AddTypeLibrary(MLIB / "data_types" / tl)
    for tl in ("raw.yml", "interm.yml", "bench.yml", "buildlib.yml"):
        inputs.AddTypeLibrary(BREF / "data_types" / tl)

    metacyc = DATA / "raw" / METACYC_DROP_IN
    if not metacyc.exists():
        raise SystemExit(
            f"the MetaCyc drop-in is not at {metacyc}.\n"
            f"  MetaCyc/BioCyc flat-files are licensed and not redistributable; place the "
            f"distribution there manually. It is the ONE given of this graph, and without "
            f"it neither ensemble has its independent member.")
    inputs.AddItem(metacyc, "raw::metacyc_flatfiles")
    n = 1

    if local_raw:
        missing = []
        for dtype, rel in LOCAL_RAW.items():
            p = DATA / "raw" / rel
            if not p.exists():
                missing.append((dtype, p))
                continue
            inputs.AddItem(p, dtype)
            n += 1
        if missing:
            # Not fatal: an absent chunk simply means its acquire transform runs, which is
            # the correct behaviour. Printed by name so a surprise download is a thing you
            # were told about rather than a thing you notice from the network graph.
            print("\n--local-raw: not staged locally, so these WILL be acquired:")
            for dtype, p in missing:
                print(f"    {dtype:34s} {p.relative_to(REPO)}")
            print()
    inputs.Save()
    print(f"inputs: {n} staged item(s)"
          + ("" if local_raw else "  (the one given; everything else is fetched)"))
    return inputs


def plan(work: Path, local_raw: bool):
    inputs = build_inputs(work, local_raw)
    resources = [
        DataInstanceLibrary.Load(MLIB / "resources" / "env"),
        DataInstanceLibrary.Load(MLIB / "resources" / "lib"),
        DataInstanceLibrary.Load(BREF / "resources" / "buildlib"),
        inputs,
    ]
    # functionalAnnotation only, from the shipped library -- NOT logistics. Its logistics
    # sibling produces ref::kofamscan_* and ref::uniref50_diamond_db, exactly the types
    # compile/ produces, and two producers for one reference means a tiebreak decides
    # provenance. This is the same restriction the gate applies, and it must stay true
    # here or the two disagree about what built a reference.
    transforms = [
        TransformInstanceLibrary.Load(MLIB / "transforms" / "functionalAnnotation"),
        TransformInstanceLibrary.Load(BREF / "transforms" / "acquire"),
        TransformInstanceLibrary.Load(BREF / "transforms" / "compile"),
        TransformInstanceLibrary.Load(BREF / "transforms" / "benchmark"),
    ]

    targets = TargetBuilder()
    for _id, dtype in TARGETS:
        targets.Add(dtype)

    # MAMBA, not a container runtime. Three of the envs this build needs
    # (build-refs-{rdkit,equilibrator,cobra}) are conda environments with no published
    # image -- rdkit+rxnmapper+localmapper, equilibrator+dGbyG, and cobra -- so the
    # runtime that can run all of them is the one that resolves an env declaration's
    # `conda:` key. See src/metasmith_libraries/resources/env/*.env.
    agent = Agent(home=Source.FromLocal(work / "agent_home"),
                  runtime=Runtime.MAMBA)
    return agent, agent.GenerateWorkflow(
        samples=list(inputs.AsSamples("raw::metacyc_flatfiles")),
        resources=resources,
        transforms=transforms,
        targets=targets,
    )


def wait_for_run(staging: Path, task_key: str, timeout_s: int) -> Path:
    """Poll the run's own log for its completion marker.

    Same shape as src/fabfos/pipeline.py::_wait_for_run. The reference build's longest
    step is the AAM ensemble (both neural mappers over ~57.5k reactions), so the default
    timeout is days rather than hours -- a driver that gives up before the work does would
    leave a half-written cache and report a failure that did not happen.
    """
    internals = staging / "agent_home" / "runs" / task_key / "_metasmith"
    deadline = time.time() + timeout_s
    while time.time() < deadline:
        if internals.exists():
            log_dirs = sorted(p for p in internals.glob("logs.*") if "latest" not in p.name)
            if log_dirs:
                last_log = log_dirs[-1] / "main.log"
                if last_log.exists():
                    text = last_log.read_text(errors="ignore")
                    if "run completed at" in text:
                        return last_log
                    if "ERROR" in text and "nextflow" in text.lower():
                        tail = "".join(text.splitlines(keepends=True)[-40:])
                        raise RuntimeError(f"workflow {task_key} failed:\n{tail}")
        time.sleep(10)
    raise TimeoutError(f"workflow {task_key} did not finish within {timeout_s}s")


def publish(results: Path, dry_run: bool) -> int:
    """Land what the run produced at the paths data/ declares.

    Deliberately a separate invocation. Publication rewrites DVC directory hashes, so it
    is something you decide to do, not something that happens because a build finished.
    `--publish-dry-run` prints every copy without making one.
    """
    if not results.exists():
        raise SystemExit(f"no results at {results}; run with --run first")

    moved = 0
    for src in sorted(results.rglob("*")):
        rel = src.relative_to(results)
        dest_rel = None
        for dtype, target in PUBLISH.items():
            if rel.parts and dtype.split("::")[1] in rel.name:
                dest_rel = target
                break
        if dest_rel is None:
            continue
        dest = DATA / dest_rel
        print(f"  {rel}  ->  data/{dest_rel}")
        if not dry_run:
            dest.parent.mkdir(parents=True, exist_ok=True)
            if src.is_dir():
                if dest.exists():
                    shutil.rmtree(dest)
                shutil.copytree(src, dest)
            else:
                shutil.copyfile(src, dest)
        moved += 1

    vec_src, vec_dest = PUBLISH_VECTOR
    print(f"  data/{vec_src}  ->  data/{vec_dest}   (R2: pinned raw, not in the DAG)")
    if not dry_run and (DATA / vec_src).exists():
        shutil.copyfile(DATA / vec_src, DATA / vec_dest)

    print("\nDVC chunks whose directory hash this changes -- run these yourself:")
    for chunk in DVC_CHUNKS:
        print(f"    dvc add data/{chunk}")
    print("\ndata/.gitignore is the only .gitignore under data/; check that still holds:")
    print("    [ $(find data -name .gitignore | wc -l) -eq 1 ]")
    return 0 if moved else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--local-raw", action="store_true",
                    help="stage data/raw/ chunks as inputs so their fetches are skipped")
    ap.add_argument("--run", action="store_true",
                    help="execute the plan; without it this prints the plan and stops")
    ap.add_argument("--publish", action="store_true",
                    help="copy a finished run's results to the paths data/ declares")
    ap.add_argument("--publish-dry-run", action="store_true",
                    help="print what --publish would copy, and copy nothing")
    ap.add_argument("--work", type=Path, default=REPO / "_build_references_run",
                    help="the run directory; kept rather than temporary, because the AAM "
                         "caches inside it are what make a killed run resumable")
    ap.add_argument("--timeout-hours", type=float, default=72.0)
    a = ap.parse_args()

    a.work.mkdir(parents=True, exist_ok=True)

    if a.publish or a.publish_dry_run:
        staging = a.work
        runs = sorted((staging / "agent_home" / "runs").glob("*"))
        if not runs:
            raise SystemExit(f"no runs under {staging}/agent_home/runs")
        return publish(runs[-1] / "results", dry_run=a.publish_dry_run)

    agent, task = plan(a.work, a.local_raw)
    if not task.ok:
        print(f"FAILED to plan:\n{task.plan}")
        return 1
    print(f"\nPlan OK -- {len(task.plan.steps)} steps")
    for step in sorted(task.plan.steps, key=lambda s: s.order):
        prods = [i.dtype_name for g in step.produces for i in g]
        print(f"  {step.order:>3}  {Path(step.transform._path).stem:<28} -> {prods}")

    if not a.run:
        print("\n(plan only -- pass --run to execute)")
        return 0

    agent.Deploy()
    agent.StageWorkflow(task, on_exist="clear")
    agent.RunWorkflow(task)
    wait_for_run(a.work, task._key, int(a.timeout_hours * 3600))
    agent.CheckWorkflow(task)
    print(f"\nresults: {a.work}/agent_home/runs/{task._key}/results")
    print("publish them with --publish (never automatic: it rewrites DVC dir hashes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
