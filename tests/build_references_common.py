"""Shared wiring for the two-stage reference build.

The build splits at the tier boundary, not at a convenient point in the code:

    stage 1  acquire/   -- every byte that comes from OUTSIDE, plus the licensed
                           drop-in's split. Fills data/raw/. Network-bound.
    stage 2  compile/ + benchmark/ + the shipped annotation lanes -- everything
                           derived HERE from what stage 1 landed. Fills
                           data/reference/ and data/benchmark/. Compute-bound.

Why split at all, when the planner is perfectly happy to resolve one graph end to end:
the two halves fail differently and are re-run on different schedules. An acquisition
fails on a dead mirror, a rate limit or a licence, and re-running it costs bandwidth. A
curation step fails on a bug in a transform, and re-running it costs a day of CPU. Fusing
them means every method fix re-enters a graph whose first third is a 15 GB download, and
means a download failure at hour 30 has nothing to show for the 30 hours.

The seam is `data/raw/`, which is exactly what the tier rule already says it is: stage 1
writes it, stage 2 only reads it. That makes the handoff inspectable on disk rather than
internal to a run directory -- you can look at what stage 1 produced, pin it with DVC,
and hand it to a stage 2 running anywhere.

STAGE 2 STILL LOADS acquire/. That looks like it defeats the split, and it is the thing
that makes the split checkable: stage 2 stages everything stage 1 was supposed to
produce, so if an acquire transform appears in stage 2's plan, stage 1 did not deliver
it. `assert_no_acquisition` turns that into a named failure. Hiding acquire/ from stage 2
would instead produce an unresolvable plan whose error message names a type, not a gap.
"""
from __future__ import annotations

import shutil
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# The PINNED engine, ahead of whatever is installed in the env. src/metasmith is on a
# feat/fabfos carrying the mamba executor; the metasmith installed in the `msm` env
# predates it and has no MAMBA runtime at all, so importing that one makes this build
# unplannable in a way whose error message points at an enum rather than at a version.
_ENGINE = REPO / "src" / "metasmith" / "src"
if (_ENGINE / "metasmith").is_dir():
    sys.path.insert(0, str(_ENGINE))

from metasmith.python_api import (                                      # noqa: E402
    Agent,
    DataInstanceLibrary,
    Source,
    TargetBuilder,
    TransformInstanceLibrary,
)

# The runtime enum was renamed ContainerRuntime -> Runtime when the mamba executor
# landed (a container runtime is no longer the only kind).
try:
    from metasmith.python_api import Runtime                            # noqa: E402
except ImportError:                                            # pragma: no cover
    from metasmith.python_api import ContainerRuntime as Runtime
if not hasattr(Runtime, "MAMBA"):                              # pragma: no cover
    import metasmith
    raise SystemExit(
        f"the metasmith at {metasmith.__file__} predates the mamba executor.\n"
        f"  This build needs Runtime.MAMBA: three of its envs (build-refs-rdkit, "
        f"-equilibrator, -cobra) are conda environments with no published image, and "
        f"metasmith picks conda: vs container: off one global runtime.\n"
        f"  The pinned engine is src/metasmith (feat/fabfos). Check the submodule is "
        f"checked out: git submodule update --init src/metasmith")
MLIB = REPO / "src" / "metasmith_libraries"
BREF = REPO / "build_references"
DATA = REPO / "data"
ARTIFACTS = REPO / "tests" / "artifacts"

# Both stages' run directories live under data/scratch/. Everything in a run directory IS
# data -- the staged input library, the agent home, the per-step work dirs, the results --
# and it is the one tier with no permanence guarantee. Reading the tier rule, `scratch/`
# means "derived, transient, safe to delete", as against raw/ (acquired), reference/ and
# benchmark/ (compiled here) and curated/ (hand-authored). Nothing in it is DVC-pinned and
# nothing outside it points in.
SCRATCH = DATA / "scratch"

# THE ONE GIVEN. Licensed and not redistributable, which is the whole reason this graph
# has a given at all; every other input is produced by a transform in the library.
METACYC_GIVEN = DATA / "raw" / "metacyc"

# ---------------------------------------------------------------------------
# The seam, written twice in opposite directions.
#
# One map, read forwards by stage 2 (stage these so they are not re-fetched) and
# backwards by stage 1's publisher (put what you fetched here). Two maps would drift, and
# the drift is silent: a publisher writing where the stager does not look means stage 2
# re-downloads something that is sitting on disk.
# ---------------------------------------------------------------------------
RAW_AT = {
    "raw::metanetx_chem_prop":     "metanetx/4.5/chem_prop.tsv",
    "raw::metanetx_chem_xref":     "metanetx/4.5/chem_xref.tsv",
    "raw::metanetx_reac_prop":     "metanetx/4.5/reac_prop.tsv",
    "raw::metanetx_reac_xref":     "metanetx/4.5/reac_xref.tsv",
    "raw::rhea2uniprot":           "rhea/rhea2uniprot.tsv",
    "raw::rhea2uniprot_trembl":    "rhea/rhea2uniprot_trembl.tsv.gz",
    "raw::kofam_profiles_archive": "kofam/profiles.tar.gz",
    "raw::kofam_ko_list":          "kofam/ko_list",
    "raw::kegg_ko_to_kegg_r":      "kegg/ko_to_kegg_r.tsv",
    "raw::uniref50_fasta":         "uniref/uniref50.fasta.gz",
    "raw::equilibrator_cache":     "equilibrator",
    "raw::laser_records":          "literature/laser",
    "raw::keio_records":           "literature/keio",
    "raw::eydallin_records":       "literature/eydallin",
    "raw::het_screen_records":     "literature/het_screen",
    # The two slices metacyc_licensed cuts out of the given. Named here so stage 2 stages
    # them rather than re-running the verify -- which would be harmless, but would put an
    # acquire transform in stage 2's plan and trip the check that proves the split held.
    "raw::metacyc_atom_mappings_smiles": "metacyc/atom-mappings-smiles.dat",
    "raw::metacyc_reactions":            "metacyc/reactions.dat",
}

# The host set fans out three ways, so its products are staged by glob rather than by
# path. host_genome emits four types and only two are consumed downstream -- orfs feed
# the four annotation lanes, the assembly IS reference R1 -- but all four are stage 1
# products and all four are staged, so stage 2's "no acquisition" check is answered by
# what is on disk rather than by which types happen to matter today.
HOSTS_DIR = DATA / "raw" / "hosts"
RAW_PER_HOST = {
    "sequences::isolate_assembly": "genome/*.fna",
    "sequences::orfs":             "genome/*.faa",
    "sequences::gbk":              "genome/*.gbk",
    "sequences::gff":              "genome/*.gff",
    "raw::host_gem":               "GEM/*.json",
}

TYPE_LIBRARIES = (
    [MLIB / "data_types" / f for f in
     ("ncbi.yml", "sequences.yml", "annotation.yml", "ref.yml", "lib.yml")]
    + [BREF / "data_types" / f for f in
       ("raw.yml", "interm.yml", "bench.yml", "buildlib.yml")]
)


def new_inputs(work: Path) -> DataInstanceLibrary:
    """An input library carrying every type library and the one given."""
    inputs = DataInstanceLibrary(work / "inputs.xgdb")
    for tl in TYPE_LIBRARIES:
        inputs.AddTypeLibrary(tl)
    if not METACYC_GIVEN.exists():
        raise SystemExit(
            f"the MetaCyc drop-in is not at {METACYC_GIVEN}.\n"
            f"  MetaCyc/BioCyc flat-files are licensed and not redistributable; place "
            f"the distribution there manually. It is the ONE given of this graph, and "
            f"without it neither ensemble has its independent member.")
    inputs.AddItem(METACYC_GIVEN, "raw::metacyc_flatfiles")
    return inputs


def make_standins(gaps: list, at: Path) -> list[tuple[Path, str]]:
    """Empty files standing in for stage 1 products that are not on disk yet.

    Planning is type-driven and never opens an input, so a zero-byte file resolves a type
    exactly as the real 12 GB does. This is what lets stage 2's curation DAG be rendered
    and inspected before stage 1 has ever run -- the same trick the plan-only gate uses
    for the MetaCyc given.

    Only ever reachable behind an explicit flag. A stand-in that reached a RUN would hand
    a transform an empty file, and most of them would produce an empty output and succeed.
    """
    at.mkdir(parents=True, exist_ok=True)
    made = []
    for i, (dtype, where) in enumerate(gaps):
        # Host products fan out three ways, so one stand-in per host, not one per type.
        n = 3 if "/hosts/" in where else 1
        for k in range(n):
            p = at / f"{dtype.replace('::', '.')}.{i}.{k}"
            p.write_bytes(b"")
            made.append((p, dtype))
    return made


def stage_raw(inputs: DataInstanceLibrary, *, quiet: bool = False) -> tuple[int, list]:
    """Add everything already sitting in data/raw/ to the input library.

    Returns (staged, gaps, staged_types). A gap -- (dtype, where) for a product not on
    disk -- is NOT fatal; it means that acquisition runs, which is the correct behaviour.
    It is named, so a surprise download is a thing you were told about rather than a thing
    you notice from the network graph an hour later.

    `staged_types` is what a caller filters its TARGETS against. The planner treats a
    target it could not put in the plan as a failure, and a target already satisfied by a
    staged input is exactly that -- so asking for one fails a build whose actual problem
    is that there is nothing left to do.

    THE HOST BLOCK IS ALL-OR-NOTHING, and that is the one non-obvious rule here. Its five
    products fan out three ways, and the planner reasons about TYPES, not paths: staging
    two hosts' assemblies and not the third's would satisfy `sequences::isolate_assembly`
    outright, and the missing host would never be fetched. So if any host is short any
    file, none of the host products are staged and host_genome/host_gem re-run for all
    three -- 61 MB, against a silently two-host reference build.
    """
    staged, gaps = 0, []
    for dtype, rel in RAW_AT.items():
        p = DATA / "raw" / rel
        if p.exists():
            inputs.AddItem(p, dtype)
            staged += 1
        else:
            gaps.append((dtype, f"raw/{rel}"))

    host_items, host_gaps = [], []
    for host_dir in sorted(p for p in HOSTS_DIR.glob("*") if p.is_dir()):
        for dtype, pattern in RAW_PER_HOST.items():
            hits = sorted(host_dir.glob(pattern))
            if hits:
                host_items += [(h, dtype) for h in hits]
            else:
                host_gaps.append((dtype, f"raw/hosts/{host_dir.name}/{pattern}"))
    if host_gaps:
        # Report EVERY host product as a gap, not only the files that were literally
        # absent. Nothing in the block was staged, so a caller reading this list to decide
        # what still has to be produced would otherwise be told that four types are
        # missing when in fact all five are.
        if not quiet:
            print(f"\nthe host set is incomplete ({len(host_gaps)} file(s) short), so "
                  f"none of it is staged and all three hosts are re-acquired:")
            for dtype, where in host_gaps:
                print(f"    short   {dtype:34s} {where}")
        gaps += [(dtype, f"raw/hosts/*/{pattern}") for dtype, pattern in RAW_PER_HOST.items()]
    else:
        for p, dtype in host_items:
            inputs.AddItem(p, dtype)
        staged += len(host_items)

    if gaps and not quiet:
        print("\nnot on disk, so the acquisition for each of these will RUN:")
        for dtype, where in gaps:
            print(f"    {dtype:36s} {where}")
        print()
    staged_types = (set(RAW_AT) | set(RAW_PER_HOST)) - {d for d, _ in gaps}
    return staged, gaps, staged_types


def resource_libraries() -> list:
    return [
        DataInstanceLibrary.Load(MLIB / "resources" / "env"),
        DataInstanceLibrary.Load(MLIB / "resources" / "lib"),
        # buildlib:: -- the ported method modules (both ensembles, the bake encoding, the
        # benchmark cohort readers). Build side only, so it is loaded here and never
        # bundled into the wheel. A resource library, like env and lib, so it adds no given.
        DataInstanceLibrary.Load(BREF / "resources" / "buildlib"),
    ]


def transform_libraries(with_curation: bool) -> list:
    """acquire/ always; the curation side only for stage 2.

    functionalAnnotation only from the shipped library -- NOT logistics. Its logistics
    sibling produces ref::kofamscan_* and ref::uniref50_diamond_db, exactly the types
    compile/ produces, and two producers for one reference means a tiebreak decides
    provenance. The plan-only gate applies the same restriction; if the two disagree they
    disagree about what built a reference.
    """
    libs = [TransformInstanceLibrary.Load(BREF / "transforms" / "acquire")]
    if with_curation:
        libs += [
            TransformInstanceLibrary.Load(MLIB / "transforms" / "functionalAnnotation"),
            TransformInstanceLibrary.Load(BREF / "transforms" / "compile"),
            TransformInstanceLibrary.Load(BREF / "transforms" / "benchmark"),
        ]
    return libs


# The conda env metasmith ITSELF runs in, as against the envs the tools run in.
#
# `Agent.container` defaults to a docker:// URI, and under MAMBA that field is read as a
# conda environment NAME -- both name "the thing this runs in", which is why they share a
# field. Leaving the default makes staging try `conda run -n docker://quay.io/...` and
# fail on the characters in the URI, several steps after the runtime was chosen.
#
# Built from metasmith's own envs/base.yml plus a .pth onto src/metasmith/src, so the
# agent runs the PINNED engine. The `msm` env on this machine resolves metasmith through
# a global symlink into a different worktree, which is the one without the mamba
# executor -- exactly the version that cannot run this.
AGENT_ENV = "msm-fabfos"


def make_agent(work: Path) -> Agent:
    """MAMBA, not a container runtime.

    Three of the envs this build needs (build-refs-{rdkit,equilibrator,cobra}) are conda
    environments with no published image. metasmith resolves an env declaration's
    `conda:` key when the runtime is MAMBA and its `container:` key otherwise, off ONE
    global setting -- so the runtime that can run all of them is the only one that can run
    any of this. Both stages use it, so a step cannot mean different things in each.
    """
    return Agent(home=Source.FromLocal(work / "agent_home"),
                 runtime=Runtime.MAMBA, container=AGENT_ENV)


def generate(work: Path, inputs: DataInstanceLibrary, targets: list[str],
             *, with_curation: bool):
    inputs.Save()
    tb = TargetBuilder()
    for dtype in targets:
        tb.Add(dtype)
    agent = make_agent(work)
    return agent, agent.GenerateWorkflow(
        samples=list(inputs.AsSamples("raw::metacyc_flatfiles")),
        resources=resource_libraries() + [inputs],
        transforms=transform_libraries(with_curation),
        targets=tb,
    )


def print_plan(task) -> set[str]:
    used = {Path(s.transform._path).stem for s in task.plan.steps}
    print(f"\nPlan OK -- {len(task.plan.steps)} steps, {len(used)} distinct transforms\n")
    for step in sorted(task.plan.steps, key=lambda s: s.order):
        prods = [i.dtype_name for g in step.produces for i in g]
        print(f"  {step.order:>3}  {Path(step.transform._path).stem:<28} -> {prods}")
    return used


def render(task, name: str) -> Path:
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    svg = ARTIFACTS / f"{name}.svg"
    # lib/env/buildlib are resources every step carries; drawing them turns the graph
    # into a hairball where every node touches every other one.
    task.plan.RenderDAG(svg, blacklist_namespaces={"lib", "env", "buildlib"})
    print(f"\nDAG -> {svg}")
    return svg


def wait_for_run(work: Path, task_key: str, timeout_s: int) -> Path:
    """Poll the run's own log for its completion marker.

    Same shape as src/fabfos/pipeline.py::_wait_for_run. Stage 2's longest step is the AAM
    ensemble (both neural mappers over ~57.5k reactions), so its default timeout is days
    rather than hours -- a driver that gives up before the work does would leave a
    half-written cache and report a failure that did not happen.
    """
    internals = work / "agent_home" / "runs" / task_key / "_metasmith"
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


def execute(agent, task, work: Path, timeout_hours: float) -> Path:
    agent.Deploy()
    agent.StageWorkflow(task, on_exist="clear")
    agent.RunWorkflow(task)
    wait_for_run(work, task._key, int(timeout_hours * 3600))
    agent.CheckWorkflow(task)
    results = work / "agent_home" / "runs" / task._key / "results"
    print(f"\nresults: {results}")
    return results


def latest_results(work: Path) -> Path:
    runs = sorted((work / "agent_home" / "runs").glob("*"))
    if not runs:
        raise SystemExit(f"no runs under {work}/agent_home/runs")
    return runs[-1] / "results"


def add_common_args(ap) -> None:
    ap.add_argument("--run", action="store_true",
                    help="execute the plan; without it this prints the plan, renders the "
                         "DAG and stops")
    ap.add_argument("--publish", action="store_true",
                    help="copy a finished run's results to the paths data/ declares")
    ap.add_argument("--publish-dry-run", action="store_true",
                    help="print what --publish would copy, and copy nothing")
    ap.add_argument("--timeout-hours", type=float, default=72.0)


def publish_by_type(results: Path, mapping: dict[str, str], dest_root: Path,
                    *, dry_run: bool) -> int:
    """Copy a run's results to the paths data/ declares, keyed by produced type.

    A result is matched to a destination by its type's NAME appearing in the result path.
    That is a convention, not a guarantee the engine makes, so every match is printed and
    an empty match set is a non-zero exit rather than a quiet success -- publishing
    nothing looks exactly like publishing everything if you do not check.
    """
    if not results.exists():
        raise SystemExit(f"no results at {results}; run with --run first")
    moved = 0
    for src in sorted(results.rglob("*")):
        rel = src.relative_to(results)
        for dtype, target in mapping.items():
            if dtype.split("::")[1] in rel.name:
                dest = dest_root / target
                print(f"  {rel}  ->  {dest.relative_to(REPO)}")
                if not dry_run:
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    if src.is_dir():
                        if dest.exists():
                            shutil.rmtree(dest)
                        shutil.copytree(src, dest)
                    else:
                        shutil.copyfile(src, dest)
                moved += 1
                break
    if not moved:
        print("  nothing matched -- the results layout does not carry type names in "
              "its paths, so this mapping needs the run's manifest instead")
    return moved
