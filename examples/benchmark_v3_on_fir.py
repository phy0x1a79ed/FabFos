#!/usr/bin/env python3
"""Run the frozen ECSPr X/Y benchmark v3 on an HPC host, then retrieve the result.

    PYTHONPATH=src/metasmith/src:src python examples/benchmark_v3_on_fir.py --user phyberos
    PYTHONPATH=src/metasmith/src:src python examples/benchmark_v3_on_fir.py --user phyberos --plan-only

THE ENGINE ON PYTHONPATH MUST BE THE PINNED ONE -- READ THIS BEFORE CHANGING IT
-------------------------------------------------------------------------------
Note the `src/metasmith/src` FIRST. The `fabfos` conda env resolves `metasmith`
to a local editable checkout (~/lib/locals/metasmith) which is ahead of the
submodule pin. That matters here in a way it does not for a local plan, because
`Deploy()` pulls the agent container by a tag derived from the engine's OWN
version: an unpublished dev version asks quay for a manifest that does not
exist, and the run dies ~40 s in with "manifest unknown" followed by a confusing
cascade about a missing relay binary. Only released versions have images.

`assert_pinned_engine()` below turns that into an immediate, named refusal.

Deploy -> Generate -> Stage -> Run -> Wait -> Retrieve. The scoring step is
NOT here: it is a sub-minute local operation on the merged table, and running
it beside the solve would put a multi-hour cluster task behind every change to
how a score is computed.

THE USERNAME IS AN ARGUMENT, AND THAT IS THE POINT
---------------------------------------------------
`main/local_mock/smoke_hpc_deploy.py` resolves it by shelling `ssh <host> echo
$USER`. On this workstation that cannot work and must not be retried: the ssh
config sets `ControlMaster no` behind a ProxyCommand guard, so a bare `ssh fir`
cannot open a fresh authenticated session -- it can only ride an existing
multiplexed one. A loop around that call is a Duo push per iteration, and a
prior ECSPr run in this project was halted by an account lockout caused exactly
that way.

So: `--user` is required, exactly one session is opened, and on failure this
exits telling the human to connect once by hand. Never delete the
ControlMaster socket and never retry the connect in a loop.

WHY THERE IS NO ARRAY, AND SO NOTHING TO DISABLE
-------------------------------------------------
The plan called for disabling the engine's job-array batching, because array
contention is the exact condition under which fir's overlay filesystem throws
bus errors (recorded in the engine's own source). That turns out to be moot by
construction rather than by configuration: this workflow is TWO tasks, a solve
and a merge. `solve_benchmark.py` loops over the 16 (facet, element) shards
in-process -- see its docstring for why the fan-out is not real yet -- so there
is no array to contend. Concurrency comes from 32 worker PROCESSES inside the
one task, which is the shape that actually scales here: the solver factorizes
with SuperLU, which is serial, so the engine pins OMP/OPENBLAS/MKL/NUMEXPR to 1
before numpy loads and forks instead. Oversubscribing measured >10x slower.

If the shards are ever split into 16 real instances, revisit this: at that
point there IS an array, and the batching should be disabled in favour of
queue-size concurrency.

WHAT IS STAGED, AND WHAT IS DELIBERATELY NOT
---------------------------------------------
Five items, all from the frozen benchmark tree, all addressed through canon so
no absolute path appears here. Note what is absent: no annotation lane, no
evidence chain, no recovery experiment. That whole chain is upstream of the
freeze. A benchmark score has to measure the SOLVER, and anything that could
re-derive X would make it measure the annotation instead.

`benchmark_universe` is 57.6 MB and is staged even though it looks like a
build-time artifact, because `base_plus_reactions` re-reads it at SOLVE time to
find the atom-transit weights of an inserted reaction. X withholds the atom
mapping on purpose, so without the universe the run does not fail fast -- it
fails on the first gain-of-function condition, deep into the job.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

from fabfos import canon  # noqa: E402
from fabfos.library import domains_for, resolve_library_root  # noqa: E402

from metasmith.python_api import (  # noqa: E402
    Agent, ContainerRuntime, DataInstanceLibrary, DataTypeLibrary, Source,
    SshSource, TargetBuilder, TransformInstanceLibrary,
)

LIB = resolve_library_root()

# The benchmark lane. `domains_for` drops ecsprNetA/ecsprNetB, which is not
# tidiness: on the benchmark, the annotation chain that feeds ecspr::base_graphs
# is upstream of the freeze and must not appear in the plan at all.
DOMAINS = domains_for(network="benchmark")

# type name -> canon symbol. All five live in the v3 tree.
STAGED: dict[str, str] = {
    "ecspr::benchmark_answer_key":  "BENCH_V3_Y",
    "ecspr::benchmark_base_graphs": "BENCH_V3_BASE_GRAPHS",
    "ecspr::benchmark_universe":    "BENCH_V3_UNIVERSE",
    "ecspr::benchmark_observations": "BENCH_V3_OBSERVATIONS",
    "ecspr::benchmark_inputs":      "BENCH_V3_X",
}

SETUP_COMMANDS = ["module load apptainer"]


def assert_pinned_engine() -> str:
    """Refuse to deploy an engine version that has no published container.

    `Deploy()` derives the agent image tag from the engine's own version, so an
    unreleased local checkout asks quay for a manifest that does not exist. The
    failure surfaces as a registry error and then an assertion about a missing
    relay binary -- neither of which names the actual cause. Checked here, on
    the version the interpreter actually imported, so the message arrives before
    a remote directory is created rather than after.

    Only enforced for a real deploy; --plan-only touches no registry.
    """
    import metasmith
    got = (Path(metasmith.__file__).parent / "version.txt").read_text().strip()
    # noqa: E501 -- see agent_container() for why the version alone is not the tag
    want = (REPO / "src/metasmith/src/metasmith/version.txt").read_text().strip()
    if got != want:
        raise SystemExit(
            f"engine version [{got}] is not the pin [{want}].\n"
            f"  imported from: {Path(metasmith.__file__).parent}\n"
            f"The agent container tag is derived from this version, and only "
            f"RELEASED versions have images on quay -- deploying [{got}] would "
            f"fail with 'manifest unknown' after creating a remote directory.\n"
            f"Re-run with the pin first on the path:\n"
            f"  PYTHONPATH=src/metasmith/src:src python {Path(__file__).name} ..."
        )
    return got


def agent_container() -> str:
    """The published agent image for the engine we actually imported.

    `Agent.container` defaults to `metasmith:{CONTAINER_TAG}`, and CONTAINER_TAG
    is `{VERSION}-{BUILD_HASH}` where BUILD_HASH is a content hash of the engine
    source tree written by `_build_hash.py` AT BUILD TIME. A source checkout --
    which is what the submodule pin is -- has no `build_hash.txt`, so the tag
    silently degrades to the bare version, and bare `0.18.8` was never pushed:
    quay carries `0.18.8-60556ca`. The run then dies on "manifest unknown".

    So the hash is COMPUTED here from the same function the build uses, rather
    than the tag being hardcoded. That is the difference between "this tag
    happens to work today" and "this image is provably built from the source on
    our PYTHONPATH" -- if the pin moves, this follows it, and if the resulting
    image was never published the pull fails loudly instead of running an engine
    that does not match the planner that produced the workflow.
    """
    from metasmith._build_hash import compute_build_hash
    from metasmith.constants import VERSION
    return f"docker://quay.io/hallamlab/metasmith:{VERSION}-{compute_build_hash()}"


def build_inputs(staging: Path) -> DataInstanceLibrary:
    xgdb = staging / "inputs.xgdb"
    if xgdb.exists():
        shutil.rmtree(xgdb)
    inputs = DataInstanceLibrary(xgdb)
    inputs.AddTypeLibrary(namespace="ecspr",
                          lib=DataTypeLibrary.Load(LIB / "data_types/ecspr.yml"))

    missing = []
    for type_name, symbol in STAGED.items():
        p = Path(getattr(canon, symbol))
        # Checked HERE rather than left to the scheduler. The planner resolves
        # on types and lineage, not on existence, so a missing tree plans
        # perfectly and then fails hours later inside a container.
        if not p.exists():
            missing.append(f"{type_name} -> canon.{symbol} -> {p}")
            continue
        inputs.AddItem(p, type_name)
    if missing:
        raise SystemExit("benchmark tree incomplete:\n  " + "\n  ".join(missing))

    inputs.Save()
    return inputs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--host", default="fir")
    ap.add_argument("--user", required=True,
                    help="remote username. REQUIRED and never auto-resolved -- "
                         "see this module's docstring on the lockout.")
    ap.add_argument("--scratch-root", default="/scratch")
    ap.add_argument("--timeout-s", type=float, default=6 * 3600)
    ap.add_argument("--poll-s", type=float, default=60.0)
    ap.add_argument("--plan-only", action="store_true",
                    help="stage and plan, render the DAG, touch no remote host")
    ap.add_argument("--out", default="transforms/build/benchmark/v3_build/run",
                    help="local directory to retrieve the merged result into")
    a = ap.parse_args()

    ts = int(time.time())
    staging = REPO / ".awm" / "data" / "runs" / f"bench_v3_{ts}"
    staging.mkdir(parents=True, exist_ok=True)

    print(f"=== staging {len(STAGED)} benchmark items ===", flush=True)
    inputs = build_inputs(staging)
    for t in sorted(STAGED):
        print(f"  lib   {t}")

    resources = [DataInstanceLibrary.Load(LIB / f"resources/{n}")
                 for n in ("containers", "lib")]
    transforms = [TransformInstanceLibrary.Load(LIB / f"transforms/{d}") for d in DOMAINS]

    if a.plan_only:
        agent = Agent(home=Source.FromLocal(staging / "agent_home"),
                      runtime=ContainerRuntime.APPTAINER)
    else:
        print(f"=== engine pin: {assert_pinned_engine()} ===", flush=True)
        agent_path = f"{a.scratch_root}/{a.user}/ecspr_bench_v3_{ts}"
        print(f"=== remote: {a.host}:{agent_path} ===", flush=True)
        container = agent_container()
        print(f"=== agent image: {container} ===", flush=True)
        agent = Agent(home=SshSource(host=a.host, path=agent_path).AsSource(),
                      runtime=ContainerRuntime.APPTAINER,
                      container=container,
                      setup_commands=SETUP_COMMANDS)

    print("=== planning ===", flush=True)
    targets = TargetBuilder()
    # The MERGED table only. Targeting the shards as well would let the planner
    # satisfy the merge from a separately-planned solve; one target, one chain.
    targets.Add("ecspr::benchmark_result")
    task = agent.GenerateWorkflow(
        samples=[inputs],
        resources=resources + [inputs],
        transforms=transforms,
        targets=targets,
    )
    if not task.ok:
        print("\nPLAN DID NOT RESOLVE. Planner hints:", file=sys.stderr)
        print(getattr(task.plan, "hints", task), file=sys.stderr)
        return 3

    print(f"resolved workflow: {len(task.plan.steps)} steps", flush=True)
    for i, step in enumerate(task.plan.steps):
        name = getattr(getattr(step, "transform", None), "name", None) or f"step{i}"
        print(f"  [{i}] {name}")

    dag = (REPO / "reports/dag/benchmark_v3").resolve()
    dag.parent.mkdir(parents=True, exist_ok=True)
    task.plan.RenderDAG(dag)
    print(f"DAG -> {dag.with_suffix('.svg')}", flush=True)

    if a.plan_only:
        print("\n--plan-only: nothing deployed, nothing run.")
        return 0

    print("=== Deploy() ===", flush=True)
    try:
        agent.Deploy()
    except subprocess.CalledProcessError as e:
        # One session, one failure, one message. NOT a retry loop -- each
        # attempt is a Duo push and a prior run here was halted by a lockout.
        print(f"\ndeploy failed ({e}). The connection is multiplexed: open ONE "
              f"session by hand (`ssh {a.host}`), leave it open, and re-run. "
              f"Do NOT delete the ControlMaster socket and do NOT retry in a "
              f"loop -- that is what causes an account lockout.", file=sys.stderr)
        return 4

    print(f"=== task key: {task.GetKey()} ===", flush=True)
    # on_exist="clear" is safe HERE and only here: agent_path carries a
    # timestamp, so it is a fresh directory every run and there is no prior
    # intermediate to destroy. Never carry this flag onto a resubmission.
    agent.StageWorkflow(task, on_exist="clear")
    agent.RunWorkflow(task)

    print(f"=== waiting (timeout {a.timeout_s / 3600:.1f}h, poll {a.poll_s:.0f}s) ===",
          flush=True)
    result = agent.WaitForWorkflow(task, timeout_s=a.timeout_s, poll_s=a.poll_s)
    print(f"=== status: {result['status']} after {result['elapsed_s'] / 60:.1f} min ===",
          flush=True)
    for line in result["tail"]:
        print(f"    {line}")
    if result["status"] != "completed":
        return 2

    src = agent.GetResultSource(task)
    out = (REPO / a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    print(f"=== result: {src.GetPath()} -> {out} ===", flush=True)
    print(f"\nnext: PYTHONPATH=src python transforms/build/benchmark/30_score_v3.py "
          f"--observations {out}/observations.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
