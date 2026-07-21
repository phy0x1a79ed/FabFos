#!/usr/bin/env python
"""The spine: one command to run a FabFos experiment.

    python run_experiment.py scadc_main --generate    # plan + DAG, no staging, no run
    python run_experiment.py scadc_main               # stage, plan, run, manifest

An experiment names its inputs and its targets (see `spec.py` and `experiments/`).
The spine stages the inputs, asks the engine for the targets, and writes a run
manifest. It is not ECSPr-specific: ECSPr significance is one target, recovery and
annotation are others.

WHAT THIS DELIBERATELY DOES NOT DO
----------------------------------
It does not run fresh annotation. The template this replaces re-ran annotation from
gated databases on every invocation, which changes the evidence, which changes the
numbers, which destroys the parity signal that makes any of this trustworthy. A
canonical driver stages an experiment's evidence as an INPUT. Fresh annotation is a
separate, later, deliberate question -- not a side effect of running the pipeline.

It also holds no paths of its own. Every path lives in the experiment spec, which
imports `canon`. The engine library holds none at all: that reverse edge -- engine
knowing about one experiment -- is what put absolute paths from each of two branches
into the other, and is the thing this whole exercise removes.

THE MANIFEST
------------
Every choice that changes a number is recorded automatically rather than remembered:
the spec digest, the engine library commit, the engine version, and a hash per input
file. If two runs disagree, diff their manifests.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
sys.stderr.reconfigure(line_buffering=True)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from fabfos import canon                                    # noqa: E402
from spec import DirInput, ExperimentSpec, Value  # noqa: E402

from metasmith.python_api import (              # noqa: E402
    Agent, ContainerRuntime, Duration, Resources, Size, Source,
    DataInstanceLibrary, DataTypeLibrary, TransformInstanceLibrary,
    TargetBuilder, METASMITH_VERSION,
)

RUNS = HERE / ".runs"


# =====================================================================
# hashing / provenance
# =====================================================================

def hash_file(path: Path, _chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(_chunk):
            h.update(chunk)
    return h.hexdigest()


def hash_path(path: Path) -> str:
    """sha256 of a file, or of a directory's (name, hash) listing."""
    path = Path(path)
    if path.is_dir():
        h = hashlib.sha256()
        for p in sorted(path.rglob("*")):
            if p.is_file():
                h.update(p.relative_to(path).as_posix().encode())
                h.update(hash_file(p).encode())
        return h.hexdigest()
    return hash_file(path)


def git_commit(repo: Path) -> str:
    try:
        r = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                           capture_output=True, text=True, check=True)
        dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain"],
                               capture_output=True, text=True).stdout.strip()
        return r.stdout.strip() + ("-dirty" if dirty else "")
    except Exception as e:
        return f"unknown ({e})"


def spec_digest(spec: ExperimentSpec) -> str:
    """Stable digest of the spec's MEANING -- targets, domains, input identities.

    Deliberately excludes mtimes and absolute staging paths: a spec that means the
    same thing must digest the same, or the manifest cannot answer "did anything
    that changes a number change?".
    """
    payload = {
        "name": spec.name,
        "targets": sorted(spec.targets),
        "domains": sorted(spec.domains),
        "namespaces": sorted(spec.namespaces),
        "per_experiment": sorted(spec.per_experiment),
        "inputs": {
            k: (
                {"value": v.content}
                if isinstance(v, Value) else
                {"dir": [p.name for p in v.files]}
                if isinstance(v, DirInput) else
                {"file": Path(v).name}
            )
            for k, v in sorted(spec.inputs.items())
        },
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode()).hexdigest()


# =====================================================================
# staging
# =====================================================================

def build_dir_input(staging: Path, type_name: str, item: DirInput) -> Path:
    """Symlink an explicit file list into a clean directory.

    Rebuilt each run so a removed entry actually disappears -- an accreting staged
    directory is the same failure as globbing a cache, one step later.
    """
    d = staging / "refs" / type_name.replace("::", "_")
    if d.exists():
        for f in d.iterdir():
            f.unlink()
    d.mkdir(parents=True, exist_ok=True)
    for f in item.files:
        if not f.exists():
            raise SystemExit(f"[{type_name}] curated file missing: {f}")
        (d / f.name).symlink_to(f)
    return d


def stage_inputs(spec: ExperimentSpec, staging: Path) -> DataInstanceLibrary:
    import shutil
    inputs_dir = staging / "inputs.xgdb"
    if inputs_dir.exists():
        shutil.rmtree(inputs_dir)
    inputs = DataInstanceLibrary(inputs_dir)
    inputs.Purge()
    for ns in spec.namespaces:
        inputs.AddTypeLibrary(
            namespace=ns,
            lib=DataTypeLibrary.Load(canon.ENGINE_LIB / f"data_types/{ns}.yml"))

    exp = inputs.AddValue(f"{spec.name}.txt", spec.name,
                          "fosmids::recovery_experiment")

    for type_name, item in spec.inputs.items():
        parents = {exp} if type_name in spec.per_experiment else None
        if isinstance(item, Value):
            inputs.AddValue(item.name, item.content, type_name, parents=parents)
        elif isinstance(item, DirInput):
            inputs.AddItem(build_dir_input(staging, type_name, item),
                           type_name, parents=parents)
        else:
            inputs.AddItem(Path(item), type_name, parents=parents)
    inputs.Save()
    return inputs


# =====================================================================
# manifest
# =====================================================================

def write_manifest(spec: ExperimentSpec, staging: Path, task, plan_steps) -> Path:
    man = {
        "experiment": spec.name,
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "spec_digest": spec_digest(spec),
        "engine": {
            "library": str(canon.ENGINE_LIB),
            "library_commit": git_commit(canon.ENGINE_LIB),
            "metasmith_version": METASMITH_VERSION,
        },
        "methods": {
            "methods_path": str(HERE),
            "methods_commit": git_commit(HERE),
            "canon_status": canon.STATUS,
            "canon_status_since": canon.STATUS_SINCE,
        },
        "targets": list(spec.targets),
        "domains": list(spec.domains),
        "plan_steps": plan_steps,
        "task_key": getattr(task, "_key", None),
        "inputs": {},
    }
    for type_name, paths in spec.resolve_paths().items():
        man["inputs"][type_name] = [
            {"path": str(p), "sha256": hash_path(p)} for p in paths
        ]
    for type_name, item in spec.inputs.items():
        if isinstance(item, Value):
            man["inputs"][type_name] = [{
                "value": item.content,
                "sha256": hashlib.sha256(item.content.encode()).hexdigest(),
            }]

    out = staging / "run_manifest.json"
    out.write_text(json.dumps(man, indent=2, sort_keys=True) + "\n")
    return out


# =====================================================================
# spine
# =====================================================================

def load_spec(name: str) -> ExperimentSpec:
    try:
        mod = importlib.import_module(f"experiments.{name}")
    except ModuleNotFoundError as e:
        avail = sorted(p.stem for p in (HERE / "experiments").glob("*.py")
                       if not p.stem.startswith("_"))
        raise SystemExit(f"no experiment [{name}] ({e}); available: {avail}")
    spec = getattr(mod, "SPEC", None)
    if not isinstance(spec, ExperimentSpec):
        raise SystemExit(f"experiments/{name}.py must define SPEC: ExperimentSpec")
    return spec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("experiment")
    ap.add_argument("--generate", action="store_true",
                    help="plan + DAG only; no staging, no run")
    ap.add_argument("--runtime", choices=[r.name for r in ContainerRuntime],
                    default=ContainerRuntime.APPTAINER.name)
    # A transform declares the resources it wants on a CLUSTER. Running the same
    # plan on a workstation needs those clamped -- nextflow refuses a step whose
    # request exceeds what the executor can see, and what it sees is the agent
    # container's limit, not the host's free memory. These flags scale the plan to
    # the machine WITHOUT touching the transform, so the cluster request stays
    # honest and the local run is not silently a different plan.
    ap.add_argument("--max-memory-gb", type=int, default=None,
                    help="clamp every step's memory (local runs)")
    ap.add_argument("--max-cpus", type=int, default=None,
                    help="clamp every step's cpus (local runs)")
    ap.add_argument("--serial", action="store_true",
                    help="one step at a time (local runs)")
    a = ap.parse_args()

    spec = load_spec(a.experiment)
    staging = RUNS / spec.name
    staging.mkdir(parents=True, exist_ok=True)

    print(f"=== experiment [{spec.name}] ===")
    print(f"canon: {canon.STATUS} since {canon.STATUS_SINCE}")
    print(f"spec digest: {spec_digest(spec)[:16]}")

    print("=== preflight ===")
    # Preflight FIRST: its checks give specific errors (and frame pending-production
    # state) before the generic missing-inputs list, which is the cue to PRODUCE the
    # canonical method's outputs, not a sign of misconfiguration.
    for check in spec.preflight:
        check()
        print(f"  ok: {getattr(check, '__name__', check)}")
    spec.check_inputs_exist()

    print("=== staging inputs ===")
    inputs = stage_inputs(spec, staging)

    print("=== loading resources & transforms ===")
    resources = [DataInstanceLibrary.Load(canon.ENGINE_LIB / f"resources/{n}")
                 for n in ("containers", "envs", "lib")]
    transforms = [TransformInstanceLibrary.Load(canon.ENGINE_LIB / f"transforms/{d}")
                  for d in spec.domains]

    print("=== generating workflow ===")
    targets = TargetBuilder()
    for t in spec.targets:
        targets.Add(t)
    agent = Agent(home=Source.FromLocal(staging / "agent_home"),
                  runtime=ContainerRuntime[a.runtime])
    task = agent.GenerateWorkflow(
        samples=list(inputs.AsSamples("fosmids::recovery_experiment")),
        resources=resources + [inputs],
        transforms=transforms,
        targets=targets,
    )
    if not task.ok:
        print(f"PLANNING FAILED: {task}")
        return 1

    plan_steps = []
    print(f"plan has {len(task.plan.steps)} steps")
    for step in task.plan.steps:
        name = Path(step.transform._path).stem
        prods = [i.dtype_name for g in step.produces for i in g]
        plan_steps.append({"order": step.order, "transform": name, "produces": prods})
        print(f"  step {step.order}: {name} -> {prods}")

    man = write_manifest(spec, staging, task, plan_steps)
    print(f"=== manifest -> {man} ===")

    dag = staging / f"{spec.name}_dag.svg"
    try:
        task.plan.RenderDAG(dag)
        print(f"DAG -> {dag}")
    except Exception as e:
        print(f"DAG render failed (non-fatal): {e}")

    if a.generate:
        print("[--generate] plan + manifest + DAG only; skipping stage/run.")
        return 0

    print("=== staging workflow ===")
    agent.Deploy()
    agent.StageWorkflow(task, on_exist="update", verify_external_paths=False)

    print("=== running ===")
    overrides = None
    if a.max_memory_gb or a.max_cpus:
        overrides = {"all": Resources(
            cpus=a.max_cpus or 2,
            memory=Size.GB(a.max_memory_gb or 6),
            duration=Duration(hours=12),
        )}
        print(f"  resource override (all steps): cpus={a.max_cpus or 2} "
              f"memory={a.max_memory_gb or 6}GB")
    params = (dict(executor=dict(queueSize=1), process=dict(maxForks=1))
              if a.serial else None)
    agent.RunWorkflow(task, params=params, resource_overrides=overrides)
    print("submitted; monitor at",
          staging / "agent_home" / "runs" / task._key / "_metasmith" / "logs.latest")
    return 0


if __name__ == "__main__":
    sys.exit(main())
