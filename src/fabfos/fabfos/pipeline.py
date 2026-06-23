"""Drive a FabFos run through metasmith.

Builds typed metasmith input libraries from the CLI arguments, asks the
planner for a workflow that produces fosmid scaffolds (and, when a vector
backbone is given, a pool-size estimate), then optionally executes it under
the chosen runtime. The default runtime is MAMBA so a conda/mamba install of
FabFos is self-contained — each tool runs via ``mamba run -n <env>`` with no
container relay.
"""
from __future__ import annotations

import json
import shutil
import time
from dataclasses import dataclass
from pathlib import Path

from metasmith.python_api import (
    Agent,
    Source,
    Runtime,
    DataTypeLibrary,
    DataInstanceLibrary,
    TransformInstanceLibrary,
    TargetBuilder,
    Resources,
    Size,
    Duration,
)

from .library import resolve_library_root, DOMAINS


@dataclass
class FabFosInputs:
    reads: Path                 # forward or interleaved reads (fastq[.gz])
    output: Path
    reverse: Path | None = None
    interleaved: bool = False
    background: Path | None = None     # host background genome
    vector: Path | None = None         # vector backbone (drives pool-size est.)
    end_forward: Path | None = None
    end_reverse: Path | None = None
    ends_facing: bool = False
    runtime: Runtime = Runtime.MAMBA
    threads: int = 8
    # ECSPr prerequisite: also build the per-fosmid bipartite metabolic graphs.
    # Requires the bundled/host reference inputs below.
    ecspr: bool = False
    base_graphs: Path | None = None        # dir of base_{C,N,S,P}.pkl
    element_bipartite: Path | None = None  # dir of mnx_bipartite_{C,N,S,P}.pkl
    reaction_db: Path | None = None        # dir of reactions.dmnd + bridge.tsv


def _parity(inp: FabFosInputs) -> str:
    return "single" if (inp.reverse is None and not inp.interleaved) else "paired"


def _write_read_metadata(staging: Path, parity: str) -> Path:
    p = staging / "read_metadata.json"
    p.write_text(json.dumps({"parity": parity, "length_class": "short"}, indent=2))
    return p


def _write_ends_table(staging: Path, inp: FabFosInputs) -> Path:
    p = staging / "ends_table.json"
    # insert ids are discovered downstream; the table primarily carries the
    # junction orientation flag the scaffold step needs.
    p.write_text(json.dumps({"insert_ids": [], "ends_facing": inp.ends_facing}, indent=2))
    return p


def build_inputs(inp: FabFosInputs, staging: Path, lib_root: Path):
    """Materialize the per-sample and per-run metasmith input libraries.

    Returns (samples_lib, resources_lib) mirroring the layout the fosmids
    transforms expect: read_metadata is the sample key, short reads are
    parented to it; host/vector/ends are per-run resources.
    """
    seq_types = DataTypeLibrary.Load(lib_root / "data_types/sequences.yml")
    fos_types = DataTypeLibrary.Load(lib_root / "data_types/fosmids.yml")

    samples_dir = staging / "samples"
    if samples_dir.exists():
        shutil.rmtree(samples_dir)
    samples = DataInstanceLibrary(samples_dir)
    samples.Purge()
    samples.AddTypeLibrary(namespace="sequences", lib=seq_types)
    meta = _write_read_metadata(staging, _parity(inp))
    meta_p = samples.AddItem(str(meta), "sequences::read_metadata")
    read_type = "sequences::short_reads_pe" if _parity(inp) == "paired" else "sequences::short_reads_se"
    samples.AddItem(str(inp.reads), read_type, parents={meta_p})
    samples.Save()

    res_dir = staging / "resources"
    if res_dir.exists():
        shutil.rmtree(res_dir)
    res = DataInstanceLibrary(res_dir)
    res.Purge()
    res.AddTypeLibrary(namespace="sequences", lib=seq_types)
    res.AddTypeLibrary(namespace="fosmids", lib=fos_types)
    if inp.background is not None:
        res.AddItem(str(inp.background), "sequences::background_genome")
    if inp.vector is not None:
        res.AddItem(str(inp.vector), "fosmids::vector_backbone")
    if inp.end_forward is not None and inp.end_reverse is not None:
        table = _write_ends_table(staging, inp)
        table_p = res.AddItem(str(table), "fosmids::end_sequences_table")
        res.AddItem(str(inp.end_forward), "fosmids::end_sequences_forward", parents={table_p})
        res.AddItem(str(inp.end_reverse), "fosmids::end_sequences_reverse", parents={table_p})

    if inp.ecspr:
        met_types = DataTypeLibrary.Load(lib_root / "data_types/metabolic.yml")
        res.AddTypeLibrary(namespace="metabolic", lib=met_types)
        if inp.base_graphs is not None:
            res.AddItem(str(inp.base_graphs), "metabolic::base_graphs")
        if inp.element_bipartite is not None:
            res.AddItem(str(inp.element_bipartite), "metabolic::element_bipartite")
        if inp.reaction_db is not None:
            res.AddItem(str(inp.reaction_db), "metabolic::reaction_reference_db")
    res.Save()
    return samples, res


def _targets(inp: FabFosInputs) -> TargetBuilder:
    targets = TargetBuilder()
    targets.Add("fosmids::scaffolds")
    if inp.vector is not None:
        targets.Add("fosmids::pool_size_estimate")
    if inp.ecspr:
        targets.Add("metabolic::fosmid_bipartite_graph")
    return targets


def generate_workflow(inp: FabFosInputs, staging: Path):
    """Plan-only: resolve the transform chain for the requested targets."""
    lib_root = resolve_library_root()
    samples, res = build_inputs(inp, staging, lib_root)
    resources = [DataInstanceLibrary.Load(lib_root / f"resources/{n}") for n in ("containers", "lib")]
    transforms = [TransformInstanceLibrary.Load(lib_root / f"transforms/{d}") for d in DOMAINS]
    agent = Agent(home=Source.FromLocal(staging / "agent_home"), runtime=inp.runtime)
    task = agent.GenerateWorkflow(
        samples=samples.AsSamples("sequences::read_metadata"),
        resources=resources + [res],
        transforms=transforms,
        targets=_targets(inp),
    )
    return agent, task


def _wait_for_run(staging: Path, task_key: str, timeout_s: int = 7200) -> Path:
    internals = staging / "agent_home" / "runs" / task_key / "_metasmith"
    deadline = time.time() + timeout_s
    last_log: Path | None = None
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
        time.sleep(5)
    raise TimeoutError(f"workflow {task_key} did not finish within {timeout_s}s")


def run_pipeline(inp: FabFosInputs) -> Path:
    """Plan, stage, execute, and collect results into ``inp.output``."""
    inp.output.mkdir(parents=True, exist_ok=True)
    staging = inp.output / "_fabfos"
    staging.mkdir(parents=True, exist_ok=True)

    agent, task = generate_workflow(inp, staging)
    assert task.ok, "FabFos workflow generation failed; see planner hints above"

    agent.Deploy()
    agent.StageWorkflow(task, on_exist="clear")
    agent.RunWorkflow(
        task,
        resource_overrides={
            "all": Resources(cpus=inp.threads, memory=Size.GB(max(2, inp.threads)), duration=Duration(hours=24)),
        },
    )
    _wait_for_run(staging, task._key)
    agent.CheckWorkflow(task)

    results = staging / "agent_home" / "runs" / task._key / "results"
    if results.exists():
        dest = inp.output / "results"
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(results, dest)
        return dest
    return results
