"""FabFos command line — a thin front end over the metasmith planner.

``fabfos`` builds a metasmith workflow that resolves fosmid inserts from
pooled reads and runs it through the chosen runtime (mamba by default, so a
conda install is self-contained). Use ``--plan-only`` to just resolve and
print the DAG without executing.
"""
import argparse
import multiprocessing
import sys
from pathlib import Path

from metasmith.python_api import Runtime

from . import __version__, NAME, SHORT_SUMMARY
from .pipeline import FabFosInputs, generate_workflow, run_pipeline


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog=NAME, description=SHORT_SUMMARY)
    io = p.add_argument_group("inputs")
    io.add_argument("-r", "--reads", metavar="FASTQ", required=True,
                    help="forward (paired) or interleaved or single-end reads, fastq[.gz]")
    io.add_argument("-2", "--reverse", metavar="FASTQ", default=None,
                    help="reverse reads for paired-end input")
    io.add_argument("-i", "--interleaved", action="store_true", default=False,
                    help="--reads holds interleaved paired-end reads")
    io.add_argument("-o", "--output", metavar="PATH", required=True,
                    help="output directory")

    fos = p.add_argument_group("fosmid pool")
    fos.add_argument("-b", "--background", metavar="FASTA", default=None,
                     help="host background genome to filter out")
    fos.add_argument("--vector", metavar="FASTA", default=None,
                     help="vector backbone fasta; enables pool-size estimation")
    fos.add_argument("--endf", metavar="FASTA", default=None,
                     help="forward-junction end sequences fasta")
    fos.add_argument("--endr", metavar="FASTA", default=None,
                     help="reverse-junction end sequences fasta")
    fos.add_argument("--ends-facing", action="store_true", default=False,
                     help="end sequences face inward across both junctions")

    ec = p.add_argument_group("ecspr (metabolic graph prerequisite)")
    ec.add_argument("--ecspr", action="store_true", default=False,
                    help="also build per-fosmid bipartite metabolic graphs (stops before the axis/conductance step)")
    ec.add_argument("--base-graphs", metavar="DIR", default=None,
                    help="reference dir of base_{C,N,S,P}.pkl host graphs")
    ec.add_argument("--element-bipartite", metavar="DIR", default=None,
                    help="reference dir of mnx_bipartite_{C,N,S,P}.pkl universe graphs")
    ec.add_argument("--reaction-db", metavar="DIR", default=None,
                    help="reference dir with reactions.dmnd + bridge.tsv")

    run = p.add_argument_group("execution")
    run.add_argument("--runtime", choices=[r.value for r in Runtime], default=Runtime.MAMBA.value,
                     help="tool runtime (default: mamba)")
    run.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(),
                     help="max threads per step")
    run.add_argument("--plan-only", action="store_true", default=False,
                     help="resolve and print the workflow DAG without executing")
    p.add_argument("-v", "--version", action="version", version=f"{NAME} {__version__}")
    return p


def _inputs_from_args(a: argparse.Namespace) -> FabFosInputs:
    return FabFosInputs(
        reads=Path(a.reads).resolve(),
        output=Path(a.output).resolve(),
        reverse=Path(a.reverse).resolve() if a.reverse else None,
        interleaved=a.interleaved,
        background=Path(a.background).resolve() if a.background else None,
        vector=Path(a.vector).resolve() if a.vector else None,
        end_forward=Path(a.endf).resolve() if a.endf else None,
        end_reverse=Path(a.endr).resolve() if a.endr else None,
        ends_facing=a.ends_facing,
        runtime=Runtime(a.runtime),
        threads=a.threads,
        ecspr=a.ecspr,
        base_graphs=Path(a.base_graphs).resolve() if a.base_graphs else None,
        element_bipartite=Path(a.element_bipartite).resolve() if a.element_bipartite else None,
        reaction_db=Path(a.reaction_db).resolve() if a.reaction_db else None,
    )


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    inp = _inputs_from_args(args)

    if args.plan_only:
        inp.output.mkdir(parents=True, exist_ok=True)
        _agent, task = generate_workflow(inp, inp.output / "_fabfos")
        if not task.ok:
            print("workflow generation FAILED; planner hints:", file=sys.stderr)
            print(task.plan.RenderHints() if hasattr(task.plan, "RenderHints") else task, file=sys.stderr)
            return 1
        print(f"resolved workflow: {len(task.plan.steps)} steps")
        for i, step in enumerate(task.plan.steps):
            name = getattr(getattr(step, "transform", None), "name", None) or f"step{i}"
            print(f"  [{i}] {name}")
        return 0

    run_pipeline(inp)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
