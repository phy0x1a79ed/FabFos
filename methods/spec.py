"""What an experiment IS, to the spine.

An experiment names its inputs and its targets. That is the whole contract. The
spine stages the inputs, asks the engine for the targets, and writes a manifest.

The spine is deliberately NOT ECSPr-specific. ECSPr significance is one target;
recovery and annotation are others. An experiment that wants a different question
answered names different targets, not a different driver.

This is a dataclass, not a config framework. There is no YAML, no loader, no
plugin registry. An experiment spec is a Python module that builds one of these
and calls it SPEC, so it can import `canon` and compute rather than transcribe --
which is the entire point (see canon.py).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable


@dataclass(frozen=True)
class DirInput:
    """A DIRECTORY-typed input, curated from an EXPLICIT list of files.

    Never a glob. Every directory this pipeline reads from -- the null cache
    above all -- holds retired artifacts beside canonical ones, so a glob does not
    mean "everything relevant", it means "whatever is lying around". The staged
    directory is built by symlinking exactly `files`, and a missing one is an
    error rather than a smaller directory.
    """
    files: tuple[Path, ...]

    def __init__(self, files):
        object.__setattr__(self, "files", tuple(Path(f) for f in files))


@dataclass(frozen=True)
class Value:
    """A small literal input, written to a file so it can be content-hashed.

    This is how a knob that changes a number enters the task hash. It is NOT
    `context.params` -- params come only from the runtime resources line and never
    reach the hash, so a params-borne knob is both silently fixed at its default
    and invisible to the cache.
    """
    name: str
    content: str


@dataclass
class ExperimentSpec:
    name: str

    # type_name -> path | DirInput | Value.
    # Inputs keyed by a type the engine declares; the spine does not interpret them.
    inputs: dict[str, Path | DirInput | Value] = field(default_factory=dict)

    # type_names staged under the experiment (parents={exp}) rather than globally.
    per_experiment: frozenset[str] = frozenset()

    # Engine targets, as "namespace::type".
    targets: tuple[str, ...] = ()

    # Transform domains to LOAD. Shrink this deliberately: a domain that is loaded
    # can be selected structurally by the planner, so leaving one in is how you get
    # a lane you did not ask for -- e.g. a fresh annotation lane silently changing
    # the evidence, which changes the numbers, which destroys the parity signal.
    domains: tuple[str, ...] = ()

    # Type-library namespaces to attach to the staged inputs.
    namespaces: tuple[str, ...] = ()

    # Preflight assertions, run BEFORE planning. Each raises on a bad basis.
    # This is where an experiment asserts what it means -- e.g. the axis count --
    # rather than trusting a filename.
    preflight: tuple[Callable[[], None], ...] = ()

    def resolve_paths(self) -> dict[str, list[Path]]:
        """Every real file this spec depends on, per input type. Manifest input."""
        out: dict[str, list[Path]] = {}
        for type_name, item in self.inputs.items():
            if isinstance(item, DirInput):
                out[type_name] = list(item.files)
            elif isinstance(item, Value):
                continue                      # hashed from content, not from disk
            else:
                out[type_name] = [Path(item)]
        return out

    def check_inputs_exist(self) -> None:
        """Fail before planning, naming every missing path at once.

        A driver that dies on the first missing file makes you rerun it once per
        missing file. This is cheap and the whole list is more useful.
        """
        missing = [
            f"  {type_name}: {p}"
            for type_name, paths in self.resolve_paths().items()
            for p in paths if not p.exists()
        ]
        if missing:
            raise SystemExit(
                f"experiment [{self.name}] is missing {len(missing)} input(s):\n"
                + "\n".join(missing))
