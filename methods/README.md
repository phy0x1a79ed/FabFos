# FabFos — the canonical methods path

All FabFos/ECSPr methods work lives here. This directory is the one place that
says how to run the method; the engine that implements it is a separate library
(`canon.ENGINE_LIB`, branch `feat/fabfos`).

## The one rule

**`canon.py` holds the basis. Prose may name a value; prose may not restate it.**

That includes this file. If you want to know how many inserts the basis has, what
the draw sizes are, which scorer is canonical, or where the axes live — import
`canon` and read the attribute. Do not write the number down here, in a docstring,
in a commit message, or in a figure caption.

This is not fussiness. The failure being corrected was four live documents each
declaring itself the canonical method and contradicting the others, because
nothing checked them. A fifth document would join that set. The only version of
this file that helps is one with nothing in it to go stale, and the rule is
enforced mechanically:

```bash
python check_no_transcribed_numbers.py
```

Any hit is a bug. The checker carves out exactly one exception, `MIGRATION.md`,
which has to name the dead.

## Layout

| Path | What it is |
|---|---|
| `canon.py` | The basis, as data, plus the two assertion helpers. The only source. |
| `check_no_transcribed_numbers.py` | Enforces the rule above. |
| `run_experiment.py` | The spine. One command runs an experiment. |
| `spec.py` | What an experiment IS: inputs + targets. A dataclass, not a framework. |
| `experiments/` | One module per experiment, each defining `SPEC`. |
| `dev_null_family.py` | Iterate a null-scoring family solo — no solver, no nextflow. |
| `parity/` | The gates. Run them; do not trust this file over them. |
| `MIGRATION.md` | What is retired, why, and what replaced it. |
| `SOCKEYE.md` | Cluster groundwork, and what is not done yet. |
| `ECSPR_METHOD_MAP.md` | Where the method's symbols live when run from a figure scope — cross-scope imports, solver, coercion, drivers. Read before re-deriving it. |

## Running things

```bash
python run_experiment.py scadc_main --generate    # plan + manifest + DAG, no run
python run_experiment.py scadc_main               # stage, plan, run, record
```

Add `--max-memory-gb`/`--max-cpus`/`--serial` on a workstation: transforms declare
what they want on a cluster, and nextflow refuses a step whose request exceeds what
the executor can see. Those flags scale the plan to the machine without editing the
transform, so the cluster request stays honest.

An experiment is one file in `experiments/`. It names inputs and targets and imports
`canon`. Nothing else.

## The gates

Run these rather than believing anything written here:

```bash
python parity/run_parity.py --all --python <env-with-sklearn>/bin/python
python parity/run_solver_parity.py          # the SOLVER vs a from-scratch rebuild
python parity/check_spine_output.py         # the whole spine vs the incumbent
python check_no_transcribed_numbers.py      # the rule at the top of this file
python pulsechase/run_pulsechase.py --all   # does ECSPr behave like a pulse-chase?
```

**These gates cover different halves, and the split is the point.** `run_parity.py` and
`check_spine_output.py` prove `(delta_obs, nulls) -> p/q` — the SCORER and the spine.
Neither runs the solver: they read `delta_obs` from `{lane}_axes_report.tsv`, and every
experiment spec stages `ecspr::{reff,ieff}_axes_report` as an INPUT (`scadc_main.py`
says so outright: "the staged, canonical solve"). `run_solver_parity.py` covers the
other half, `graph + evidence -> delta_obs`, which nothing checked before.

The first four are REPRODUCTION gates — they prove the engine reproduces the incumbent,
and every one stays green through the star's clique artifact, the MIN scorer, and the
refusal-gap. `pulsechase/run_pulsechase.py` is the missing BEHAVIOURAL gate: it asks
whether ECSPr behaves like the radio-labelled pulse-chase it claims to be, encoding the
two-stage contract (graph construction, then the `(network, source, sink) -> Ieff` solve)
as checks that are each able to fail. It runs the star against the atom graph and
deliberately contains checks the star fails *as designed* (cofactor-theft, input-leak,
role-respecting) — an affirmative artifact magnitude, not an absence of a pass.

The suite was first built to a red/green baseline with four checks *expected-red*, each
naming a correction; those corrections have since landed and the checks are now green
gates: `I3` (ambiguous atom correspondences are DILUTED by fanout, not refused), `II2` and
`II9` (Ieff is the all-paths super-node solve — combining parallel routes and edge-
preserving — replacing the MIN-over-pairs scorer), and `I9` (`Grid(prune=False)` retains
every component, so no reaction's atoms are silently dropped). Each stays *able to fail*:
break the property and it goes red. Its exit code keys ONLY on a green-expected atom check
regressing; the star-fails are the named baseline, not failures. Needs the full ECSPr
stack incl. rdkit — it gates its env and names any missing dependency
(`mamba run -n scadc-metabolic-model …`).

It scores against a **from-scratch rebuild**, not against the incumbent table — the
incumbent is the less accurate artifact (~1e-7 vs the rebuild, where the engine is
~1e-13). Gating the solver against it would pin a correct implementation to an
incorrect referent and call it parity. The incumbent's error is reported as
information; `--strict` opts into failing on it.

The `--python` is not optional and not a convenience. The env that runs this spine
carries metasmith, not scikit-learn, and the mixture fitter **degrades instead of
failing** without it — a single log-normal, flagged `no-sklearn`, diffed against an
incumbent that was fitted properly. That reports `FAIL` underneath a banner telling
you a real divergence is a finding about the incumbent, which would be exactly the
wrong conclusion. The gate now refuses that interpreter by name rather than
answering it.

Each is built to fail: the parity gate refuses a retired table by schema, the spine
gate fails when it finds nothing to compare, and the tripwire was probed with
deliberate violations. A gate that cannot fail is not a gate — which is the failure
mode this whole directory exists to remove.

## What none of these gates covered: direction

The reproduction gates and the pulse-chase GROUPS I–III all exercise the **undirected**
contract. Graph construction is symmetric, and the all-paths super-node solve is symmetric
in `source ↔ sink`: a reaction that runs one way in the cell and the other way only against
a steep thermodynamic gradient is scored identically in both directions. Directionality is
a correctness axis every one of those checks is structurally blind to — which is why
acceptance criteria that stop there are underspecified, not merely incomplete.

**The undirected solve is scientifically WRONG, not just incomplete** (see
`MIGRATION.md` 2026-07-20). Because it is symmetric, flux runs backward through
irreversible reactions against the thermodynamic gradient — connectivity the cell cannot
realise — so it systematically over-connects the graph. The **directed (diode) solve is
THE method**; the undirected path survives only as the *mathematical* symmetric-limit
target of the reduction gate below (`IV2`), never as a valid model. Never cite an
undirected result as a finding.

`pulsechase` GROUP IV closes that gap. It runs the engine's real rectified-network solver
(`ecspr_directed`, a committed primitive: signed incidence + per-edge forward/backward
conductance → `C_eff` by semismooth Newton) and pins the directed contract from both sides:
`IV1` — an irreversible reaction conducts forward but is throttled in reverse (≥ a committed
factor on a toy diode chain); `IV2` — at the symmetric limit `ratio == 1` the directed solve
reproduces the undirected `1/R_eff` exactly, so where direction is *unknown* the directed
model **is** the undirected model and any divergence in a comparison is real chemistry, not
solver machinery. The direction itself comes from the directionality ensemble's per-reaction
`g_rev/g_fwd = exp(ΔG'/RT)` ratio (`direction/combine.py`, staged content-hashed), wired onto
the base graph by `directed/directed_network.py`; a reaction with no evidence gets ratio 1.0,
a *proven* no-op.

**The behavioural gate still cannot see biology.** A green pulse-chase suite proves ECSPr
behaves like a pulse-chase — atoms conserved, all-paths measured, direction throttled — but
it says nothing about whether the *numbers track real metabolism*. That is a different
instrument: the X/Y validation benchmark (`fig-model`, `main/ecspr/benchmark/`) scores an
implementation's LOF/GOF observations against measured biology by directional-specificity
AUC. The benchmark can rank two behaviourally-valid solvers against ground truth; the gate
cannot. They are complementary halves — the gate proves the method is *self-consistent*, the
benchmark proves it is *right* — and neither substitutes for the other.

## Trying a new null

Null *generation* (style, K, N, seed) and null *scoring* (how a p-value comes off the
draws) are orthogonal. The scoring half is swappable today: add a family to `FAMILIES`
in the engine's `resources/lib/ecspr_null_study.py`, then

```bash
python dev_null_family.py --family <name>              # solo: no solver, no nextflow
FABFOS_NULL_FAMILY=<name> python run_experiment.py scadc_null_study   # same file, in the DAG
```

`RunTransform` routes through the same `ExecuteStep` the DAG uses, so what you
iterate on solo is what runs. The family is a staged, content-hashed input, so two
families cannot collide on one cache entry.

Alternative families produce `*_significance_study` tables. The canonical
significance type keeps exactly one producer and no staged input can change what it
means — that asymmetry is deliberate, and `MIGRATION.md` says why.

Null *generation* is **not** yet a transform: the draws are still a reused opaque
artifact rather than a cached instance of a declared computation.

## Which way the dependency runs

One way only:

```
main/fabfos/  ──imports──>  canon.py
     │
     └──stages inputs into──>  engine library (metasmith-libraries/fabfos)
```

The engine **must not** import `canon.py`, and must not contain a SCADC path. It
receives everything it needs as staged, content-hashed inputs. An engine that
knows about one experiment is an engine that cannot run a second one — that
reverse edge is the disease, and it is what put absolute paths from one branch
into the other and back.

## Using the assertion helpers

A figure generator, or anything else that reads a scored table, should refuse a
table that is not the canonical one rather than quietly plotting it:

```python
import canon
sig = canon.assert_canonical_significance(path)   # returns the DataFrame
axes = canon.assert_canonical_axes(canon.AXES_TSV)
```

These raise `canon.CanonError` — loudly, never a warning. A figure that asserted
the opposite of its own annotations has already shipped into the publish tree
once; the point of these helpers is that the same mistake becomes impossible to
make silently.

## Status

**"Canonical" is the current standard method to use — not a certificate that its
outputs have been validated.** `canon.py` declares the standard; declaring it is what
licenses validators to produce and pin its outputs. Validation (parity gates, produced
tables) *verifies* the standard — it never confers canonicity, and its absence never
withholds it. This ordering breaks the deadlock where implementers won't call a change
canonical without validation and validators won't validate what isn't canonical.

`canon.STATUS` records where the standard is; `canon.production_status()` reports whether
a given frozen output has been produced yet or is *pending production* (a legitimate state
— the run produces it). The ECSPr standard is the DIRECTED solve as of 2026-07-18
(`canon.CANONICAL_ORIENTATION`). Its directed reference solve + matched null were
**PRODUCED and frozen on 2026-07-20** with the fixed diode solver at the current draw grid
(`canon.production_status()` reports `pending=False`); the directed significance was scored
from them (`assert_canonical_significance` green, both lanes). See `MIGRATION.md` 2026-07-20.
