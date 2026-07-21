# Migration — what is retired, and what replaced it

The tombstone. This is the one file under `main/fabfos/` allowed to write a
retired value as a literal, because a tombstone has to name its dead and a dead
value cannot drift. Live basis values still may not appear here — import them
from `canon.py` (see `README.md`).

Nothing on this list may be quietly resurrected. Each entry says what dies, why,
and what to use instead. Where an artifact still exists on disk, it exists as a
**frozen referent** — readable, never authoritative.

## Status of this migration — the handoff happened

**2026-07-14: the parity gate went green and the entries below took effect.**

The engine reproduces the incumbent on both lanes to ~1e-16 — machine epsilon —
on every one of the joined rows, with exact `at_floor` agreement and the split
contigs present in the output. Re-provable in seconds, on CPU:

```bash
python parity/run_parity.py --all
```

That handoff is the point of the whole exercise. Without it we would not have
replaced the competing claims; we would have added a fifth. It is dated here and
in `canon.STATUS`, and it means two things at once:

| Handoff | State |
|---|---|
| Parity gate green, both lanes | **done 2026-07-14** |
| `canon.STATUS` → CANONICAL | **done 2026-07-14** |
| `04_reaction_network` → frozen referent | **in effect 2026-07-14** |

The gate proves **reproduction, not truth**. It shows the engine reproduces
today's incumbent. It does not revisit whether the incumbent is right.

The gate is also demonstrably able to fail, which is the only thing that makes a
green one worth anything. Pointed at a retired sibling table it refuses by schema;
with the retired ORF-id regex it loses exactly the split contigs — see below.

---

## 2026-07-17: the canonical graph moved to the honest reference (T1)

**The solver graph is no longer the incumbent star cache. It is the frozen
MetaNetX-4.5 reference graph** (`reference.GRAPH_DIR`; the reference-fed solve at
`reference.REF_ROOT/solve`). The star mapped ~27k carbon transits that no atom map
supports (electron/photon acceptors given atom edges, name-trust fabrications,
pseudo-reactions); the reference refuses them. A prior pass measured the swap
AUC-neutral, so this is a **correctness / honesty promotion** — the honest graph is
the one we publish on — not a score change. The star run stays on disk as a frozen
referent.

What moved, in `canon.py`:

| pin | was | now |
|---|---|---|
| `BIPARTITE_DIR` | `INCUMBENT_CACHE` | `REFERENCE_GRAPH_DIR` |
| `SOLVE_BASE_DIR` | `INCUMBENT_CACHE` | `REFERENCE_SOLVE_DIR` |
| `AXES_TESTABLE_JSON` | incumbent cache | reference `solve/axes_testable.json` |
| staged solve resolver | `incumbent_axes_report()` | `reference_axes_report()` (new) |

`incumbent_axes_report()` / `incumbent_sig_table()` are **kept** as frozen historical
referents: the scorer-parity gate still regresses the frozen incumbent inputs to the
frozen incumbent outputs (a check of the mixture-SF port, independent of which graph
is canonical), so that gate never went red.

**Testability, not the axis set, shrank.** The axis *set* is unchanged (assert against
`canon.AXES_JSON`, the graph-independent definition). On the honest graph one carbon
axis (a phospholipid endpoint the star reached only through a fabricated transit) is no
longer in the LCC, so `AXES_TESTABLE_JSON` is now a strict subset of the set — asserted
by count against `AXES_JSON`, never against the testable subset.

**Gates that proved the move (all green before the commit):**

- `parity/run_solver_parity.py` — the engine solve on the reference graph reproduces a
  from-scratch dense rebuild (shares no math) to ~1e-14 on all four elements.
- `pulsechase/run_pulsechase.py --all` — 27 checks, 0 regressions on the reference base,
  including the directed symmetric-limit (IV2) and the reference-honesty check (V1: zero
  refused reactions carry a fabricated transit pair).
- `parity/run_parity.py --all` — the scorer regression, still green (untouched inputs).
- `check_no_transcribed_numbers.py` + `assert_canonical_reference()` — green.

The staged **experiment** (`scadc_main`) is repointed onto the reference solve + a
reference-graph null **together** in the next step, so `delta_obs` and its null are
never on different graphs.

---

## 2026-07-18: the null moved to the reference graph, and the axis set was re-scored (T2)

The null now runs on the **same honest reference graph** as the observed solve. The draws
themselves are graph-independent ORF lists, reused byte-for-byte from the frozen cache;
only the SMW pass was re-run, pointed at the reference universe graph
(`reference/build_reference_null.py`, the redirect mirroring `_ref_patch`). It ran on
**Sockeye** (V100, esmc.sif CUDA build) and was proven a faithful port of the local
generator: Sockeye vs local agree to ~1e-13 absolute, row-for-row, on every shared draw
size. `canon.REFERENCE_NULL_DIR` pins the result; `experiments/scadc_main.py` now stages
`reference_axes_report` + that null **together**, so `delta_obs` and its null are on one
graph.

**The canonical axis set (`canon.AXIS_SET`) was re-scored on the reference graph**
(`reference/finish_t2_significance.py`,
canonical mixture scorer). The promotion is validated **by effect size**, which is how
ECSPr significance is read — the mixture `survives` q<0.05 flag saturates on the heavy
tail and is a diagnostic, not the ranking:

| lane | axis-level effect-size Spearman | top-20 axes shared |
|---|---|---|
| reff | 0.98 | 19/20 |
| ieff | 0.97 | 18/20 |

The honest graph refuses the fabricated carbon transits yet **preserves the axis-level
ranking the paper reads** — the honesty promotion is signal-neutral, as the earlier
AUC-neutral pass predicted. The saturating survivor-flag Jaccard is near zero, which is
the flag reshuffling a near-threshold set (survivorship is axis-level; all fosmids on an
axis share one `delta_obs`), **not** a divergence in the result. Rank by effect size, not
the flag — the same rule the negbin/`survives` retirements below rest on.

---

## 2026-07-18: the canonical ECSPr solve IS now DIRECTED (T3)

**The directed (diode) solve is the current standard method.** "Canonical" here is the
standard method to use, not a certificate that its frozen outputs have been produced (see
canon.py's header). The directed solve is declared canonical because the honest atom graph
IS directed; producing its frozen outputs is validation OF that standard, licensed BY the
declaration, and is tracked below as pending — it does not gate the declaration.

**The canonical solve gains an ORIENTATION, and the honest one is DIRECTED.** The
undirected solve treats every reaction as a symmetric two-way conductance; the atom
graph it runs on is in fact DIRECTED, and the direction annotator (`canon.DIR_TABLE`)
already gives every base-graph reaction a conductance ratio `g_rev/g_fwd = exp(dG'/RT)`.
The directed (diode) solve orients the base per reaction roles + that ratio and solves
each `(fosmid, axis)` as a rectified-network `C_eff`. It emits the SAME reff/ieff schema,
so the directed `delta_ieff` flows through `canon.CANONICAL_LANE` and the significance
scorer UNCHANGED — no new lane, no new column.

**This is a correctness / completeness promotion, NOT a score chase.** Measured on the
assembled directed prototype: the swap is net-neutral (ΔAUC −0.005) though it moves
~17.5% of cells, and the AAM/direction ensemble agrees with the curated reference on
93.9% of reactions. As with the honest-graph promotion (T1), the honest model preserves
the ranking the paper reads; we adopt it because it is the honest model, not because it
scores better.

**The symmetric limit proves reduction to the undirected answer.** `solve-directed
--force-noop` forces every ratio to 1.0, degrading the diode model to the undirected
one; it must reproduce the undirected reference `axes_report` to solver tolerance. That
parity gate is what licenses calling the directed solve a strict generalization rather
than a different method.

**The undirected solve is RETAINED, but only as a mathematical referent — it is
scientifically WRONG as a model** (see the 2026-07-20 clarification below).
`reference_axes_report(orientation="undirected")` and `transforms/ecspr/solve.py` stay
intact and readable *for the symmetric-limit reduction gate only*; the directed tables are
siblings beside them (`_directed`-suffixed). An undirected result is never a valid finding.

What landed (the whole directed method + its wiring):

| piece | where |
|---|---|
| engine domain `ecsprDirected` (its own `solve_directed.py`); undirected `solve.py` split into its own `ecsprUndirected` domain | separate domains so the planner never coin-flips between two producers of `ecspr::{reff,ieff}_axes_report` — the ecsprNetA precedent |
| generic input types `ecspr::reaction_roles` + `ecspr::direction_ratios` | `data_types/ecspr.yml`; the two network-agnostic inputs `solve-directed` needs beyond the undirected solve |
| `canon.CANONICAL_ORIENTATION` (+ `ORIENTATIONS`) | dispatches `reference_axes_report()` and `reference_null_files()`; **set to `"directed"` — the standard** |
| `canon.reference_axes_report_directed()` / `reference_null_files()` / `production_status()` | name the directed frozen outputs + matched null, and report whether they are produced yet |
| Network A flipped to directed via native GEM flux bounds | `scadc_netA` + `transforms/ecsprNetA/_gem_direction.py` (`canon.NETA_DIR_TABLE`) — reversible=1.0, irreversible=diode |

**The directed solver defect is FIXED (2026-07-18).** The engine solver
(`resources/lib/ecspr_directed.py`) carried a non-reproducible defect on backflow axes:
warm≠cold, machine-dependent (e.g. `cofactor__MNXM738068__MNXM256__S` gave 0.092 vs 0.070
across starts), passing only the toy/symmetric gates that never exercise the real backflow
path. Root cause: active-set CHATTERING at the diode kink — edges sitting at `x≈0` flip
their branch stiffness by ×1e9 (`g+` vs the `1e-9·g+` floor) between Newton iterates, so the
hard-diode Newton cycles and stalls at a start-dependent point. Fix: a CONSISTENT
softplus-smoothed diode (width `DIODE_SMOOTH_DELTA=1e-6`) makes the conductance vary smoothly
and the energy strictly convex + C², so Newton reaches the UNIQUE minimiser
start-independently; globalised by an energy-Armijo line search with an energy-stagnation
convergence flag; SPD systems factored by CHOLMOD (symbolic reuse) with an `splu`+ridge
fallback. At `g-=g+` the smoothing is an exact no-op, so the symmetric-limit parity with the
undirected `R_eff` is preserved to machine precision (integration gate `directed_parity.py`:
4e-13 on the real C/N/S/P graphs). Validated on the real element-S reference: the worst
defect axis went 36% → 3.7e-7 reproducible; all 7 S axes `conv=True`; CHOLMOD==splu 4.4e-13.

**The frozen directed outputs are PENDING PRODUCTION (not a gate on the standard).**
`CANONICAL_ORIENTATION = "directed"` is set; `canon.production_status()` reports whether the
directed outputs are on disk yet. Two must be produced with the fixed solver: (a) the directed
OBSERVED reference `axes_report` (cheap — minutes locally) and (b) a matched directed NULL. The
previously-computed directed null (`REF/null_directed/`, artifact #860) was generated by the
DEFECTIVE solver and is CONTAMINATED on exactly the backflow axes the fix repairs, so it cannot
be reused — the null is re-run (HPC, ~1h on Sockeye; the frozen graph/solver-independent draws
are reused verbatim). Until both are frozen+hash-pinned, a consumer that stages them (e.g.
`scadc_main`) fails loudly listing the missing files — that is the cue to PRODUCE them, not a
sign the standard is wrong. `delta_obs` and its null are already on the directed orientation in
lockstep via the resolvers; producing the tables is the remaining validation work.

---

## 2026-07-20: undirected is SCIENTIFICALLY WRONG; the null draw grid moved to the CLEAN-rebuild ORF-percentile sizes

Two changes, recorded together because the directed production run that closes the
pending-production state (T3 above) carries both.

**1. The undirected solve is scientifically WRONG, not merely superseded.** The undirected
model treats every reaction as a symmetric two-way conductance, so flux runs backward
through irreversible reactions against the thermodynamic gradient — connectivity and
conductance the cell cannot realise. Metabolism is directional; letting an irreversible
reaction conduct in reverse is not a rougher approximation, it is the **wrong answer**, and
it systematically **over-connects** the graph. Earlier notes (T1, T3) framed the swap to the
directed model as an *honesty / completeness* promotion measured "net-neutral" on the
benchmark; that undersold it. The directed (diode) solve is **THE** method. The undirected
path survives for exactly ONE non-scientific role: the **mathematical symmetric-limit target**
of the reduction gate (`directed/directed_parity.py`), where directed at ratio→1 must
reproduce `1/R_eff` exactly — a check that the directed solver *generalises the resistor*,
not a claim that the resistor models biology. `canon.ORIENTATIONS` /
`reference_axes_report(orientation="undirected")` stay reachable for that gate; never cite an
undirected number as a finding.

**2. The null draw grid moved to the current CLEAN-rebuild ORF-percentile sizes** — the sizes
the DIRECTED production null is generated at. `canon.DRAW_SIZES` now holds them (import the
name; the value is not restated here).

**Retired grids (the tombstone names its dead):**
- `(14, 28, 34, 42, 51)` — the prior length-anchored grid. The directed production null moved
  off it; it is superseded.
- `(21, 34, 56)` — the older incumbent null GENERATOR grid, already retired above.

Both now live in `canon.RETIRED_DRAW_SIZES`.

**Consequence, recorded honestly.** The scorer/spine parity gates that regress the **incumbent**
tables (`parity/run_parity.py`, `parity/check_spine_output.py`) join incumbent nulls at
`canon.DRAW_SIZES`, and no incumbent null exists at the new grid — those gates cannot join at
the new grid and are **no longer green under it**. This is expected and
acceptable: the incumbent (undirected, old-grid) path is the *wrong, superseded* one, so a
green regression against it is not a property worth preserving. The gate that matters for the
directed method is `directed/directed_parity.py` (the directed symmetric-limit reduction),
which is grid-independent. The incumbent-parity gates are frozen historical checks of a
retired scorer port; if a future task needs them green, pin them to `RETIRED_DRAW_SIZES`
explicitly rather than reading `DRAW_SIZES`.

**In progress at this migration.** The directed OBSERVED solve (the full fosmid basis × the
canonical axis set) and the directed NULL are being produced on Sockeye at the new grid with
the fixed diode solver
(engine `ecspr_directed.py` sha `e16c8584…`, commit 35ea3ee). Until frozen + hash-pinned into
the `_directed` canon paths, `canon.production_status()` reports pending — the cue to finish
the run, not a defect in the standard.

## Retired: the negbin scorer

`04_reaction_network/reff/05r_reff_significance.py`.

**Why.** Wrong distributional family for power-law-tailed data; measured 1.9–2.7×
anti-conservative. Not fixable by tuning — the family is wrong.

**Instead.** `canon.SCORER`. See `canon.incumbent_sig_table()`.

**Note.** The published `sig_negbin_*` tables predate the background-selection
guard (built before it landed) and today's code reproduces them on only ~90% of
cells. Any number quoted from them is unreproducible. Do not cite them.

## Retired: `sig_emp_*` and `sig_negbin_*` as separate tables

**Why.** These were never competing tables. The canonical scorer emits the
empirical-CCDF p-value as a **column on every row** (`p_emp`), fed the same raw
pool. A separate empirical table adds a second thing to accidentally read, which
is exactly how three figure generators ended up on three different null tables
across two days of scorer swaps, one of them shipping a title that asserted the
opposite of its own annotations.

**Instead.** One table, `canon.SCORER`; cross-check a survivor against the
`p_emp` column on its own row. Both are diagnostics of the same row, not rival
answers. `canon.assert_canonical_significance()` refuses the retired tables by
name.

## Retired: draw sizes {21, 34, 56}

**Why.** Superseded by `canon.DRAW_SIZES` (length-anchored, from the CLEAN
rebuild). This one is the sharpest illustration of the disease: the retired sizes
are still declared by the null **generator** while the **scorer** in the same
directory had already moved to the new set. A constant forked inside a single
directory, with nothing failing.

**Instead.** `canon.DRAW_SIZES`, and no `DRAW_SIZES` constant in the engine at
all — the engine derives sizes from its staged nulls directory, which the driver
builds from an explicit list (`canon.FROZEN_NULL_FILES`).

**Hazard.** The live cache still holds the retired sizes beside the canonical
ones, plus an un-suffixed draws file left by the overwrite incident. **Never glob
the cache.** Even the hand-made byte-copy backup (`canon.INCUMBENT_K1000`) holds
retired sizes — it is a safer source, not a curated one.

## Retired: the set2cat axes

**Why.** Superseded by `canon.AXIS_SET`. The retired set is what the library's
data-type description and the incumbent E2E driver still name.

**Instead.** `canon.AXES_TSV` / `canon.AXES_JSON`, asserted by count and
per-element split via `canon.assert_canonical_axes()` — never by filename, since
the filename is what drifted.

## Retired: both chain drivers

`04_reaction_network/run_reff_chain.sh` and `run_set2_t3_chain.sh`.

**Why.** Both invoke a superseded scorer (`run_reff_chain.sh` → the negbin one;
`run_set2_t3_chain.sh` → the background-partition scorer on the betweenness
lane), and both operate on the retired axis set. Reproducing today's published
state from either is not possible.

**Caveat, recorded honestly.** Neither carries an in-repo retirement marker, and
at least one external context doc still references them as current. They are
retired by this document, and that is the first time it has been written down.

**Instead.** The spine: one experiment driver, `run_experiment.py` (T3).

## Retired: the `matched_N` column

**Why.** It belonged to the nearest-size scorer. The canonical scorer interpolates
survival functions across flanking anchors, so there is no single matched size;
the bracket columns replace it. Its presence in a table is positive evidence the
table is stale, and `canon.assert_canonical_significance()` treats it that way.

**Instead.** `canon.SIG_COLUMNS`.

## Retired: the `\w`-based ORF-id regex

The engine's ORF counter used `^>(\w+?)_\d+`.

**Why.** `\w` excludes the dot, so an ORF on a split contig matches nothing — and
the counter tested the match before using it, so those ORFs were **dropped rather
than erroring**. The affected fosmids then had no ORF count, were filtered out of
the observation set, and vanished from the scored table. No exception, no warning,
no row.

**Measured.** On the canonical basis the retired regex counts 195 contigs where
the `rsplit` form counts `canon.FOSMID_BASIS`. The missing four are exactly
`canon.SPLIT_CONTIGS`.

**Instead.** The `rsplit("_", 1)` form, which the canonical scorer already used
and already carried a comment warning about precisely this. The parity gate now
asserts `canon.SPLIT_CONTIGS` are present in the output **positively** — a silent
absence cannot be noticed, so it has to be checked rather than watched for.

## Retired: `context.params` as a carrier of science

**Why.** `params` is populated solely from the runtime resources metadata line
(`cpus`/`memory`/`attempt`). A transform reading `context.params.get("device")`
therefore *always* gets the default — the GPU path was unreachable through the
driver, which is why the last end-to-end run bypassed the engine with raw sbatch.
Worse, params do not enter the task hash, so two runs intending different knobs
would collide on one cache entry. That is the overwrite bug, relocated into the
engine.

**Instead.** The rule the design rests on: **anything that changes a number is a
staged, content-hashed input.** A GPU run and a CPU run of the same step must
produce different cache keys. In the engine, device/dtype are now the staged
`ecspr::compute_profile` type; the ECSPr transforms carry no `context.params`.

**The bug class is wider than ECSPr — not fixed here.** Found while doing the
above: at least eight other transforms in the engine carry *science* knobs on
`context.params`, so each is silently pinned at its default and invisible to the
cache key. Legitimate uses (`cpus`, `memory`, `attempt`) are excluded — these are
not those:

| Transform | Knob(s) borne by params |
|---|---|
| `functionalAnnotation/uniref_lane.py` | `max_target_seqs`, `evalue` |
| `functionalAnnotation/dl_ec_lane.py` | `device` |
| `functionalAnnotation/proteinbert_embed.py` | `device` |
| `functionalAnnotation/compile_evidence.py` | `source` |
| `fosmids/cluster_contigs.py` | `min_contig_length`, `cluster_identity` |
| `fosmids/select_reference_inserts.py` | `min_insert_length` |
| `assembly/filter_contigs_by_length.py` | `min_length`, `max_pident` |

Each is the same defect with the same two consequences. Deliberately left alone:
they are outside the ECSPr path, none is on the canonical chain, and fixing them
blind — without a parity gate on their lanes of the kind T2 built for this one —
would be the same mistake at a different address. Fix them when a driver needs
them, with a gate.

---

## Frozen referents — read, never write

`04_reaction_network/` becomes read-only once the gate is green. It is **not**
deleted: the parity gate joins against it, so it must stay intact. It is the
thing being reproduced, not a thing to run.

The gate proves reproduction, not truth. It shows the engine reproduces today's
incumbent. It does not revisit whether the incumbent is right.

---

## Known inversions, not yet fixed

Recorded so they are not rediscovered as surprises.

- **The axis table lives in the publish *output* tree while being read as a
  pipeline *input*** (`canon.AXES_TSV`), by many scripts. It has no writer in the
  repo — it is user-pointed, so publish is currently its only home.
- **Result tables are tracked in git while load-bearing scripts are not** —
  including the lane that produced the current basis.
- **The engine's own premises have drifted** (found while grounding this work,
  detailed in the T2 notes): its vendored mixture fitter is a pre-guard copy of
  the incumbent's, and its E2E driver imports a runtime that no longer exists in
  the engine it targets.
