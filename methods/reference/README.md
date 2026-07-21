# The frozen ECSPr reference: atom-mapping + directionality, pre-baked

Atom-atom mapping (AAM) and directionality are **not dynamic features of a FabFos
run.** Every reaction in any FabFos network is a MetaNetX reaction id; its
substrates and products come from MetaNetX `reac_prop`; annotation lanes only ever
*project into* an MNXR, never mint a reaction outside the index. So `AAM(mnxr)` and
`direction(mnxr)` are **static functions of the reaction id**, and the reaction
universe is **closed and finite** — bounded by one MetaNetX release. What is
experiment-dynamic is only *which* MNXRs are present, their evidence weight, and
reinforcement copies. The chemistry never is.

This directory pre-bakes both, **once, for the whole universe**, as a frozen,
version-pinned, MNXR-keyed reference that any consumer reads as static input —
decoupled from any FabFos experiment. The goal is **100% adjudication**: every
reaction carries exactly one AAM verdict and one direction verdict, and none is
silently absent.

## Layout

| Path | What it is |
|---|---|
| `reference.py` | The pin (MetaNetX release), the paths, the schema, the verdict vocabulary and reason taxonomy, and `assert_reference` — the one place these live. |
| `build_ledger.py` | The closure ledger: enumerates every reaction, assigns each a verdict, writes the ledger + the manifest pin + the T4 worklist. |

The frozen asset itself lives in the **data tree**, not here (code and data do not
mix): `reference.REF_ROOT` → `.awm/data/ecspr_reference/mnxref-<version>/`, holding
`atom_pairs.parquet`, `direction.parquet`, `closure_ledger.parquet`, `worklist.tsv`,
and `MANIFEST.json`.

## The one rule, again

Like `canon.py`, this is a **data module**: prose may name a value defined here,
prose may not restate it. The live numbers — how many reactions are adjudicated, how
the refusals break down by class — live in `MANIFEST.json` and the ledger parquet, and
are read, never transcribed. `check_no_transcribed_numbers.py` governs this file too.

## The verdict model

Every reaction gets one verdict-state per lane, from `reference.VERDICT_STATES`:

- **resolved** — a confident map (AAM) / a ratio with evidence (direction).
- **diluted-ambiguous** — known but spread: fanout, symmetry, or member disagreement
  for AAM; ratio 1.0 / no-evidence for direction (the shrink-to-no-op verdict).
- **refused** — an explicit, reasoned refusal, carrying one **reason class**.

Plus one non-verdict the ledger exists to drive to zero: **pending** — "a layer is
still expected to act on this" (a mappable reaction not yet mapped, a carrier not yet
rescued, a direction not yet computed). **A reasoned refusal is a verdict; only
`pending` is un-adjudicated.** Closure is `pending == 0`.

### The reason taxonomy is the engine's

The refusal reason classes are `ecspr_aam_rescue`'s — electron carrier, acyl carrier,
non-molecule, SEED lump, unreconciled stub — *imported*, plus the structural boundary
classes a universe ledger needs (pseudo-reaction, generic R-group, unmapped-computable,
no-transfer). The ledger classifies a blocking participant with the **same**
`placeholder_for` / `REFUSE` the rescue path uses, so the ledger and the mapper cannot
drift: a reaction the ledger calls `electron_carrier` is exactly one the mapper would
stand a placeholder in for.

A reason class is either **terminal** (a reasoned refusal that closes the reaction —
pseudo-reaction, non-molecule, SEED lump, no-transfer) or **still-resolvable / pending**
(a class a T2 automatic pass or a T4 class rule is expected to convert). See
`reference.REASON_TERMINAL` / `AAM_PENDING_REASONS`.

## Running it

```bash
mamba run -n scadc-metabolic-model python build_ledger.py \
    --worklist $REF_ROOT/worklist.tsv
```

Reads whatever AAM + direction tables it is pointed at — the **frozen reference tables
once they exist, the incumbent single-source tables until then** — and measures the gap.
Run it before T2/T3/T4 to get the worklist; run it with `--strict` after, as the closure
gate (exit non-zero unless 0 un-adjudicated). Same script, same universe, both times.

- `--worklist` writes, per pending reason class, the **frequency-ranked structureless
  blocker names, keyed by MNXM id** — the T4 authoring aid. One class rule on the
  most-blocking carrier reclaims many reactions at once (this is what keeps
  full-universe manual adjudication finite).

## Reachability

The target is full coverage, but effort is prioritised: the ledger marks each reaction
`in_evidence` (any host/fosmid ORF projects onto it — the reachable set) and
`in_base_graph` (the strict in-network subset). Pseudo-reactions (transport / exchange /
degenerate) are classified **out** explicitly — non-chemistry, a terminal refusal — never
left ambiguous.

## Freezing (T5)

`build_ledger.py` refreshes `MANIFEST.json` every run: the MetaNetX release, the
`reac_prop` content hash (a moved universe fails `assert_reference` loudly), the closure
state, and whether the frozen tables exist. `assert_reference(require_closed=True)` — the
consumer's guard, mirroring `canon.assert_canonical_*` — passes only when the tables are
present, the universe hash still matches, and `n_pending == 0`.
