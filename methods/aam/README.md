# Ensemble atom-atom mapping (AAM)

Makes the atom lane's correspondence universe an **ensemble** instead of a single mapper.
Today `atom_pairs_universe.parquet` is RXNMapper-only with no record of who mapped what;
this directory fuses multiple mappers into one correspondence table with provenance,
mirroring the directionality ensemble (`../direction/`) member-for-member.

Scope now spans the two neural members **plus the independent curated member (MetaCyc)** —
the combiner, the provenance schema, three wired members, and the pairwise agreement between
every pair measured. The full universe-scale LocalMapper run and promotion of the fused table
into production are named below as the remaining steps, not done here.

## The members (method families, not databases)

| # | member | family | routing | state |
|---|--------|--------|---------|-------|
| 1 | RXNMapper | neural transformer (attention) | mapped reaction SMILES | full universe cache (57,522 rxns) |
| 2 | LocalMapper | neural transformer (template-augmented) | mapped reaction SMILES | cache covers 403 rxns; universe run pending |
| 3 | MetaCyc curated | database atom maps (`atom-mappings-smiles.dat`) | `metacyc.reaction:`→MNXR crosswalk | **wired** (`metacyc_member.py`, 16,526 MNXR) |

Members 1 and 2 are both transformers on reaction SMILES → **correlated**. Member 3 is a
database, not a model → **independent**. That split is the whole design, exactly as
eQuilibrator+dGbyG (correlated) vs BioCyc (independent) is for direction. MetaCyc is
deliberately **not** in `NEURAL_MEMBERS`, so a MetaCyc∪neural consensus fuses undiscounted
(full weight) and MetaCyc breaks ties in disagreements — an independent vote moves the
ensemble where a correlated one cannot.

### The curated member (`metacyc_member.py`)

MetaCyc ships `atom-mappings-smiles.dat`: mapped reaction SMILES in the **exact** `[C:n]…>>…`
form the engine's `pairs_from_mapped` already consumes, with R-groups embedded in the SMILES.
So there is no index decoder to re-derive and no `classes.dat` carrier-formula blocker — the
prior round's "absent decoder" worry was about the *other* route (`atom-mapping.dat` INDICES
against compound formulae), which this supersedes. The loader joins each MetaCyc reaction id to
MNXR through the `metacyc.reaction:` prefix of `reac_xref` (99.6% of lines join; the ~292
stereo/direction frame variants that collapse to a shared MNXR are deduped deterministically,
keep-first), and returns the same `{mnxr → (mapped_smiles, confidence)}` shape a neural cache
does. A loud floor (`JOIN_FLOOR`) fires if the crosswalk breaks. Curated confidence is a named
constant (`CURATED_CONFIDENCE`): a curated map is asserted, not scored.

## Common identity space

Every member's correspondence is expressed in the engine's canonical-rank identity:
`(metabolite, CanonicalRankAtoms(breakTies=True))` on the sanitized per-metabolite molecule
**is** an atom, invariant to how a reaction happened to write it. Because the identity is a
property of the metabolite and not the mapper, two mappers that map the same physical atom
name the same node. Each member's correspondence is obtained by running the engine's own
`pairs_from_mapped` on that member's mapped SMILES — one extractor, so the only thing that
varies between members is the mapping, never the identity bookkeeping. (This is why the older
`20_build_aam_cache.py`'s `MolFromSmiles`-position scheme is **not** reused: it is a different,
per-reaction-unstable identity; the canonical-rank scheme is the audited one.)

## The combiner (`combine_aam.py`)

A per-atom vote in correspondence space, the discrete analogue of direction's ΔG' vote:

1. **Consensus.** Members that map a substrate atom to the *same* product atom corroborate
   it. Neural-neural agreement is **discounted** by `NEURAL_SHARED_FLOOR` — two transformers
   agreeing is partly shared architectural bias, not independent confirmation, so they cannot
   count as two independent votes (the TAU_SHARED analogue). A neural∪curated consensus is
   undiscounted.
2. **Single member.** An atom only one member maps is kept at `SINGLE_MEMBER_CREDIT` —
   present, not corroborated. Coverage the ensemble gains over any one member.
3. **Disagreement → dilution.** Where members map an atom to *different* products, the
   correspondence is **not** decided by majority or confidence — it is spread across the
   disputed products by member weight, so a contested atom dilutes its transfer rather than
   committing. This is the fanout dilution the atom graph already applies to name-ambiguous
   pairs (`ambiguous_diluted`), extended across mappers, and it is the AAM analogue of
   direction's shrink-to-no-op: no member maps an atom → no pair → ratio-1.0-equivalent.
4. **Provenance, not selection.** Every emitted pair carries `method` / `source` /
   `confidence`. Like `dir_tier`, this is a record of which members spoke and how they
   agreed — **not** a gate on usability. Every pair is used; the graph reads all of them.

Output schema grafts the provenance columns onto the universe's pair schema:
`mnxr, element, substrate, product, sub_idx, prod_idx, pair_w, method, source, confidence`
(`sub_idx`/`prod_idx` are canonical ranks; `method` ∈ {`consensus`, `<member>_only`,
`disagree_diluted`}). This mirrors the tier/method/confidence schema the prior
`20_build_aam_cache.py` already proved out — reused as a schema, not as code.

## The measurements this round exists to make

**Neural-neural agreement.** On the 403 reactions both neural members cover, fused into 14,224
provenance-carrying atom-pairs: **the two neural mappers agree on 93.9% of the atoms both
address** (reproduce with `--no-metacyc --overlap-only`). A second *neural* member confirms the
first ~94% of the time — useful for robustness and a little coverage, but strongly correlated,
exactly as `NEURAL_SHARED_FLOOR` assumes. That is the empirical basis for wiring an
*independent* third member rather than a second neural one.

**Curated-vs-neural agreement (does MetaCyc argue?).** Over the full universe MetaCyc and
RXNMapper co-map ~24.7k C/N/S/P atoms; **they agree on 92.4% at unique-atom identity.** That is
just below the ~95% sanity bar, so the shortfall was diagnosed rather than reported bare:
collapsing **symmetry-equivalent** product atoms (`CanonicalRankAtoms(breakTies=False)`)
lifts agreement to **95.0%** — 41% of the disagreement is the two mappers picking tie-broken
*different but physically equivalent* atoms, a canonicalization artifact, not chemistry. So the
chemically-meaningful curated-vs-neural agreement sits **at the bar**; the residual ~5% is
genuine, and since MetaCyc is curated it is where the neural mapper is the more likely one to
be wrong. The ensemble **records** that as provenance (`disagree_diluted`), never resolves it —
a curated disagreement is a finding, exactly as the design intends.

**Coverage and its honest loss.** MetaCyc's member is 16,526 MNXR, but only a fraction fuse
pairs: the re-extraction status splits ~`stripped`/`unparseable`/`no_pairs`/`ok`. The dominant
losses are *safe* and *measured, not hidden*:
- `unparseable` (~28%) — MetaCyc curates **generic reaction classes** whose SMILES carry
  R-group tokens (`[R:n]`); RDKit's SMILES parser rejects `R` as an atom. Real property of the
  source, not a wiring bug; no parse → no pairs → never a wrong pair.
- `stripped` (~56%) — MetaCyc and MetaNetX write a reaction's participant set differently
  (protons/water/generic frames), so the template count ≠ the MNXR equation's. Cross-namespace,
  safe (`pairs_from_mapped` refuses rather than mis-pairs), lowers coverage only.

**Base-graph reach.** Of the 4,675 reactions in the ECSPr base graphs (C/N/S/P union), MetaCyc
covers **1,702 (36.4%)** by MNXR — the number that decides downstream usefulness: an
independent vote available on roughly a third of the base graph.

## What is deferred, and why

- **Universe-scale LocalMapper.** The on-disk LocalMapper cache covers 403 reactions (39 at
  τ≥0.9), not the 57,522-reaction universe. Both mappers are installed (`scadc-aam` env:
  rxnmapper 0.4.3 + localmapper 0.1.5), so this is a compute run over
  `aam_rxnmapper_universe.tsv`'s SMILES, not a blocker — just not this round's demonstration.
- **Promotion to the universe.** This combiner writes a demonstration/validation table. Making
  the fused, provenance-carrying table the *production* `atom_pairs_universe.parquet` is a
  separate step: the engine's `ecspr_atom_pairs.py` would consume per-member mapped-SMILES
  caches and emit the ensemble columns, and `ecspr_atom_graph.py` would read
  `confidence`/`source`. The tunables would then graduate into `canon.py`. Named, not done —
  this round wires MetaCyc and measures agreement, it does not promote or re-freeze the referent.
- **R-group coverage.** ~28% of MetaCyc's maps are generic (R-group) reaction classes RDKit
  cannot parse as SMILES. They are the curated database's reach into reaction *classes* rather
  than concrete reactions; recovering them would need a wildcard-substitution path that the
  concrete-participant `match_mols` would then mostly drop anyway. Left as documented coverage.

## Tunables

`NEURAL_SHARED_FLOOR`, `SINGLE_MEMBER_CREDIT`, `NEURAL_MEMBERS` are committed in
`combine_aam.py` before the run (the same discipline as the pulse-chase `thresholds.py`).
They are provisional: when the ensemble becomes canonical they graduate into `canon.py` the
way `DIR_TAU_SHARED` did, and stop being restated anywhere else.

## Run

```
# three-member ensemble over the full universe (MetaCyc's reach visible)
mamba run -n scadc-metabolic-model python combine_aam.py

# materialise the curated member as an inspectable cached TSV
mamba run -n scadc-metabolic-model python metacyc_member.py

# acceptance self-test (crosswalk floor, collision floor, curated-vs-neural spot-check)
mamba run -n scadc-metabolic-model python selftest_metacyc.py

# reproduce the original two-member neural demonstration (93.9%)
mamba run -n scadc-metabolic-model python combine_aam.py --no-metacyc --overlap-only
```

Pure re-extraction over cached / curated mapped SMILES — runs anywhere rdkit + the engine
import; no mapper is invoked. `--overlap-only` restricts to reactions both *neural* members
cover; `--no-metacyc` drops the curated member.
