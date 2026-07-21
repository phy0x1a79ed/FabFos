# ECSPr execution-surface map

Where the ECSPr method actually lives when you run it from a figure scope. This is
the map that kept getting re-derived by throwaway explore agents; it is now durable.
It names **symbols and file:line**, never basis values — for any number, import it
from `canon.py` (the `check_no_transcribed_numbers.py` rule governs this file too).

Verified by grep at every cited line on 2026-07-20. Line numbers drift; the symbol
names and the import topology are the load-bearing part.

## The one thing that surprises everyone: cross-scope imports

A figure scope (`fig-diagrams`, `fig-model`, …) has **no local `_val.py` and no local
`_gof_specs.py`**. The validation library is imported *across scopes* by hard-coded
`sys.path` inserts:

- `main/ecspr/validation/_directed_panel.py:40-53` inserts
  `metabolic-modelling/main/metabolic-modelling/04_reaction_network` and its
  `validation/` onto `sys.path`, then (`:50`) appends `canon.ENGINE_LIB/resources/lib`.
- `main/ecspr/validation/_panelspec.py:68-73` does the same.

So `_val`, `_gof_specs`, `_smw`, `_graphs`, `_reff` **resolve to the metabolic-modelling
copies**, while `_panelspec`, `_directed_panel`, `canon` are the local figure-scope
copies. Any code you write in a figure scope against `coerce_target` / `gene_lookup` /
`gene_add_map` binds to the metabolic-modelling `_val.py` at import time. Edit T0/T2
code there, not here.

Absolute roots (this machine):
- **MM** = `/home/tony/agentic_workspace/projects/scadc/metabolic-modelling/main/metabolic-modelling/04_reaction_network`
- **FD** = `/home/tony/agentic_workspace/projects/scadc/fig-diagrams/main`
- **ENGINE** = `canon.ENGINE_LIB` = `/home/tony/agentic_workspace/projects/metasmith-libraries/fabfos` (branch `feat/fabfos`); diode kernel under `ENGINE/resources/lib/`.

## A. Coercion path (near-product / stereo-sibling aliasing)

| Symbol | Location | Note |
|---|---|---|
| `gof_axis_endpoints()` | `MM/validation/73_gen_coercion_alias.py:37` | endpoint source = every `(s,d)` in `GOF_SPECS` (imported `:28`) |
| block-1 sibling search | `73_*:81-98` | `b1_to_conn` = {InChIKey-block1 → [(mnxm, base.degree)]}; picks max-degree connected sibling, skips ids already in LCC or lacking an InChIKey |
| writes alias TSV | `73_*:102` | `MM/validation/ground_truth/gof_coercion_alias.tsv` |
| alias TSV schema | `73_*:17-18` | `network, from_mnxm, to_mnxm, from_name, to_name, block1, to_deg, n_designs` |
| `NETWORKS` list | `73_*:33` | `A_iML1515, A_iECDH10B, B_iML1515, B_iECDH10B` — **no `dh10b`** |
| `coerce_target(mnxm, net=None)` | `MM/validation/_val.py:319` | returns `_load_coerce().get((net or _NET, mnxm), mnxm)` |
| `_load_coerce()` | `_val.py:308` | reads alias TSV → `{(network, from_mnxm): to_mnxm}` |
| `panel_products(specs=None)` | `FD/ecspr/validation/_panelspec.py:233` | plain union of GOF+LOF products |

**The ban you will hit at T0.** Coercion is documented as GOF-target-only and
**forbidden on the canonical axis set / baseline / LOF endpoints** — stated twice:
`_val.py:322` (coerce_target docstring) and `_panelspec.py:237-256` (with a measured
"all four intersections empty" proof). That proof is *circular for our purpose*: the
alias table is built only from `gof_axis_endpoints()`, so canonical-axis endpoints were
never offered to the aliaser — the intersection is empty by construction, not because
no canonical endpoint is coercible. `_panelspec.py:83-111` already documents the live
defect this blocks: the alpha-G6P id `MNXM1105855` is ABSENT from both net-A carbon
LCCs and should resolve to `MNXM1364111` (same InChIKey block-1 `NBSCHQHZLSJFNQ`).

## B. Directed solver + facet (the scoring kernel)

| Symbol | Location | Returns / note |
|---|---|---|
| `class DirectedFacet` | `FD/ecspr/validation/_directed_panel.py:71` | `__init__(facet, X, roles, ratios, *, force_noop=False)` |
| `score_add_many(addmaps)` | `_directed_panel.py:155` | `{cid: {G_base, G_pert, p, n_added, usable=(p>=2), touched}}`; one `directed_solve_element(warm=True)` sweep per axis |
| `score_add(addmap)` | `_directed_panel.py:178` | wraps `score_add_many({"_":addmap})` |
| `score_delete(mnxrs)` | `_directed_panel.py:183` | LOF |
| `score_perturbation(addmap, delete)` | `_directed_panel.py:205` | combined; **forfeits factor reuse** (one full solve per condition) |
| diode kernel | `ENGINE/resources/lib/ecspr_directed.py` | `orient_and_weight`, `build_ar2m_directed`, `directed_solve_element`, `_softmax0`, `DIODE_SMOOTH_DELTA` |
| network helpers | `ENGINE/resources/lib/ecspr_network.py` | `load_reaction_roles`, `load_direction_ratios`, `_base_lcc` |

**Warm-start granularity:** `score_add_many` reuses the base solve **per axis**, not one
global factorization — the softplus diode has no Woodbury collapse, so each cell is an
independent Newton solve. Budget cluster time as (axes × conditions) full solves, not a
low-rank update. `score_perturbation` reuses nothing.

## C. Canon + axes

| Symbol | Location | Note |
|---|---|---|
| `AXES_JSON` | `FD/fabfos/canon.py:225` | `…/04_reaction_network/cache/biomass_dag_axes_<AXIS_SET>.json`; keyed by axis-name string → `{source:[mnxm], sink:[mnxm], category, …}` |
| `AXES_TSV` | `canon.py:221` | `…/figures/publish/03_model/ecspr_scadc/biomass_edges_<AXIS_SET>.tsv` — **path does not exist**; read axes from `AXES_JSON` |
| `AXES_TESTABLE_JSON` | `canon.py:232` | graph-dependent strict subset (one carbon phospholipid axis drops on the honest reference graph) |
| `assert_canonical_axes(table)` | `canon.py:580` | asserts count == `AXES_N` and per-element == `AXES_PER_ELEMENT`; refuses retired sets by filename |
| `AXES_N`, `AXES_PER_ELEMENT` | `canon.py:233-234` | import; do not transcribe |
| `REFERENCE_REAC_PROP_SHA256` | `canon.py:294` | reac_prop gate hash |
| `ENGINE_LIB` | `canon.py:69` | engine root |

**Axis-name format:** `{category}__{SRC_MNXM}__{DST_MNXM}__{ELEM}` (e.g.
`catabolic__MNXM1364061__MNXM46__C`). SRC/DST are raw MNXM ids; element is the trailing
token. Adding a key to `AXES_JSON` breaks `assert_canonical_axes` — declare any new axis
(e.g. a glycogen validation axis) in a **sibling** registry, never here.

**Curated transplant table** (distinct from axis definitions): `gof_*.tsv` under
`MM/validation/ground_truth/` (e.g. `gof_forsberg.tsv`). Schema:
`dataset  role  gene  ec  mnxr  in_base  src_mnxm  src_name  sink_mnxm  sink_name  element  expected_dir  citation  note`. Read by `_panelspec.curated_addmap` (`:153`, filtering `role=="add"`).

## C-bis. Gene → reaction (GEM/GPR)

All in `MM/validation/_val.py`:

| Symbol | Location | Note |
|---|---|---|
| `_MODEL_ORACLE` | `_val.py:64-67` | `{iML1515:{gem,crosswalk}, iECDH10B:{gem,crosswalk}}` |
| `NETWORKS` registry | `_val.py:68-74` | keys `dh10b, A_iML1515, A_iECDH10B, B_iML1515, B_iECDH10B` |
| `_NET` default | `_val.py:75` | `dh10b` — **always `set_network` explicitly** before scoring |
| `set_network(net)` | `_val.py:78` | resets cobra/crosswalk/evidence singletons |
| `gene_lookup()` | `_val.py:375` | `{locus or mnemonic: cobra.Gene}` — keyed both ways (the join for name-only screens) |
| `ko_dead_mnxrs(gene_key)` | `_val.py:388` | correct AND/OR GPR via `r.gpr.eval([g.id])` |
| `gene_add_map(gene_key, weight=1.0)` | `_val.py:411` | `{mnxr: weight}` add-map (reverses a KO) |

iECDH10B GEM = `data/scadc/metabolic_modelling/networks/iECDH10B_1368.json`; iML1515 GEM =
`fig-model/main/metabolic_model/models/iML1515.json`; crosswalks
`gem_rxn_to_mnxr_{iML1515,iECDH10B}.parquet`.

## D. Cluster drivers — two families, do not confuse them

Both use `DirectedFacet` + `score_add_many`. Pick by cluster; recipes differ.

- **fir (Alliance)** — `FD/fabfos/directed/fir/`: `fir_null.py`, `fir_observed.py`,
  `run_null.sbatch`, `run_obs.sbatch`, `manifest.txt`. sbatch: account `def-shallam_cpu`
  (batch alt `rrg-shallam-ab_cpu`), `StdEnv/2023 python/3.12` + a **venv** (np2.4/sp1.17/
  pd3.0/nx3.6), **no container**, `splu` path (no CHOLMOD). CLI:
  `--stage-dir --element --draw-size --workers --out-reff --out-ieff --resume --force-noop`.
- **sockeye (UBC ARC)** — staging `/scratch/st-shallam-1/txyliu/scadc_dir_null`, account
  `st-shallam-1`, partition `skylake` 32c/node, container `esmc.sif`, `gcc/9.4.0` +
  `apptainer/1.3.1` as **two separate module loads**, BLAS pinned 1 thread/worker.
  Recorded run: one job per draw-size, balanced shards by element. Recipe of record:
  `data/scadc/plans/09-directed-canonical-fir-run.md` (sockeye RESTART-2 section).

**GOF/LOF panel cluster twin** (uses the add-map path, the T4 template):
`FD/ecspr/validation/fir_panel.py` (`DirectedFacet` copy, `score_add_many`), local
counterpart `20d_run_directed_panel.py`.

Standalone drivers are **canon-free and parameterised by stage-dir** on purpose:
`canon.py` resolves local paths absent on the cluster. Clear `__pycache__` on scratch
before launch (stale bytecode has bitten this lane). Do **not** reuse `reff/sockeye/` —
that is the retired undirected Woodbury/GPU path, never fired.

## E. Data conventions

- References root: `data/scadc/references/` (= `.awm/data/references/`). Datasets:
  `laser/`, `keio/`, `eydallin/`, `metanetx/`, `bigg/`, … Shape: immutable source docs
  under `raw/` (or `inputs/`) beside a parsed top-level TSV. LASER is the richest exemplar;
  `laser/inputs/EcAccession.txt` (`Gene Synonym Location Strand Length PID Code COG Product`)
  is the name→b-number fallback.
- Directed outputs of record live under `data/scadc/ecspr_reference/mnxref-4_5/{solve_directed,null_directed}/`
  (reference graph) — netA/netB caches under `MM/cache/<facet>/base_<ELEM>.pkl`.
- The canonical-run recipes and outcomes: `data/scadc/plans/08-*`, `09-*` (directed).

## graphify status

`graphify build` currently refuses this corpus: every relevant subtree carries at least
one doc/paper file and no LLM API key is set for semantic extraction (a code-only corpus
would need none). Until a key is available or the docs are excluded, verify against source
directly (grep at the file:line above) rather than the graph. Re-attempt with
`graphify build target=<code-only-subtree>` once a key is present.
