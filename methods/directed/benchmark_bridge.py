"""Benchmark bridge -- score the DIRECTED prototype on the X/Y validation benchmark,
head-to-head against the undirected incumbent, on identical cells.

The undirected arm already exists: the incumbent validation panel (``panel_{lof,gof}.parquet``)
holds ``G_base``/``G_pert`` per (facet, condition, anchor, product) cell, scored as the
benchmark baseline (directional-specificity AUC ~0.83 carbon). This bridge recomputes only
the **directed** conductance on the *same* cells, so any AUC difference is attributable to
direction alone -- not to a reimplementation of the panel, the perturbations, or the join.

For each incumbent cell it:
  * reconstructs the perturbation from ``condition_id`` + ``unit`` -- LOF: the gene's dead
    reactions via ``_val.ko_dead_mnxrs`` -> ``base_minus_reactions``; GOF: the added mnxr(s)
    from ``unit`` -> ``base_plus_reactions`` -- reusing the incumbent's own perturbation
    code, never a copy;
  * solves the two-terminal conductance with the engine's rectified solver, forward
    weights ``gp = w`` and backward ``gm = ratio * w`` from the directionality ensemble
    (``DirectedNet``), on the perturbed graph and on the base graph;
  * emits ``effect = log(C_pert / C_base)`` -- the same quantity the incumbent emits as
    ``s = log(G_pert/G_base)`` -- through the incumbent baseline adapter's exact formatting
    (``FACET_MAP``, ``canon_id``, the GOF host-split), so the observation join into ``Y`` is
    bit-identical to the baseline's.

``--validate`` first proves the pipe: on a sample of conditions it solves with the ratios
forced to 1.0 (the symmetric limit) and checks the result reproduces the incumbent's ``s``
per cell. Only once that holds is the directed arm (``--run``) trustworthy.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from directed_network import (DirectedNet, load_roles, load_ratios,          # noqa: E402
                              sha256_of, DIRECTION_ANNOTATION)
sys.path.insert(0, str(HERE.parent))
from fabfos import canon                                                                 # noqa: E402

# Incumbent validation lane (perturbation maps + the panel parquets) and the benchmark
# scorer/adapter -- both read-only.
VAL_LANE = Path("/home/tony/agentic_workspace/projects/scadc/metabolic-modelling/main/"
                "metabolic-modelling/04_reaction_network/validation")
FIGMODEL_VAL = Path("/home/tony/agentic_workspace/projects/scadc/fig-model/main/ecspr/validation")
PANEL_OUT = FIGMODEL_VAL / "out"
BENCH_DIR = Path("/home/tony/agentic_workspace/projects/scadc/fig-model/main/ecspr/benchmark")
Y_DIR = canon.DATA / "ecspr" / "benchmark" / "v1" / "Y"
sys.path.insert(0, str(VAL_LANE))
sys.path.insert(0, str(FIGMODEL_VAL))
sys.path.insert(0, str(BENCH_DIR))
import _val                                                                  # noqa: E402
import _panelspec as ps                                                      # noqa: E402
from _anchors import canon_id                                               # noqa: E402

# Incumbent baseline adapter constants (mirrored from 40_baseline.py so the observation
# join is identical). facet id "A_iECDH10B" -> benchmark facet "netA_iECDH10B".
FACET_MAP = {"A_iML1515": "netA_iML1515", "A_iECDH10B": "netA_iECDH10B",
             "B_iML1515": "netB_iML1515", "B_iECDH10B": "netB_iECDH10B"}
GEM_SUFFIX = {"A_iML1515": "_ML", "A_iECDH10B": "_DH",
              "B_iML1515": "_ML", "B_iECDH10B": "_DH"}

ANCHOR_MNXM = {k: m for k, m, *_ in ps.ANCHORS}     # anchor key -> mnxm


def log_ratio(cp, cb):
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.log(cp / cb)
    return s if np.isfinite(s) else np.nan


def perturbed_graph(base, arm, condition_id, unit):
    """Reconstruct the incumbent condition's perturbed graph from its id + unit."""
    if arm == "lof":
        _bigg, mnxrs = _val.ko_dead_mnxrs(unit)          # unit == gene
        return _val.base_minus_reactions(base, mnxrs), len(mnxrs)
    # GOF: unit lists the added mnxr(s); condition_id gof_ext_<mnxr> is single-reaction
    mnxrs = [m for m in str(unit).split(";") if m.startswith("MNXR")]
    addmap = {m: 1.0 for m in mnxrs}
    return _val.base_plus_reactions(base, "C", addmap), len(mnxrs)


def _ceff(net, s_node, t_node, tol=None):
    """Directed conductance; 0.0 if either terminal left the perturbed LCC."""
    if s_node not in net.idx or t_node not in net.idx:
        return 0.0
    return net.ceff(s_node, t_node, tol=tol)


def _subsample(grp, offtarget, seed_base):
    """Per condition, keep ALL target cells (is_diagonal) + a SEEDED random sample of
    `offtarget` off-target cells. Returns the kept sub-DataFrame. This bounds wall-clock
    without silently shrinking coverage -- the specificity AUC is target-vs-off-target
    within a condition, so all targets + a fixed off-target sample keeps the metric valid
    with a stated (smaller) negative set. Seeded per condition for reproducibility."""
    if offtarget is None:
        return grp
    tgt = grp[grp.is_diagonal]
    off = grp[~grp.is_diagonal]
    if len(off) > offtarget:
        # deterministic seed from the condition id so the sample is reproducible
        seed = (seed_base + abs(hash(str(grp.condition_id.iloc[0]))) % 100000) & 0x7fffffff
        off = off.sample(n=offtarget, random_state=seed)
    return pd.concat([tgt, off])


def solve_panel(facet, roles, ratios, force_noop, limit_conditions=None,
                arms=("lof", "gof"), progress_every=10, tol=None, offtarget=None, seed=0):
    """Return the incumbent panel rows for `facet` with an added `effect` column
    (directed, or symmetric-limit if force_noop). Base conductance per (anchor,product)
    is solved once and reused as the denominator across all conditions.

    `tol` loosens the Newton tolerance (safe for the rank-only metric); `offtarget` caps
    the off-target cells per condition (all targets always kept) to bound wall-clock."""
    _val.set_network(facet)
    base = _val.load_base("C")
    Dbase = DirectedNet.from_graph(base, "C", roles, ratios, force_noop=force_noop)
    inc_facet = facet  # parquet uses the short facet id

    frames = []
    kept_cells = total_cells = 0
    for arm in arms:
        df = pd.read_parquet(PANEL_OUT / f"panel_{arm}.parquet")
        df = df[(df.facet == inc_facet) & df.usable & ~df.anchor_incident].copy()
        if df.empty:
            continue
        # base conductance cache over the (anchor, product) pairs this arm needs
        cbase = {}
        conds = list(df.groupby("condition_id"))
        if limit_conditions:
            conds = conds[:limit_conditions]
        rows_out = []
        for ci, (cid, grp) in enumerate(conds):
            total_cells += len(grp)
            grp = _subsample(grp, offtarget, seed)
            kept_cells += len(grp)
            unit = grp.unit.iloc[0]
            try:
                Gp, npert = perturbed_graph(base, arm, cid, unit)
                Dp = DirectedNet.from_graph(Gp, "C", roles, ratios, force_noop=force_noop)
            except Exception as e:                       # a broken perturbation -> skip cond
                print(f"  [{arm}] {cid}: perturbation failed ({e}); skipping", flush=True)
                continue
            for r in grp.itertuples():
                s_node = ("met", ANCHOR_MNXM.get(r.anchor))
                t_node = ("met", r.product)
                if s_node not in cbase:
                    cbase[s_node] = {}
                if t_node not in cbase[s_node]:
                    cbase[s_node][t_node] = _ceff(Dbase, s_node, t_node, tol=tol)
                cb = cbase[s_node][t_node]
                cp = _ceff(Dp, s_node, t_node, tol=tol)
                eff = 0.0 if (cb == cp) else log_ratio(cp, cb)
                rows_out.append((r.Index, eff, cb, cp))
            if progress_every and (ci + 1) % progress_every == 0:
                print(f"  [{arm}] {ci + 1}/{len(conds)} conditions ({kept_cells} cells)",
                      flush=True)
        if rows_out:
            idx, eff, cb, cp = zip(*rows_out)
            sub = df.loc[list(idx)].copy()
            sub["effect"] = eff
            sub["c_base"] = cb
            sub["c_pert"] = cp
            frames.append(sub)
    if offtarget is not None:
        print(f"[subsample] off-target cap={offtarget}, seed={seed}: kept {kept_cells:,} of "
              f"{total_cells:,} cells ({kept_cells / max(total_cells, 1):.1%}); all target "
              f"cells retained", flush=True)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def adapt(df, facet):
    """Incumbent baseline adapter (mirrored): panel rows -> observations schema, with
    `effect` taken from our directed column, not from `s`."""
    known = set(pd.read_csv(Y_DIR / "conditions.tsv", sep="\t", dtype=str).condition_id)
    edge_id = "C__" + df.anchor + "__" + df["product"].map(canon_id)

    def cid(row):
        c = row.condition_id
        if row.arm == "gof" and c.startswith("gof_ext_"):
            cand = c + GEM_SUFFIX[facet]
            return cand if cand in known else c
        return c

    return pd.DataFrame(dict(
        facet=df.facet.map(FACET_MAP),
        condition_id=df.apply(cid, axis=1),
        edge_id=edge_id,
        effect=df.effect.astype(float),
    ))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--facet", default="A_iECDH10B", choices=list(FACET_MAP))
    ap.add_argument("--out", type=Path, default=HERE / "bench_out")
    ap.add_argument("--validate", type=int, default=0,
                    help="N: symmetric-limit reproduction check on N conditions/arm")
    ap.add_argument("--run", action="store_true", help="full directed arm -> observations.tsv")
    ap.add_argument("--tol", type=float, default=1e-6,
                    help="Newton tolerance (rank-only metric; 1e-6 is safe and faster)")
    ap.add_argument("--offtarget", type=int, default=None,
                    help="cap off-target cells per condition (all targets kept); bounds wall-clock")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    roles, ratios = load_roles(), load_ratios()
    print(f"facet={args.facet}  direction sha256={sha256_of(DIRECTION_ANNOTATION)[:16]}…", flush=True)

    if args.validate:
        # Symmetric limit (ratio==1) must reproduce the incumbent's undirected s per cell.
        df = solve_panel(args.facet, roles, ratios, force_noop=True,
                         limit_conditions=args.validate)
        m = df.dropna(subset=["s", "effect"])
        m = m[np.isfinite(m.s) & np.isfinite(m.effect)]
        err = (m.effect - m.s).abs()
        print(f"\n[validate] symmetric-limit reproduction on {df.condition_id.nunique()} "
              f"conditions, {len(m)} finite cells:")
        print(f"  |effect_ratio1 - incumbent_s|: max={err.max():.2e} "
              f"mean={err.mean():.2e} median={err.median():.2e}")
        print(f"  cells with |d|>1e-6: {(err > 1e-6).sum()} / {len(m)}")
        return 0

    if args.run:
        df = solve_panel(args.facet, roles, ratios, force_noop=False,
                         tol=args.tol, offtarget=args.offtarget, seed=args.seed)
        obs = adapt(df, args.facet)
        out = args.out / f"observations_directed_{args.facet}.tsv"
        obs.to_csv(out, sep="\t", index=False)
        # directed fraction actually exercised
        dfrac = float((df.c_base != df.c_pert).mean())
        print(f"\n[run] wrote {len(obs):,} directed observations -> {out}")
        print(f"  conditions={df.condition_id.nunique()} cells_moved={(df.c_base != df.c_pert).sum():,} "
              f"({dfrac:.1%})")
        print(f"  score with:\n    python {BENCH_DIR}/30_score.py --obs {out} --mode signed "
              f"--out {args.out}/scored_directed_{args.facet}")
        return 0

    ap.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
