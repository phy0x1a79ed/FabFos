"""LASER gain-of-function arm, measured under three ground constructions, on both lanes.

The question is not "is there an effect" -- it is whether the media->biomass probe produces
numbers that are (a) SENSIBLE against each condition's stated expectation and (b)
DISTINGUISHABLE between conditions, and whether the choice of ground changes either.

One row per (lane, mode, element, condition, precursor). Deltas are against that
lane/mode/element's OWN baseline solve, so a mode is never compared to another mode's zero
point.

LANES
    ``netB`` -- evidence-weighted (epi300, E_full), built here from the atom-pair table.
    ``netA`` -- curated GEM (iECDH10B) at uniform E=1.0, PRE-BUILT and loaded from .npz
      because `ecspr_build._load_model` needs cobra, which the pinned ecspr container does
      not carry. The graph is therefore identical to what the container would have built.

Both lanes take the perturbation through `ground.merge_reactions`, verified equal to a full
rebuild at machine precision -- applying it by two different mechanisms would confound a
netA-vs-netB difference with a difference in how it was perturbed.

Incremental and resumable: every condition appends and flushes, and a restart skips
(lane, mode, element, condition) keys already present. Assume this dies at 60%.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

import ecspr_build as B                                             # noqa: E402
import ecspr_probe as P                                             # noqa: E402
import ground as G                                                  # noqa: E402
from ecspr_graph import AtomGraph                                   # noqa: E402
from conditions import gof_conditions                               # noqa: E402

COLS = ["lane", "mode", "element", "media", "condition_id", "n_rxns", "rxns", "metabolite",
        "draw_base", "draw_cond", "delta_draw", "rel_draw",
        "total_base", "total_cond", "delta_total", "rel_total",
        "prec_share_base", "prec_share_cond", "leak_frac_base", "leak_frac_cond",
        "is_novel", "n_pairs_element", "n_edges_added",
        "n_edges_touching_source_component", "reachable",
        "converged", "port", "leak", "add_weight", "expected"]


def _done(path: Path) -> set:
    if not path.exists():
        return set()
    out = set()
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out.add((r["lane"], r["mode"], r["element"], r["condition_id"]))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--netA-dir", default=None, help="pre-built netA .npz directory")
    ap.add_argument("--out", required=True)
    ap.add_argument("--lanes", nargs="*", default=["netB", "netA"])
    ap.add_argument("--elements", nargs="*", default=["C", "N", "S", "P"])
    ap.add_argument("--modes", nargs="*", default=list(G.MODES))
    ap.add_argument("--port-mult", type=float, default=1e3)
    ap.add_argument("--leak-mult", type=float, default=1e-7)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args(argv)

    ref = Path(a.ref)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    done, new = _done(out), not out.exists()

    axes = json.loads((ref / "derived/axes/biomass_dag_axes_set4.json").read_text())
    ratios = B.load_direction_ratios(ref / "derived/mnxref-4_5/direction.parquet")
    weights = B.load_evidence_weights(ref / "derived/evidence/evidence_weights.parquet")
    w_median = float(np.median(list(weights.values())))

    conds = [c for c in gof_conditions()[0] if c["n_rxns"]]
    if a.limit:
        conds = conds[: a.limit]
    mine = [c for i, c in enumerate(conds) if i % a.nshards == a.shard]
    print(f"[shard {a.shard}/{a.nshards}] {len(mine)}/{len(conds)} conditions "
          f"lanes={a.lanes} elements={a.elements} modes={a.modes}", flush=True)

    fh = open(out, "a", buffering=1, newline="")
    w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
    if new:
        w.writeheader()

    for X in a.elements:
        sets = P.axis_node_sets(axes, X)
        if not sets["precursors"]:
            print(f"  [{X}] no precursors in the axis set -- skipped", flush=True)
            continue
        media = (sets["catabolic_sources"] or sets["media"])[0]
        pairs = B.load_pairs(ref / "derived/mnxref-4_5/atom_pairs.parquet", X)
        npair = pairs.groupby("mnxr").size().to_dict()

        for lane in a.lanes:
            if lane == "netB":
                g0, w_add = B.graph_from_pairs(pairs, X, weights, ratios), w_median
            else:
                f = Path(a.netA_dir) / f"atom_graph_{X}.npz"
                if not f.exists():
                    print(f"  [{X}/{lane}] {f} missing -- skipped", flush=True)
                    continue
                # uniform E=1.0: a curated model asserts presence, not evidence strength
                g0, w_add = AtomGraph.load(f), 1.0
            if not g0.m:
                print(f"  [{X}/{lane}] empty graph -- skipped", flush=True)
                continue
            med = float(np.median(g0.gp))
            port, leak = a.port_mult * med, a.leak_mult * med
            comp = G.source_component(g0, media)
            print(f"  [{X}/{lane}] media={media} {g0.n:,}n/{g0.m:,}e med_gp={med:.6g} "
                  f"src-component={len(comp):,}n add_w={w_add:.6g}", flush=True)

            for mode in a.modes:
                base = G.measure(g0, media, sets, mode, port=port, leak=leak)
                print(f"    [{X}/{lane}/{mode}] base total={base['total']:.6g} "
                      f"prec_share={base['prec_share']:.6f}", flush=True)
                t0 = time.time()
                for k, c in enumerate(mine):
                    if (lane, mode, X, c["condition_id"]) in done:
                        continue
                    try:
                        gc, nadd = G.merge_reactions(g0, pairs, X, c["rxns"], w_add, ratios)
                        cur = G.measure(gc, media, sets, mode, port=port, leak=leak)
                        rep = G.added_edge_report(g0, gc, comp, media)
                    except Exception as e:      # one condition must not kill the shard
                        print(f"      !! {c['condition_id']} [{X}/{lane}/{mode}] "
                              f"{type(e).__name__}: {e}", flush=True)
                        continue
                    bt, ct = base["total"], cur["total"]
                    for m, dv in cur["draw"].items():
                        bv = base["draw"].get(m, 0.0)
                        w.writerow(dict(
                            lane=lane, mode=mode, element=X, media=media,
                            condition_id=c["condition_id"], n_rxns=c["n_rxns"],
                            rxns=";".join(c["rxns"]), metabolite=m,
                            draw_base=bv, draw_cond=dv, delta_draw=dv - bv,
                            rel_draw=((dv - bv) / bv) if abs(bv) > 0 else float("nan"),
                            total_base=bt, total_cond=ct, delta_total=ct - bt,
                            rel_total=((ct - bt) / bt) if abs(bt) > 0 else float("nan"),
                            prec_share_base=base["prec_share"],
                            prec_share_cond=cur["prec_share"],
                            leak_frac_base=base["leak_frac"],
                            leak_frac_cond=cur["leak_frac"],
                            is_novel=int(not any(r in weights for r in c["rxns"])),
                            n_pairs_element=sum(npair.get(r, 0) for r in c["rxns"]),
                            converged=cur["converged"], port=port, leak=leak,
                            add_weight=w_add, expected=c["expected"][:200], **rep))
                print(f"    [{X}/{lane}/{mode}] {len(mine)} conditions in "
                      f"{time.time() - t0:.0f}s", flush=True)
    fh.close()
    print(f"[shard {a.shard}] DONE -> {out}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
