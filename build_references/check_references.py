#!/usr/bin/env python3
"""Gate: the compiled reference tables must agree with their sources and each other.

    mamba run -n fabfos python build_references/check_references.py

Exits 1 on any failure. Findings that are true, permanent and NOT failures -- coverage
gaps, a lane that cannot be built, pseudogenes the model carries -- are printed under
NOTE, because a gap you can see is a different thing from a gap you cannot.

The load-bearing check is the equivalence one: build all four per-element atom graphs
from the compiled tables AND from ``ecspr_build.graph_from_pairs`` on the frozen string
tables, and require the same nodes, the same edge count and the same conductances. A
mismatch there means the encoding is wrong in a way no smaller test catches -- and the
dangerous direction of wrong (a too-narrow rank field merging two atoms onto one node)
RAISES conductance, so it reads as an improvement.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "src" / "metasmith_libraries" / "resources" / "lib"))
import hosts  # noqa: E402
import refs  # noqa: E402

failures: list[str] = []
notes: list[str] = []

# The epi300 evidence weights, used only as a realistic fixture for the equivalence
# check. Any weight set would prove the encoding; a real one also exercises the fanout.
WEIGHTS = refs.REPO / "data" / "_old" / "derived" / "evidence" / "evidence_weights.parquet"
DEPLOYED_KO_BRIDGE = Path(
    "/home/tony/agentic_workspace/projects/scadc/metabolic-modelling/main/"
    "metabolic-modelling/_reference_try1/betweenness/cache/ko_to_mnxr.tsv")


def check(label: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""),
          flush=True)
    if not ok:
        failures.append(f"{label}: {detail}")


def note(msg: str) -> None:
    print(f"  [NOTE] {msg}", flush=True)
    notes.append(msg)


# =====================================================================

def check_bake() -> dict:
    print("\nmetabolism -- the compiled tables are one artifact in three files")
    try:
        ident = refs.assert_same_bake()
        check("bake identity is byte-identical across all three files", True,
              f"bake {ident['vocab_sha256'][:16]}")
    except Exception as e:
        check("bake identity is byte-identical across all three files", False, str(e))
        return {}

    V = refs.load_vocab()
    check("vocabulary content hash matches the identity block",
          refs.vocab_sha256(V.df) == ident["vocab_sha256"])
    check("sources unchanged since the bake",
          refs.sha256_file(refs.SRC_ATOM_PAIRS) == ident["src_atom_pairs_sha256"]
          and refs.sha256_file(refs.SRC_DIRECTION) == ident["src_direction_sha256"])
    check("bit widths are measured, not assumed, and fit the budget",
          ident["max_atom_rank"] < (1 << ident["rank_bits"])
          and ident["n_met"] <= (1 << ident["met_bits"])
          and 2 * (ident["met_bits"] + ident["rank_bits"]) <= refs.NODE_KEY_BUDGET,
          f"met {ident['met_bits']} + rank {ident['rank_bits']} bits, "
          f"max rank {ident['max_atom_rank']} < {1 << ident['rank_bits']}, "
          f"edge key {2*(ident['met_bits']+ident['rank_bits'])} of {refs.NODE_KEY_BUDGET}")

    # The reaction vocabulary must be the UNION. Coding against the atom-pair set alone
    # would drop every reaction with direction evidence and no atom pairs, and each one
    # dropped falls back to ratio 1.0 -- fully reversible, i.e. MORE conductance than the
    # evidence supports. Silent, and in the flattering direction.
    d = refs.load_direction()
    ap = pd.read_parquet(refs.ATOM_PAIRS, columns=["rxn"])
    n_dir_only = len(set(d["rxn"]) - set(ap["rxn"]))
    check("reaction vocabulary is the union of both sources",
          V.size("rxn") == ident["n_rxn"] and n_dir_only > 0,
          f"{V.size('rxn'):,} reactions, {n_dir_only:,} of them carry direction "
          f"evidence but no atom pairs")

    # ratio > 1 FLIPS the edge, it does not amplify it. Pinning the count means a bake
    # that quietly pre-flipped, or narrowed the column, fails here.
    n_flip = int((d["ratio"].to_numpy() > 1.0).sum())
    check("edge-flipping ratio count is preserved", n_flip == 10485,
          f"{n_flip:,} of {len(d):,} ratios exceed 1")
    check("orientation is stored as written", ident["orientation"] == refs.ORIENTATION,
          ident["orientation"])

    total = sum(p.stat().st_size for p in refs.BAKED)
    src = refs.SRC_ATOM_PAIRS.resolve().stat().st_size + \
        refs.SRC_DIRECTION.resolve().stat().st_size
    check("compiled tables are smaller than their sources", total < src,
          f"{total/1e6:.2f} MB against {src/1e6:.2f} MB")
    return ident


def check_equivalence() -> None:
    print("\nequivalence -- the compiled tables build the same graph as the reference "
          "builder")
    if not WEIGHTS.exists():
        check("evidence-weight fixture is available", False, str(WEIGHTS))
        return
    from ecspr_build import graph_from_pairs

    w = pd.read_parquet(WEIGHTS)
    w = w[w["source"] == "epi300"]
    weights = {str(r): float(e) for r, e in zip(w["mnxr"], w["E_full"]) if float(e) > 0}
    src = pd.read_parquet(refs.SRC_ATOM_PAIRS)
    dsrc = pd.read_parquet(refs.SRC_DIRECTION, columns=["mnxr", "ratio"])
    ratios = {str(r): float(v) for r, v in zip(dsrc["mnxr"], dsrc["ratio"])
              if str(r) != "EMPTY"}

    V = refs.load_vocab()
    lut = refs.ratio_by_code(V)
    t0 = time.perf_counter()
    allp = refs.load_atom_pairs()
    t_load = time.perf_counter() - t0

    t_ref = t_bake = 0.0
    for X in refs.ELEMENT_ORDER:
        t0 = time.perf_counter()
        g = graph_from_pairs(src, X, weights, ratios)
        t_ref += time.perf_counter() - t0
        t0 = time.perf_counter()
        b = refs.compile_atom_graph(X, weights, vocab=V, pairs=allp, ratio_lut=lut)
        t_bake += time.perf_counter() - t0
        same = (g.n == b.n and g.m == b.m and set(g.nodes) == set(b.nodes)
                and np.allclose(np.sort(g.gp), np.sort(b.gp), rtol=1e-6, atol=0)
                and np.allclose(np.sort(g.gm), np.sort(b.gm), rtol=1e-6, atol=0))
        check(f"element {X}: same nodes, edges and conductances", same,
              f"{g.n:,} nodes / {g.m:,} edges from "
              f"{g.meta['n_reactions_used']:,} reactions")
    note(f"compile: {t_bake*1000:.0f} ms + {t_load*1000:.0f} ms load, against the "
         f"reference builder's {t_ref*1000:.0f} ms")


def check_gem_tables() -> None:
    print("\nGEM lane -- one table per host, in the shared schema")
    frames = {}
    for host in sorted(hosts.HOSTS):
        p = hosts.out_dir(host) / "gpr_gem.parquet"
        if not p.exists():
            check(f"{host}: gpr_gem.parquet exists", False, str(p))
            continue
        d = pd.read_parquet(p)
        frames[host] = d
        check(f"{host}: schema is the shared GPR schema",
              tuple(d.columns) == hosts.GPR_COLS,
              f"{len(d):,} rows, {d['mnxr'].nunique():,} MNXR")
        check(f"{host}: in_atom_universe is carried on every row",
              d["in_atom_universe"].notna().all() and d["mnxr"].notna().all())
        gap = d.loc[~d["in_atom_universe"], "mnxr"].nunique()
        note(f"{host} GEM: {gap:,} of {d['mnxr'].nunique():,} reactions have no "
             f"atom-pair coverage -- kept and flagged, never inner-joined away")

    # EPI300 has no model of its own; its table is the DH10B model under another tag.
    if {"e_coli_dh10b", "e_coli_epi300"} <= set(frames):
        cols = [c for c in hosts.GPR_COLS if c != "host"]
        check("e_coli_epi300's GEM table equals e_coli_dh10b's apart from the host tag",
              frames["e_coli_epi300"][cols].equals(frames["e_coli_dh10b"][cols]),
              "the edit list is empty; see gem_identity_report.md")
    rep = hosts.out_dir("e_coli_epi300") / "gem_identity_report.md"
    check("the EPI300 identity report exists", rep.exists(),
          str(rep.relative_to(refs.REPO)) if rep.exists() else "run "
          "check_epi300_identity.py")


def check_pseudogenes() -> None:
    print("\nannotation vs model -- genes the genome calls broken that the model keeps")
    gbk = hosts.genome_gbk("e_coli_dh10b")
    if not gbk.exists():
        note("dh10b genome not checked out; skipping the pseudogene report")
        return
    feat = hosts.parse_genbank_features(gbk)
    pseudo = set(feat.loc[feat["pseudo"], "gene"].dropna())
    d = pd.read_parquet(hosts.out_dir("e_coli_dh10b") / "gpr_gem.parquet")
    d = d[d["feature_kind"] == "gem_gene"]
    hit = d[d["feature_name"].isin(pseudo)]
    if hit.empty:
        note("no model gene is flagged /pseudo in the annotation")
        return
    for g, grp in hit.groupby("feature_name"):
        # every reaction they touch has an `or` alternative, so knocking them out changes
        # no reaction's liveness -- but a knockout study needs to know they are there
        solo = grp[grp["gpr_rule"] == grp["feature_id"].iloc[0]]
        note(f"{g} ({grp['feature_id'].iloc[0]}) is /pseudo in the annotation but is a "
             f"model gene on {grp['mnxr'].nunique():,} reactions "
             f"({len(solo)} where it is the only gene). Kept, so the table stays a "
             f"faithful reading of the model file.")


def check_denovo_tables() -> None:
    print("\nde-novo lane -- declared lane sets, and what could not be built")
    cfg = yaml.safe_load((HERE / "lanes.yml").read_text())
    unavailable = {n: s for n, s in cfg["bridges"].items()
                   if s.get("available") is False or not s.get("path")}
    for n, s in unavailable.items():
        note(f"bridge {n!r} is unavailable: {' '.join((s.get('reason') or '').split())}")

    for host in sorted(cfg["hosts"]):
        p = hosts.out_dir(host) / "gpr_denovo.parquet"
        if not p.exists():
            check(f"{host}: gpr_denovo.parquet exists", False, str(p))
            continue
        d = pd.read_parquet(p)
        check(f"{host}: schema is the shared GPR schema",
              tuple(d.columns) == hosts.GPR_COLS,
              f"{len(d):,} rows, {d['feature_id'].nunique():,} ORFs, "
              f"{d['mnxr'].nunique():,} MNXR")

        declared = set(cfg["hosts"][host]["lanes"])
        expected_blocked = {l for l, s in cfg["hosts"][host]["lanes"].items()
                            if (s.get("requires") or cfg["lanes"][l].get("bridge"))
                            in unavailable}
        built = set(d["channel"].unique())
        # A host that quietly produced a smaller table is the failure this asserts
        # against: the lane set is part of the file's identity, not an incidental.
        check(f"{host}: lane set is exactly the declared, buildable one",
              built == declared - expected_blocked,
              f"built {sorted(built)}; blocked {sorted(expected_blocked)}")
        gap = d.loc[~d["in_atom_universe"], "mnxr"].nunique()
        note(f"{host} de-novo: {gap:,} of {d['mnxr'].nunique():,} reactions have no "
             f"atom-pair coverage")

    # raw_score is per-lane and per-host in scale; asserting it is comparable would be
    # the error. Assert instead that every row declares its channel.
    for host in sorted(cfg["hosts"]):
        p = hosts.out_dir(host) / "gpr_denovo.parquet"
        if p.exists():
            d = pd.read_parquet(p, columns=["channel", "projection_via", "raw_score"])
            check(f"{host}: every row names its channel and projection",
                  d["channel"].notna().all() and d["projection_via"].notna().all())


def check_ko_bridge() -> None:
    print("\nKO bridge -- the ported builder must agree with the deployed one")
    p = refs.REPO / "data/reference/functional_annotation/bridges/ko_to_mnxr.tsv"
    if not p.exists():
        check("ko_to_mnxr exists", False, str(p))
        return
    new = pd.read_csv(p, sep="\t")
    note(f"ko_to_mnxr: {len(new):,} rows over {new['ko'].nunique():,} KOs, "
         f"{new['mnxr'].nunique():,} reactions "
         f"(median {new.groupby('ko')['mnxr'].nunique().median():.0f} reactions per KO)")
    if not DEPLOYED_KO_BRIDGE.exists():
        note("the deployed bridge is not readable here; skipping the agreement check")
        return
    dep = pd.read_csv(DEPLOYED_KO_BRIDGE, sep="\t")
    shared = set(dep["ko"]) & set(new["ko"])
    dp = {(k, m) for k, m in zip(dep["ko"], dep["mnx_r"]) if k in shared}
    np_ = {(k, m) for k, m in zip(new["ko"], new["mnxr"]) if k in shared}
    check("agrees with the deployed bridge on every shared KO", dp == np_,
          f"{len(shared):,} shared KOs, {len(dp & np_):,} pairs agree, "
          f"{len(dp - np_):,} deployed-only, {len(np_ - dp):,} new-only")
    note(f"the deployed bridge was scoped to the fosmid KOs ({dep['ko'].nunique():,}); "
         f"this one to the host KOs ({new['ko'].nunique():,}), which is why it is larger")


def main() -> int:
    print("== compiled reference tables ==")
    check_bake()
    check_equivalence()
    check_gem_tables()
    check_pseudogenes()
    check_denovo_tables()
    check_ko_bridge()

    print()
    if failures:
        print(f"{len(failures)} FAILED:")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"all checks passed ({len(notes)} notes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
