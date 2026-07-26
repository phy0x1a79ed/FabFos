"""The reference gate: check what the build produced, against what it claims.

    PATH="/home/tony/lib/miniforge3/envs/msm/bin:$PATH" \\
        python build_references/check_references.py --results <run>/results

Ported from the pre-library ``build_references/check_references.py``, with its path
registry replaced by explicit arguments -- the tables now come out of a metasmith run
directory, not out of a fixed layout under ``data/``.

NOTE vs FAILURE, which is the distinction the whole script is built around: a gap you can
SEE is a different thing from a gap you cannot. A reference that covers less than you
hoped is a NOTE; a reference that disagrees with itself is a FAILURE. Only the second
exits non-zero.

The load-bearing check is the equivalence one. ``bake_metabolism``'s own selftest proves
the encoding round-trips row by row; this proves the encoded tables BUILD THE SAME GRAPH
as the reference builder does from the string tables -- same nodes, same edge count, same
conductances, per element. Those are different claims: a bake can round-trip perfectly and
still be joined wrongly at compile time, and the resulting graph is plausible rather than
broken.

The constants check exists for a subtler reason. ``buildlib::dir_canon`` is a copy of
``src/fabfos/canon.py``'s ``DIR_*`` block, because the direction ensemble runs in a conda
env that has no import path to the fabfos package. Two copies of a constant drift, and
this pair drifts silently: one COMPUTES a ratio and the other VALIDATES a table carrying
one, so a divergence produces a table that passes its own validator while meaning
something else.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
MLIB = REPO / "src" / "metasmith_libraries"
BREF = REPO / "build_references"
sys.path.insert(0, str(BREF / "resources" / "buildlib"))
sys.path.insert(0, str(MLIB / "resources" / "lib"))

import refs_encoding as refs  # noqa: E402

FAILURES: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""),
          flush=True)
    if not ok:
        FAILURES.append(label)


def note(msg: str) -> None:
    print(f"  [NOTE] {msg}", flush=True)


def find(results: Path, name: str) -> Path | None:
    hits = sorted(results.rglob(name))
    return hits[0] if hits else None


# =====================================================================

def check_bake(vocab_p: Path, pairs_p: Path, dir_p: Path) -> dict | None:
    print("\nbake -- the trio is one artifact")
    try:
        ident = refs.assert_same_bake(vocab_p, pairs_p, dir_p)
    except Exception as e:
        check("all three files carry the same bake identity", False, str(e))
        return None
    check("all three files carry the same bake identity", True,
          f"bake {ident['vocab_sha256'][:16]}")
    check("rank field is wide enough for the observed maximum",
          ident["max_atom_rank"] < (1 << ident["rank_bits"]),
          f"{ident['max_atom_rank']} < {1 << ident['rank_bits']}")
    check("metabolite field is wide enough for the vocabulary",
          ident["n_met"] <= (1 << ident["met_bits"]),
          f"{ident['n_met']:,} <= {1 << ident['met_bits']:,}")
    check("edge key fits in int64",
          2 * (ident["met_bits"] + ident["rank_bits"]) <= refs.NODE_KEY_BUDGET)
    note(f"{ident['atom_pairs_rows']:,} atom pairs over {ident['n_rxn']:,} reactions, "
         f"{ident['n_met']:,} metabolites")
    return ident


def check_equivalence(ident: dict, vocab_p: Path, pairs_p: Path, dir_p: Path,
                      src_pairs: Path, src_dir: Path) -> None:
    """The compiled tables must build the same graph as the reference builder.

    ``ecspr_build.graph_from_pairs`` works on the STRING tables and is the definition;
    ``refs_encoding.compile_atom_graph`` works on the compiled ones and is the thing being
    checked. Node ORDER differs -- the baked table is sorted, so first-seen order differs
    -- which is why nodes are compared as sets and conductances sorted before comparison.
    """
    print("\nequivalence -- the compiled tables build the same graph as the builder")
    try:
        from ecspr_build import graph_from_pairs
    except ImportError as e:
        note(f"ecspr_build not importable ({e}); equivalence not checked")
        return
    if not (src_pairs.exists() and src_dir.exists()):
        note("the ensemble intermediates are not in this run's results, so the string "
             "tables the builder needs are unavailable; equivalence not checked")
        return

    src = pd.read_parquet(src_pairs)
    dsrc = pd.read_parquet(src_dir, columns=["mnxr", "ratio"])
    ratios = {str(r): float(v) for r, v in zip(dsrc["mnxr"], dsrc["ratio"])
              if str(r) != "EMPTY"}
    # Uniform weights over every reaction the pairs cover. The deployed gate used a
    # staged evidence-weight fixture; uniform E is the same test with one fewer input --
    # what is being compared is the join, the flip and the edge factorisation, none of
    # which depends on the weights being interesting.
    weights = {str(r): 1.0 for r in src["mnxr"].unique()}

    V = refs.load_vocab(vocab_p)
    direction = refs.load_direction(dir_p)
    lut = refs.ratio_by_code(V, direction)
    t0 = time.perf_counter()
    allp = refs.load_atom_pairs(pairs_p)
    t_load = time.perf_counter() - t0

    t_ref = t_bake = 0.0
    for X in refs.ELEMENT_ORDER:
        t0 = time.perf_counter()
        g = graph_from_pairs(src, X, weights, ratios)
        t_ref += time.perf_counter() - t0
        t0 = time.perf_counter()
        b = refs.compile_atom_graph(X, weights, ident=ident, vocab=V, pairs=allp,
                                    ratio_lut=lut)
        t_bake += time.perf_counter() - t0
        same = (g.n == b.n and g.m == b.m and set(g.nodes) == set(b.nodes)
                and np.allclose(np.sort(g.gp), np.sort(b.gp), rtol=1e-6, atol=0)
                and np.allclose(np.sort(g.gm), np.sort(b.gm), rtol=1e-6, atol=0))
        check(f"element {X}: same nodes, edges and conductances", same,
              f"{g.n:,} nodes / {g.m:,} edges from "
              f"{g.meta['n_reactions_used']:,} reactions")
    note(f"compile: {t_bake*1000:.0f} ms + {t_load*1000:.0f} ms load, against the "
         f"reference builder's {t_ref*1000:.0f} ms")


def check_direction_constants() -> None:
    """The build-side copy of the DIR_* block must equal the run-side one.

    See the module docstring: these are the same numbers used at two different times, and
    a divergence is invisible from either side alone.
    """
    print("\nconstants -- the direction ensemble's two copies agree")
    try:
        import dir_canon
        sys.path.insert(0, str(REPO / "src"))
        from fabfos import canon
    except Exception as e:                                       # pragma: no cover
        note(f"could not import both copies ({e}); constants not checked")
        return
    names = [n for n in dir(dir_canon) if n.startswith("DIR_")]
    bad = [n for n in names
           if not hasattr(canon, n) or getattr(canon, n) != getattr(dir_canon, n)]
    check(f"all {len(names)} DIR_* constants match src/fabfos/canon.py", not bad,
          f"diverged: {bad}" if bad else "")


def check_bridge(bridge_p: Path) -> None:
    print("\nbridge -- one table, three id spaces")
    if not bridge_p.exists():
        note("mnxr_lookup.parquet is not in this run's results; bridge not checked")
        return
    b = pd.read_parquet(bridge_p, columns=["id", "id_source", "mnxr", "evidence_quality"])
    per_id = b.groupby("id")["id_source"].nunique()
    clashes = int((per_id > 1).sum())
    # `id_source` is a LABEL, not a disambiguator -- the consumer slices on it and joins on
    # `id` alone, so a collision would mix two namespaces' claims into one lane.
    check("no id appears in more than one id_source", clashes == 0,
          f"{clashes:,} colliding ids" if clashes else f"{b['id'].nunique():,} distinct ids")
    dupes = int(b.duplicated(subset=["id", "mnxr"]).sum())
    check("deduped to distinct (id, mnxr)", dupes == 0, f"{dupes:,} duplicate pairs")
    note(f"{len(b):,} rows, {b['mnxr'].nunique():,} MNXR, "
         f"sources {b['id_source'].value_counts().to_dict()}")
    note(f"evidence_quality {b['evidence_quality'].value_counts().to_dict()}")


def check_gem_tables(results: Path) -> None:
    print("\nGEM GPR -- ruleless rows kept, and EPI300 == DH10B modulo host")
    tables = sorted(results.rglob("*gpr_table_gem*")) + sorted(results.rglob("gpr_gem*"))
    tables = [t for t in tables if t.suffix == ".parquet"]
    if not tables:
        note("no GEM GPR tables in this run's results; not checked")
        return
    frames = {}
    for t in tables:
        d = pd.read_parquet(t)
        host = d["host"].iat[0] if len(d) else t.stem
        frames[host] = d
        n_ruleless = int((d["feature_kind"] == "ruleless").sum())
        # Dropping ruleless reactions makes every gene set look like starvation, because
        # exchanges and spontaneous chemistry are live in every condition.
        check(f"{host}: ruleless reactions have rows", n_ruleless > 0,
              f"{n_ruleless:,} of {len(d):,}")
        uniform = bool((d["raw_score"] == 1.0).all())
        check(f"{host}: raw_score is uniform 1.0", uniform,
              "a curated model asserts presence, not evidence strength")

    if "e_coli_epi300" in frames and "e_coli_dh10b" in frames:
        cols = [c for c in frames["e_coli_epi300"].columns if c != "host"]
        same = frames["e_coli_epi300"][cols].reset_index(drop=True).equals(
            frames["e_coli_dh10b"][cols].reset_index(drop=True))
        # They share iECDH10B_1368 and the measured edit list is EMPTY, so identical is
        # the correct outcome; a divergence means one of them silently used another model.
        check("EPI300 and DH10B tables are identical apart from the host tag", same)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", type=Path, required=True,
                    help="a metasmith run's results directory")
    a = ap.parse_args()
    results = a.results
    if not results.exists():
        raise SystemExit(f"no results at {results}")

    print(f"== reference gate over {results} ==")
    vocab_p = find(results, "*metabolism_vocab*") or find(results, "vocab.parquet")
    pairs_p = find(results, "*atom_pairs*") or find(results, "atom_pairs.parquet")
    dir_p = find(results, "*direction_ratios*") or find(results, "direction.parquet")
    src_pairs = find(results, "*aam_pairs*")
    src_dir = find(results, "*direction_annotation*")
    bridge_p = find(results, "*mnxr_lookup*")

    if vocab_p and pairs_p and dir_p:
        ident = check_bake(vocab_p, pairs_p, dir_p)
        if ident and src_pairs and src_dir:
            check_equivalence(ident, vocab_p, pairs_p, dir_p, src_pairs, src_dir)
        elif ident:
            note("the ensemble intermediates are transient and were not collected, so "
                 "the equivalence check has no string tables to compare against")
    else:
        note("the metabolism trio is not in this run's results; bake not checked")

    check_direction_constants()
    if bridge_p:
        check_bridge(bridge_p)
    else:
        note("mnxr_lookup.parquet is not in this run's results; bridge not checked")
    check_gem_tables(results)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} FAILED: " + "; ".join(FAILURES))
        return 1
    print("all checks passed (see NOTEs above for what was not checked)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
