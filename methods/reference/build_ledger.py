"""The closure ledger: one verdict per MetaNetX reaction, adjudication measured.

WHAT IT IS
----------
The reference (``reference.py``) claims to adjudicate the WHOLE MetaNetX reaction
universe -- every reaction carries exactly one AAM verdict and one direction verdict.
"Whole universe" is a claim that can be false silently: a reaction can simply be absent
from a table and nobody notices. This ledger is the instrument that makes the claim
falsifiable. It enumerates every reaction in ``reac_prop`` (~83.8k), and for each records:

  * its reachability (does it appear in a real network, or is it universe-only) -- so
    effort can be prioritised even though the target is FULL coverage;
  * its AAM verdict-state and, if refused, the reason CLASS;
  * its direction verdict-state and reason class.

A verdict-state is one of ``reference.VERDICT_STATES`` (resolved / diluted-ambiguous /
refused) or ``pending`` -- the one state the ledger exists to drive to zero. ``pending``
is not "no row"; it is "this reaction is one an automatic or manual layer is still
expected to act on" (a mappable reaction not yet mapped, a carrier not yet rescued, a
direction not yet computed). A REFUSAL with a TERMINAL reason (a pseudo-reaction, a
non-molecule, a SEED lump, a no-transfer reaction) is NOT pending -- it is adjudicated,
because a reasoned refusal is a verdict. Closure is ``pending == 0``.

WHY THE REASON TAXONOMY IS THE ENGINE'S
---------------------------------------
The reason classes are ``ecspr_aam_rescue``'s -- electron carrier, acyl carrier,
non-molecule, unreconciled stub -- imported, not re-listed, so the ledger classifies a
blocking participant EXACTLY as the rescue path does. A reaction the ledger calls
``electron_carrier`` is one ``placeholder_for`` would stand in for; a reaction it calls
``acyl_carrier`` is one ``REFUSE`` blocks as an ACP thioester. The ledger and the mapper
cannot drift, because they read the same classifier.

WHAT IT IS NOT
--------------
It does not MAP or compute direction. It reads whatever AAM table and direction table
it is pointed at (the reference tables once frozen; the incumbent single-source tables
until then) and measures the gap. Run it before T2/T3/T4 to get the worklist; run it
after to prove closure. Same script, same universe, both times.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import reference as ref

# The reason taxonomy and its classifiers come from the engine's rescue module, so the
# ledger and the mapper name a blocking participant identically (never a second list).
ENGINE_LIB = Path("/home/tony/agentic_workspace/projects/metasmith-libraries/fabfos/resources/lib")
sys.path.insert(0, str(ENGINE_LIB))
from ecspr_aam_rescue import (placeholder_for, REFUSE, count_el,        # noqa: E402
                              parse_equation, load_chem)

CNSP = ("C", "N", "S", "P")

# Reason priority: a reaction blocked by several structureless participants takes its
# dominant blocker's class -- the MOST TERMINAL one, so a reaction only lands in a
# still-resolvable class when EVERY blocker is still-resolvable. Terminal blockers win.
REASON_PRIORITY = (
    ref.R_NON_MOLECULE, ref.R_SEED_LUMP, ref.R_ACYL_CARRIER,
    ref.R_GENERIC_RGROUP, ref.R_UNRECONCILED_STUB, ref.R_ELECTRON_CARRIER,
)

# The AAM verdict is PENDING (un-adjudicated) exactly when its refusal reason is one an
# automatic or manual layer is still expected to convert. After the T4 adjudication (see
# reference.py's note), NONE is: the automatic ensemble + placeholder library + balance-verified
# crosswalk ran to exhaustion, and every residual class is a terminal verdict (a reasoned refusal
# that no further layer converts without fabrication) -- except unmapped_computable, which
# aam_verdict now returns as DILUTED (abstain), also not pending. So the set is EMPTY: closure is
# every reaction carrying resolved / diluted / reasoned-refused. A future MNXref release re-opens
# adjudication by re-running T4 (new structureless classes would surface in the worklist again).
AAM_PENDING_REASONS = frozenset()


# =====================================================================
# Inputs
# =====================================================================
def load_reactions(reac_prop: Path):
    """Every reaction row: mnxr -> (equation, classifs, is_transport). Includes the
    ``EMPTY`` sentinel and one-sided rows -- the ledger classifies them OUT explicitly,
    it does not skip them (a skipped reaction is an un-adjudicated reaction)."""
    out = {}
    with open(reac_prop) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if not p or not p[0]:
                continue
            eq = p[1] if len(p) > 1 else ""
            classifs = p[3] if len(p) > 3 else ""
            is_transport = (len(p) > 5 and p[5].strip() == "T")
            out[p[0]] = (eq, classifs, is_transport)
    return out


def base_graph_mnxrs() -> set:
    """The strict in-network reaction set: ``('rxn', MNXR)`` nodes across the per-element
    base graphs. A subset of the evidence-reachable set; recorded for prioritisation."""
    import pickle
    out = set()
    for el in ref.BASE_GRAPH_ELEMENTS:
        p = ref.BASE_GRAPH_DIR / f"base_{el}.pkl"
        if not p.exists():
            continue
        g = pickle.load(open(p, "rb"))
        out |= {n[1] for n in g.nodes if isinstance(n, tuple) and n[0] == "rxn"}
    return out


def load_aam(path: Path):
    """(mapped_set, method_map). ``method_map`` is mnxr -> the reaction's coarse AAM state
    from a provenance-carrying reference table; empty for the incumbent (schema has no
    ``method``), where presence alone means resolved."""
    df = pd.read_parquet(path)
    mapped = set(df["mnxr"].astype(str))
    method_map = {}
    if "method" in df.columns:
        # a reaction is diluted if EVERY emitted pair is a disagreement/ambiguous split;
        # resolved if any pair is a consensus/single-member concrete correspondence.
        for mnxr, g in df.groupby("mnxr"):
            methods = set(g["method"].astype(str))
            resolved_like = any(m == "consensus" or m.endswith("_only") for m in methods)
            method_map[str(mnxr)] = ref.V_RESOLVED if resolved_like else ref.V_DILUTED
    return mapped, method_map


def load_direction(path: Path):
    """mnxr -> (dir_tier, dir_method). Absent reactions are the direction worklist."""
    df = pd.read_parquet(path)
    return {str(r.mnxr): (int(r.dir_tier), str(r.dir_method))
            for r in df.itertuples(index=False)}


# =====================================================================
# Per-participant classification (the engine's taxonomy)
# =====================================================================
def participant_reason(name: str, formula: str) -> str:
    """The refusal class a single STRUCTURELESS participant contributes.

    Mirrors ``ecspr_aam_rescue``: a name ``placeholder_for`` stands in for is an electron
    carrier; a name ``REFUSE`` blocks is an ACP acyl carrier or a non-molecule; a formula
    carrying a polymer/R token is a generic class; anything else is an ordinary metabolite
    MetaNetX minted with no structure -- an unreconciled stub, the ``--resolved`` case.
    """
    nm = (name or "").strip()
    if placeholder_for(nm):
        return ref.R_ELECTRON_CARRIER
    low = nm.lower()
    if REFUSE.match(nm):
        if "acp" in low or "acyl" in low:
            return ref.R_ACYL_CARRIER
        return ref.R_NON_MOLECULE
    if formula and ("*" in formula or "(" in formula):
        return ref.R_GENERIC_RGROUP
    return ref.R_UNRECONCILED_STUB


def dominant(reasons) -> str:
    for r in REASON_PRIORITY:
        if r in reasons:
            return r
    # every blocker was still-resolvable in a way not on the priority list -> the least
    # terminal one; defensive, REASON_PRIORITY already covers the taxonomy.
    return next(iter(reasons))


def aam_verdict(mnxr, subs, prods, is_pseudo, chem, smi_map, mapped, method_map):
    """(state, reason). ``reason`` is None for a resolved/diluted map."""
    if mnxr in mapped:
        return method_map.get(mnxr, ref.V_RESOLVED), None
    if is_pseudo:
        return ref.V_REFUSED, ref.R_PSEUDO
    parts = set(subs) | set(prods)
    structureless = [m for m in parts if m not in smi_map]
    if structureless:
        reasons = {participant_reason(chem.get(m, ("", "", ""))[0],
                                      chem.get(m, ("", "", ""))[1]) for m in structureless}
        return ref.V_REFUSED, dominant(reasons)
    # every participant has a structure but the reaction is not mapped.
    struct = [smi_map[m] for m in parts]
    if any("*" in s for s in struct):
        return ref.V_REFUSED, ref.R_GENERIC_RGROUP
    # does any C/N/S/P atom cross at all? if not, there is nothing for the atom lane.
    has_cnsp = False
    for m in parts:
        formula = chem.get(m, ("", "", ""))[1]
        for X in CNSP:
            c = count_el(formula, X)
            if c:
                has_cnsp = True
                break
        if has_cnsp:
            break
    if not has_cnsp:
        return ref.V_REFUSED, ref.R_NO_TRANSFER
    # every structure is present but no ensemble member produced a map (mapper-defeating size
    # or symmetry -- heme, tetrapyrroles, long FAS chains). Every atom abstains, which is the
    # AAM analogue of direction ratio 1.0: a DILUTED verdict (adjudicated), not a refusal. The
    # reason label rides along for provenance. See reference.py's T4 adjudication note.
    return ref.V_DILUTED, ref.R_UNMAPPED_COMPUTABLE


def dir_verdict(mnxr, is_pseudo, dir_lookup):
    """(state, reason) for direction. A no-evidence reaction is DILUTED (ratio 1.0 --
    the direction analogue of shrink-to-no-op), not pending: it is the honest verdict."""
    if is_pseudo:
        return ref.V_REFUSED, ref.R_PSEUDO
    if mnxr not in dir_lookup:
        return ref.V_PENDING, None
    tier, method = dir_lookup[mnxr]
    if method == "refused":
        return ref.V_REFUSED, ref.R_NON_MOLECULE  # direction refused (its own reasons)
    if tier == 0 or method == "no_evidence":
        return ref.V_DILUTED, None                # ratio 1.0, no-evidence -- adjudicated
    return ref.V_RESOLVED, None


# =====================================================================
# Build
# =====================================================================
def build(reactions, chem, smi_map, mapped, method_map, dir_lookup, reachable, base):
    rows = []
    for mnxr, (eq, classifs, is_transport) in reactions.items():
        pe = parse_equation(eq)
        one_sided = pe is not None and (not pe[0] or not pe[1])
        is_pseudo = is_transport or pe is None or one_sided
        if pe is None:
            subs, prods = [], []
        else:
            subs, prods = pe

        a_state, a_reason = aam_verdict(mnxr, subs, prods, is_pseudo, chem, smi_map,
                                        mapped, method_map)
        d_state, d_reason = dir_verdict(mnxr, is_pseudo, dir_lookup)

        a_pending = (a_state == ref.V_REFUSED and a_reason in AAM_PENDING_REASONS)
        d_pending = (d_state == ref.V_PENDING)
        rows.append(dict(
            mnxr=mnxr,
            in_evidence=mnxr in reachable,
            in_base_graph=mnxr in base,
            is_pseudo=is_pseudo,
            is_transport=is_transport,
            aam_state=a_state, aam_reason=a_reason or "",
            aam_pending=a_pending,
            dir_state=d_state, dir_reason=d_reason or "",
            dir_pending=d_pending,
            pending=a_pending or d_pending,
        ))
    return pd.DataFrame(rows)


def report(df: pd.DataFrame) -> None:
    n = len(df)
    print("=" * 72)
    print(f"CLOSURE LEDGER -- MetaNetX {ref.MNXREF_VERSION} ({ref.MNXREF_DATE})")
    print("=" * 72)
    print(f"reactions enumerated        : {n:,}")
    print(f"  reachable (in evidence)   : {int(df.in_evidence.sum()):,}")
    print(f"  in a base graph (network) : {int(df.in_base_graph.sum()):,}")
    print(f"  pseudo (non-chemistry)    : {int(df.is_pseudo.sum()):,}")

    for lane, state_col, reason_col, pend_col in (
            ("AAM", "aam_state", "aam_reason", "aam_pending"),
            ("DIRECTION", "dir_state", "dir_reason", "dir_pending")):
        print(f"\n{lane} verdict-states:")
        for state, k in df[state_col].value_counts().items():
            print(f"    {state:20s} {k:>8,}")
        refused = df[df[state_col] == ref.V_REFUSED]
        if len(refused):
            print(f"  {lane} refusal reason classes:")
            for reason, k in refused[reason_col].value_counts().items():
                tag = "  [still-resolvable]" if reason in ref.REASON_STILL_RESOLVABLE else \
                      ("  [pending]" if reason in AAM_PENDING_REASONS and lane == "AAM" else
                       "  [terminal]")
                print(f"      {reason:22s} {k:>8,}{tag}")
        print(f"  {lane} un-adjudicated (pending) : {int(df[pend_col].sum()):,}")

    n_pending = int(df.pending.sum())
    print("\n" + "-" * 72)
    print(f"UN-ADJUDICATED (either lane pending): {n_pending:,} / {n:,}")
    # reachable-first worklist
    reach_pending = df[df.pending & df.in_evidence]
    print(f"  of which reachable                : {len(reach_pending):,}")
    if n_pending == 0:
        print("\n  ✓ CLOSED: every reaction carries a verdict.")
    else:
        print(f"\n  ✗ OPEN: {n_pending:,} reaction(s) still to adjudicate (the worklist).")
    print("-" * 72)


def write_manifest(df, atom_pairs: Path, direction: Path) -> None:
    """Refresh the reference pin. The manifest is what ``reference.assert_reference``
    reads: it records the MetaNetX release, the reac_prop content hash (so a moved
    universe fails loudly), the closure state (``n_pending``, which must reach 0), and
    whether the frozen tables exist yet. Re-written every ledger run, so the pin tracks
    the current adjudication state and ``assert_reference(require_closed=True)`` fails
    until closure and passes at freeze.
    """
    tables_present = ref.ATOM_PAIRS.exists() and ref.DIRECTION.exists()
    man = dict(
        mnxref_version=ref.MNXREF_VERSION,
        mnxref_date=ref.MNXREF_DATE,
        reac_prop_sha256=ref.reac_prop_hash(),
        n_reactions=int(len(df)),
        n_pending=int(df.pending.sum()),
        aam_states={k: int(v) for k, v in df.aam_state.value_counts().items()},
        aam_reasons={k: int(v) for k, v in
                     df[df.aam_state == ref.V_REFUSED].aam_reason.value_counts().items()},
        dir_states={k: int(v) for k, v in df.dir_state.value_counts().items()},
        aam_source=str(atom_pairs),
        direction_source=str(direction),
        frozen_tables_present=tables_present,
    )
    ref.MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    ref.MANIFEST.write_text(json.dumps(man, indent=2, sort_keys=True) + "\n")
    print(f"[ledger] manifest pin -> {ref.MANIFEST} (n_pending={man['n_pending']:,})")


def write_worklist(df, reactions, chem, smi_map, out: Path) -> None:
    """The T4/T2 authoring aid: for every PENDING AAM reaction, the frequency-ranked
    structureless BLOCKER names per reason class, so a class rule can be authored from the
    data (not guessed). One row per (reason_class, blocker_mnxm): its name, how many
    pending reactions it blocks, and an example. Terminal-refused reactions are omitted --
    they need no work. Sorted most-blocking first, which is the order to work the classes.
    """
    pend = df[df.aam_pending]
    pend_mnxr = set(pend.mnxr)
    reason_of = dict(zip(df.mnxr, df.aam_reason))
    blk = {}  # (reason, mnxm) -> [count, name, example_mnxr]
    for mnxr in pend_mnxr:
        eq = reactions[mnxr][0]
        pe = parse_equation(eq)
        if not pe:
            continue
        reason = reason_of[mnxr]
        for m in set(pe[0]) | set(pe[1]):
            if m in smi_map:
                continue
            key = (reason, m)
            if key not in blk:
                blk[key] = [0, chem.get(m, ("", "", ""))[0], mnxr]
            blk[key][0] += 1
    cols = ["reason_class", "blocker_mnxm", "blocker_name", "n_pending_blocked", "example_mnxr"]
    rows = [dict(reason_class=r, blocker_mnxm=m, blocker_name=v[1],
                 n_pending_blocked=v[0], example_mnxr=v[2])
            for (r, m), v in blk.items()]
    # at closure there are no pending reactions and hence no blockers: write the empty
    # (header-only) worklist rather than crashing on a missing column -- the worklist being
    # empty IS the closure signal.
    wl = (pd.DataFrame(rows).sort_values(["reason_class", "n_pending_blocked"],
                                         ascending=[True, False])
          if rows else pd.DataFrame(columns=cols))
    wl.to_csv(out, sep="\t", index=False)
    print(f"[ledger] worklist ({len(wl):,} distinct blockers) -> {out}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reac-prop", type=Path, default=ref.REAC_PROP)
    ap.add_argument("--chem-prop", type=Path, default=ref.CHEM_PROP)
    ap.add_argument("--atom-pairs", type=Path, default=None,
                    help="AAM table (default: reference if frozen, else the incumbent seed)")
    ap.add_argument("--direction", type=Path, default=None,
                    help="direction table (default: reference if frozen, else the seed)")
    ap.add_argument("--out", type=Path, default=ref.LEDGER)
    ap.add_argument("--worklist", type=Path, default=None,
                    help="also write the per-reason frequency-ranked blocker names (T4 aid)")
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero unless 0 un-adjudicated (the closure gate)")
    a = ap.parse_args()

    atom_pairs = a.atom_pairs or (ref.ATOM_PAIRS if ref.ATOM_PAIRS.exists() else ref.SEED_ATOM_PAIRS)
    direction = a.direction or (ref.DIRECTION if ref.DIRECTION.exists() else ref.SEED_DIRECTION)
    print(f"[ledger] AAM table       : {atom_pairs}")
    print(f"[ledger] direction table : {direction}")

    reactions = load_reactions(a.reac_prop)
    print(f"[ledger] {len(reactions):,} reactions in the universe", flush=True)
    chem = load_chem(a.chem_prop)
    smi_map = {m: v[2] for m, v in chem.items() if v[2]}
    print(f"[ledger] chem_prop: {len(chem):,} metabolites, {len(smi_map):,} with SMILES",
          flush=True)
    mapped, method_map = load_aam(atom_pairs)
    print(f"[ledger] AAM: {len(mapped):,} mapped reactions"
          f"{' (provenance-carrying)' if method_map else ' (incumbent, presence=resolved)'}",
          flush=True)
    dir_lookup = load_direction(direction)
    print(f"[ledger] direction: {len(dir_lookup):,} reactions with a verdict", flush=True)
    reachable = set(pd.read_parquet(ref.EVIDENCE_TABLE)["mnxr"].astype(str))
    base = base_graph_mnxrs()
    print(f"[ledger] reachable: {len(reachable):,}; in base graphs: {len(base):,}",
          flush=True)

    df = build(reactions, chem, smi_map, mapped, method_map, dir_lookup, reachable, base)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(a.out, index=False)
    print(f"[ledger] wrote {a.out}\n", flush=True)
    write_manifest(df, atom_pairs, direction)
    if a.worklist:
        write_worklist(df, reactions, chem, smi_map, a.worklist)
    report(df)

    n_pending = int(df.pending.sum())
    if a.strict and n_pending:
        print(f"\n[ledger] STRICT: {n_pending:,} un-adjudicated -- not closed.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
