"""Structural MCS member: Indigo ``automap`` over the MetaNetX universe.

The FOURTH AAM ensemble member, and the second NON-NEURAL one. RXNMapper and
LocalMapper are transformers; MetaCyc is a curated database. This member is neither --
it is a STRUCTURAL mapper: it reads the two molecules and computes the maximum common
substructure correspondence directly, the way RDT does, with no training corpus and no
lookup table. So like MetaCyc it is a genuinely INDEPENDENT vote (it is not in
``combine_aam.NEURAL_MEMBERS``, so a consensus that includes it is undiscounted), and
unlike MetaCyc it exists for EVERY reaction whose molecules have structure, not only the
ones a curator wrote down.

WHY A STRUCTURAL MEMBER EARNS ITS PLACE
---------------------------------------
Two gaps the SMILES-based members leave, that a structural mapper closes:

  * THE 512-TOKEN CASUALTIES (``unmapped_computable``). RXNMapper's transformer accepts
    at most 512 tokens, and the long reactions -- the ones with the most cofactors -- blow
    that limit and come back empty. ``forced_pairs`` rescues only the n==1 single-atom
    cases. Indigo ``automap`` has NO context window: it maps the long reactions the neural
    universe never reached, at ~0.5 ms each.
  * ATOMS THAT FLOW THROUGH A CARRIER (``acyl_carrier`` / ``generic_rgroup``). Where a
    reaction routes atoms through an acyl-[ACP] thioester or a generic R-group, the carrier
    body is drawn as a ``*`` super-node -- an atom the mapping must CONSERVE, never invent a
    bond into. ``automap`` maps the concrete atoms straight THROUGH the ``*`` (verified: the
    acetyl methyl of acetyl-[ACP] lands on acetoacetyl-[ACP] while the ``*`` bodies stay
    unmapped), and the engine's ``pairs_from_mapped`` then suppresses the ``*`` for free
    (``el not in ELEMENTS``). A naive placeholder cannot do this -- it has no atoms to route
    through -- which is exactly why ``ecspr_aam_rescue`` REFUSES the acyl-carrier class.

WHAT IT MAPS, AND WHAT IT LEAVES
--------------------------------
The member builds a reaction SMILES for every reaction whose participants ALL have a
structure, drawing on three structure sources, in order:

  1. MetaNetX ``chem_prop`` -- the concrete metabolites (including the ones written with a
     ``*`` R-group: automap conserves the ``*`` and maps the concrete remainder).
  2. ``ecspr_aam_rescue.placeholder_for`` -- deterministic redox/carrier STAND-INS for the
     structureless electron carriers (ferredoxin, quinone, thioredoxin, ...). Their atoms
     are scaffolding, suppressed downstream, exactly as in the rescue path.
  3. ``--resolved`` -- an identifier-keyed curated crosswalk (the ``ecspr_aam_rescue``
     schema and its gates), for the acyl carriers and unreconciled stubs MetaNetX left
     structureless. A resolved carrier is drawn with a ``*`` body so its concrete atoms are
     KEPT (they are the metabolite's own) while the unknown body cancels across the reaction.

A reaction with a participant none of these three can supply is simply absent from this
member's output -- the honest "this member did not map it", the AAM analogue of ratio 1.0,
left for another member or for T4.

THE STAND-IN CLAIM IS STILL TESTED, NOT TRUSTED. As in the rescue path, standing a
placeholder or a ``*``-bodied carrier into a reaction asserts the carrier is CONSERVED for
the scored element. That is checked per reaction and per element by
``ecspr_aam_rescue.concrete_balance`` and written to ``--out-balance``; the fused
extractor applies it as a filter, so a carrier that actually donated an atom is refused for
that element while remaining usable for the elements it is genuinely inert to. This member
EMITS the placeholder map and the per-element balance verdicts (the same two side tables the
rescue path emits) so the fusion step can suppress and gate exactly as it does for the
neural universe -- one balance instrument, shared across members.

DROP-IN. ``load_member`` returns ``{mnxr -> (mapped_rxn_smiles, confidence)}`` -- the shape
``combine_aam`` consumes -- with ``confidence = MCS_CONFIDENCE`` for every reaction. The
mapped SMILES flow through the engine's own ``pairs_from_mapped`` unchanged, so the only
thing that varies between this member and the neural ones is the mapping, never the identity
bookkeeping.

Runs under ``scadc-metabolic-model`` (rdkit + pandas + ``epam-indigo``). Indigo is the only
new dependency; the structure sources and the balance gate are all ``ecspr_aam_rescue``'s.
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

# The structure sources, the placeholder taxonomy, and the balance gate are the engine's --
# never a second copy, so this member and the rescue path stand in for a carrier identically.
ENGINE_LIB = Path("/home/tony/agentic_workspace/projects/metasmith-libraries/fabfos/resources/lib")
sys.path.insert(0, str(ENGINE_LIB))
from ecspr_aam_rescue import (placeholder_for, load_chem, load_equations,   # noqa: E402
                              parse_equation, load_resolved, gate_bodies_cancel,
                              build_rxn_smiles, concrete_balance, ELEMENTS,
                              SMILES_LEN_LIMIT)

# A structural MCS map is a deterministic INFERENCE from the two molecules -- not a neural
# model's scored guess, and not a curator's asserted fact. It carries full weight as an
# independent vote (it is not neural, so it is never subject to NEURAL_SHARED_FLOOR), and the
# combiner's disagreement-dilution is what tempers it where a symmetric molecule leaves the
# atom correspondence genuinely ambiguous. Kept a named constant so the combiner and the
# materialised TSV agree on one value, the way CURATED_CONFIDENCE does for MetaCyc.
MCS_CONFIDENCE = 1.0

# automap mode: DISCARD any incoming maps and compute fresh. Our built SMILES carry none,
# so this only guards against a structure source that smuggled one in.
AUTOMAP_MODE = "discard"

# Per-reaction MCS budget. Indigo's embedding enumerator is combinatorial: on a large or
# symmetric molecule the maximum-common-substructure search can blow up, and a few percent
# of MetaNetX reactions would otherwise hang the whole universe run. ``aam-timeout`` (ms)
# makes Indigo RAISE rather than spin, so a pathological reaction costs one budget and is
# skipped (``automap_failed``), never a stall. A reaction that times out is one whose atom
# correspondence is genuinely ambiguous anyway (symmetry) -- the graph would dilute it -- so
# losing it to this member is not losing a confident map. 1 s trades a little coverage for a
# bounded universe wall-clock.
AUTOMAP_TIMEOUT_MS = 1000  # canon-ok: Indigo automap timeout in ms, not canon.DRAW_K


def _indigo(timeout_ms: int = AUTOMAP_TIMEOUT_MS):
    """Lazily construct one Indigo session with the per-reaction MCS budget set. Imported
    here, not at module load, so the member's join/build logic can be exercised (and this
    file imported) under an env without ``epam-indigo`` -- only the mapping needs it."""
    from indigo import Indigo
    ind = Indigo()
    ind.setOption("aam-timeout", int(timeout_ms))
    return ind


def build_structure_map(chem: dict, resolved: dict) -> dict:
    """mnxm -> SMILES: the concrete metabolites, plus any curated structures ``--resolved``
    supplied. Placeholders are added per-reaction (they depend on which generics a reaction
    actually contains), so they are NOT folded in here."""
    smi = {m: v[2] for m, v in chem.items() if v[2]}
    smi.update(resolved)      # resolved rows are additive (load_resolved refuses overrides)
    return smi


def map_universe(reac_prop: Path, chem_prop: Path, resolved_path: Path | None,
                 targets: set | None = None, shard: int = 0, nshards: int = 1,
                 timeout_ms: int = AUTOMAP_TIMEOUT_MS,
                 out_fh=None, bal_fh=None, attempted_fh=None, done: set | None = None):
    """Map every buildable reaction in this shard with Indigo ``automap``.

    ``shard``/``nshards`` partition the universe DETERMINISTICALLY by sorted-index modulo,
    so K worker processes cover disjoint slices and the union is the whole universe -- the
    embedding search is single-threaded per Indigo session, so cores are used by running
    shards side by side, not by threading one session.

    RESUMABLE AGAINST A C-LEVEL HANG. Indigo's ``aam-timeout`` bounds the embedding
    enumerator, but a few reactions blow up in a DIFFERENT C path it does not cover and spin
    a full core indefinitely -- and a Python signal cannot interrupt a C loop that never
    returns to the interpreter. So the guard is external: this streams each result to
    ``out_fh``/``bal_fh`` (append + flush) as it is produced, and records every reaction it
    is ABOUT to hand to Indigo in ``attempted_fh`` (flush) BEFORE the call. A watchdog that
    kills a stalled shard and relaunches it with ``done`` = the attempted set therefore skips
    the one reaction that hung and never loses the thousands already mapped. Without a
    watchdog (``out_fh`` None) it behaves as before, returning everything in memory.

    Returns ``(member, placeholders, balance, tally)`` where
      * ``member``      : ``{mnxr -> (mapped_rxn_smiles, MCS_CONFIDENCE)}``
      * ``placeholders``: ``{mnxm -> (smiles, tag, name)}`` -- generics stood in for
      * ``balance``     : list of ``{mnxr, element, balanced}`` -- concrete-balance verdicts
      * ``tally``       : Counter of per-reaction outcomes (why each reaction landed where)
    """
    done = done or set()
    chem = load_chem(chem_prop)
    resolved = load_resolved(resolved_path, chem) if resolved_path else {}
    smi = build_structure_map(chem, resolved)
    print(f"[mcs] chem_prop: {len(chem):,} metabolites, {len(smi):,} with structure "
          f"({len(resolved):,} curated); {len(done):,} already done (resume)", flush=True)

    eqs = load_equations(reac_prop, targets)
    if nshards > 1:
        keys = [k for i, k in enumerate(sorted(eqs)) if i % nshards == shard]
        eqs = {k: eqs[k] for k in keys}
        print(f"[mcs] shard {shard}/{nshards}: {len(eqs):,} equations", flush=True)
    else:
        print(f"[mcs] {len(eqs):,} equations to consider", flush=True)

    ind = _indigo(timeout_ms)
    member, ph, balance_rows, tally = {}, {}, [], Counter()
    for mnxr, eq in eqs.items():
        if mnxr in done:
            tally["resumed_skip"] += 1
            continue
        pe = parse_equation(eq)
        if not pe:
            tally["unparseable_equation"] += 1
            continue
        subs, prods = pe
        parts = set(subs) | set(prods)

        # generics: participants with no concrete structure. Try a placeholder stand-in.
        gens = [m for m in parts if m not in smi]
        gen_ph = {m: placeholder_for(chem.get(m, ("", "", ""))[0]) for m in gens}
        if any(v is None for v in gen_ph.values()):
            tally["blocked_no_structure"] += 1     # acyl/stub/non-molecule -> another layer
            continue

        # a curated `*` body only cancels for balance when it stands on both sides
        if resolved and not gate_bodies_cancel(subs, prods, resolved):
            tally["refused_bodies_dont_cancel"] += 1
            continue

        merged = dict(smi)
        for m, (s, t) in gen_ph.items():
            merged[m] = s
        rxn_smi, bad = build_rxn_smiles(subs, prods, merged)
        if not rxn_smi:
            tally["build_failed"] += 1
            continue
        if len(rxn_smi) > SMILES_LEN_LIMIT:
            tally["too_long"] += 1
            continue

        # mark ATTEMPTED before the call that can hang, so a kill+resume skips exactly it.
        if attempted_fh is not None:
            attempted_fh.write(mnxr + "\n")
            attempted_fh.flush()
        try:
            rxn = ind.loadReaction(rxn_smi)
            rxn.automap(AUTOMAP_MODE)
            mapped = rxn.smiles()
        except Exception:
            tally["automap_failed"] += 1
            continue

        member[mnxr] = (mapped, MCS_CONFIDENCE)
        if out_fh is not None:
            out_fh.write(f"{mnxr}\t{mapped}\t{MCS_CONFIDENCE}\n")
            out_fh.flush()
        for m, (s, t) in gen_ph.items():
            ph[m] = (s, t, chem.get(m, ("", "", ""))[0])
        # the conservation claim, tested per element (the stand-ins are the ph set here)
        ph_set = set(gen_ph)
        for X in ELEMENTS:
            b = concrete_balance(subs, prods, chem, ph_set, X, resolved)
            if b is not None:
                balance_rows.append(dict(mnxr=mnxr, element=X, balanced=bool(b)))
                if bal_fh is not None:
                    bal_fh.write(f"{mnxr}\t{X}\t{bool(b)}\n")
                    bal_fh.flush()
        tally["mapped"] += 1
        if tally["mapped"] % 2000 == 0:
            print(f"[mcs]   mapped {tally['mapped']:,} (attempted {sum(tally.values()):,})",
                  flush=True)

    return member, ph, balance_rows, tally


def load_member(reac_prop: Path, chem_prop: Path, resolved: Path | None = None,
                targets: set | None = None) -> dict:
    """``{mnxr -> (mapped_rxn_smiles, confidence)}`` -- drop-in for ``combine_aam``'s member
    shape, so the combiner treats the MCS mapper as just another member."""
    member, _ph, _bal, _t = map_universe(Path(reac_prop), Path(chem_prop),
                                          Path(resolved) if resolved else None, targets)
    return member


def _default_paths():
    data = Path("/home/tony/agentic_workspace/data/scadc")
    return (data / "references/metanetx/reac_prop.tsv",
            data / "references/metanetx/chem_prop.tsv")


def main():
    reac_default, chem_default = _default_paths()
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reac-prop", type=Path, default=reac_default)
    ap.add_argument("--chem-prop", type=Path, default=chem_default)
    ap.add_argument("--resolved", type=Path, default=None,
                    help="curated structure crosswalk (ecspr_aam_rescue schema): acyl "
                         "carriers and unreconciled stubs MetaNetX left structureless. "
                         "Each row is an identifier-keyed ASSERTION, gated by its citation.")
    ap.add_argument("--targets", type=Path, default=None,
                    help="file of MNXR ids, one per line (default: the whole universe)")
    ap.add_argument("--nshards", type=int, default=1,
                    help="partition the universe into N deterministic shards for parallel "
                         "workers; the union covers everything (Indigo is single-threaded "
                         "per session, so cores are used by running shards side by side)")
    ap.add_argument("--shard", type=int, default=0, help="which shard (0..nshards-1)")
    ap.add_argument("--timeout-ms", type=int, default=AUTOMAP_TIMEOUT_MS,
                    help="per-reaction Indigo MCS budget (ms); a reaction that exceeds it is "
                         "skipped, not stalled")
    ap.add_argument("--out", type=Path,
                    default=Path(__file__).resolve().parent / "mcs_aam.tsv")
    ap.add_argument("--out-placeholders", type=Path, default=None,
                    help="placeholder map (mnxm, smiles, tag, name) -- the generics stood "
                         "in for; the extractor names then SUPPRESSES these. "
                         "Default: alongside --out.")
    ap.add_argument("--out-balance", type=Path, default=None,
                    help="per-(reaction, element) concrete-balance verdicts -- the extractor "
                         "drops unbalanced elements. Default: alongside --out.")
    ap.add_argument("--resume", action="store_true",
                    help="skip reactions already ATTEMPTED (in the .attempted sidecar) and "
                         "append -- so a watchdog can kill a stalled shard and relaunch it "
                         "past the reaction that hung, keeping everything already mapped.")
    a = ap.parse_args()

    targets = None
    if a.targets:
        targets = {l.strip() for l in open(a.targets) if l.strip()}

    # per-shard output suffix so K workers never collide; a merge step concatenates them.
    suf = f".shard{a.shard}" if a.nshards > 1 else ""
    out = a.out.with_name(a.out.stem + suf + a.out.suffix)
    out_ph = a.out_placeholders or a.out.with_name(a.out.stem + suf + "_placeholders.tsv")
    out_bal = a.out_balance or a.out.with_name(a.out.stem + suf + "_balance.tsv")
    out_att = out.with_suffix(out.suffix + ".attempted")
    out_done = out.with_suffix(out.suffix + ".done")

    # resume = skip everything already ATTEMPTED (a superset of what was mapped), so the one
    # reaction a prior run hung on -- recorded attempted just before it hung -- is skipped.
    done = set()
    fresh = not (a.resume and out.exists())
    if a.resume and out_att.exists():
        done = {l.strip() for l in open(out_att) if l.strip()}
    mode = "w" if fresh else "a"

    out_fh = open(out, mode)
    if fresh:
        out_fh.write("mnxr\tmapped_rxn_smiles\tconfidence\n")
        out_fh.flush()
    bal_fh = open(out_bal, mode)
    if fresh:
        bal_fh.write("mnxr\telement\tbalanced\n")
        bal_fh.flush()
    att_fh = open(out_att, mode)

    member, ph, balance_rows, tally = map_universe(
        a.reac_prop, a.chem_prop, a.resolved, targets,
        shard=a.shard, nshards=a.nshards, timeout_ms=a.timeout_ms,
        out_fh=out_fh, bal_fh=bal_fh, attempted_fh=att_fh, done=done)
    out_fh.close(); bal_fh.close(); att_fh.close()

    # placeholders are a static function of chem (every generic with an admissible stand-in);
    # written whole at the end. A kill loses only this small table, and the merge step
    # regenerates it globally, so it is never the resumability-critical artifact.
    pd.DataFrame([dict(mnxm=m, smiles=s, tag=t, name=n)
                  for m, (s, t, n) in ph.items()]).to_csv(out_ph, sep="\t", index=False)
    out_done.write_text("ok\n")   # clean-completion marker for the watchdog

    print("\n" + "=" * 64)
    print("MCS structural AAM member (Indigo automap) -- build report")
    print("=" * 64)
    for k, n in tally.most_common():
        print(f"    {k:<28} {n:>8,}")
    print(f"\n  mapped this run                : {tally.get('mapped', 0):,}")
    print(f"  distinct generics stood in for : {len(ph):,}")
    nbad = sum(1 for r in balance_rows if not r["balanced"])
    print(f"  balance verdicts               : {len(balance_rows):,} "
          f"({nbad:,} unbalanced -> refused for that element)")
    print(f"\nwrote {out}")
    print(f"      {out_ph}")
    print(f"      {out_bal}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
