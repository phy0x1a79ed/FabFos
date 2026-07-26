"""The two neural AAM members: build reaction SMILES, map them, cache the result.

Ported from the deployed pair ``19b_rxnmapper_universe.py`` (RXNMapper over the whole
MNX reaction universe) and ``18_localmapper_inference.py`` (LocalMapper), which shared
this shape and differed only in the mapper call. Merged into one module because the two
members must see the SAME reaction SMILES: they are the CORRELATED pair whose agreement
the combiner discounts, and that discount only means anything if the disagreement it
measures is between mappers rather than between two slightly different SMILES builders.

Reaction SMILES are built from ``reac_prop``'s MNXR equation and ``chem_prop``'s
per-metabolite SMILES, stoichiometry expanded (a coefficient of 2 writes the metabolite
twice), so an atom map's indices line up with the participants the equation names. A
reaction whose participants are not all structurally known is not mapped at all -- a
partial reaction SMILES would map the atoms it does have onto the wrong destinations.

RESUMABLE, and it has to be: this is the longest single step in the reference build
(~57.5k reactions, hours per member). The cache is an append-only TSV keyed on ``mnxr``;
a restart re-reads it and does only what is missing, so a killed run costs the reaction
it was on and nothing else.

Output is one TSV per member, in the schema ``aam_combine.load_mapper`` reads:

    mnxr | rxn_smiles | mapped_rxn_smiles | confidence

No tau cut is applied here. The combiner is what decides how much a member's confidence
is worth, and a low-confidence map is still a vote to dilute against -- dropping it here
would remove the disagreement the ensemble is built to record.
"""
from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import pandas as pd

# The equation term grammar, verbatim from the deployed builder. NOTE it matches only
# `MNXM...@compartment` terms -- specials like `WATER@MNXD1` and `BIOMASS@MNXD1` do NOT
# match, which is deliberate here (they have no chem_prop SMILES) and is also the exact
# mismatch that made the extractor refuse every water-bearing reaction when the same
# regex was reused downstream. See ecspr_atom_pairs.EQ_TERM for that side of it.
EQ_TERM = re.compile(r"(\d+(?:\.\d+)?)\s+(MNXM\w+)@\w+")

# A reaction SMILES past this length is a polymer/macromolecule the transformers cannot
# usefully attend over, and attempting it costs minutes for a map nothing will trust.
SMILES_LEN_LIMIT = 8000

COLUMNS = ("mnxr", "rxn_smiles", "mapped_rxn_smiles", "confidence")


def parse_equation(eq: str):
    if "=" not in eq:
        return None
    lhs, rhs = eq.split("=", 1)

    def side(s):
        return [(float(c), m) for c, m in EQ_TERM.findall(s)]

    L, R = side(lhs), side(rhs)
    if not L or not R:
        return None
    return L, R


def load_smiles(chem_prop: Path, mnxm_set: set[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(chem_prop) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            mnxm = parts[0]
            if mnxm not in mnxm_set:
                continue
            smi = parts[8].strip()
            if smi:
                out[mnxm] = smi
    return out


def load_all_equations(reac_prop: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(reac_prop) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            mnxr = parts[0]
            if mnxr == "EMPTY":
                continue
            out[mnxr] = parts[1]
    return out


def build_rxn_smiles(eq: str, smi_map: dict[str, str]):
    parsed = parse_equation(eq)
    if parsed is None:
        return None, "unparseable_eq"
    L, R = parsed

    def side(side_list):
        out = []
        for coef, m in side_list:
            s = smi_map.get(m)
            if not s:
                # One unknown structure poisons the whole reaction: mapping the rest
                # would hand the missing participant's atoms to whatever is left.
                return None, f"no_smiles:{m}"
            n = max(1, int(round(coef)))
            for _ in range(n):
                out.append(s)
        return out, None

    Ls, err = side(L)
    if err:
        return None, err
    Rs, err = side(R)
    if err:
        return None, err
    return ".".join(Ls) + ">>" + ".".join(Rs), None


def build_universe(reac_prop: Path, chem_prop: Path) -> dict[str, str]:
    """``{mnxr -> reaction SMILES}`` for every reaction both members will be asked to map."""
    equations = load_all_equations(reac_prop)
    needed = set()
    for eq in equations.values():
        for _, m in EQ_TERM.findall(eq):
            needed.add(m)
    smi_map = load_smiles(chem_prop, needed)
    print(f"[aam] {len(equations):,} equations, {len(needed):,} MNXM needed, "
          f"{len(smi_map):,} with SMILES", flush=True)

    rxn_smiles: dict[str, str] = {}
    n_unbuild = 0
    for mnxr, eq in equations.items():
        s, err = build_rxn_smiles(eq, smi_map)
        if err or not s or len(s) > SMILES_LEN_LIMIT:
            n_unbuild += 1
            continue
        rxn_smiles[mnxr] = s
    print(f"[aam] buildable reaction SMILES: {len(rxn_smiles):,}  "
          f"unbuildable/too-large: {n_unbuild:,}", flush=True)
    return rxn_smiles


def _resume(out_tsv: Path) -> set[str]:
    if not out_tsv.exists() or out_tsv.stat().st_size == 0:
        return set()
    try:
        prev = pd.read_csv(out_tsv, sep="\t")
        done = set(prev["mnxr"].astype(str))
        print(f"[aam] resume: {len(done):,} reactions already cached", flush=True)
        return done
    except Exception as e:                                   # a truncated final line
        print(f"[aam] resume failed ({e}); starting over", file=sys.stderr, flush=True)
        return set()


def _writer(out_tsv: Path):
    write_header = not out_tsv.exists() or out_tsv.stat().st_size == 0
    fh = open(out_tsv, "a", buffering=1)                     # line-buffered: a kill -9
    if write_header:                                          # loses at most one row
        fh.write("\t".join(COLUMNS) + "\n")

    def emit(mnxr, smi, mapped, conf):
        vals = []
        for v in (mnxr, smi, mapped, conf):
            s = "" if v is None else str(v)
            vals.append(s.replace("\t", " ").replace("\n", " ").replace("\r", " "))
        fh.write("\t".join(vals) + "\n")

    return fh, emit


def run_rxnmapper(todo: list[tuple[str, str]], emit, chunk_size: int = 4):
    from rxnmapper import RXNMapper
    m = RXNMapper()
    t0 = time.time()
    n_ok = n_fail = 0
    i = 0
    while i < len(todo):
        j = min(i + chunk_size, len(todo))
        chunk = todo[i:j]
        smis = [s for _, s in chunk]
        try:
            res = m.get_attention_guided_atom_maps(smis)
        except Exception:
            # Per-reaction fallback so one bad reaction does not kill a whole chunk.
            res = []
            for one in smis:
                try:
                    r = m.get_attention_guided_atom_maps([one])
                    res.append(r[0] if r else {})
                except Exception:
                    res.append({"mapped_rxn": "", "confidence": float("nan")})
        for (mnxr, smi), r in zip(chunk, res):
            mapped = r.get("mapped_rxn", "")
            n_ok, n_fail = (n_ok + 1, n_fail) if mapped else (n_ok, n_fail + 1)
            emit(mnxr, smi, mapped, float(r.get("confidence", float("nan"))))
        i = j
        if i % 400 == 0 or i == len(todo):
            dt = time.time() - t0
            rate = i / dt if dt > 0 else 0
            print(f"      {i}/{len(todo)}  ok={n_ok} fail={n_fail}  "
                  f"({dt:.0f}s, {rate:.1f} rxn/s, eta {(len(todo)-i)/rate if rate else 0:.0f}s)",
                  flush=True)


def run_localmapper(todo: list[tuple[str, str]], emit):
    from localmapper import localmapper
    mapper = localmapper(device="cpu")
    t0 = time.time()
    n_ok = n_fail = 0
    for i, (mnxr, smi) in enumerate(todo, 1):
        try:
            r = mapper.get_atom_map(smi, return_dict=True)
            mapped = r.get("mapped_rxn", "") if isinstance(r, dict) else (r or "")
            # LocalMapper reports a template-match confidence under one of two names
            # depending on version; absent means it did not score, not that it scored 0.
            conf = float(r.get("confident", r.get("confidence", float("nan")))) \
                if isinstance(r, dict) else float("nan")
        except Exception:
            mapped, conf = "", float("nan")
        n_ok, n_fail = (n_ok + 1, n_fail) if mapped else (n_ok, n_fail + 1)
        emit(mnxr, smi, mapped, conf)
        if i % 200 == 0 or i == len(todo):
            dt = time.time() - t0
            rate = i / dt if dt > 0 else 0
            print(f"      {i}/{len(todo)}  ok={n_ok} fail={n_fail}  "
                  f"({dt:.0f}s, {rate:.1f} rxn/s, eta {(len(todo)-i)/rate if rate else 0:.0f}s)",
                  flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--member", required=True, choices=["rxnmapper", "localmapper"])
    ap.add_argument("--reac-prop", type=Path, required=True)
    ap.add_argument("--chem-prop", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True, help="the member's cached TSV")
    ap.add_argument("--limit", type=int, default=None,
                    help="cap the number of reactions mapped this invocation; the cache "
                         "is resumable, so this splits one long run rather than "
                         "shrinking the member")
    a = ap.parse_args()

    universe = build_universe(a.reac_prop, a.chem_prop)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    done = _resume(a.out)
    todo = [(m, s) for m, s in universe.items() if m not in done]
    if a.limit:
        todo = todo[:a.limit]
    print(f"[aam:{a.member}] todo: {len(todo):,}", flush=True)
    if not todo:
        print(f"[aam:{a.member}] nothing to do", flush=True)
        return 0

    fh, emit = _writer(a.out)
    try:
        if a.member == "rxnmapper":
            run_rxnmapper(todo, emit)
        else:
            run_localmapper(todo, emit)
    finally:
        fh.close()

    df = pd.read_csv(a.out, sep="\t")
    n_mapped = int((df["mapped_rxn_smiles"].astype(str).str.len() > 0).sum())
    print(f"[aam:{a.member}] {len(df):,} rows cached, {n_mapped:,} carry a map "
          f"-> {a.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
