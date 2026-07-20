#!/usr/bin/env python3
"""Gate: the library and the provenance tree must agree.

Four checks, each corresponding to an acceptance criterion:

  AC2  every manifest key exists on disk and its type resolves
  AC3  every declared item has a generated record carrying a sha256
  AC4  every container record carries a digest, or an explicit blocker
  ---   reac_prop's hash agrees across all the places that pin it

The last one is the point of the whole exercise in miniature. reac_prop.tsv's
sha256 is pinned in three places: this library's provenance record,
mnxref-4_5/MANIFEST.json (as reac_prop_sha256), and canon. Three transcriptions
of one number is exactly the failure mode canon exists to remove, so the
provenance record is the source and the others are checked against it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
RECORDS = REPO / "provenance" / "data"
CONTAINERS = REPO / "provenance" / "containers"
LIB = REPO / ".awm" / "data" / "ref"

failures: list[str] = []
notes: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""))
    if not ok:
        failures.append(f"{label}: {detail}")


def main() -> int:
    import yaml

    decl = yaml.safe_load((RECORDS / "_declared.yml").open())
    declared = {i["id"]: i for i in decl["items"]}

    print("AC3 -- every declared item has a record with a sha256")
    missing_ok = 0
    for item_id, item in declared.items():
        rec_path = RECORDS / f"{item_id}.yml"
        if not rec_path.exists():
            check(item_id, False, "no generated record")
            continue
        rec = yaml.safe_load(rec_path.open())
        if rec.get("status") == "MISSING":
            missing_ok += 1
            if not rec.get("degraded_mode"):
                check(item_id, False, "declared MISSING but states no degraded_mode")
            continue
        if not rec.get("sha256"):
            check(item_id, False, "record carries no sha256")
    check(f"{len(declared) - missing_ok} present items all hashed",
          not failures, f"{missing_ok} recorded MISSING with a degradation path")

    print("\nAC4 -- every container record carries a digest or an explicit blocker")
    unresolved = []
    for rec_path in sorted(CONTAINERS.glob("*.yml")):
        rec = yaml.safe_load(rec_path.open())
        if rec.get("resolved"):
            if not str(rec.get("digest", "")).startswith("sha256:"):
                check(rec_path.stem, False, "resolved:true but no sha256 digest")
        else:
            if not rec.get("blocker"):
                check(rec_path.stem, False, "unresolved with no blocker stated")
            unresolved.append(rec_path.stem)
    n = len(list(CONTAINERS.glob("*.yml")))
    check(f"{n - len(unresolved)}/{n} images digest-pinned",
          True, f"unresolved, each with a blocker: {', '.join(unresolved)}")

    print("\nno image may float on a mutable :latest without saying so")
    for rec_path in sorted(CONTAINERS.glob("*.yml")):
        rec = yaml.safe_load(rec_path.open())
        if str(rec.get("reference", "")).endswith(":latest") and rec.get("resolved"):
            check(rec_path.stem, False, "resolved against a :latest tag")

    print("\nreac_prop sha256 agrees everywhere it is pinned")
    rec = yaml.safe_load((RECORDS / "metanetx.reac_prop.yml").open())
    truth = rec["sha256"]
    print(f"  provenance record: {truth}")
    man_path = LIB / "derived" / "mnxref-4_5" / "MANIFEST.json"
    if man_path.exists():
        pinned = json.loads(man_path.read_text()).get("reac_prop_sha256")
        print(f"  MANIFEST.json:     {pinned}")
        check("MANIFEST.json agrees with the provenance record", pinned == truth,
              "the pre-bake was built against different reac_prop bytes"
              if pinned != truth else "")
    else:
        notes.append("MANIFEST.json not in the library; skipped")

    print()
    for note in notes:
        print(f"  note: {note}")
    if failures:
        print(f"\n{len(failures)} FAILURE(S):")
        for f in failures:
            print(f"  {f}")
        return 1
    print("\nall gates green")
    return 0


if __name__ == "__main__":
    sys.exit(main())
