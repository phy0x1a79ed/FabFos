"""Read MetaCyc's ``reactions.dat`` attribute-value flat file.

Why this module exists. The deployed curated member read a **sqlite pgdb** -- a staged
per-reaction JSON table in a sibling project -- and pulled ``REACTION-DIRECTION``,
``LEFT`` and ``RIGHT`` out of it with a ``SELECT``. That pgdb is a derived artifact of a
tree this repo does not have, and staging it would put a derived file in the acquisition
tier. The licensed drop-in ships ``reactions.dat`` itself, so the member reads the source
directly and this module is the only new code the switch costs.

The format (documented at ``bioinformatics.ai.sri.com/ptools/flatfile-format.html``):

  * ``#`` at column 0 is a comment; the header block is all comments.
  * ``KEY - VALUE`` is one attribute. A key may repeat -- ``LEFT`` appears once per
    left-hand compound -- so every attribute is collected as a LIST, never overwritten.
    Overwriting is how a two-substrate reaction silently becomes a one-substrate one.
  * ``^SUBKEY - VALUE`` annotates the attribute line above it (``^COMPARTMENT - CCO-IN``
    under a ``LEFT``). Sub-slots are attached to their parent rather than dropped,
    because the compartment is what distinguishes the two sides of a transport reaction.
  * ``//`` on its own line ends a record.
  * A value may continue onto the next line with a leading ``/`` (rare, and only in free
    text such as COMMENT); joined rather than treated as a new attribute.

Compound ids are returned VERBATIM, bare -- ``GLT``, ``PROTON``, ``WATER`` -- which is
the form ``dir_refdata.load_metacyc_compound_to_mnxm`` keys on.

Encoding: the file is latin-1 in practice (the copyright line carries a non-UTF-8 (c)),
so it is opened that way rather than with errors="replace", which would corrupt a
compound id rather than a copyright symbol if one ever appeared.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterator

ENCODING = "latin-1"


def iter_records(path: Path) -> Iterator[dict]:
    """Yield one ``{ATTRIBUTE: [values...]}`` dict per record.

    Sub-slot lines are yielded under the key ``"^<PARENT>.<SUBKEY>"`` positionally
    aligned with the parent list, so ``rec["LEFT"][i]`` and
    ``rec["^LEFT.COMPARTMENT"][i]`` describe the same participant. A participant with no
    sub-slot gets ``None`` in that position rather than shifting the list.
    """
    rec: dict[str, list] = {}
    last_key: str | None = None

    def flush():
        nonlocal rec, last_key
        if rec:
            yield_rec, rec, last_key = rec, {}, None
            return yield_rec
        rec, last_key = {}, None
        return None

    with open(path, encoding=ENCODING) as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if not line or line.startswith("#"):
                continue
            if line == "//":
                out = flush()
                if out is not None:
                    yield out
                continue
            if line.startswith("^"):
                if last_key is None:
                    continue                      # a sub-slot with no parent: nothing to attach to
                sub, _, val = line[1:].partition(" - ")
                key = f"^{last_key}.{sub.strip()}"
                slot = rec.setdefault(key, [])
                # pad to the parent's current length - 1, so index i lines up with parent i
                target = len(rec[last_key]) - 1
                while len(slot) < target:
                    slot.append(None)
                if len(slot) == target:
                    slot.append(val.strip())
                else:                              # a second sub-slot of the same name
                    slot[target] = val.strip()
                continue
            if line.startswith("/"):
                # continuation of the previous value (free text only)
                if last_key and rec.get(last_key):
                    rec[last_key][-1] += " " + line[1:].strip()
                continue
            key, sep, val = line.partition(" - ")
            if not sep:
                continue
            key = key.strip()
            rec.setdefault(key, []).append(val.strip())
            last_key = key
        out = flush()
        if out is not None:
            yield out


def load_reactions(path: Path) -> list[dict]:
    """Every record that has a UNIQUE-ID, as ``{unique_id, direction, left, right}``.

    ``direction`` is the single ``REACTION-DIRECTION`` value or ``None``. It is asserted
    single-valued: MetaCyc writes at most one, and quietly taking the first of several
    would be picking a direction, which is the one thing this whole chain exists not to do.
    """
    out = []
    for rec in iter_records(path):
        ids = rec.get("UNIQUE-ID")
        if not ids:
            continue
        d = rec.get("REACTION-DIRECTION") or []
        if len(d) > 1:
            raise ValueError(
                f"{ids[0]}: {len(d)} REACTION-DIRECTION values {d} -- MetaCyc writes at "
                f"most one, and choosing between them is not this reader's call")
        out.append(dict(
            unique_id=ids[0],
            direction=(d[0] if d else None),
            left=list(rec.get("LEFT") or []),
            right=list(rec.get("RIGHT") or []),
        ))
    return out
