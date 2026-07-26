"""KO -> KEGG reaction, from each KO flat-file's REACTION block.

An ACQUISITION, not a compile step: every row comes from KEGG REST, one
`get/<ko>` per KO at ~3 requests/second. Classifying it here is what retires the
dependency on the sibling project's `kegg_requests.db` -- that cache was a
licensed artifact being read at build time, which is a licensing question rather
than a path fix.

Takes the KO list as input because that is the KO universe to fetch; it does not
read its EC tags. Routing KO through its `[EC:...]` tag instead of through KEGG
reactions is easy and wrong -- measured on the three hosts it takes the kofam lane
from 646 to 8,548 reactions and makes 91% of them reactions the EC lane already
reaches, so the two lanes stop being independent evidence.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image    = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
ko_list  = model.AddRequirement(lib.GetType("raw::kofam_ko_list"))
out      = model.AddProduct(lib.GetType("raw::kegg_ko_to_kegg_r"))

KEGG_URL   = "https://rest.kegg.jp/get/{ko}"
KEGG_PAUSE = 0.4

# The whole KO list, not a host-derived subset. The deployed builder scoped the fetch to
# the KOs the declared hosts actually called (via a lanes.yml that no longer exists), and
# a bridge scoped that way silently narrows whenever a host is added. This transform's
# only input is the KO list, so its scope is the KO universe -- ~28k requests, ~3 h at
# the 0.4 s spacing KEGG asks for, inside the declared 8 h. It is also the more honest
# scope: the bridge stops depending on which hosts happen to be declared.
DRIVER = r'''
import re, sys, time, urllib.request
from pathlib import Path

import pandas as pd

KEGG_URL = "{kegg_url}"
PAUSE = {pause}
OUT = Path("{out}")

# The REACTION block of a KEGG K-entry flat-file, and the R-numbers on it. Both lifted
# VERBATIM from the deployed builder (scadc 05_build_ko_to_mnxr.py, via
# archive/_build_references_backup/build_bridges.py) so this bridge cannot drift from
# the one that produced the deployed lane.
_REACTION_BLOCK = re.compile(r"^REACTION\s+(.+?)(?=\n[A-Z]+\s|\Z)", re.S | re.M)
_R_ID = re.compile(r"\bR\d{{5}}\b")


def reactions_of(text):
    if not text:
        return set()
    m = _REACTION_BLOCK.search(text)
    return set(_R_ID.findall(m.group(1))) if m else set()


kos = sorted(set(pd.read_csv("{ko_list}", sep="\t", dtype=str)["knum"].dropna()))
print(f"[kegg] {{len(kos):,}} KOs in the list", flush=True)

# Resume off the partial output. The fetch is hours long and rate-limited, so a restart
# that began again from KO 1 would be a second three-hour run, and a run that gave up
# would leave a bridge that is silently narrower than the deployed one.
done = set()
if OUT.exists() and OUT.stat().st_size > 0:
    prev = pd.read_csv(OUT, sep="\t", dtype=str)
    done = set(prev["ko"].dropna())
    print(f"[kegg] resume: {{len(done):,}} KOs already resolved", flush=True)

# A KO with no REACTION block emits no row, so the output alone cannot distinguish
# "not fetched" from "fetched, nothing to say". The sidecar records the second, which is
# what makes the resume exact instead of re-fetching ~20k silent KOs every restart.
SIDECAR = OUT.with_suffix(".no_reaction.txt")
if SIDECAR.exists():
    silent = {{ln.strip() for ln in SIDECAR.read_text().splitlines() if ln.strip()}}
    done |= silent
    print(f"[kegg] resume: {{len(silent):,}} KOs known to carry no REACTION block", flush=True)

todo = [k for k in kos if k not in done]
print(f"[kegg] fetching {{len(todo):,}} KOs (~{{len(todo)*PAUSE/60:.0f}} min at {{PAUSE}}s spacing)",
      flush=True)

write_header = not OUT.exists() or OUT.stat().st_size == 0
out_fh = open(OUT, "a", buffering=1)
if write_header:
    out_fh.write("ko\tkegg_r\n")
sc_fh = open(SIDECAR, "a", buffering=1)

n_rows = n_silent = n_failed = 0
for i, ko in enumerate(todo, 1):
    try:
        with urllib.request.urlopen(KEGG_URL.format(ko=ko), timeout=20) as r:
            text = r.read().decode()
        time.sleep(PAUSE)
    except Exception as e:
        n_failed += 1
        print(f"    [warn] {{ko}}: {{type(e).__name__}}: {{e}}", flush=True)
        continue
    rs = reactions_of(text)
    if not rs:
        n_silent += 1
        sc_fh.write(ko + "\n")
        continue
    for rid in sorted(rs):
        out_fh.write(f"{{ko}}\t{{rid}}\n")
        n_rows += 1
    if i % 500 == 0:
        print(f"    {{i:,}}/{{len(todo):,}}  rows={{n_rows:,}} silent={{n_silent:,}} "
              f"failed={{n_failed:,}}", flush=True)
out_fh.close()
sc_fh.close()

df = pd.read_csv(OUT, sep="\t", dtype=str).drop_duplicates().sort_values(["ko", "kegg_r"])
df.to_csv(OUT, sep="\t", index=False)
print(f"[kegg] {{len(df):,}} rows, {{df['ko'].nunique():,}} KOs carry a reaction; "
      f"{{n_silent:,}} carry no REACTION block, {{n_failed:,}} failed", flush=True)
if n_failed:
    # A partial bridge projects the kofam lane through fewer reactions than it should,
    # which reads downstream as a host simply not having those functions.
    raise SystemExit(f"{{n_failed}} KO fetches failed; re-run to resume rather than "
                     f"building the bridge from an incomplete fetch")
'''


def protocol(context: ExecutionContext):
    iko = context.Input(ko_list)
    iout = context.Output(out)
    driver = DRIVER.format(kegg_url=KEGG_URL, pause=KEGG_PAUSE,
                           ko_list=iko.container, out=iout.container)
    context.LocalShell("cat > _kegg_ko_reactions.py << 'PYEOF'\n" + driver + "\nPYEOF\n")
    context.ExecWithContainer(image=image, cmd="python3 _kegg_ko_reactions.py")
    return ExecutionResult(
        manifest=[{out: iout.local}],
        success=iout.local.exists() and iout.local.stat().st_size > 0,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(hours=8)),
)
