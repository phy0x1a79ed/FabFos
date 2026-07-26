"""KO -> KEGG reaction, from KEGG REST's bulk link endpoint.

An ACQUISITION, not a compile step. One request returns the entire mapping:

    https://rest.kegg.jp/link/reaction/ko     ~245 KB, ~2 s, 12,276 pairs

Classifying it here is what retires the dependency on the sibling project's
`kegg_requests.db` -- that cache was a licensed artifact being read at build time, which
is a licensing question rather than a path fix.

WHY THIS IS NOT THE DEPLOYED METHOD, which matters because everything else in this
library is a faithful port. The deployed builder (scadc `05_build_ko_to_mnxr.py`, via
`archive/_build_references_backup/build_bridges.py`) fetched one `get/<ko>` flat-file per
KO and parsed its REACTION block. It never ran that over the KO universe: it read a
185 MB sqlite cache of 3,855 previously-fetched flat-files and scoped the remainder to
the KOs the declared hosts actually called, because ~28k live requests for KOs nothing
annotates is three hours of rate-limited fetching to no purpose. Its own docstring says
so. The result shipped at 2,738 rows over 1,356 KOs.

That method cannot be reproduced here -- the cache is licensed and not redistributable,
and the `lanes.yml` that defined the host scope is gone. Rather than port the crawl and
widen it to the full KO list (~3 h, and ~22k of those requests return nothing, since only
22% of KOfam KOs have any reaction at all), this uses the bulk endpoint, which is a
strict improvement on both:

                              pairs      KOs     MNXR
    deployed cache + crawl    2,634    1,312    1,077
    one bulk link call       11,723    5,950    6,970

MEASURED EQUIVALENT, not assumed. Parsing the REACTION block out of all 3,854 cached
flat-files and comparing per KO against the bulk result: 3,838 identical (99.6%), 12 where
the bulk endpoint has more, 4 where the cache has more -- consistent with KEGG updates
since the cache was written, not with a parsing difference.

The consequence to know downstream: the kofam lane now reaches ~6.5x the reactions it did
in the deployed build, so its numbers are NOT comparable to the archived ones.

Takes the KO list as input to intersect against, and does not read its EC tags. Routing KO
through its `[EC:...]` tag instead of through KEGG reactions is easy and wrong -- measured
on the three hosts it takes the kofam lane from 646 to 8,548 reactions and makes 91% of
them reactions the EC lane already reaches, so the two lanes stop being independent
evidence. The same argument is why the EC route is NOT moved to KEGG's
`link/reaction/enzyme` even though that is one more cheap call: the EC lane is built from
MetaNetX's own classifications, and sourcing both lanes from KEGG would correlate them.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image    = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
ko_list  = model.AddRequirement(lib.GetType("raw::kofam_ko_list"))
out      = model.AddProduct(lib.GetType("raw::kegg_ko_to_kegg_r"))

KEGG_LINK_URL = "https://rest.kegg.jp/link/reaction/ko"

# A floor on the response, not a target. KEGG answers a malformed or unsupported query
# with HTTP 200 and an empty body, so "succeeded" and "returned nothing" are the same
# status code -- and a bridge built from an empty response is not an error anywhere
# downstream, it is a kofam lane that annotates nothing. Set well below the 12,276
# observed on 2026-07-25 so ordinary growth or shrinkage does not trip it.
MIN_PAIRS = 5000

DRIVER = r'''
import re, sys, urllib.request
from pathlib import Path

import pandas as pd

URL = "{url}"
OUT = Path("{out}")
MIN_PAIRS = {min_pairs}

# ko:K00001\trn:R00623 -- the only shape this endpoint emits.
LINE = re.compile(r"^ko:(K\d{{5}})\t+rn:(R\d{{5}})$")

print(f"[kegg] {{URL}}", flush=True)
with urllib.request.urlopen(URL, timeout=120) as r:
    if r.status != 200:
        raise SystemExit(f"KEGG returned HTTP {{r.status}}")
    text = r.read().decode()
print(f"[kegg] {{len(text):,}} bytes", flush=True)

rows, malformed = [], 0
for ln in text.splitlines():
    if not ln.strip():
        continue
    m = LINE.match(ln)
    if m:
        rows.append(m.groups())
    else:
        malformed += 1
if malformed:
    # Not fatal on its own, but it means the response is not what this parser was
    # written against, and the count below is then measuring the wrong thing.
    print(f"    [warn] {{malformed:,}} line(s) did not match ko:K#####\\trn:R#####",
          flush=True)
if len(rows) < MIN_PAIRS:
    raise SystemExit(
        f"KEGG returned {{len(rows):,}} pairs, under the {{MIN_PAIRS:,}} floor. "
        f"An empty or truncated response comes back as HTTP 200, so this is checked "
        f"rather than trusted: a bridge built from it would silently annotate nothing.")

df = pd.DataFrame(rows, columns=["ko", "kegg_r"]).drop_duplicates()

# Intersect with the KOfam profile set, the same `& known` the deployed builder applied.
# A KO with no profile can never come out of kofamscan, so its rows could only ever be
# dead weight -- but the count is printed rather than dropped in silence, because a large
# number here would mean the two KEGG releases have drifted apart.
known = set(pd.read_csv("{ko_list}", sep="\t", dtype=str)["knum"].dropna())
before = df["ko"].nunique()
df = df[df["ko"].isin(known)]
dropped = before - df["ko"].nunique()
print(f"[kegg] {{before:,}} KOs carry a reaction; {{dropped:,}} have no KOfam profile "
      f"and were dropped", flush=True)

df = df.sort_values(["ko", "kegg_r"]).reset_index(drop=True)
df.to_csv(OUT, sep="\t", index=False)
print(f"[kegg] {{len(df):,}} rows over {{df['ko'].nunique():,}} KOs and "
      f"{{df['kegg_r'].nunique():,}} KEGG reactions -> {{OUT}}", flush=True)
'''


def protocol(context: ExecutionContext):
    iko = context.Input(ko_list)
    iout = context.Output(out)
    driver = DRIVER.format(url=KEGG_LINK_URL, min_pairs=MIN_PAIRS,
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
    # One request. The 8 h this used to declare was for a 28k-KO crawl that no longer
    # happens; the whole step is now seconds, and the headroom is for KEGG being slow.
    resources=Resources(cpus=1, memory=Size.GB(2), duration=Duration(minutes=15)),
)
