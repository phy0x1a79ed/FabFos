"""eQuilibrator compound cache.

Pinned explicitly rather than left to the package's first-use fetch. Otherwise a
build silently depends on a network call whose result is not recorded anywhere,
and two builds a year apart can disagree with nothing in the provenance to show it.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image = model.AddRequirement(lib.GetType("env::equilibrator.env"))
cache = model.AddProduct(lib.GetType("raw::equilibrator_cache"))

# The package fetches its compound cache into a platformdirs location on first use.
# Pointing that at the product directory is what turns an implicit network call into a
# recorded acquisition: the same env var is set when the direction ensemble runs, so the
# member reads exactly the cache this step pinned.
CACHE_ENV = "EQUILIBRATOR_CACHE_DIR"


def protocol(context: ExecutionContext):
    icache = context.Output(cache)
    icache.local.mkdir(parents=True, exist_ok=True)
    # Instantiating ComponentContribution is what triggers the download; there is no
    # public "just fetch" entry point, and calling one reaction through it is the
    # cheapest way to prove the cache is complete rather than half-written.
    context.ExecWithContainer(image=image, cmd=f"""
        export {CACHE_ENV}={icache.container}
        export XDG_CACHE_HOME={icache.container}
        python3 - <<'PY'
from equilibrator_api import ComponentContribution
cc = ComponentContribution()
w = cc.get_compound("kegg:C00001")
print("[equilibrator] cache primed;", "water resolved" if w is not None else "water MISSING")
raise SystemExit(0 if w is not None else 1)
PY
    """)
    ok = icache.local.exists() and any(icache.local.rglob("*"))
    return ExecutionResult(
        manifest=[{cache: icache.local}],
        success=ok,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(hours=2)),
)
