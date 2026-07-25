"""Benchmark condition sources -- four publications' supplementary tables.

One transform, four products, because they are acquired together as the benchmark's
citation set and there is no case where one is wanted without the others.

  laser     engineering records; an upstream GIT CHECKOUT, pinned by url+sha256
            rather than hashed as data bytes -- hashing a tree of git repos turns
            its .git into thousands of cache objects and it stops being a checkout
  keio      single-gene knockout collection (Baba 2006) -> the `lof` cohort
  eydallin  screen hits (Eydallin 2010)                 -> the `eydallin` cohort
  het       heterologous / metagenomic screening tables -> the `gof_het` cohort

PROVISIONAL. The exact file inventory is confirmed while porting benchmark_v4,
whose `_v4common.py` reaches into absolute paths in the sibling tree that have to
be resolved to chunks one at a time.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
laser     = model.AddProduct(lib.GetType("raw::laser_records"))
keio      = model.AddProduct(lib.GetType("raw::keio_records"))
eydallin  = model.AddProduct(lib.GetType("raw::eydallin_records"))
het       = model.AddProduct(lib.GetType("raw::het_screen_records"))

LASER_URL = "https://bitbucket.org/jdwinkler/laser_release"

def protocol(context: ExecutionContext):
    raise NotImplementedError(
        "contract sketch only -- acquire/literature.py declares what it consumes and "
        "produces so the planner can resolve the DAG; the fetch is not written yet."
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(hours=1)),
)
