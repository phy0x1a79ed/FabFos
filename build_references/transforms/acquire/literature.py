"""Benchmark condition sources that can be FETCHED -- LASER and Keio.

Split from the two that cannot. This transform's two products come off the network; its
sibling `literature_extractions.py` carries the two that are hand extractions from paper
supplements and can only be verified.

The split is not tidiness. All four used to be one transform, and that transform could
never succeed: it verified the two extracted tables inside its own OUTPUT directory,
which nothing but itself creates, so the drop-in it asked for had nowhere to be dropped.
Worse, it coupled them -- LASER being absent scheduled a step that then refused on
Eydallin, so a missing git checkout reported itself as a missing PDF extraction.

Separated, each product's absence means what it says, and the extracted pair is satisfied
the way a drop-in should be: by being on disk under `data/raw/literature/` before the
build runs, where the stage 1 driver stages it and never schedules a producer at all.

  laser  engineering records; an upstream GIT CHECKOUT, pinned by url+rev rather than
         hashed as data bytes -- hashing a tree of git repos turns its .git into
         thousands of cache objects and it stops being a checkout
  keio   single-gene knockout collection (Baba 2006) -> the `lof` cohort
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image  = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
laser  = model.AddProduct(lib.GetType("raw::laser_records"))
keio   = model.AddProduct(lib.GetType("raw::keio_records"))

LASER_URL = "https://bitbucket.org/jdwinkler/laser_release"

# Keio ships as one supplementary zip off the paper's landing page; the four tables the
# `lof` cohort reads are inside it. Fetched as the zip and unpacked here rather than
# file-by-file, so the acquisition is one citable object.
KEIO_URL = ("https://www.embopress.org/action/downloadSupplement"
            "?doi=10.1038%2Fmsb4100050&file=msb4100050-sup-0001.zip")


def protocol(context: ExecutionContext):
    ilaser = context.Output(laser)
    ikeio  = context.Output(keio)

    _cmd = f"""
        git clone --depth 1 {LASER_URL}.git {ilaser.container}
        rev=$(git -C {ilaser.container} rev-parse HEAD)
        echo "$rev" > {ilaser.container}/.PINNED_REV
        echo "[literature] laser at $rev"

        mkdir -p {ikeio.container}
        wget -q "{KEIO_URL}" -O {ikeio.container}/keio_baba2006_supplementary.zip
        ( cd {ikeio.container} && unzip -o -q keio_baba2006_supplementary.zip )
    """
    context.ExecWithEnv() \
        .ifContainerDo(env=image, cmd=_cmd) \
        .ifVirtualEnvDo(env=image, cmd=_cmd)

    ok = (ilaser.local / ".PINNED_REV").exists() and any(ikeio.local.glob("*.xls"))
    return ExecutionResult(
        manifest=[{laser: ilaser.local, keio: ikeio.local}],
        success=ok,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(hours=1)),
)
