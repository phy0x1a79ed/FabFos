"""KOfam profiles + KO list.

The archive is acquired as shipped and untarred by compile/kofam_ref.py, so the
acquisition stays a byte-for-byte copy of upstream and the unpacking is a separate,
re-runnable step.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image     = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
profiles  = model.AddProduct(lib.GetType("raw::kofam_profiles_archive"))
ko_list   = model.AddProduct(lib.GetType("raw::kofam_ko_list"))

PROFILES_URL = "ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz"
KO_LIST_URL  = "ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz"

def protocol(context: ExecutionContext):
    iprofiles = context.Output(profiles)
    iko_list = context.Output(ko_list)
    # The archive is moved, not unpacked: compile/kofam_ref.py untars it. Keeping the
    # acquisition a byte-for-byte copy of upstream is what lets the unpack be re-run
    # without re-fetching 400 MB, and is why this is not the shipped
    # logistics/downloadKofamDB.py (which fuses the two steps).
    _cmd = f"""
        wget -q {PROFILES_URL} -O profiles.tar.gz
        wget -q {KO_LIST_URL} -O ko_list.gz
        gunzip -c ko_list.gz > {iko_list.container}
        mv profiles.tar.gz {iprofiles.container}
    """
    context.ExecWithEnv() \
        .ifContainerDo(env=image, cmd=_cmd) \
        .ifVirtualEnvDo(env=image, cmd=_cmd)
    return ExecutionResult(
        manifest=[{profiles: iprofiles.local, ko_list: iko_list.local}],
        success=iprofiles.local.exists() and iko_list.local.exists(),
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(hours=4)),
)
