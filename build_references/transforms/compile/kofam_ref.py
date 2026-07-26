"""R3 -- unpack the KOfam profiles and carry the KO list through.

The ko_list passes through unchanged. It is re-emitted under the run-side type
rather than consumed straight from raw:: so that every input a run tool sees comes
from the compiled tier -- a run tool reaching into raw/ is how a pipeline quietly
acquires an undeclared dependency on the acquisition layout.
"""
from metasmith.python_api import *

lib   = TransformInstanceLibrary.ResolveParentLibrary(__file__)
model = Transform()

image      = model.AddRequirement(lib.GetType("env::python_for_data_science.env"))
archive    = model.AddRequirement(lib.GetType("raw::kofam_profiles_archive"))
raw_kolist = model.AddRequirement(lib.GetType("raw::kofam_ko_list"))
profiles   = model.AddProduct(lib.GetType("ref::kofamscan_profiles"))
ko_list    = model.AddProduct(lib.GetType("ref::kofamscan_ko_list"))

def protocol(context: ExecutionContext):
    iarchive = context.Input(archive)
    ikolist  = context.Input(raw_kolist)
    iprof    = context.Output(profiles)
    iko      = context.Output(ko_list)

    # The tarball unpacks to a `profiles/` directory; --strip-components=1 puts the .hmm
    # files directly under the product, because kofamscan is handed a profile DIRECTORY
    # and a nested extra level makes it find nothing while raising nothing.
    _cmd = f"""
        mkdir -p {iprof.container}
        tar -xzf {iarchive.container} -C {iprof.container} --strip-components=1
        cp {ikolist.container} {iko.container}
    """
    context.ExecWithEnv() \
        .ifContainerDo(env=image, cmd=_cmd) \
        .ifVirtualEnvDo(env=image, cmd=_cmd)

    n_hmm = len(list(iprof.local.glob("*.hmm"))) if iprof.local.is_dir() else 0
    Log.Info(f"unpacked {n_hmm:,} HMM profiles")
    return ExecutionResult(
        manifest=[{profiles: iprof.local, ko_list: iko.local}],
        success=n_hmm > 0 and iko.local.exists() and iko.local.stat().st_size > 0,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(8), duration=Duration(hours=1)),
)
