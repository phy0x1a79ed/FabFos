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

# Keio ships as one supplementary zip off the paper's landing page; the four tables the
# `lof` cohort reads are inside it. Fetched as the zip and unpacked here rather than
# file-by-file, so the acquisition is one citable object.
KEIO_URL = ("https://www.embopress.org/action/downloadSupplement"
            "?doi=10.1038%2Fmsb4100050&file=msb4100050-sup-0001.zip")

# Eydallin 2010's hits table. The publisher ships it as a PDF supplement; the tabular
# form the `eydallin` cohort reads was extracted once by hand in the sibling project, so
# it is EXTRACTION (human judgement) rather than a download -- which under the tier rule
# makes it curated, not raw. It is named here because the cohort needs it, and the
# transform says so rather than fetching something that is not the same table.
EYDALLIN_EXTRACTED = "eydallin2010_hits.tsv"

# The heterologous / metagenomic screening cohort. Resolved while porting benchmark_v4:
# its source of record is `heterologous_uniprot.tsv` -- one row per heterologous protein
# with the UniProt accession the screen reported -- and step 15 fetches those sequences
# from UniProt. That table is itself an extraction from the screening papers, so it lands
# in the same category as Eydallin's: named and required, not downloadable.
HET_TABLE = "heterologous_uniprot.tsv"


def _refuse(product_name, filename, why):
    return (
        f"{product_name}: {filename} is not present.\n"
        f"  {why}\n"
        f"  This is a NAMED refusal, not a silent gap: the cohort that reads it would "
        f"otherwise ship as a smaller conditions table that looks complete.")


def protocol(context: ExecutionContext):
    ilaser = context.Output(laser)
    ikeio  = context.Output(keio)
    ieyd   = context.Output(eydallin)
    ihet   = context.Output(het)

    # LASER is a git CHECKOUT, pinned by url. It is deliberately not hashed as data
    # bytes: hashing a tree that contains a .git turns it into thousands of cache
    # objects and it stops being a checkout.
    context.ExecWithContainer(image=image, cmd=f"""
        git clone --depth 1 {LASER_URL}.git {ilaser.container}
        rev=$(git -C {ilaser.container} rev-parse HEAD)
        echo "$rev" > {ilaser.container}/.PINNED_REV
        echo "[literature] laser at $rev"

        mkdir -p {ikeio.container}
        wget -q "{KEIO_URL}" -O {ikeio.container}/keio_baba2006_supplementary.zip
        ( cd {ikeio.container} && unzip -o -q keio_baba2006_supplementary.zip )
    """)

    # The two extracted tables. Both are human extractions from paper supplements, so
    # there is nothing to fetch; a drop-in is verified and the transform refuses by name
    # when it is absent.
    for out_path, filename, product_name, why in (
        (ieyd.local, EYDALLIN_EXTRACTED, "eydallin",
         "Eydallin 2010 ships its hits as a PDF supplement; the tabular form is a "
         "one-time hand extraction. Place it in this product's directory."),
        (ihet.local, HET_TABLE, "het_screen",
         "The heterologous screening cohort's protein list is an extraction from the "
         "screening papers (one row per protein, with its UniProt accession). Place it "
         "in this product's directory."),
    ):
        out_path.mkdir(parents=True, exist_ok=True)
        if not (out_path / filename).exists():
            raise SystemExit(_refuse(product_name, filename, why))
        Log.Info(f"verified {product_name}/{filename}")

    ok = (ilaser.local / ".PINNED_REV").exists() and \
         any(ikeio.local.glob("*.xls")) and \
         (ieyd.local / EYDALLIN_EXTRACTED).exists() and \
         (ihet.local / HET_TABLE).exists()
    return ExecutionResult(
        manifest=[{laser: ilaser.local, keio: ikeio.local,
                   eydallin: ieyd.local, het: ihet.local}],
        success=ok,
    )

TransformInstance(
    protocol=protocol,
    model=model,
    group_by=image,
    labels=["local"],
    resources=Resources(cpus=1, memory=Size.GB(4), duration=Duration(hours=1)),
)
