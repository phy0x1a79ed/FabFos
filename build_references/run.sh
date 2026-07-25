#!/usr/bin/env bash
# Rebuild every compiled reference table, end to end.
#
# Four envs, because they cannot co-reside: `fabfos` for the pandas/pyarrow steps,
# `scadc-metabolic-model` for anything that opens a curated model through cobra, and
# `dvc` for the pins. Each step prints its own banner so a failure names its step.
#
# Inputs are the DVC-pinned, immutable bundles under data/raw/ and the MetaNetX and
# Rhea trees under data/reference/. Outputs land beside them:
#
#   data/reference/metabolism/{vocab,atom_pairs,direction}.parquet
#   data/reference/functional_annotation/bridges/*
#   data/reference/hosts/<host>/gpr_{gem,denovo}.parquet
#
# Verify runs last and is the gate; it exits non-zero on any failure.
set -euo pipefail
cd "$(dirname "$0")/.."
REPO="$(pwd)"

FORCE="${1:-}"          # pass --force to overwrite existing tables

echo "== 0. inputs (dvc) =="
mamba run --no-capture-output -n dvc dvc checkout \
    data/raw/hosts/e_coli_dh10b.dvc \
    data/raw/hosts/e_coli_k12.dvc \
    data/raw/hosts/e_coli_epi300.dvc

echo
echo "== 1. metabolism bake (fabfos) =="
# frozen atom mapping + directionality -> three int-coded parquets. Self-tests on the
# way out: every one of the 2,455,235 rows round-trips, or this step fails.
mamba run --no-capture-output -n fabfos python build_references/bake_metabolism.py $FORCE

echo
echo "== 2. intermediate -> MNXR bridges (fabfos) =="
# ec/ko/uniprot. The KO step reads the scadc KEGG cache once and fetches the remainder
# from KEGG REST; after the first run ko_to_kegg_r.tsv is the only input.
mamba run --no-capture-output -n fabfos python build_references/build_bridges.py $FORCE

echo
echo "== 3. EPI300 identity (fabfos) =="
# re-measures the EPI300<->DH10B delta. Its passing is what makes step 4 correct to
# build EPI300's GEM table by tagging the DH10B model.
mamba run --no-capture-output -n fabfos python build_references/check_epi300_identity.py

echo
echo "== 4. GEM GPR tables (scadc-metabolic-model) =="
# cobra lives only in this env; it is needed for the annotation dictionaries the
# crosswalk reads.
mamba run --no-capture-output -n scadc-metabolic-model \
    python build_references/build_host_gem_gpr.py $FORCE

echo
echo "== 5. de-novo GPR tables (fabfos) =="
mamba run --no-capture-output -n fabfos python build_references/build_host_denovo_gpr.py $FORCE

echo
echo "== 6. verify (fabfos) =="
mamba run --no-capture-output -n fabfos python build_references/check_references.py

echo
echo "== 7. declare + pin =="
echo "  not run automatically -- staging rewrites the library and re-hashes DVC dirs."
echo "  when the tables above are the ones you want to keep:"
echo "    python transforms/build/_stage/build_ref_library.py --stage"
echo "    python transforms/build/_stage/build_ref_library.py --hash"
echo "    mamba run -n dvc dvc add data/reference/metabolism data/reference/hosts \\"
echo "                             data/reference/functional_annotation"
