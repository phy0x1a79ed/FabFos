#!/usr/bin/env bash
# Regenerate the per-reaction direction annotator end to end.
#
# Three members spanning method families, fused into one directional conductance
# ratio per base-graph reaction (no evidence -> ratio 1.0 -> reversible). Runs
# across three envs (they cannot co-reside): p312 for pure-pandas steps,
# `equilibrator` (py3.11) for eQuilibrator, `dgbyg` (py3.13) for the GNN.
#
# Paths mirror canon.py (canon.DIR_*). Outputs land in data/scadc/direction/.
set -euo pipefail
cd "$(dirname "$0")"
D="$(pwd)"

DATA=/home/tony/agentic_workspace/data/scadc
REF=$DATA/references/metanetx
STAGING=/home/tony/agentic_workspace/projects/self-improvement/main/staging
DEST=$DATA/direction
mkdir -p "$DEST"

# dGbyG is vendored (bundled weights, gitignored). Provision it once.
VENDOR="$D/../../../lib/vendor/dGbyG"
if [ ! -d "$VENDOR/.git" ]; then
  echo "== cloning dGbyG (~369MB, weights bundled) =="
  git clone --depth 1 https://github.com/f-wc/dGbyG.git "$VENDOR"
  mamba run -n dgbyg pip install -e "$VENDOR"
fi

# dgbyg's env rdkit needs the env's newer libstdc++, and the module dir on PYTHONPATH.
dgbyg() { mamba run -n dgbyg bash -c "cd $D && LD_LIBRARY_PATH=\$CONDA_PREFIX/lib PYTHONPATH=$D python $*"; }

echo "== base-graph MNXR universe (union over C,N,S,P) =="
mamba run -n p312 python - "$DEST/base_graph_mnxrs.json" <<'PY'
import sys, json, pickle
sys.path.insert(0, "..")            # main/fabfos on path
import canon
mnxrs = set()
for f in canon.BASE_GRAPH_FILES:
    with open(canon.SOLVE_BASE_DIR / f, "rb") as fh:
        G = pickle.load(fh)
    mnxrs |= {n[1] for n in G.nodes if n[0] == "rxn"}
json.dump(sorted(mnxrs), open(sys.argv[1], "w"))
print(f"  {len(mnxrs)} reactions")
PY

echo "== member 3: curated BioCyc direction, orientation-aligned (p312) =="
mamba run -n p312 python curated.py \
  --metacyc-pgdb "$STAGING/metacyc26.pgdb" --ecocyc-pgdb "$STAGING/ecocyc26.pgdb" \
  --reac-xref "$REF/reac_xref.tsv" --reac-prop "$REF/reac_prop.tsv" --chem-xref "$REF/chem_xref.tsv" \
  --out "$DEST/curated_per_mnxr.parquet" --out-per-reaction "$DEST/curated_per_reaction.parquet"

echo "== T3 calibration: category -> dG' on the eQ measured arm (equilibrator) =="
mamba run -n equilibrator python calibrate.py \
  --curated "$DEST/curated_per_mnxr.parquet" \
  --reac-prop "$REF/reac_prop.tsv" --chem-prop "$REF/chem_prop.tsv" \
  --out-calibration "$DEST/calibration.parquet" --out-points "$DEST/calib_points.parquet"
# NB: SIGMA_0 = robust marginal spread of the measured arm; committed in canon.DIR_SIGMA_0.
# If the calibration changes, recompute and re-commit it (must stay in DIR_SIGMA_0_BAND).

echo "== member 1: eQuilibrator over the base graph (equilibrator) =="
mamba run -n equilibrator env PYTHONPATH="$D" python eval_members.py --member eq \
  --mnxrs "$DEST/base_graph_mnxrs.json" --reac-prop "$REF/reac_prop.tsv" \
  --chem-prop "$REF/chem_prop.tsv" --out "$DEST/member_eq_base.parquet"

echo "== member 2: dGbyG over the base graph (dgbyg) =="
dgbyg "eval_members.py --member dgbyg --mnxrs $DEST/base_graph_mnxrs.json \
  --reac-prop $REF/reac_prop.tsv --chem-prop $REF/chem_prop.tsv --out $DEST/member_dgbyg_base.parquet"

echo "== T4 combine into the annotator table (p312) =="
mamba run -n p312 python combine.py \
  --base-mnxrs "$DEST/base_graph_mnxrs.json" \
  --eq "$DEST/member_eq_base.parquet" --dgbyg "$DEST/member_dgbyg_base.parquet" \
  --curated "$DEST/curated_per_mnxr.parquet" --calibration "$DEST/calibration.parquet" \
  --out "$DEST/direction_annotation.parquet"

echo "== verify =="
mamba run -n p312 python selftest.py \
  --table "$DEST/direction_annotation.parquet" --base-mnxrs "$DEST/base_graph_mnxrs.json"
