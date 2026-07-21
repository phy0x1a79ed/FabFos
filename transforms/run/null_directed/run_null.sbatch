#!/bin/bash
#SBATCH --account=st-shallam-1
#SBATCH --job-name=dir_null
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --array=0-19
#SBATCH --output=/scratch/st-shallam-1/txyliu/scadc_dir/logs/null-%A_%a.out
#
# The directed null, one (element, draw-size) shard per array task, on SOCKEYE.
#
# Ported from fir. Every #SBATCH line above differs from the fir original, and
# the reasons are worth stating because three of them are traps:
#
#   account      fir used the Alliance-style `def-shallam_cpu`; sockeye
#                allocations are `st-<pi>-1`.
#   cpus/mem     fir asked 192 cores / 256 GB -- a fir full-node shape. Sockeye's
#                skylake and cascade nodes are 32 and 40 cores at 190 GB, so a
#                256 GB request NEVER SCHEDULES here; it queues forever rather
#                than failing, which is the worst way for this to go wrong.
#
#                96 G is MEASURED, not guessed. Smoke on the worst cell (C/N43,
#                the largest LCC at n=12,884) at 100 and 300 draws gave a peak
#                RSS of 548 and 558 MB -- a slope of ~0.05 MB/draw, so the full
#                4,000-draw `additions` dict adds only ~195 MB and the parent
#                peaks near 0.75 GB. Workers fork (CoW), so a pessimistic
#                1.5 GB/worker x 32 is ~48 GB. 96 G is 2x that.
#
#                Take the two-point measurement seriously if you re-tune this:
#                `additions` is built for ALL draws before the pool starts, so
#                it is the dominant term AND it is exactly what --smoke-draws
#                truncates. A single small smoke underestimates the real
#                high-water mark ~20x. Measure a slope, not a point.
#
#                If a high-water mark ever does exceed this, move to
#                `--partition=big_mem` (385 GB) and raise it -- do NOT quietly
#                drop draws to fit.
#   time         fir capped at 1 h because fir capped at 1 h. Sockeye allows 7
#                days on every batch partition, so walltime is essentially free
#                and there is no reason to run this tight. The coherent base
#                graph is ~1.9x fir's, and per-cell Newton cost grows with the
#                active set, so a 1 h limit would have died at the wall.
#   partition    fir auto-routed full-node core counts; sockeye does not, but the
#                default (skylake) is correct for this shape.
#
# The array range MUST equal the manifest line count minus one. stage_null.py
# generates manifest.txt from canon.DRAW_SIZES, so 4 elements x 5 sizes = 20
# lines = --array=0-19. If you narrow the grid, narrow this too: a task that
# reads past the end of the manifest gets an empty X/N and fails obscurely.
set -e

SC=/scratch/st-shallam-1/$USER/scadc_dir
mkdir -p "$SC/logs" "$SC/out"

# Run inside the SAME ecspr image the benchmark solve used, rather than a pip
# venv. This is not a style preference: the whole sockeye path for this project
# is apptainer-based (examples/benchmark_v3_on_hpc.py pre-pulls on the login
# node because compute nodes have no outbound network), and building a venv
# here would have introduced a second, unpinned dependency surface for the one
# step whose output feeds every significance call.
#
# Sockeye's Lmod needs these as SEPARATE load commands. Combined into one
# `module load gcc/9.4.0 apptainer/1.3.1`, the second resolves against the
# pre-gcc module tree, silently finds nothing, and you get
# `apptainer: command not found` with no hint a module was skipped.
module load gcc/9.4.0
module load apptainer/1.3.1

SIF=/scratch/st-shallam-1/$USER/apptainer_cache/docker..quay.io_hallamlab_external_ecspr..2026.07.14.sif
if [[ ! -f "$SIF" ]]; then
    echo "FATAL: ecspr image missing at $SIF -- pre-pull it on the LOGIN node" >&2
    exit 1
fi

# The image carries networkx 3.6.1 while the base/bipartite pickles were written
# locally under networkx 3.3, and mnxref.graph's declaration warns these graphs
# are version-coupled to networkx. That coupling was CHECKED here, not assumed:
# base_C.pkl and mnx_bipartite_C.pkl both unpickle in this image and report the
# expected 13,524 nodes / 362,482 skeleton edges. Re-check it if the image moves.

# One BLAS thread per process: the parallelism here is the fork pool over
# (axis x draw-chunk), and letting BLAS also thread oversubscribes the node
# badly enough to run slower than serial.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

read -r X N < <(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SC/manifest.txt")
if [[ -z "$X" || -z "$N" ]]; then
    echo "FATAL: manifest.txt has no line $((SLURM_ARRAY_TASK_ID + 1))" >&2
    exit 1
fi
echo "shard $SLURM_ARRAY_TASK_ID -> element=$X draw-size=$N on $(hostname)"

apptainer exec -B /scratch:/scratch "$SIF" python "$SC/null_shard.py" \
    --stage-dir "$SC" \
    --element "$X" \
    --draw-size "$N" \
    --workers "${SLURM_CPUS_PER_TASK}" \
    --out-reff "$SC/out/reff_${X}_N${N}.tsv" \
    --out-ieff "$SC/out/ieff_${X}_N${N}.tsv"

# Deliberately NO --resume. Its check is expected-vs-actual ROW COUNT, and the
# row count barely moves under a universe change (the testable set went 78 -> 79
# here, +1.3%). So on a re-run against new base graphs it would judge every stale
# shard "complete" and skip the entire sweep, leaving the previous universe's
# output sitting in out/ looking finished. Re-running a shard costs compute;
# silently keeping the wrong null costs the result.
echo "shard $SLURM_ARRAY_TASK_ID done: $(wc -l < "$SC/out/ieff_${X}_N${N}.tsv") lines"
