#!/bin/bash
#SBATCH --account=st-shallam-1
#SBATCH --job-name=dir_null
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
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
#                than failing, which is the worst way for this to go wrong. 180 G
#                fits a standard node. If the measured high-water mark exceeds
#                it, move to `--partition=big_mem` (385 GB) and raise this --
#                do NOT quietly drop draws to fit.
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

# Sockeye's Lmod needs these as SEPARATE load commands. Combined into one
# `module load gcc/9.4.0 python/3.12`, the second resolves against the pre-gcc
# module tree, silently finds nothing, and you get the system python.
module load gcc/9.4.0
module load python/3.12

source "$SC/venv/bin/activate"

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

python "$SC/null_shard.py" \
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
