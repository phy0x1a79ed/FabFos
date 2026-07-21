#!/bin/bash
#SBATCH --account=def-shallam_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=256G
#SBATCH --time=1:00:00
#SBATCH --array=0-19
#SBATCH --job-name=dir_null
#SBATCH --output=/scratch/phyberos/scadc_dir/logs/null-%A_%a.out
set -e
SC=/scratch/$USER/scadc_dir
module load StdEnv/2023 python/3.12
source $SC/venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
read X N < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $SC/manifest.txt)
echo "task $SLURM_ARRAY_TASK_ID -> $X N$N on $(hostname) nproc=$(nproc) $(date)"
python $SC/fir_null.py --stage-dir $SC --element $X --draw-size $N \
  --workers ${SLURM_CPUS_PER_TASK} --resume \
  --out-reff $SC/out/reff_${X}_N${N}.tsv --out-ieff $SC/out/ieff_${X}_N${N}.tsv
echo "done $X N$N $(date)"
