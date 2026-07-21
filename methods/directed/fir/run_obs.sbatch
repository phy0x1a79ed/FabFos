#!/bin/bash
#SBATCH --account=def-shallam_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem=256G
#SBATCH --time=0:30:00
#SBATCH --job-name=dir_obs
#SBATCH --output=/scratch/phyberos/scadc_dir/logs/obs-%j.out
set -e
SC=/scratch/$USER/scadc_dir
module load StdEnv/2023 python/3.12
source $SC/venv/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p $SC/obs
echo "observed directed solve on $(hostname) nproc=$(nproc) $(date)"
python $SC/fir_observed.py --stage-dir $SC --workers ${SLURM_CPUS_PER_TASK} \
  --out-reff $SC/obs/reff_axes_report.tsv --out-ieff $SC/obs/ieff_axes_report.tsv
echo "observed done $(date)"; wc -l $SC/obs/reff_axes_report.tsv
