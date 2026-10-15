#!/bin/bash
#SBATCH --job-name=sdoh_v2
#SBATCH --account=PCON0023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=12:00:00
#SBATCH --output=sdoh_v2_%j.log
#SBATCH --mail-type=END,FAIL

# Quick test:  sbatch --export=ALL,QUICK=TRUE run_sdoh_v2.slurm
# Full run:    sbatch run_sdoh_v2.slurm
cd "/users/PAS2982/summertea1109/Random Forest ML/NHANES_BMD"
module load gcc/12.3.0 R/4.5.2
export QUICK=${QUICK:-FALSE}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript sdoh_v2_prediction.R
