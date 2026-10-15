#!/bin/bash
#SBATCH --job-name=sdoh_v2_sens
#SBATCH --account=PCON0023
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=12:00:00
#SBATCH --array=0-5
#SBATCH --output=sdoh_v2_sens_%A_%a.log
#SBATCH --mail-type=END,FAIL
# Full run (all six scenarios in parallel):  sbatch run_sdoh_v2_sens.slurm
# Quick test first:                          sbatch --export=ALL,QUICK=TRUE run_sdoh_v2_sens.slurm
# One scenario only, e.g. hip (index 5):     sbatch --array=5 run_sdoh_v2_sens.slurm
SCENARIOS=(main_extra complete mnar cyc1318 nomed hip)
export SCENARIO=${SCENARIOS[$SLURM_ARRAY_TASK_ID]}
cd "/users/PAS2982/summertea1109/Random Forest ML/NHANES_BMD"
module load gcc/12.3.0 R/4.5.2
export QUICK=${QUICK:-FALSE}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript sdoh_v2_sensitivity.R
