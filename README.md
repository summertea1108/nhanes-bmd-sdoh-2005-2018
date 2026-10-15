# Race, social determinants and bone density prediction in US adults (NHANES 2005–2018)

Analysis code for:

> Liu J, Jung J, Li Q, Wu L, Wu Q. *Race, social determinants and bone density prediction in US adults: implications for race-neutral risk tools.* Submitted to the *Journal of Bone and Mineral Research*.

All data are public NHANES files from the National Center for Health Statistics (https://wwwn.cdc.gov/nchs/nhanes/). The scripts download the required `.xpt` files into `nhanes_cache/` if they are not already there. No data are included in this repository.

## Scripts and what they produce

Run the scripts in this order. Each later script reads `results/imputations.rds`, which the first script creates.

| Step | Script | Produces (manuscript items) |
|---|---|---|
| 0 | `install_packages.R` | Installs the required R packages. |
| 1 | `bmd_sdoh_reanalysis.R` | Downloads and harmonizes NHANES 2005–2010, 2013–2014 and 2017–2018; builds the analytic sample (adults ≥50 years with a valid femoral neck DXA scan, n = 10,981); multiple imputation by sex (`results/imputations.rds`); descriptive statistics (Table 1); sample selection (Fig. S1); PIR missingness by education (Supplementary File 3). Variable definitions are listed in Supplementary File 1. The cross-validation outputs of this script are not used in the paper. |
| 2 | `oi_revision.R` | Associations of education and income with femoral neck T-score, using 20 imputations that include T-score (Table 3); model with sex interactions and joint Wald test (Supplementary File 2); included vs excluded participants (Table S2). |
| 3 | `sdoh_bundle.R` | Eight-item social determinants of health (SDoH) bundle: harmonization check and item prevalence (Table 2, Table S4), item correlations (Table S5), associations of the disadvantage score with T-score and osteoporosis (Table 4, Fig. 1), sex and race-and-ethnicity interactions (Table S9), and race-and-ethnicity coefficients with and without SDoH and body size (Table 6). Output files `00_` to `13_` are used in the paper; its cross-validation outputs (`20_` to `23_`) are not. |
| 4 | `sdoh_v2_prediction.R` | Main prediction analysis: four strategies (A clinical; B + SDoH block; C + race and ethnicity; D + both), secondary strategies and the OST benchmark; imputation inside each training fold using only the variables each strategy may use; 10 × 10-fold cross-validation; Rao–Wu survey-design bootstrap (500 replicates). Produces Table 5, Fig. 2 and Tables S6–S8. Slurm job: `run_sdoh_v2.slurm`. |
| 5 | `sdoh_v2_sensitivity.R` | Same procedure as step 4 for the secondary and sensitivity analyses: education and income added, unemployment added, complete cases, imputed PIR shifted by −0.5 and −1.0 (missing not at random), 2013–2018 cycles only, osteoporosis-medication users excluded, and total hip T-score. Produces Supplementary File 4, Table S1 and Table S10. Slurm array job (six scenarios in parallel): `run_sdoh_v2_sens.slurm`. |
| — | `figures/make_fig1.py`, `figures/make_fig2.py` | Draw Figs. 1 and 2 (PDF, EPS and 600-dpi TIFF) from the CSV outputs of steps 3 and 4. |

`oi_addons.R` is kept from an earlier version of this repository; its results are not used in the paper.

## How to run

```r
# R 4.5.2 (packages: survey 4.5, mice 3.19.0, mitools 2.7, ranger 0.18.0, dplyr, haven and related)
setwd("path/to/this/folder")
source("install_packages.R")
source("bmd_sdoh_reanalysis.R")   # several hours; checkpointed in results/checkpoints/
source("oi_revision.R")
source("sdoh_bundle.R")
```

Steps 4 and 5 were run on the Ohio Supercomputer Center (Slurm, 28 cores per job):

```bash
sbatch run_sdoh_v2.slurm                                   # main prediction analysis (~75 min)
sbatch run_sdoh_v2_sens.slurm                              # six sensitivity scenarios (~1.5–2 h each, in parallel)
sbatch --export=ALL,QUICK=TRUE run_sdoh_v2_sens.slurm      # quick test first (a few minutes)
```

Edit the `cd` path and `--account` in the `.slurm` files for your own system, or set `NHANES_BMD_DIR` to the project folder. Every script has a quick-test mode (`QUICK=TRUE`) and resumes from checkpoints if interrupted. Fold seeds are fixed, so strategies A–D in the full-sample sensitivity scenario reproduce the main results of step 4 exactly.

Figures (Python 3 with pandas and matplotlib):

```bash
python figures/make_fig1.py results/sdoh_bundle/10_aim1_score_associations.csv fig1.tiff
python figures/make_fig2.py results/sdoh_v2/32_v2_paired_contrasts.csv fig2.pdf
```

## License

MIT (see `LICENSE`).
