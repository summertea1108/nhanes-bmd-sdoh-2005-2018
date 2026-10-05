# Education, income and bone mineral density in US adults aged ≥50 (NHANES 2005–2018)

Analysis code for:

> Liu J, Jung J, Li Q, Wu L, Wu Q. *Do education and income improve identification of low bone mineral density beyond clinical risk factors? A cross-sectional analysis of 10,981 US adults aged 50 years and older (NHANES 2005–2018).* Submitted to *Osteoporosis International*.

## What the code does

| Script | Purpose |
|---|---|
| `install_packages.R` | Installs the required R packages. |
| `bmd_sdoh_reanalysis.R` | Main analysis: downloads and harmonizes NHANES 2005–2010, 2013–2014 and 2017–2018; builds the analytic sample (adults ≥50 years with femoral neck DXA); multiple imputation by sex; 10 × 10-fold cross-validation comparing survey-weighted random forest and spline regression with and without education and income; permutation importance; accumulated local effects; sensitivity analyses; survey-weighted regression with sex interactions. |
| `oi_addons.R` | Additional analyses: crude and progressively adjusted associations of education and income with T-score (Table 4), and total hip BMD as a secondary outcome. Requires `results/imputations.rds` from the main script. |
| `oi_revision.R` | Revision analyses: association analyses re-run with a separate imputation that includes T-score (Table 4 and the sex-interaction model), the medication-free sensitivity analysis with both methods, and the comparison of included and excluded participants. Requires `results/imputations.rds` from the main script. |

## Data

All data are public NHANES files from the National Center for Health Statistics (https://wwwn.cdc.gov/nchs/nhanes/). The main script downloads the required `.xpt` files into `nhanes_cache/` if they are not already present. No data are included in this repository.

## How to run

```r
# R >= 4.4
setwd("path/to/this/folder")
source("install_packages.R")
source("bmd_sdoh_reanalysis.R")   # several hours; progress is checkpointed in results/checkpoints/
source("oi_addons.R")             # about 2 hours
source("oi_revision.R")           # about 1-2 hours
```

Set `Sys.setenv(QUICK = "1")` before running for a fast test with fewer cross-validation repeats. Results (tables and figures as CSV/PNG) are written to `results/`. The run resumes from checkpoints if interrupted.

Software used: R 4.4.0 with ranger, survey, mice, mitools, iml, splines, dplyr and related packages (see `install_packages.R`).

## License

MIT (see `LICENSE`).
