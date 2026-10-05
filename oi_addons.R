# =============================================================================
#  Add-on analyses for the Osteoporosis International version
#   Part A (minutes): crude vs progressively adjusted survey-weighted
#          associations of education and income with femoral neck T-score,
#          by sex, pooled over the saved imputations  -> explains the null
#   Part B (~1.5-2 h): total hip BMD as a secondary outcome - repeats the
#          primary contrast (full vs base) with random forest and spline
#          regression, 10 x 10-fold CV, 5 imputations, with checkpoints
#
#  Needs (already on OSC from the main run):
#    results/imputations.rds   and   nhanes_cache/DXXFEM_*.xpt
#  Run in RStudio:  setwd(".../NHANES_BMD"); source("oi_addons.R")
#  Re-running resumes from saved folds.
# =============================================================================
suppressPackageStartupMessages({
  library(haven); library(dplyr); library(tidyr); library(purrr); library(readr)
  library(survey); library(mitools); library(mice); library(ranger); library(splines)
})
options(survey.lonely.psu = "adjust", dplyr.summarise.inform = FALSE)
PROJECT_DIR <- Sys.getenv("NHANES_BMD_DIR", getwd())  # folder containing nhanes_cache/; defaults to the working directory
if (dir.exists(PROJECT_DIR)) setwd(PROJECT_DIR)
out <- function(f) file.path("results", f)

detect_threads <- function() {
  for (v in c("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE", "SLURM_JOB_CPUS_PER_NODE")) {
    n <- suppressWarnings(as.integer(sub("\\(.*", "", Sys.getenv(v))))
    if (!is.na(n) && n > 0) return(n)
  }
  n <- tryCatch(length(parallel::mcaffinity()), error = function(e) NA_integer_)
  if (!is.na(n) && n > 0) return(n)
  max(1L, min(4L, parallel::detectCores() - 1L))
}
CFG <- list(seed = 20261001, k_folds = 10, n_repeats = 10, num_trees = 500, tune_trees = 300,
            grid = expand.grid(mtry = c(2, 4, 6), min.node.size = c(10, 50)),
            threads = detect_threads(),
            # NHANES III total femur reference, non-Hispanic White women 20-29 y
            # (Looker et al., Osteoporos Int 1998)                          # VERIFY
            hip_mean = 0.942, hip_sd = 0.122)
message("Cores used: ", CFG$threads)

imps <- readRDS(out("imputations.rds"))
BASE_COMMON <- c("age", "weight", "height", "race", "cycle", "smoking", "alcohol",
                 "phys_active", "thyroid", "ra", "calcium_supp", "vitd_supp",
                 "prior_fx", "glucocorticoid", "parent_hipfx", "osteo_med")
WOMEN_EXTRA <- c("postmeno", "ever_ht", "ever_preg")
keep <- c("SEQN", "tscore", "osteoporosis", "wt_pool", "SDMVPSU", "SDMVSTRA", "sex")
imputed_list <- function(sx) {
  lapply(seq_len(imps[[sx]]$imp$m), function(i)
    bind_cols(imps[[sx]]$data[, keep], complete(imps[[sx]]$imp, i)) %>%
      mutate(pir_cat = cut(pir, c(-Inf, 1.3, 3.5, Inf), right = FALSE,
                           labels = c("Low", "Middle", "High"))))
}

# ---- Part A: crude vs adjusted associations ---------------------------------
adj_sets <- list(
  "M0 crude (education + income only)"      = character(0),
  "M1 + age, race and ethnicity, cycle"      = c("age", "race", "cycle"),
  "M2 + weight, height"                      = c("age", "race", "cycle", "weight", "height"),
  "M3 + smoking, alcohol, physical activity" = c("age", "race", "cycle", "weight", "height",
                                                 "smoking", "alcohol", "phys_active"),
  "M4 fully adjusted (all base predictors)"  = "ALL")
partA <- list()
for (sx in c("Male", "Female")) {
  dl <- imputed_list(sx)
  des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~wt_pool, nest = TRUE,
                   data = imputationList(dl))
  for (mn in names(adj_sets)) {
    cov <- adj_sets[[mn]]
    if (identical(cov, "ALL")) cov <- c(BASE_COMMON, if (sx == "Female") WOMEN_EXTRA)
    f <- as.formula(paste("tscore ~ education + pir_cat", if (length(cov)) paste("+", paste(cov, collapse = " + ")) else ""))
    pooled <- MIcombine(with(des, svyglm(f)))
    b <- coef(pooled); V <- vcov(pooled)
    terms <- grep("^education|^pir_cat", names(b), value = TRUE)
    wald <- function(tt) { w <- as.numeric(t(b[tt]) %*% solve(V[tt, tt]) %*% b[tt]); pchisq(w, length(tt), lower.tail = FALSE) }
    partA[[length(partA) + 1]] <- tibble(sex = sx, model = mn, term = terms, est = b[terms],
      se = sqrt(diag(V))[terms]) %>%
      mutate(lo = est - 1.96 * se, hi = est + 1.96 * se, p = 2 * pnorm(-abs(est / se)),
             p_joint_education = wald(grep("^education", terms, value = TRUE)),
             p_joint_income = wald(grep("^pir_cat", terms, value = TRUE)))
  }
  message("Part A done: ", sx)
}
write_csv(bind_rows(partA), out("table_crude_vs_adjusted.csv"))

# ---- Part B: total hip BMD --------------------------------------------------
hip <- map_dfr(c("D", "E", "F", "H", "J"), function(s)
  read_xpt(file.path("nhanes_cache", paste0("DXXFEM_", s, ".xpt"))) %>% select(SEQN, DXXOFBMD)) %>%
  zap_labels() %>% filter(!is.na(DXXOFBMD)) %>%
  mutate(hip_t = (DXXOFBMD - CFG$hip_mean) / CFG$hip_sd, hip_op = as.integer(hip_t <= -2.5))

wmetrics <- function(y, yhat, w) { w <- w / sum(w); e <- y - yhat; ybar <- sum(w * y); mse <- sum(w * e^2)
  c(MSE = mse, RMSE = sqrt(mse), MAE = sum(w * abs(e)), R2 = 1 - mse / sum(w * (y - ybar)^2)) }
wauc <- function(risk, case, w) {
  ok <- !is.na(risk) & !is.na(case); risk <- risk[ok]; case <- case[ok]; w <- w[ok]
  if (length(unique(case)) < 2) return(NA_real_)
  u <- sort(unique(risk)); g <- factor(match(risk, u), levels = seq_along(u))
  wn <- tapply(w * (case == 0), g, sum); wn[is.na(wn)] <- 0
  wp <- tapply(w * (case == 1), g, sum); wp[is.na(wp)] <- 0
  sum(wp * (cumsum(wn) - wn + 0.5 * wn)) / (sum(wp) * sum(wn)) }
nb_ci <- function(x, ntr, nts) { x <- x[!is.na(x)]; J <- length(x); m <- mean(x)
  se <- sqrt((1 / J + nts / ntr) * var(x)); tq <- qt(0.975, J - 1)
  c(est = m, lo = m - tq * se, hi = m + tq * se, p = 2 * pt(-abs(m / se), J - 1)) }
tune_rf <- function(tr, preds, w, seed) {
  f <- reformulate(preds, "hip_t")
  errs <- vapply(seq_len(nrow(CFG$grid)), function(i) ranger(f, data = tr, num.trees = CFG$tune_trees,
    mtry = min(CFG$grid$mtry[i], length(preds)), min.node.size = CFG$grid$min.node.size[i],
    case.weights = w, respect.unordered.factors = "order", num.threads = CFG$threads,
    seed = seed)$prediction.error, numeric(1))
  as.list(CFG$grid[which.min(errs), ]) }
fit_rf <- function(tr, preds, w, seed, p) ranger(reformulate(preds, "hip_t"), data = tr,
  num.trees = CFG$num_trees, mtry = min(p$mtry, length(preds)), min.node.size = p$min.node.size,
  case.weights = w, respect.unordered.factors = "order", num.threads = CFG$threads, seed = seed)
lm_formula <- function(preds) as.formula(paste("hip_t ~", paste(vapply(preds, function(p)
  if (p %in% c("age", "weight", "height", "pir")) sprintf("ns(%s, df = 3)", p) else p, ""), collapse = " + ")))

perf <- list()
for (sx in c("Male", "Female")) {
  dl <- lapply(imputed_list(sx), function(x) inner_join(x, hip, by = "SEQN"))
  base <- c(BASE_COMMON, if (sx == "Female") WOMEN_EXTRA); full <- c(base, "education", "pir")
  set.seed(CFG$seed + sum(utf8ToInt(paste("hip", sx))))
  n <- nrow(dl[[1]]); folds <- list()
  for (r in seq_len(CFG$n_repeats)) { id <- sample(rep(seq_len(CFG$k_folds), length.out = n))
    for (k in seq_len(CFG$k_folds)) folds[[length(folds) + 1]] <- list(rep = r, fold = k, test = which(id == k)) }
  ck <- out(file.path("checkpoints", paste0("hip_", sx))); dir.create(ck, recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time(); nnew <- 0
  for (fi in seq_along(folds)) {
    f_ck <- file.path(ck, sprintf("fold_%03d.rds", fi))
    if (file.exists(f_ck)) { perf[[length(perf) + 1]] <- readRDS(f_ck); next }
    te <- folds[[fi]]$test; tr1 <- dl[[1]][-te, ]
    params <- tune_rf(tr1, full, tr1$wt_pool / mean(tr1$wt_pool), CFG$seed + fi)
    res <- list()
    for (mi in seq_along(dl)) {
      tr <- dl[[mi]][-te, ]; ts <- dl[[mi]][te, ]; w_tr <- tr$wt_pool / mean(tr$wt_pool); tr$.w <- w_tr
      for (mn in c("base", "full")) {
        preds <- if (mn == "base") base else full
        p <- predict(fit_rf(tr, preds, w_tr, CFG$seed + fi, params), ts, num.threads = CFG$threads)$predictions
        pl <- tryCatch(predict(lm(lm_formula(preds), data = tr, weights = .w), ts), error = function(e) rep(NA, nrow(ts)))
        for (L in list(list("Random forest", p), list("Linear (splines)", pl)))
          res[[length(res) + 1]] <- tibble(sex = sx, rep = folds[[fi]]$rep, fold = folds[[fi]]$fold, imp = mi,
            learner = L[[1]], model = mn, !!!as.list(wmetrics(ts$hip_t, L[[2]], ts$wt_pool)),
            AUC_osteoporosis = wauc(-L[[2]], ts$hip_op, ts$wt_pool), n_train = nrow(tr), n_test = nrow(ts))
      }
    }
    r <- bind_rows(res); saveRDS(r, f_ck); perf[[length(perf) + 1]] <- r; nnew <- nnew + 1
    mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    message(sprintf("  hip | %s: fold %d / %d done (%.1f min, about %.0f min left)", sx, fi,
                    length(folds), mins, (length(folds) - fi) * mins / nnew))
  }
}
perf <- bind_rows(perf); write_csv(perf, out("cv_raw_performance_hip.csv"))
pf <- perf %>% group_by(sex, learner, model, rep, fold) %>%
  summarise(across(c(R2, RMSE, AUC_osteoporosis, n_train, n_test), mean), .groups = "drop")
summ <- pf %>% pivot_longer(c(R2, RMSE, AUC_osteoporosis), names_to = "metric") %>%
  select(-n_train, -n_test) %>% left_join(pf %>% select(sex, learner, model, rep, fold, n_train, n_test),
                                          by = c("sex", "learner", "model", "rep", "fold")) %>%
  pivot_wider(names_from = model, values_from = value) %>%
  group_by(sex, learner, metric) %>%
  summarise(ci = list(nb_ci(full - base, mean(n_train), mean(n_test))),  # paired differences first
            base = mean(base, na.rm = TRUE), full_mean = mean(full, na.rm = TRUE), .groups = "drop") %>%
  mutate(delta = map_dbl(ci, ~ .x[["est"]]), lo = map_dbl(ci, ~ .x[["lo"]]),
         hi = map_dbl(ci, ~ .x[["hi"]]), p = map_dbl(ci, ~ .x[["p"]])) %>% select(-ci)
write_csv(summ, out("table_hip_deltas.csv"))
message("All add-ons done. New files: table_crude_vs_adjusted.csv, table_hip_deltas.csv, cv_raw_performance_hip.csv")
