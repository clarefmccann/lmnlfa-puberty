## lmnlfa_growth_staged.R
## Staged longitudinal MNLFA: extends the validated growth + informant-DIF
## model (lmnlfa_growth_informant.R) to also test race/ethnicity and
## waist-to-height ratio (WtHR) as growth-factor IMPACT predictors, and
## age/race/WtHR/informant as DIF ("background variable") candidates --
## the extension deferred from that file's original header comment.
##
## WHY STAGED (not one saturated model): the cross-sectional analog of this
## exact setup (age/race/WtHR/informant as both impact and DIF candidates,
## estimated jointly) produced Rhat 1.37-1.54 and ESS 7-9/4000 despite zero
## divergences -- a between-chain multimodality problem, because a
## covariate's effect on item responses is only weakly separable between
## "real impact on the trait" and "DIF" when both are free at once. The fix
## there (mnlfa_crosssectional_staged.R) was 3 stages, replicated here:
##   A. Impact only: race + WtHR (baseline) on the growth intercept and
##      slope means. No DIF at all (not even informant) -- matches
##      mnlfa-crosssectional-impact.stan's precedent of a clean, DIF-free
##      impact stage.
##      -> scripts/stan/lmnlfa-growth-impact.stan
##   B. DIF screening: impact FIXED at Stage A's posterior means. Every
##      item x DIF-covariate (age, race, WtHR, informant) loading/intercept
##      DIF term estimated freely.
##      -> scripts/stan/lmnlfa-growth-difscreen.stan
##   C. Final: impact free again, DIF restricted to the item x covariate
##      cells retained by Stage B's screening (BH-FDR + magnitude floor,
##      same rule as the cross-sectional model).
##      -> scripts/stan/lmnlfa-growth-final.stan
##
## AGE IS NOT AN IMPACT COVARIATE HERE, unlike the cross-sectional model:
## age is already this model's growth axis (mu_slp *is* the age effect),
## so "age impact on the growth factor" isn't a separable concept the way
## it is cross-sectionally. Age-varying DIF -- item parameters shifting
## with age independent of the latent trajectory -- is the correct
## longitudinal analog, and IS tested (in kdif, alongside race/WtHR/
## informant).
##
## WtHR is dual-role, like age: the growth-factor IMPACT covariate uses
## each person's BASELINE (earliest available occasion) WtHR -- a
## between-person question about adiposity predicting the trajectory
## overall -- while the DIF covariate uses each OCCASION's own measured
## WtHR -- a within-occasion question about whether that moment's body
## composition affects how that moment's item response is being
## interpreted. Race is time-invariant, so it's the same value in both
## roles.
##
## COST NOTE: Stage B here has up to p x kdif x 2 = 4 x 6 x 2 = 48 DIF
## parameters, vs. only 8 in the informant-only model (l_dif_informant +
## n_dif_informant) -- a real increase in per-iteration cost on top of the
## already-expensive per-person growth-factor parameters. Use n_subsample
## below for a first diagnostic pass at each stage before committing to a
## full-sample HPC run, same caution that applied to
## lmnlfa_growth_informant.R.
##
## IMPORTANT: this script does NOT run any of the three fits automatically
## -- each stage's fit call is gated behind `run_stage_X <- TRUE/FALSE`
## flags below. Set them and run yourself (locally or on HPC).
##
## Usage:
##   Rscript lmnlfa_growth_staged.R <sex>
##   sex: female | male
##
## Outputs: written to <OUT_DIR>/lmnlfa_growth_staged/

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(posterior)
})

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop(
    "cmdstanr not found. Run: install.packages('cmdstanr', repos=c('https://mc-stan.org/r-packages/', getOption('repos')))"
  )
}
library(cmdstanr)

cmdstan_env <- Sys.getenv("CMDSTAN")
if (nzchar(cmdstan_env) && dir.exists(cmdstan_env)) {
  set_cmdstan_path(cmdstan_env)
}
cat("CmdStan path:", cmdstan_path(), "\n")
options(mc.cores = as.integer(Sys.getenv("NSLOTS", unset = "4")))
set.seed(90025)

# ---------------------------------------------------------------------------
# STAGE SWITCHES + SUBSAMPLE -- review before running
# ---------------------------------------------------------------------------
run_stage_A <- FALSE
run_stage_B <- FALSE
run_stage_C <- TRUE
n_subsample <- 1500 # e.g. 1500 for a diagnostic pass; NULL = full sample

# ---------------------------------------------------------------------------
# ARGUMENTS
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript lmnlfa_growth_staged.R <female|male>")
}
sx <- args[1]
if (!sx %in% c("female", "male")) {
  stop("sex must be 'female' or 'male'")
}
cat("Sex:", sx, "\n")

run_tag <- paste0(
  sx,
  if (is.null(n_subsample)) "" else paste0("_n", n_subsample)
)
run_label <- paste0(
  sx,
  if (is.null(n_subsample)) "" else paste0(" (n=", n_subsample, " subsample)")
)

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
root_path <- Sys.getenv("HOME_DIR")
if (!nzchar(root_path)) {
  root_path <- Sys.getenv("HOME")
}

data_dir <- Sys.getenv("DATA_DIR")
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  data_dir <- file.path(
    root_path,
    "projects/abcd-projs/dissertation/study1/outputs"
  )
}
if (!dir.exists(data_dir)) {
  stop("Cannot locate data directory: ", data_dir)
}

out_base <- Sys.getenv("OUT_DIR")
if (!nzchar(out_base)) {
  out_base <- file.path(
    root_path,
    "projects/abcd-projs/dissertation/study1/outputs"
  )
}
out_dir <- file.path(out_base, "lmnlfa_growth_staged")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# outputs are split into subfolders by type, not dumped flat into out_dir
fits_dir <- file.path(out_dir, "fits")
dif_selection_dir <- file.path(out_dir, "dif-selection")
pubertal_estimates_dir <- file.path(out_dir, "pubertal-estimates")
trajectories_dir <- file.path(out_dir, "trajectories")
impact_comparisons_dir <- file.path(out_dir, "impact-comparisons")
dif_illustrations_dir <- file.path(out_dir, "dif-illustrations")
for (d in c(
  fits_dir, dif_selection_dir, pubertal_estimates_dir,
  trajectories_dir, impact_comparisons_dir, dif_illustrations_dir
)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

script_dir <- Sys.getenv("SGE_O_WORKDIR")
if (!nzchar(script_dir)) {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
  script_dir <- if (length(file_arg) > 0) {
    dirname(normalizePath(file_arg[1]))
  } else {
    "scripts"
  }
}

palette_file <- file.path(script_dir, "color_palette.R")
if (!file.exists(palette_file)) {
  palette_file <- file.path("scripts", "color_palette.R")
}
source(palette_file)

find_stan <- function(name) {
  f <- file.path(script_dir, "stan", name)
  if (!file.exists(f)) {
    f <- file.path("scripts", "stan", name)
  }
  if (!file.exists(f)) {
    stop("Cannot find ", name, ": ", f)
  }
  f
}
stan_file_A <- find_stan("lmnlfa-growth-impact.stan")
stan_file_B <- find_stan("lmnlfa-growth-difscreen.stan")
stan_file_C <- find_stan("lmnlfa-growth-final.stan")

cat("Data dir:  ", data_dir, "\n")
cat("Output dir:", out_dir, "\n")

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------
parent_df <- read.csv(file.path(data_dir, paste0(sx, "_parent_long.csv")))
youth_df <- read.csv(file.path(data_dir, paste0(sx, "_youth_long.csv")))

if (!"whtr" %in% names(parent_df)) {
  stop(
    "whtr column not found -- re-run 00_data_foundation.R first ",
    "(it now pulls waist_in and derives whtr = waist_in / height_in)."
  )
}

wave_order <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")

# ---------------------------------------------------------------------------
# BUILD LONGITUDINAL DATA: 4 shared items, informant/age/race/WtHR as DIF
# candidates, race + baseline WtHR as growth-factor impact predictors
# ---------------------------------------------------------------------------
build_lmnlfa_data_staged <- function(
  parent_df,
  youth_df,
  sex_label,
  ordinal_items = c("peta", "petb", "petc", "petd"),
  n_subsample = NULL
) {
  cat(
    "\n=== Building longitudinal data (staged impact + DIF) |",
    sex_label,
    "===\n"
  )

  whtr_min <- 0.25
  whtr_max <- 0.75

  clean_reporter <- function(df, reporter_label, informant_val) {
    df %>%
      select(id, wave, age, race, whtr, all_of(ordinal_items)) %>%
      filter(
        !is.na(age),
        !is.na(race),
        !is.na(whtr),
        whtr >= whtr_min,
        whtr <= whtr_max
      ) %>%
      filter(if_all(
        all_of(ordinal_items),
        ~ !is.na(.) & as.integer(.) %in% 1:4
      )) %>%
      mutate(
        wave = factor(wave, levels = wave_order),
        reporter = reporter_label,
        informant_c = informant_val
      )
  }

  dat <- bind_rows(
    clean_reporter(parent_df, "parent", 1),
    clean_reporter(youth_df, "youth", -1)
  ) %>%
    arrange(id, wave, reporter)

  if (nrow(dat) < 500) {
    stop("Insufficient data for ", sex_label)
  }

  # optional person-level subsample -- see lmnlfa_growth_informant.R for
  # rationale (cuts both dominant cost drivers, per-person parameter count
  # and observation count, roughly proportionally)
  all_ids <- sort(unique(dat$id))
  if (!is.null(n_subsample) && n_subsample < length(all_ids)) {
    sub_ids <- sort(sample(all_ids, n_subsample))
    dat <- dat %>% filter(id %in% sub_ids)
    cat(
      "  Subsampled to",
      n_subsample,
      "of",
      length(all_ids),
      "people (",
      round(100 * n_subsample / length(all_ids), 1),
      "% )\n"
    )
  }

  dat <- dat %>%
    mutate(
      race_grp = case_when(
        race == 1 ~ "Hispanic",
        race == 2 ~ "White",
        race == 3 ~ "Black",
        race %in% c(7, 11, 12, 13) ~ "Other",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(race_grp)) %>%
    mutate(
      race_grp = factor(
        race_grp,
        levels = c("Hispanic", "White", "Black", "Other")
      )
    )

  ids <- sort(unique(dat$id))
  age_mean <- mean(dat$age, na.rm = TRUE)
  age_sd <- sd(dat$age, na.rm = TRUE)
  # pooled across all occasions, so baseline-WtHR and occasion-WtHR share
  # one centering/scaling and stay directly comparable
  whtr_mean <- mean(dat$whtr, na.rm = TRUE)
  whtr_sd <- sd(dat$whtr, na.rm = TRUE)

  race_contr <- contr.sum(4)
  colnames(race_contr) <- paste0("race_c", 1:3)
  race_mat <- race_contr[as.integer(dat$race_grp), , drop = FALSE]

  dat <- dat %>%
    mutate(
      person_idx = match(id, ids),
      time_idx = as.integer(wave),
      age_c = (age - age_mean) / age_sd,
      whtr_c = (whtr - whtr_mean) / whtr_sd,
      race_c1 = race_mat[, 1],
      race_c2 = race_mat[, 2],
      race_c3 = race_mat[, 3]
    )

  # baseline (earliest available occasion) WtHR per person -- the
  # growth-factor impact covariate
  baseline_whtr <- dat %>%
    arrange(person_idx, time_idx) %>%
    distinct(person_idx, .keep_all = TRUE) %>%
    select(person_idx, whtr_baseline_c = whtr_c)

  person_covars <- dat %>%
    distinct(person_idx, race_c1, race_c2, race_c3) %>%
    left_join(baseline_whtr, by = "person_idx") %>%
    arrange(person_idx)

  # raw (uncoded) covariates, for interpretable plotting downstream
  person_raw <- dat %>%
    distinct(person_idx, id, race_grp) %>%
    left_join(
      dat %>%
        arrange(person_idx, time_idx) %>%
        distinct(person_idx, .keep_all = TRUE) %>%
        select(person_idx, whtr_baseline = whtr),
      by = "person_idx"
    ) %>%
    arrange(person_idx)

  dat_long <- dat %>%
    select(
      id,
      person_idx,
      time_idx,
      age_c,
      whtr_c,
      race_c1,
      race_c2,
      race_c3,
      informant_c,
      all_of(ordinal_items)
    ) %>%
    pivot_longer(
      cols = all_of(ordinal_items),
      names_to = "item",
      values_to = "y_raw"
    ) %>%
    mutate(item_idx = match(item, ordinal_items), y_int = as.integer(y_raw)) %>%
    arrange(person_idx, time_idx, item_idx)

  cat(
    "  n persons:",
    length(ids),
    "| n items:",
    length(ordinal_items),
    "| n obs:",
    nrow(dat_long),
    "\n"
  )
  cat("  Race group counts:\n")
  print(table(dat$race_grp[!duplicated(dat$person_idx)]))

  list(
    dat_long = dat_long,
    person_covars = person_covars,
    person_raw = person_raw,
    item_names = ordinal_items,
    ids = ids,
    age_mean = age_mean,
    age_sd = age_sd,
    whtr_mean = whtr_mean,
    whtr_sd = whtr_sd,
    ni = length(ids),
    d = 7L,
    p = length(ordinal_items),
    nobs = nrow(dat_long),
    is_binary = rep(0L, length(ordinal_items)),
    k_items = rep(4L, length(ordinal_items)),
    k_max = 4L
  )
}

prep <- build_lmnlfa_data_staged(
  parent_df,
  youth_df,
  sx,
  n_subsample = n_subsample
)

ximp <- prep$person_covars %>%
  select(race_c1, race_c2, race_c3, whtr_baseline_c) %>%
  as.matrix()
dif_covar_names <- c(
  "age",
  "race_c1",
  "race_c2",
  "race_c3",
  "whtr",
  "informant"
)
xdif <- prep$dat_long %>%
  select(age_c, race_c1, race_c2, race_c3, whtr_c, informant_c) %>%
  as.matrix()

kimp <- ncol(ximp)
kdif <- ncol(xdif)
cat("\nkimp =", kimp, " kdif =", kdif, "\n")

base_stan_data <- list(
  nobs = prep$nobs,
  p = prep$p,
  ni = prep$ni,
  d = prep$d,
  person = prep$dat_long$person_idx,
  itm = prep$dat_long$item_idx,
  time = prep$dat_long$time_idx,
  age_c = prep$dat_long$age_c,
  y = prep$dat_long$y_int,
  is_binary = prep$is_binary,
  k_item = prep$k_items,
  k_max = prep$k_max,
  sigma_l = 1.5,
  sigma_nu = 1.5,
  sigma_cor = 1.0,
  sigma_f = 1.5,
  sigma_di = 0.5
)

fit_settings <- list(
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 400,
  iter_sampling = 400,
  adapt_delta = 0.95,
  init = 0.1,
  refresh = 100,
  show_messages = TRUE,
  seed = 90025
)

bin_dir <- file.path(out_dir, "bin")
dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# STAGE A: GROWTH-FACTOR IMPACT ONLY
# ---------------------------------------------------------------------------
rds_A <- file.path(fits_dir, paste0("fitA_impact_", run_tag, ".rds"))

if (file.exists(rds_A)) {
  cat("\nLoading cached Stage A fit:", rds_A, "\n")
  fitA <- readRDS(rds_A)
} else if (run_stage_A) {
  cat("\n--- Stage A: growth-factor impact only ---\n")
  stan_model_A <- cmdstan_model(
    stan_file_A,
    exe_file = file.path(bin_dir, "lmnlfa-growth-impact-local"),
    force_recompile = TRUE
  )
  stan_data_A <- c(base_stan_data, list(kimp = kimp, ximp = ximp))
  fitA <- do.call(
    stan_model_A$sample,
    c(list(data = stan_data_A), fit_settings)
  )
  fitA$save_object(rds_A)
  diagA <- fitA$diagnostic_summary(quiet = TRUE)
  cat(
    "Stage A divergences:",
    sum(diagA$num_divergent),
    " max treedepth hits:",
    sum(diagA$num_max_treedepth),
    "\n"
  )
  print(
    fitA$summary(
      variables = c("mu_slp", "b_mu_int", "b_mu_slp", "phi_int", "phi_slp")
    ),
    digits = 3
  )
} else {
  stop(
    "Stage A not yet run and no cached fit found. Set run_stage_A <- TRUE ",
    "at the top of this script and re-run when you're ready to sample."
  )
}

b_mu_int_fixed <- fitA$summary(variables = "b_mu_int")$mean
b_mu_slp_fixed <- fitA$summary(variables = "b_mu_slp")$mean
cat("\nStage A impact estimates (fixed going into Stage B):\n")
print(setNames(round(b_mu_int_fixed, 3), paste0("b_mu_int_", colnames(ximp))))
print(setNames(round(b_mu_slp_fixed, 3), paste0("b_mu_slp_", colnames(ximp))))

# ---------------------------------------------------------------------------
# STAGE B: DIF SCREENING (impact fixed)
# ---------------------------------------------------------------------------
rds_B <- file.path(fits_dir, paste0("fitB_difscreen_", run_tag, ".rds"))

if (file.exists(rds_B)) {
  cat("\nLoading cached Stage B fit:", rds_B, "\n")
  fitB <- readRDS(rds_B)
} else if (run_stage_B) {
  cat("\n--- Stage B: DIF screening (impact fixed at Stage A estimates) ---\n")
  stan_model_B <- cmdstan_model(
    stan_file_B,
    exe_file = file.path(bin_dir, "lmnlfa-growth-difscreen-local"),
    force_recompile = TRUE
  )
  stan_data_B <- c(
    base_stan_data,
    list(
      kimp = kimp,
      kdif = kdif,
      ximp = ximp,
      xdif = xdif,
      b_mu_int_fixed = b_mu_int_fixed,
      b_mu_slp_fixed = b_mu_slp_fixed
    )
  )
  fitB <- do.call(
    stan_model_B$sample,
    c(list(data = stan_data_B), fit_settings)
  )
  fitB$save_object(rds_B)
  diagB <- fitB$diagnostic_summary(quiet = TRUE)
  cat(
    "Stage B divergences:",
    sum(diagB$num_divergent),
    " max treedepth hits:",
    sum(diagB$num_max_treedepth),
    "\n"
  )
} else {
  stop(
    "Stage B not yet run and no cached fit found. Set run_stage_B <- TRUE ",
    "at the top of this script and re-run when you're ready to sample."
  )
}

# ---------------------------------------------------------------------------
# DIF SELECTION: Benjamini-Hochberg FDR correction on posterior tail
# probabilities, matching mnlfa_crosssectional_staged.R's convention
# (Gottfredson et al. 2019 / aMNLFA): BH-correct loading-DIF tests first
# (family of p x kdif tests, q = .05); any item x covariate cell with
# significant loading DIF automatically retains intercept DIF too;
# BH-correct intercept DIF separately for the remaining cells. On top of
# BH, also require |posterior mean| > mag_floor, since large N clears BH
# on statistical power alone even for trivial effects.
# ---------------------------------------------------------------------------
select_dif <- function(
  fitB,
  p,
  kdif,
  covar_names,
  fdr_q = 0.05,
  mag_floor = 0.05
) {
  l_draws <- fitB$draws(variables = "l_dif", format = "matrix")
  n_draws <- fitB$draws(variables = "n_dif", format = "matrix")

  tail_prob <- function(draws, varname, p, kdif) {
    m <- matrix(NA_real_, p, kdif)
    for (i in 1:p) {
      for (k in 1:kdif) {
        col <- draws[, paste0(varname, "[", i, ",", k, "]")]
        m[i, k] <- 2 * min(mean(col > 0), mean(col < 0))
      }
    }
    m
  }
  post_mean <- function(draws, varname, p, kdif) {
    m <- matrix(NA_real_, p, kdif)
    for (i in 1:p) {
      for (k in 1:kdif) {
        m[i, k] <- mean(draws[, paste0(varname, "[", i, ",", k, "]")])
      }
    }
    m
  }

  l_p <- tail_prob(l_draws, "l_dif", p, kdif)
  n_p <- tail_prob(n_draws, "n_dif", p, kdif)
  l_mean <- post_mean(l_draws, "l_dif", p, kdif)
  n_mean <- post_mean(n_draws, "n_dif", p, kdif)

  l_adj <- matrix(
    p.adjust(as.vector(l_p), method = "BH"),
    nrow = p,
    ncol = kdif
  )
  l_pattern <- ((l_adj < fdr_q) & (abs(l_mean) > mag_floor)) * 1

  remaining <- l_pattern == 0
  n_adj_remaining <- p.adjust(n_p[remaining], method = "BH")
  n_pattern <- l_pattern
  n_pattern[remaining] <- ((n_adj_remaining < fdr_q) &
    (abs(n_mean[remaining]) > mag_floor)) *
    1

  rownames(l_pattern) <- rownames(n_pattern) <- prep$item_names
  colnames(l_pattern) <- colnames(n_pattern) <- covar_names

  list(
    l_pattern = l_pattern,
    n_pattern = n_pattern,
    l_mean = l_mean,
    n_mean = n_mean
  )
}

dif_sel <- select_dif(fitB, prep$p, kdif, dif_covar_names)
cat(
  "\nDIF selection (Benjamini-Hochberg FDR-corrected, q = .05, |effect| > .05):\n"
)
cat("Loading DIF retained:\n")
print(dif_sel$l_pattern)
cat("Intercept DIF retained:\n")
print(dif_sel$n_pattern)

write.csv(
  as.data.frame(dif_sel$l_pattern) %>% rownames_to_column("item"),
  file.path(dif_selection_dir, paste0("dif_selection_loading_", run_tag, ".csv")),
  row.names = FALSE
)
write.csv(
  as.data.frame(dif_sel$n_pattern) %>% rownames_to_column("item"),
  file.path(dif_selection_dir, paste0("dif_selection_intercept_", run_tag, ".csv")),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# STAGE C: FINAL COMBINED FIT (impact free, screened DIF only)
# ---------------------------------------------------------------------------
rds_C <- file.path(fits_dir, paste0("fitC_final_", run_tag, ".rds"))

if (file.exists(rds_C)) {
  cat("\nLoading cached Stage C fit:", rds_C, "\n")
  fitC <- readRDS(rds_C)
} else if (run_stage_C) {
  cat("\n--- Stage C: final combined fit ---\n")
  stan_model_C <- cmdstan_model(
    stan_file_C,
    exe_file = file.path(bin_dir, "lmnlfa-growth-final-local"),
    force_recompile = TRUE
  )
  stan_data_C <- c(
    base_stan_data,
    list(
      kimp = kimp,
      kdif = kdif,
      ximp = ximp,
      xdif = xdif,
      l_pattern = dif_sel$l_pattern,
      n_pattern = dif_sel$n_pattern,
      ml = sum(dif_sel$l_pattern),
      mn = sum(dif_sel$n_pattern)
    )
  )
  fitC <- do.call(
    stan_model_C$sample,
    c(list(data = stan_data_C), fit_settings)
  )
  fitC$save_object(rds_C)
  diagC <- fitC$diagnostic_summary(quiet = TRUE)
  cat(
    "Stage C divergences:",
    sum(diagC$num_divergent),
    " max treedepth hits:",
    sum(diagC$num_max_treedepth),
    "\n"
  )
  print(
    fitC$summary(
      variables = c("mu_slp", "b_mu_int", "b_mu_slp", "phi_int", "phi_slp")
    ),
    digits = 3
  )
} else {
  cat(
    "\nStage C not yet run. Set run_stage_C <- TRUE at the top of this ",
    "script and re-run when you're ready to sample the final model.\n"
  )
}

# ---------------------------------------------------------------------------
# FACTOR SCORES + INTERPRETABLE VISUALIZATIONS (Stage C)
# fac_gr is saved directly (transformed parameters) in every stage's Stan
# file, so no reconstruction is needed the way the cross-sectional model's
# local eta required.
# ---------------------------------------------------------------------------
if (exists("fitC")) {
  covar_labels <- c("Hispanic", "White", "Black", "WtHR (baseline)")

  growth_summ <- fitC$summary(
    variables = c("mu_slp", "phi_int", "phi_slp", "rho", "eti_sd")
  )
  cat("\nGrowth parameters (Stage C):\n")
  print(
    growth_summ[, c("variable", "mean", "sd", "q5", "q95", "rhat", "ess_bulk")],
    digits = 3
  )
  write.csv(
    growth_summ,
    file.path(pubertal_estimates_dir, paste0("growth_params_", run_tag, ".csv")),
    row.names = FALSE
  )

  fac_gr_summ <- fitC$summary(variables = "fac_gr")
  intercepts <- fac_gr_summ$mean[grepl("^fac_gr\\[1,", fac_gr_summ$variable)]
  slopes <- fac_gr_summ$mean[grepl("^fac_gr\\[2,", fac_gr_summ$variable)]

  scores <- prep$person_raw %>%
    mutate(intercept = intercepts, slope = slopes)
  write.csv(
    scores,
    file.path(pubertal_estimates_dir, paste0("growth_factor_scores_", run_tag, ".csv")),
    row.names = FALSE
  )
  cat("\nGrowth factor scores written:", nrow(scores), "rows\n")

  # (1) individual trajectories by race, sampled people (excludes wobble)
  age_grid <- seq(-3, 3, by = 0.1)
  set.seed(90025)
  samp_idx <- sample(seq_len(prep$ni), min(200, prep$ni))
  traj_df <- map_dfr(samp_idx, function(k) {
    tibble(
      person_idx = k,
      race_grp = scores$race_grp[k],
      age_c = age_grid,
      age = age_grid * prep$age_sd + prep$age_mean,
      eta = intercepts[k] + slopes[k] * age_grid
    )
  })

  race_levels <- c("Hispanic", "White", "Black", "Other")
  p_spaghetti <- ggplot(
    traj_df,
    aes(x = age, y = eta, group = person_idx, colour = race_grp)
  ) +
    geom_line(alpha = 0.2, linewidth = 0.3) +
    scale_colour_manual(values = setNames(pal_chains, race_levels)) +
    labs(
      title = paste0("Individual puberty trajectories by race - ", run_label),
      subtitle = paste0(
        "Predicted from growth factors (n = ",
        length(samp_idx),
        " sampled); excludes occasion-specific wobble"
      ),
      x = "Age (years)",
      y = "Latent puberty (eta)",
      colour = "Race/ethnicity"
    ) +
    theme_minimal(base_size = 13)
  ggsave(
    file.path(trajectories_dir, paste0("trajectories_by_race_", run_tag, ".png")),
    p_spaghetti,
    width = 8.5,
    height = 5,
    dpi = 180
  )

  # population mean trajectory + 90% CI, from posterior draws directly --
  # at ximp = 0 (race averaged across groups -- these are effect/
  # contr.sum codes, so ximp = 0 is the unweighted average of the 4
  # race-group means, NOT the sample's actual (unbalanced) grand mean;
  # WtHR at its own sample mean), b_mu_int/b_mu_slp's contribution is
  # exactly 0 by construction, so this is mu_slp alone, same computation
  # as the informant-only model's mean_trajectory plot
  mu_slp_draws <- fitC$draws(variables = "mu_slp", format = "matrix")
  mean_traj <- map_dfr(age_grid, function(a) {
    vals <- as.vector(mu_slp_draws) * a
    tibble(
      age_c = a,
      eta_med = median(vals),
      eta_lo = quantile(vals, 0.05),
      eta_hi = quantile(vals, 0.95)
    )
  }) %>%
    mutate(age = age_c * prep$age_sd + prep$age_mean)

  p_meantraj <- ggplot(mean_traj, aes(x = age)) +
    geom_ribbon(
      aes(ymin = eta_lo, ymax = eta_hi),
      fill = pal_primary_fill,
      alpha = 0.35
    ) +
    geom_line(aes(y = eta_med), colour = pal_primary, linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    labs(
      title = paste0("Mean puberty growth trajectory - ", run_label),
      subtitle = "Posterior median + 90% CI; race averaged across groups, WtHR at sample mean",
      x = "Age (years)",
      y = "Latent puberty (eta)"
    ) +
    theme_minimal(base_size = 13)
  ggsave(
    file.path(trajectories_dir, paste0("mean_trajectory_", run_tag, ".png")),
    p_meantraj,
    width = 7,
    height = 5,
    dpi = 180
  )

  # (2) Stage A vs Stage C impact comparison
  a_int <- fitA$summary(variables = "b_mu_int")
  a_slp <- fitA$summary(variables = "b_mu_slp")
  c_int <- fitC$summary(variables = "b_mu_int")
  c_slp <- fitC$summary(variables = "b_mu_slp")
  cmp <- bind_rows(
    tibble(
      covariate = covar_labels,
      growth_param = "Intercept",
      mean = a_int$mean,
      q5 = a_int$q5,
      q95 = a_int$q95,
      stage = "A: Impact only"
    ),
    tibble(
      covariate = covar_labels,
      growth_param = "Slope",
      mean = a_slp$mean,
      q5 = a_slp$q5,
      q95 = a_slp$q95,
      stage = "A: Impact only"
    ),
    tibble(
      covariate = covar_labels,
      growth_param = "Intercept",
      mean = c_int$mean,
      q5 = c_int$q5,
      q95 = c_int$q95,
      stage = "C: Impact + DIF"
    ),
    tibble(
      covariate = covar_labels,
      growth_param = "Slope",
      mean = c_slp$mean,
      q5 = c_slp$q5,
      q95 = c_slp$q95,
      stage = "C: Impact + DIF"
    )
  ) %>%
    mutate(covariate = factor(covariate, levels = rev(covar_labels)))

  p_cmp <- ggplot(
    cmp,
    aes(x = mean, y = covariate, colour = stage, shape = stage)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(
      aes(xmin = q5, xmax = q95),
      position = position_dodge(width = 0.5),
      size = 0.5
    ) +
    facet_wrap(~growth_param) +
    scale_colour_manual(
      values = c("A: Impact only" = pal_two[2], "C: Impact + DIF" = pal_two[1])
    ) +
    scale_shape_manual(
      values = c(
        "A: Impact only" = pal_shapes_two[2],
        "C: Impact + DIF" = pal_shapes_two[1]
      )
    ) +
    labs(
      title = paste0(
        "Growth-factor impact, before vs after separating DIF - ",
        run_label
      ),
      subtitle = "Effect-coded deviations from grand mean (race); WtHR is a 1-SD effect; posterior mean + 90% CI",
      x = "Effect on growth factor",
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
  ggsave(
    file.path(impact_comparisons_dir, paste0("impact_comparison_A_vs_C_", run_tag, ".png")),
    p_cmp,
    width = 9,
    height = 5,
    dpi = 180
  )

  # (3) DIF illustration: item characteristic curve for the item/covariate
  # cell with the largest retained loading DIF (whichever covariate that
  # turns out to be -- age, race, WtHR, or informant)
  if (sum(dif_sel$l_pattern) > 0) {
    l_dif_summ_all <- fitC$summary(variables = "l_dif")
    retained_vals <- which(dif_sel$l_pattern == 1, arr.ind = TRUE)
    retained_names <- paste0(
      "l_dif[",
      retained_vals[, 1],
      ",",
      retained_vals[, 2],
      "]"
    )
    retained_summ <- l_dif_summ_all %>% filter(variable %in% retained_names)
    target_row <- retained_summ[which.max(abs(retained_summ$mean)), ]
    idx <- as.integer(regmatches(
      target_row$variable,
      regexec("\\[(\\d+),(\\d+)\\]", target_row$variable)
    )[[1]][2:3])
    it_idx <- idx[1]
    cov_idx <- idx[2]
    target_item <- prep$item_names[it_idx]
    target_covar <- dif_covar_names[cov_idx]

    lp_mean <- fitC$summary(variables = "lp")$mean[it_idx]
    np_mean <- fitC$summary(variables = "np")$mean[it_idx]
    tau_summ <- fitC$summary(variables = "tau")
    tau1 <- tau_summ$mean[tau_summ$variable == paste0("tau[", it_idx, ",1]")]
    l_dif_val <- fitC$summary(variables = "l_dif")$mean[
      fitC$summary(variables = "l_dif")$variable ==
        paste0("l_dif[", it_idx, ",", cov_idx, "]")
    ]
    n_dif_val <- fitC$summary(variables = "n_dif")$mean[
      fitC$summary(variables = "n_dif")$variable ==
        paste0("n_dif[", it_idx, ",", cov_idx, "]")
    ]

    eta_grid2 <- seq(-3, 3, by = 0.1)
    curve_df <- map_dfr(
      c("Low (-1 SD)" = -1, "Grand mean" = 0, "High (+1 SD)" = 1),
      function(x_val) {
        nu <- np_mean + n_dif_val * x_val
        lam <- lp_mean * exp(l_dif_val * x_val)
        tibble(eta = eta_grid2, prob = plogis(nu + lam * eta_grid2 - tau1))
      },
      .id = "level"
    )
    curve_df$level <- factor(
      curve_df$level,
      levels = c("Low (-1 SD)", "Grand mean", "High (+1 SD)")
    )

    p_dif <- ggplot(
      curve_df,
      aes(x = eta, y = prob, colour = level, linetype = level)
    ) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(
        values = setNames(pal_three_ref, levels(curve_df$level))
      ) +
      scale_linetype_manual(
        values = setNames(pal_linetypes_three, levels(curve_df$level))
      ) +
      labs(
        title = paste0(
          "Item characteristic curve for '",
          target_item,
          "' by ",
          target_covar,
          " - ",
          run_label
        ),
        subtitle = "P(response above lowest category) vs. latent puberty; other DIF covariates held at reference",
        x = "Latent puberty (eta)",
        y = "P(response > lowest category)",
        colour = target_covar,
        linetype = target_covar
      ) +
      theme_minimal(base_size = 13)
    ggsave(
      file.path(
        dif_illustrations_dir,
        paste0(
          "dif_illustration_",
          target_item,
          "_",
          target_covar,
          "_",
          run_tag,
          ".png"
        )
      ),
      p_dif,
      width = 8.5,
      height = 5.5,
      dpi = 180
    )
    cat(
      "\nLargest retained DIF: item '",
      target_item,
      "' x covariate '",
      target_covar,
      "' | loading DIF:",
      round(l_dif_val, 3),
      "| intercept DIF:",
      round(n_dif_val, 3),
      "\n"
    )
  } else {
    cat("\nNo DIF terms retained -- skipping DIF illustration plot.\n")
  }

  cat("\nVisualizations written to:", out_dir, "\n")
}

cat("\nAll outputs written to:", out_dir, "\n")
cat("Done:", sx, "\n")
