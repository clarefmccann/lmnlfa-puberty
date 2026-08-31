## mnlfa_crosssectional_staged.R
## Staged version of mnlfa_crosssectional.R, fixing the impact-vs-DIF
## confound diagnosed in the fully-saturated single-stage fit (Rhat
## 1.37-1.54, ESS 7-9/4000 despite 0 divergences/max-treedepth -- a
## between-chain multimodality problem: age's effect on item responses is
## only weakly separable between "real impact" and "DIF" when both are
## estimated freely together with no screening).
##
## Three stages, each cached to an RDS so a re-run resumes rather than
## refitting (same pattern as lmnlfa_hpc.R's Step 1 / Step 2):
##   A. Impact only (age, race, WHtR on factor mean/variance), no DIF.
##      -> scripts/stan/mnlfa-crosssectional-impact.stan
##   B. DIF screening: impact FIXED at Stage A's posterior means, every
##      item x DIF-covariate (age, race, WHtR, informant) loading/intercept
##      DIF term estimated freely. Fixing impact removes its ability to
##      compete with DIF for the same variance.
##      -> scripts/stan/mnlfa-crosssectional-difscreen.stan
##   C. Final: impact free again, DIF restricted to only the item x
##      covariate cells whose 90% CI excludes 0 in Stage B (with the usual
##      MNLFA rule: if a loading-DIF term is retained, the matching
##      intercept-DIF term is retained too).
##      -> scripts/stan/mnlfa-crosssectional-final.stan
##
## IMPORTANT: this script does NOT run any of the three fits automatically
## -- each stage's fit call is present but gated behind `run_stage_X <-
## TRUE/FALSE` flags below so you can review before committing compute.
## Set the flags and run yourself (locally or on HPC).
##
## Usage:
##   Rscript mnlfa_crosssectional_staged.R <sex>
##   sex: female | male
##
## Outputs: written to <OUT_DIR>/mnlfa_crosssectional_staged/

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
# STAGE SWITCHES -- review before running
# ---------------------------------------------------------------------------
run_stage_A <- FALSE
run_stage_B <- FALSE
run_stage_C <- TRUE

# ---------------------------------------------------------------------------
# ARGUMENTS
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript mnlfa_crosssectional_staged.R <female|male>")
}
sx <- args[1]
if (!sx %in% c("female", "male")) {
  stop("sex must be 'female' or 'male'")
}
cat("Sex:", sx, "\n")

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
out_dir <- file.path(out_base, "mnlfa_crosssectional_staged")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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
if (!file.exists(palette_file)) palette_file <- file.path("scripts", "color_palette.R")
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
stan_file_A <- find_stan("mnlfa-crosssectional-impact.stan")
stan_file_B <- find_stan("mnlfa-crosssectional-difscreen.stan")
stan_file_C <- find_stan("mnlfa-crosssectional-final.stan")

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

# ---------------------------------------------------------------------------
# BUILD CROSS-SECTIONAL DATA (identical to mnlfa_crosssectional.R)
# ---------------------------------------------------------------------------
build_crosssectional_data <- function(
  parent_df,
  youth_df,
  sex_label,
  ordinal_items = c("peta", "petb", "petc", "petd")
) {
  cat("\n=== Building cross-sectional MNLFA data |", sex_label, "===\n")

  covars <- c("id", "wave", "age", "race", "whtr")

  whtr_min <- 0.25
  whtr_max <- 0.75

  clean_reporter <- function(df, reporter_label) {
    df %>%
      select(all_of(covars), all_of(ordinal_items)) %>%
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
      mutate(reporter = reporter_label)
  }

  n_whtr_implausible <- sum(
    bind_rows(parent_df, youth_df)$whtr < whtr_min |
      bind_rows(parent_df, youth_df)$whtr > whtr_max,
    na.rm = TRUE
  )
  cat(
    "  Rows excluded for implausible WHtR (outside",
    whtr_min,
    "-",
    whtr_max,
    "):",
    n_whtr_implausible,
    "\n"
  )

  candidates <- bind_rows(
    clean_reporter(parent_df, "parent"),
    clean_reporter(youth_df, "youth")
  )

  chosen_waves <- candidates %>%
    distinct(id, wave) %>%
    group_by(id) %>%
    slice_sample(n = 1) %>%
    ungroup()

  dat <- candidates %>%
    inner_join(chosen_waves, by = c("id", "wave"))

  if (nrow(dat) < 200) {
    stop("Insufficient data for ", sex_label)
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
  whtr_mean <- mean(dat$whtr, na.rm = TRUE)
  whtr_sd <- sd(dat$whtr, na.rm = TRUE)

  race_contr <- contr.sum(4)
  colnames(race_contr) <- paste0("race_c", 1:3)
  race_mat <- race_contr[as.integer(dat$race_grp), , drop = FALSE]

  dat <- dat %>%
    mutate(
      person_idx = match(id, ids),
      age_c = (age - age_mean) / age_sd,
      whtr_c = (whtr - whtr_mean) / whtr_sd,
      informant_c = if_else(reporter == "parent", 1, -1),
      race_c1 = race_mat[, 1],
      race_c2 = race_mat[, 2],
      race_c3 = race_mat[, 3]
    ) %>%
    arrange(person_idx, reporter)

  person_covars <- dat %>%
    distinct(person_idx, age_c, race_c1, race_c2, race_c3, whtr_c) %>%
    arrange(person_idx)

  # raw (uncoded) covariates, for interpretable plotting downstream
  person_raw <- dat %>%
    distinct(person_idx, id, age, race_grp, whtr) %>%
    arrange(person_idx)

  dat_long <- dat %>%
    select(
      id,
      person_idx,
      reporter,
      age_c,
      race_c1,
      race_c2,
      race_c3,
      whtr_c,
      informant_c,
      all_of(ordinal_items)
    ) %>%
    pivot_longer(
      cols = all_of(ordinal_items),
      names_to = "item",
      values_to = "y_raw"
    ) %>%
    mutate(item_idx = match(item, ordinal_items), y_int = as.integer(y_raw)) %>%
    arrange(person_idx, reporter, item_idx)

  cat(
    "  n persons:",
    length(ids),
    "| n rows (person x reporter):",
    nrow(dat),
    "| n obs (x item):",
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
    ni = length(ids),
    p = length(ordinal_items),
    nobs = nrow(dat_long),
    k_items = rep(4L, length(ordinal_items)),
    k_max = 4L
  )
}

prep <- build_crosssectional_data(parent_df, youth_df, sx)

ximp <- prep$person_covars %>%
  select(age_c, race_c1, race_c2, race_c3, whtr_c) %>%
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

# ---------------------------------------------------------------------------
# STAGE A: IMPACT ONLY
# ---------------------------------------------------------------------------
rds_A <- file.path(out_dir, paste0("fitA_impact_", sx, ".rds"))

if (file.exists(rds_A)) {
  cat("\nLoading cached Stage A fit:", rds_A, "\n")
  fitA <- readRDS(rds_A)
} else if (run_stage_A) {
  cat("\n--- Stage A: impact only ---\n")
  bin_dir <- file.path(out_dir, "bin")
  dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)
  stan_model_A <- cmdstan_model(
    stan_file_A,
    exe_file = file.path(bin_dir, "mnlfa-crosssectional-impact-local"),
    force_recompile = TRUE
  )
  stan_data_A <- list(
    nobs = prep$nobs,
    p = prep$p,
    ni = prep$ni,
    person = prep$dat_long$person_idx,
    itm = prep$dat_long$item_idx,
    y = prep$dat_long$y_int,
    k_item = prep$k_items,
    k_max = prep$k_max,
    kimp = kimp,
    ximp = ximp,
    sigma_l = 1.5,
    sigma_nu = 1.5,
    sigma_f = 1.5
  )
  fitA <- stan_model_A$sample(
    data = stan_data_A,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.9,
    init = 0.1,
    refresh = 100,
    show_messages = TRUE,
    seed = 90025
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
  print(fitA$summary(variables = c("b_mu", "b_phi", "phi0")), digits = 3)
} else {
  stop(
    "Stage A not yet run and no cached fit found. Set run_stage_A <- TRUE ",
    "at the top of this script and re-run when you're ready to sample."
  )
}

b_mu_fixed <- fitA$summary(variables = "b_mu")$mean
b_phi_fixed <- fitA$summary(variables = "b_phi")$mean
cat("\nStage A impact estimates (fixed going into Stage B):\n")
print(setNames(round(b_mu_fixed, 3), paste0("b_mu_", colnames(ximp))))
print(setNames(round(b_phi_fixed, 3), paste0("b_phi_", colnames(ximp))))

# ---------------------------------------------------------------------------
# STAGE B: DIF SCREENING (impact fixed)
# ---------------------------------------------------------------------------
rds_B <- file.path(out_dir, paste0("fitB_difscreen_", sx, ".rds"))

if (file.exists(rds_B)) {
  cat("\nLoading cached Stage B fit:", rds_B, "\n")
  fitB <- readRDS(rds_B)
} else if (run_stage_B) {
  cat("\n--- Stage B: DIF screening (impact fixed at Stage A estimates) ---\n")
  bin_dir <- file.path(out_dir, "bin")
  dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)
  stan_model_B <- cmdstan_model(
    stan_file_B,
    exe_file = file.path(bin_dir, "mnlfa-crosssectional-difscreen-local"),
    force_recompile = TRUE
  )
  stan_data_B <- list(
    nobs = prep$nobs,
    p = prep$p,
    ni = prep$ni,
    person = prep$dat_long$person_idx,
    itm = prep$dat_long$item_idx,
    y = prep$dat_long$y_int,
    k_item = prep$k_items,
    k_max = prep$k_max,
    kimp = kimp,
    kdif = kdif,
    ximp = ximp,
    xdif = xdif,
    b_mu_fixed = b_mu_fixed,
    b_phi_fixed = b_phi_fixed,
    sigma_l = 1.5,
    sigma_nu = 1.5,
    sigma_f = 1.5,
    sigma_di = 0.5
  )
  fitB <- stan_model_B$sample(
    data = stan_data_B,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.9,
    init = 0.1,
    refresh = 100,
    show_messages = TRUE,
    seed = 90025
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
# probabilities, matching the aMNLFA convention (Gottfredson et al. 2019):
# BH-correct loading-DIF tests first (family of p x kdif tests, q = .05);
# any item x covariate cell with significant loading DIF necessarily
# retains intercept DIF too (predictors in the loading equation must also
# appear in the intercept equation); BH-correct intercept DIF separately,
# but only for the remaining cells not already retained via loading.
#
# On top of BH, also require |posterior mean| > mag_floor: with ~9,800
# person-reporter rows and ~39,000 item-level observations, even trivial
# effects (2-5% loading shifts, <0.03-logit intercept shifts) clear BH
# comfortably on statistical power alone. The floor targets that specific
# big-N-only tail without touching the substantial effects (loading shifts
# of 15-35%; intercept shifts comparable to or exceeding items' own
# baseline intercepts). The floor is NOT applied to intercept cells forced
# in via the loading rule -- that inclusion is a model-specification
# requirement, not a significance decision, so it stands regardless of the
# forced cell's own intercept magnitude.
#
# Replaces the earlier 90%-CI-with-no-correction screen, which retained an
# implausibly dense 43/48 cells -- far more than the ~5 false positives
# expected by chance alone at that threshold across 48 simultaneous tests.
# This reuses Stage B's already-converged draws; no new sampling needed.
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
        # two-sided posterior tail probability: the Bayesian analog of a
        # p-value -- proportion of the posterior on the "wrong side" of 0
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

  # intercept DIF: cells retained via loading are automatically included
  # (regardless of the forced cell's own intercept magnitude); BH-correct
  # + magnitude-floor only the remaining cells' own intercept tests
  remaining <- l_pattern == 0
  n_adj_remaining <- p.adjust(n_p[remaining], method = "BH")
  n_pattern <- l_pattern # start from loading-forced inclusions
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
  file.path(out_dir, paste0("dif_selection_loading_", sx, ".csv")),
  row.names = FALSE
)
write.csv(
  as.data.frame(dif_sel$n_pattern) %>% rownames_to_column("item"),
  file.path(out_dir, paste0("dif_selection_intercept_", sx, ".csv")),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# STAGE C: FINAL COMBINED FIT (impact free, screened DIF only)
# ---------------------------------------------------------------------------
rds_C <- file.path(out_dir, paste0("fitC_final_", sx, ".rds"))

if (file.exists(rds_C)) {
  cat("\nLoading cached Stage C fit:", rds_C, "\n")
  fitC <- readRDS(rds_C)
} else if (run_stage_C) {
  cat("\n--- Stage C: final combined fit ---\n")
  bin_dir <- file.path(out_dir, "bin")
  dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)
  stan_model_C <- cmdstan_model(
    stan_file_C,
    exe_file = file.path(bin_dir, "mnlfa-crosssectional-final-local"),
    force_recompile = TRUE
  )
  stan_data_C <- list(
    nobs = prep$nobs,
    p = prep$p,
    ni = prep$ni,
    person = prep$dat_long$person_idx,
    itm = prep$dat_long$item_idx,
    y = prep$dat_long$y_int,
    k_item = prep$k_items,
    k_max = prep$k_max,
    kimp = kimp,
    kdif = kdif,
    ximp = ximp,
    xdif = xdif,
    l_pattern = dif_sel$l_pattern,
    n_pattern = dif_sel$n_pattern,
    ml = sum(dif_sel$l_pattern),
    mn = sum(dif_sel$n_pattern),
    sigma_l = 1.5,
    sigma_nu = 1.5,
    sigma_f = 1.5,
    sigma_di = 0.5
  )
  fitC <- stan_model_C$sample(
    data = stan_data_C,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.9,
    init = 0.1,
    refresh = 100,
    show_messages = TRUE,
    seed = 90025
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
  print(fitC$summary(variables = c("b_mu", "b_phi", "phi0")), digits = 3)
} else {
  cat(
    "\nStage C not yet run. Set run_stage_C <- TRUE at the top of this ",
    "script and re-run when you're ready to sample the final model.\n"
  )
}

# ---------------------------------------------------------------------------
# FACTOR SCORES + INTERPRETABLE VISUALIZATIONS (Stage C)
# ---------------------------------------------------------------------------
if (exists("fitC")) {
  covar_labels <- c("Age", "Hispanic", "White", "Black", "WHtR")

  # eta is a local variable inside model{} in the Stan file, so it was never
  # saved (only parameters/transformed parameters/generated quantities are).
  # Reconstruct it per posterior draw from what WAS saved (eta_raw, b_mu,
  # b_phi, phi0) rather than re-running Stan -- same non-centered formula
  # the model itself uses: eta = ximp*b_mu + phi0*exp(ximp*b_phi)*eta_raw.
  draws_eta_raw <- fitC$draws(variables = "eta_raw", format = "matrix")
  draws_b_mu <- fitC$draws(variables = "b_mu", format = "matrix")
  draws_b_phi <- fitC$draws(variables = "b_phi", format = "matrix")
  draws_phi0 <- fitC$draws(variables = "phi0", format = "matrix")

  mu_eta_draws <- draws_b_mu %*% t(ximp)
  sd_eta_draws <- as.vector(draws_phi0) * exp(draws_b_phi %*% t(ximp))
  eta_draws <- mu_eta_draws + sd_eta_draws * draws_eta_raw

  scores <- prep$person_raw %>%
    mutate(eta_mean = colMeans(eta_draws), eta_sd = apply(eta_draws, 2, sd))
  write.csv(
    scores,
    file.path(out_dir, paste0("factor_scores_", sx, ".csv")),
    row.names = FALSE
  )
  cat("\nFactor scores written:", nrow(scores), "rows\n")

  # (1) cross-sectional "trajectory": factor score vs age by race
  race_levels <- c("Hispanic", "White", "Black", "Other")
  p_traj <- ggplot(
    scores,
    aes(x = age, y = eta_mean, colour = race_grp, shape = race_grp, linetype = race_grp)
  ) +
    geom_point(alpha = 0.10, size = 0.6) +
    geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
    scale_colour_manual(values = setNames(pal_chains, race_levels)) +
    scale_shape_manual(values = setNames(pal_shapes_chains, race_levels)) +
    scale_linetype_manual(values = setNames(pal_linetypes_chains, race_levels)) +
    labs(
      title = paste0("Cross-sectional latent puberty by age and race - ", sx),
      subtitle = "One randomly selected wave per person; posterior mean factor scores (Stage C)",
      x = "Age (years)",
      y = "Latent puberty (eta)",
      colour = "Race/ethnicity",
      shape = "Race/ethnicity",
      linetype = "Race/ethnicity"
    ) +
    theme_minimal(base_size = 13)
  ggsave(
    file.path(out_dir, paste0("factor_scores_by_age_race_", sx, ".png")),
    p_traj,
    width = 8,
    height = 5,
    dpi = 180
  )

  # (2) factor score vs WHtR
  p_whtr <- ggplot(scores, aes(x = whtr, y = eta_mean)) +
    geom_point(alpha = 0.10, size = 0.6, colour = pal_primary) +
    geom_smooth(method = "loess", se = TRUE, colour = "black") +
    labs(
      title = paste0("Latent puberty by waist-to-height ratio - ", sx),
      x = "Waist-to-height ratio",
      y = "Latent puberty (eta)"
    ) +
    theme_minimal(base_size = 13)
  ggsave(
    file.path(out_dir, paste0("factor_scores_by_whtr_", sx, ".png")),
    p_whtr,
    width = 7,
    height = 5,
    dpi = 180
  )

  # (3) Stage A vs Stage C impact comparison -- visualizes the impact-vs-DIF
  # confound directly: how much did each covariate's mean-impact estimate
  # move once DIF was properly separated out?
  a_mu <- fitA$summary(variables = "b_mu")
  c_mu <- fitC$summary(variables = "b_mu")
  cmp <- bind_rows(
    tibble(
      covariate = covar_labels,
      mean = a_mu$mean,
      q5 = a_mu$q5,
      q95 = a_mu$q95,
      stage = "A: Impact only"
    ),
    tibble(
      covariate = covar_labels,
      mean = c_mu$mean,
      q5 = c_mu$q5,
      q95 = c_mu$q95,
      stage = "C: Impact + DIF"
    )
  ) %>%
    mutate(covariate = factor(covariate, levels = rev(covar_labels)))

  p_cmp <- ggplot(cmp, aes(x = mean, y = covariate, colour = stage, shape = stage)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(
      aes(xmin = q5, xmax = q95),
      position = position_dodge(width = 0.5),
      size = 0.5
    ) +
    scale_colour_manual(
      values = c("A: Impact only" = pal_two[2], "C: Impact + DIF" = pal_two[1])
    ) +
    scale_shape_manual(
      values = c("A: Impact only" = pal_shapes_two[2], "C: Impact + DIF" = pal_shapes_two[1])
    ) +
    labs(
      title = paste0(
        "Mean impact on latent puberty, before vs after separating DIF - ",
        sx
      ),
      subtitle = "Effect-coded deviations from grand mean (race); posterior mean + 90% CI",
      x = "Effect on latent puberty (mean impact)",
      y = NULL,
      colour = NULL,
      shape = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
  ggsave(
    file.path(out_dir, paste0("impact_comparison_A_vs_C_", sx, ".png")),
    p_cmp,
    width = 7.5,
    height = 5,
    dpi = 180
  )

  # (4) DIF illustration: item characteristic curve for the item with the
  # largest race DIF, showing measurement bias directly -- at the SAME true
  # latent puberty level, how differently does this item get endorsed by
  # race?
  l_dif_summ_all <- fitC$summary(variables = "l_dif")
  target_row <- l_dif_summ_all[which.max(abs(l_dif_summ_all$mean)), ]
  target_idx <- as.integer(regmatches(
    target_row$variable,
    regexec("\\[(\\d+),(\\d+)\\]", target_row$variable)
  )[[1]][2])
  target_item <- prep$item_names[target_idx]
  it_idx <- target_idx

  lp_mean <- fitC$summary(variables = "lp")$mean[it_idx]
  np_mean <- fitC$summary(variables = "np")$mean[it_idx]
  tau_summ <- fitC$summary(variables = "tau")
  tau1 <- tau_summ$mean[tau_summ$variable == paste0("tau[", it_idx, ",1]")]

  l_dif_summ <- fitC$summary(variables = "l_dif")
  n_dif_summ <- fitC$summary(variables = "n_dif")
  l_dif_cell <- function(k) {
    l_dif_summ$mean[
      l_dif_summ$variable == paste0("l_dif[", it_idx, ",", k, "]")
    ]
  }
  n_dif_cell <- function(k) {
    n_dif_summ$mean[
      n_dif_summ$variable == paste0("n_dif[", it_idx, ",", k, "]")
    ]
  }

  race_conditions <- list(
    "Grand mean" = c(0, 0, 0),
    "White" = c(0, 1, 0),
    "Black" = c(0, 0, 1)
  )
  eta_grid <- seq(-3, 3, by = 0.1)
  curve_df <- map_dfr(names(race_conditions), function(rc_label) {
    race_vals <- race_conditions[[rc_label]]
    xdif_vec <- c(0, race_vals, 0, 0) # age, race_c1-3, whtr, informant, all at reference except race
    n_dif_val <- sum(sapply(1:6, n_dif_cell) * xdif_vec)
    l_dif_val <- sum(sapply(1:6, l_dif_cell) * xdif_vec)
    nu <- np_mean + n_dif_val
    lam <- lp_mean * exp(l_dif_val)
    eta_lin <- nu + lam * eta_grid
    tibble(eta = eta_grid, prob = plogis(eta_lin - tau1), race = rc_label)
  })

  p_dif <- ggplot(curve_df, aes(x = eta, y = prob, colour = race, linetype = race)) +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(
      values = c(
        "Grand mean" = pal_three_ref[1],
        "White" = pal_three_ref[2],
        "Black" = pal_three_ref[3]
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Grand mean" = pal_linetypes_three[1],
        "White" = pal_linetypes_three[2],
        "Black" = pal_linetypes_three[3]
      )
    ) +
    labs(
      title = paste0(
        "Item characteristic curve for '",
        target_item,
        "' by race - ",
        sx
      ),
      subtitle = "P(response above lowest category) vs. latent puberty\n(age/WHtR/informant held at reference)",
      x = "Latent puberty (eta)",
      y = "P(response > lowest category)",
      colour = "Race/ethnicity",
      linetype = "Race/ethnicity"
    ) +
    theme_minimal(base_size = 13)
  ggsave(
    file.path(
      out_dir,
      paste0("dif_illustration_", target_item, "_", sx, ".png")
    ),
    p_dif,
    width = 8.5,
    height = 5.5,
    dpi = 180
  )

  cat("Visualizations written to:", out_dir, "\n")
}

cat("\nAll outputs written to:", out_dir, "\n")
cat("Done:", sx, "\n")
