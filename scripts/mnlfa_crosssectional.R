## mnlfa_crosssectional.R
## Cross-sectional MNLFA: randomly selects ONE wave per person, keeps
## whichever reporter(s) (parent/youth) have complete data at that wave, and
## tests impact (age, race, WHtR) and DIF (age, race, WHtR, informant) on the
## 4 sex-specific PDS items.
##
## IMPORTANT: this script only builds the data and compiles the model. It
## deliberately does NOT run the sampler -- run it yourself (locally or on
## HPC) once you've reviewed the data assembly. See the bottom of this file
## for the commented-out fit call.
##
## Requires 00_data_foundation.R to have been re-run with the waist/WHtR
## addition (waist_in, whtr now in covariate_cols) before this will find
## those columns.
##
## Usage:
##   Rscript mnlfa_crosssectional.R <sex>
##   sex: female | male
##
## Stan:    scripts/stan/mnlfa-crosssectional.stan
## Outputs: written to <OUT_DIR>/mnlfa_crosssectional/

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

set.seed(90025)

# ---------------------------------------------------------------------------
# ARGUMENTS
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript mnlfa_crosssectional.R <female|male>")
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
out_dir <- file.path(out_base, "mnlfa_crosssectional")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# SGE copies the job script, so use SGE_O_WORKDIR there. For a direct
# `Rscript path/to/this.R` invocation (e.g. running interactively from
# inside scripts/), sys.frames()[[1]]$ofile is unset -- that trick only
# works when a script is source()'d, not when Rscript runs it directly --
# so fall back to parsing the --file= argument Rscript always passes, which
# resolves correctly regardless of the current working directory.
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
stan_file <- file.path(script_dir, "stan", "mnlfa-crosssectional.stan")
if (!file.exists(stan_file)) {
  stan_file <- file.path("scripts", "stan", "mnlfa-crosssectional.stan")
}
if (!file.exists(stan_file)) {
  stop("Cannot find mnlfa-crosssectional.stan: ", stan_file)
}

cat("Data dir:  ", data_dir, "\n")
cat("Output dir:", out_dir, "\n")
cat("Stan file: ", stan_file, "\n")

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
# BUILD CROSS-SECTIONAL DATA: ONE RANDOM WAVE PER PERSON,
# KEEPING WHICHEVER REPORTER(S) ARE COMPLETE AT THAT WAVE
# ---------------------------------------------------------------------------
build_crosssectional_data <- function(
  parent_df,
  youth_df,
  sex_label,
  ordinal_items = c("peta", "petb", "petc", "petd")
) {
  cat("\n=== Building cross-sectional MNLFA data |", sex_label, "===\n")

  covars <- c("id", "wave", "age", "race", "whtr")

  # Plausibility filter for WHtR: the raw waist circumference variable has a
  # small number of clearly erroneous values (e.g. a 14-year-old with a
  # 150-inch waist, or others under 16 inches) that would otherwise badly
  # distort the mean/SD used to center this covariate. 0.25-0.75 comfortably
  # covers the realistic adolescent range with margin.
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
    "  Rows excluded for implausible WHtR (outside", whtr_min, "-", whtr_max,
    "):", n_whtr_implausible, "\n"
  )

  candidates <- bind_rows(
    clean_reporter(parent_df, "parent"),
    clean_reporter(youth_df, "youth")
  )

  # one randomly chosen wave per person, among waves with >=1 complete reporter
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

  ids <- sort(unique(dat$id))
  age_mean <- mean(dat$age, na.rm = TRUE)
  age_sd <- sd(dat$age, na.rm = TRUE)
  whtr_mean <- mean(dat$whtr, na.rm = TRUE)
  whtr_sd <- sd(dat$whtr, na.rm = TRUE)

  # race: collapse to 4 groups (Hispanic, White, Black, Other), effect-coded
  # (contr.sum) so item intercepts reflect the grand mean across groups.
  # Codes per NBDCtools (ab_g_stc__cohort_ethnrace__mhisp):
  #   1 Hispanic | 2 White | 3 Black | 7/11/12/13 -> Other
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
    mutate(race_grp = factor(race_grp, levels = c("Hispanic", "White", "Black", "Other")))

  # re-derive ids/means after the race_grp filter potentially drops rows
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

  # person-level impact covariates (one row per person; take the first row
  # since age/race/whtr are identical across a person's reporter rows)
  person_covars <- dat %>%
    distinct(person_idx, age_c, race_c1, race_c2, race_c3, whtr_c) %>%
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
    "  n persons:", length(ids),
    "| n rows (person x reporter):", nrow(dat),
    "| n obs (x item):", nrow(dat_long), "\n"
  )
  cat(
    "  Reporter counts -- parent:", sum(dat$reporter == "parent"),
    "| youth:", sum(dat$reporter == "youth"), "\n"
  )
  cat("  Race group counts:\n")
  print(table(dat$race_grp[!duplicated(dat$person_idx)]))
  cat(
    "  Age range:", round(range(dat$age), 1),
    "| mean:", round(age_mean, 2), "| sd:", round(age_sd, 2), "\n"
  )
  cat(
    "  WHtR range:", round(range(dat$whtr, na.rm = TRUE), 3),
    "| mean:", round(whtr_mean, 3), "| sd:", round(whtr_sd, 3), "\n"
  )

  list(
    dat_long = dat_long,
    person_covars = person_covars,
    item_names = ordinal_items,
    ids = ids,
    age_mean = age_mean,
    age_sd = age_sd,
    whtr_mean = whtr_mean,
    whtr_sd = whtr_sd,
    ni = length(ids),
    p = length(ordinal_items),
    nobs = nrow(dat_long),
    k_items = rep(4L, length(ordinal_items)),
    k_max = 4L
  )
}

prep <- build_crosssectional_data(parent_df, youth_df, sx)

# ---------------------------------------------------------------------------
# ASSEMBLE STAN DATA
# ---------------------------------------------------------------------------
ximp <- prep$person_covars %>%
  select(age_c, race_c1, race_c2, race_c3, whtr_c) %>%
  as.matrix()

xdif <- prep$dat_long %>%
  select(age_c, race_c1, race_c2, race_c3, whtr_c, informant_c) %>%
  as.matrix()

stan_data <- list(
  nobs = prep$nobs,
  p = prep$p,
  ni = prep$ni,
  person = prep$dat_long$person_idx,
  itm = prep$dat_long$item_idx,
  y = prep$dat_long$y_int,
  k_item = prep$k_items,
  k_max = prep$k_max,
  kimp = ncol(ximp),
  kdif = ncol(xdif),
  ximp = ximp,
  xdif = xdif,
  sigma_l = 1.5,
  sigma_nu = 1.5,
  sigma_f = 1.5,
  sigma_di = 0.5
)

cat("\nstan_data assembled: kimp =", stan_data$kimp, " kdif =", stan_data$kdif, "\n")

# ---------------------------------------------------------------------------
# COMPILE STAN MODEL (isolated exe -- see lmnlfa_growth_only_tanhcor.R for why)
# ---------------------------------------------------------------------------
bin_dir <- file.path(out_dir, "bin")
dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)
cat("\nCompiling Stan model (local exe)...\n")
stan_model <- cmdstan_model(
  stan_file,
  exe_file = file.path(bin_dir, "mnlfa-crosssectional-local"),
  force_recompile = TRUE
)
cat("Compiled.\n")

saveRDS(stan_data, file.path(out_dir, paste0("stan_data_", sx, ".rds")))
cat("\nstan_data saved to:", file.path(out_dir, paste0("stan_data_", sx, ".rds")), "\n")

# ---------------------------------------------------------------------------
# FIT
# ---------------------------------------------------------------------------
fit <- stan_model$sample(
  data = stan_data,
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
fit$save_object(file.path(out_dir, paste0("fit_crosssectional_", sx, ".rds")))

diag <- fit$diagnostic_summary(quiet = TRUE)
cat("Divergences:", diag$num_divergent, " Max treedepth hits:", diag$num_max_treedepth, "\n")
print(fit$summary(variables = c("b_mu", "b_phi", "phi0")), digits = 3)

cat("\nAll outputs written to:", out_dir, "\n")
cat("Done:", sx, "\n")
