## lmnlfa_growth_sigmoid_staged_ssp.R
## Alternative to lmnlfa_growth_sigmoid_staged.R's DIF-selection step: uses
## Bayesian spike-and-slab regularization (SSP; Chen & Bauer, 2024,
## "Modeling Construct Change Over Time Amidst Potential Changes in
## Construct Measurement: A Longitudinal Moderated Factor Analysis
## Approach") instead of freely estimating every DIF term with a diffuse
## prior and screening post hoc with Benjamini-Hochberg FDR + a magnitude
## floor. Built so the two methods can be run and compared side by side on
## the SAME data, SAME growth curve, and SAME Stage A impact estimates.
##
## WHAT'S SHARED WITH THE BH-FDR SCRIPT (lmnlfa_growth_sigmoid_staged.R),
## AND WHY THIS SCRIPT DOESN'T REFIT THEM:
##   - Stage A (impact only) doesn't involve DIF at all, so it's identical
##     regardless of DIF-selection method. This script LOADS the existing
##     cached Stage A fit from lmnlfa_growth_sigmoid_staged/fits/ rather
##     than refitting it -- if that cache doesn't exist yet, run
##     lmnlfa_growth_sigmoid_staged.R first (with run_stage_A <- TRUE).
##   - Stage C (final) only needs a 0/1 retained-cell pattern as input --
##     its own mechanics don't care HOW that pattern was derived. So this
##     script reuses lmnlfa-growth-sigmoid-final.stan UNCHANGED, just
##     feeding it the SSP-derived pattern instead of the BH-FDR one. This
##     means the comparison between methods is isolated to Stage B's
##     selection logic, not confounded by any other difference.
##
## WHAT'S DIFFERENT:
##   - Stage B here uses scripts/stan/lmnlfa-growth-sigmoid-difscreen-ssp.stan
##     (spike-and-slab DIF parameterization) instead of -difscreen.stan.
##   - DIF is selected by thresholding the posterior mean of the shared
##     inclusion parameter r_incl (Chen & Bauer used 0.7-0.8; default here
##     is 0.8) instead of BH-FDR + magnitude floor.
##   - All outputs go to a SEPARATE lmnlfa_growth_sigmoid_staged_ssp/
##     directory (same fits/trajectories/etc. subfolder structure) so
##     nothing here overwrites the BH-FDR run's results.
##
## SSP HYPERPARAMETERS (ssp_scale, phi_shape, phi_rate below) are NOT
## copied from Chen & Bauer's own paper -- theirs were calibrated for
## binary items and different covariate scaling than ours. The defaults
## here target a lasso-prior SD of roughly 0.3 on the raw DIF magnitude
## (tighter than the BH-FDR model's diffuse sigma_di = 0.5, looser than
## Chen & Bauer's own ~0.05-0.07), as a reasoned starting point -- watch
## r_incl's posterior distribution after the first run (it should look
## clearly bimodal-ish, not stuck uniformly near 0.5) and retune if not.
##
## Usage:
##   Rscript lmnlfa_growth_sigmoid_staged_ssp.R <sex>
##   sex: female | male
##
## Outputs: written to <OUT_DIR>/lmnlfa_growth_sigmoid_staged_ssp/

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

options(error = function() {
  cat("\n\n=== FATAL ERROR TRACEBACK ===\n", file = stderr())
  calls <- sys.calls()
  for (i in rev(seq_along(calls))) {
    cat(i, ": ", paste(deparse(calls[[i]]), collapse = " "), "\n", sep = "", file = stderr())
  }
  cat("=== END TRACEBACK ===\n\n", file = stderr())
  flush(stdout())
  flush(stderr())
  if (!interactive()) quit(save = "no", status = 1, runLast = FALSE)
})

# ---------------------------------------------------------------------------
# STAGE SWITCHES + SSP SETTINGS -- review before running
# ---------------------------------------------------------------------------
run_stage_B <- FALSE
run_stage_C <- FALSE
n_subsample <- 1500 # MUST match the BH-FDR run's n_subsample -- same sample

ssp_scale <- 1 # "u" in Chen & Bauer Eq. 14 -- known SD-like reference scale
phi_shape <- 200 # gamma hyperprior shape for lasso penalties (per covariate)
phi_rate <- 42 # implies mean penalty ~4.76 -> lasso-prior SD ~0.3 on raw DIF
ssp_threshold <- 0.8 # inclusion-probability threshold for "retained" DIF

# ---------------------------------------------------------------------------
# ARGUMENTS
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript lmnlfa_growth_sigmoid_staged_ssp.R <female|male>")
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
sshfs_data_dir <- "/private/tmp/sshfs/projects/abcd-projs/dissertation/study1/outputs"
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  if (dir.exists(sshfs_data_dir)) {
    data_dir <- sshfs_data_dir
  } else {
    data_dir <- file.path(root_path, "projects/abcd-projs/dissertation/study1/outputs")
  }
}
if (!dir.exists(data_dir)) {
  stop("Cannot locate data directory: ", data_dir)
}

out_base <- Sys.getenv("OUT_DIR")
if (!nzchar(out_base)) {
  out_base <- data_dir
}

# Stage A is SHARED with the BH-FDR run -- read-only reference, never
# written to by this script
bhfdr_dir <- file.path(out_base, "lmnlfa_growth_sigmoid_staged")
bhfdr_fits_dir <- file.path(bhfdr_dir, "fits")

# this script's OWN outputs -- entirely separate tree, same subfolder
# taxonomy as the BH-FDR run for direct side-by-side comparison
out_dir <- file.path(out_base, "lmnlfa_growth_sigmoid_staged_ssp")
fits_dir <- file.path(out_dir, "fits")
dif_selection_dir <- file.path(out_dir, "dif-selection")
pubertal_estimates_dir <- file.path(out_dir, "pubertal-estimates")
trajectories_dir <- file.path(out_dir, "trajectories")
impact_comparisons_dir <- file.path(out_dir, "impact-comparisons")
dif_illustrations_dir <- file.path(out_dir, "dif-illustrations")
for (d in c(
  out_dir, fits_dir, dif_selection_dir, pubertal_estimates_dir,
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
stan_file_B_ssp <- find_stan("lmnlfa-growth-sigmoid-difscreen-ssp.stan")
stan_file_C <- find_stan("lmnlfa-growth-sigmoid-final.stan") # SHARED, unchanged

cat("Data dir:  ", data_dir, "\n")
cat("Stage A (shared, read-only):", bhfdr_fits_dir, "\n")
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
# BUILD LONGITUDINAL DATA -- identical to lmnlfa_growth_sigmoid_staged.R.
# Same seed + same build order = same n=1500 sample, so Stage A's cached
# fixed impact estimates apply here without re-deriving them.
# ---------------------------------------------------------------------------
build_lmnlfa_data_staged <- function(
  parent_df,
  youth_df,
  sex_label,
  ordinal_items = c("peta", "petb", "petc", "petd"),
  n_subsample = NULL
) {
  cat(
    "\n=== Building longitudinal data (staged impact + DIF, sigmoid, SSP) |",
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

  baseline_whtr <- dat %>%
    arrange(person_idx, time_idx) %>%
    distinct(person_idx, .keep_all = TRUE) %>%
    select(person_idx, whtr_baseline_c = whtr_c)

  person_covars <- dat %>%
    distinct(person_idx, race_c1, race_c2, race_c3) %>%
    left_join(baseline_whtr, by = "person_idx") %>%
    arrange(person_idx)

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
dif_covar_names <- c("age", "race_c1", "race_c2", "race_c3", "whtr", "informant")
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
  sigma_k = 0.5
)

fit_settings <- list(
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 400,
  iter_sampling = 400,
  adapt_delta = 0.97, # slightly higher than the BH-FDR run's 0.95 -- the
  # spike-and-slab mixture prior's near-discontinuous geometry (Laplace x
  # Beta(.5,.5)) is more prone to divergences than a plain diffuse prior
  init = 0.1,
  refresh = 100,
  show_messages = TRUE,
  seed = 90025
)

bin_dir <- file.path(out_dir, "bin")
dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# STAGE A: LOAD SHARED FIT (never refit here -- see header)
# ---------------------------------------------------------------------------
rds_A <- file.path(bhfdr_fits_dir, paste0("fitA_impact_", run_tag, ".rds"))
if (!file.exists(rds_A)) {
  stop(
    "Stage A fit not found at ", rds_A, ". Run lmnlfa_growth_sigmoid_staged.R ",
    "(with run_stage_A <- TRUE) for this sex FIRST -- Stage A is shared ",
    "between the BH-FDR and SSP pipelines and is not refit here."
  )
}
cat("\nLoading SHARED Stage A fit:", rds_A, "\n")
fitA <- readRDS(rds_A)

b_mu_logk_fixed <- fitA$summary(variables = "b_mu_logk")$mean
b_mu_alpha_fixed <- fitA$summary(variables = "b_mu_alpha")$mean
cat("\nStage A impact estimates (fixed going into Stage B):\n")
print(setNames(round(b_mu_logk_fixed, 3), paste0("b_mu_logk_", colnames(ximp))))
print(setNames(round(b_mu_alpha_fixed, 3), paste0("b_mu_alpha_", colnames(ximp))))

# ---------------------------------------------------------------------------
# STAGE B (SSP): DIF SCREENING VIA SPIKE-AND-SLAB REGULARIZATION
# ---------------------------------------------------------------------------
rds_B <- file.path(fits_dir, paste0("fitB_ssp_", run_tag, ".rds"))

if (file.exists(rds_B)) {
  cat("\nLoading cached Stage B (SSP) fit:", rds_B, "\n")
  fitB <- readRDS(rds_B)
} else if (run_stage_B) {
  cat("\n--- Stage B (SSP): DIF screening via spike-and-slab regularization ---\n")
  stan_model_B <- cmdstan_model(
    stan_file_B_ssp,
    exe_file = file.path(bin_dir, "lmnlfa-growth-sigmoid-difscreen-ssp-local"),
    force_recompile = TRUE
  )
  stan_data_B <- c(
    base_stan_data,
    list(
      kimp = kimp,
      kdif = kdif,
      ximp = ximp,
      xdif = xdif,
      b_mu_logk_fixed = b_mu_logk_fixed,
      b_mu_alpha_fixed = b_mu_alpha_fixed,
      ssp_scale = ssp_scale,
      phi_shape = phi_shape,
      phi_rate = phi_rate
    )
  )
  fitB <- do.call(
    stan_model_B$sample,
    c(list(data = stan_data_B), fit_settings)
  )
  fitB$save_object(rds_B)
  diagB <- fitB$diagnostic_summary(quiet = TRUE)
  cat(
    "Stage B (SSP) divergences:",
    sum(diagB$num_divergent),
    " max treedepth hits:",
    sum(diagB$num_max_treedepth),
    "\n"
  )
  r_summ <- fitB$summary(variables = "r_incl")
  cat(
    "\nr_incl posterior mean distribution -- min:", round(min(r_summ$mean), 3),
    "| max:", round(max(r_summ$mean), 3),
    "| n > ", ssp_threshold, ":", sum(r_summ$mean > ssp_threshold),
    "of", nrow(r_summ), "\n"
  )
} else {
  stop(
    "Stage B (SSP) not yet run and no cached fit found. Set run_stage_B <- TRUE ",
    "at the top of this script and re-run when you're ready to sample."
  )
}

# ---------------------------------------------------------------------------
# DIF SELECTION: threshold the shared inclusion parameter r_incl, rather
# than BH-FDR + magnitude floor. Chen & Bauer used 0.7-0.8; default here is
# 0.8. Both loading and intercept DIF for a cell are gated by the SAME
# r_incl (they share one inclusion parameter -- see the .stan file header),
# so l_pattern and n_pattern end up identical here, unlike the BH-FDR
# version where they can differ.
# ---------------------------------------------------------------------------
select_dif_ssp <- function(fitB, p, kdif, covar_names, threshold = 0.8) {
  r_draws <- fitB$draws(variables = "r_incl", format = "matrix")
  r_mean <- matrix(NA_real_, p, kdif)
  for (i in 1:p) {
    for (k in 1:kdif) {
      r_mean[i, k] <- mean(r_draws[, paste0("r_incl[", i, ",", k, "]")])
    }
  }
  pattern <- (r_mean > threshold) * 1
  rownames(pattern) <- prep$item_names
  colnames(pattern) <- covar_names
  list(l_pattern = pattern, n_pattern = pattern, r_mean = r_mean)
}

dif_sel <- select_dif_ssp(fitB, prep$p, kdif, dif_covar_names, threshold = ssp_threshold)
cat("\nDIF selection (spike-and-slab, inclusion threshold =", ssp_threshold, "):\n")
cat("Retained (shared loading + intercept pattern):\n")
print(dif_sel$l_pattern)

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
# DIF SELECTION HEATMAP -- same style as the BH-FDR script's, for direct
# visual comparison
# ---------------------------------------------------------------------------
covar_display <- setNames(
  c("Age", "Hispanic", "White", "Black", "WtHR", "Informant"),
  c("age", "race_c1", "race_c2", "race_c3", "whtr", "informant")
)[dif_covar_names]
item_display <- c(
  peta = "Growth spurt",
  petb = "Body hair",
  petc = "Skin changes",
  petd = if (sx == "male") "Voice changes" else "Breast growth"
)[prep$item_names]

l_heat <- as.data.frame(dif_sel$l_pattern) %>%
  rownames_to_column("item") %>%
  pivot_longer(-item, names_to = "covariate", values_to = "retained") %>%
  mutate(half = "Loading")
n_heat <- as.data.frame(dif_sel$n_pattern) %>%
  rownames_to_column("item") %>%
  pivot_longer(-item, names_to = "covariate", values_to = "retained") %>%
  mutate(half = "Intercept")

heat_df <- bind_rows(l_heat, n_heat) %>%
  mutate(
    status = case_when(
      half == "Loading" & retained == 1 ~ "Loading DIF",
      half == "Intercept" & retained == 1 ~ "Intercept DIF",
      TRUE ~ "Not retained"
    ),
    item_label = factor(item_display[item], levels = rev(unname(item_display))),
    covar_label = factor(covar_display[covariate], levels = unname(covar_display)),
    half = factor(half, levels = c("Loading", "Intercept"))
  )

universal_covars <- dif_covar_names[sapply(dif_covar_names, function(cv) {
  all(dif_sel$l_pattern[, cv] == 1 | dif_sel$n_pattern[, cv] == 1)
})]
heatmap_caption <- if (length(universal_covars) > 0) {
  covar_list <- paste(covar_display[universal_covars], collapse = " and ")
  raw_caption <- paste0(
    covar_list,
    " DIF retained on every item -> unadjusted PDS sums are not directly ",
    "comparable across levels of ", tolower(covar_list), "."
  )
  paste(strwrap(raw_caption, width = 90), collapse = "\n")
} else {
  NULL
}

p_dif_heat <- ggplot(heat_df, aes(x = half, y = item_label, fill = status)) +
  geom_tile(colour = "white", linewidth = 1) +
  facet_grid(~covar_label, switch = "x") +
  scale_fill_manual(
    values = c(
      "Loading DIF" = unname(sunset_anchors["slate"]),
      "Intercept DIF" = unname(sunset_anchors["bordeaux"]),
      "Not retained" = unname(sunset_anchors["cream"])
    ),
    breaks = c("Loading DIF", "Intercept DIF", "Not retained")
  ) +
  labs(
    title = paste0("Retained DIF terms (spike-and-slab), ", run_label),
    subtitle = paste0("Inclusion threshold = ", ssp_threshold, "; shared loading + intercept pattern"),
    caption = heatmap_caption,
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.spacing = unit(0.3, "lines"),
    legend.position = "bottom",
    plot.caption = element_text(hjust = 0, size = 11)
  )
ggsave(
  file.path(dif_selection_dir, paste0("dif_selection_heatmap_", run_tag, ".png")),
  p_dif_heat,
  width = 11,
  height = 6,
  dpi = 220
)

# ---------------------------------------------------------------------------
# STAGE C: FINAL COMBINED FIT -- reuses the SHARED, unchanged final.stan
# ---------------------------------------------------------------------------
rds_C <- file.path(fits_dir, paste0("fitC_ssp_", run_tag, ".rds"))

if (file.exists(rds_C)) {
  cat("\nLoading cached Stage C fit:", rds_C, "\n")
  fitC <- readRDS(rds_C)
} else if (run_stage_C) {
  cat("\n--- Stage C: final combined fit (SSP-derived DIF pattern) ---\n")
  stan_model_C <- cmdstan_model(
    stan_file_C,
    exe_file = file.path(bin_dir, "lmnlfa-growth-sigmoid-final-local"),
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
      variables = c("mu_logk", "b_mu_logk", "phi_logk", "mu_alpha", "b_mu_alpha", "phi_alpha")
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
# FACTOR SCORES + INTERPRETABLE VISUALIZATIONS (Stage C) -- identical to
# lmnlfa_growth_sigmoid_staged.R from here down, just writing into the SSP
# output directory
# ---------------------------------------------------------------------------
if (exists("fitC")) {
  covar_labels <- c("Hispanic", "White", "Black", "WtHR (baseline)")

  growth_summ <- fitC$summary(
    variables = c("mu_logk", "phi_logk", "mu_alpha", "phi_alpha", "rho", "eti_sd")
  )
  cat("\nGrowth parameters (Stage C, SSP):\n")
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
  log_rates <- fac_gr_summ$mean[grepl("^fac_gr\\[1,", fac_gr_summ$variable)]
  infl_ages_c <- fac_gr_summ$mean[grepl("^fac_gr\\[2,", fac_gr_summ$variable)]

  scores <- prep$person_raw %>%
    mutate(
      log_rate = log_rates,
      rate = exp(log_rates),
      inflection_age_c = infl_ages_c,
      inflection_age = infl_ages_c * prep$age_sd + prep$age_mean
    )
  write.csv(
    scores,
    file.path(pubertal_estimates_dir, paste0("growth_factor_scores_", run_tag, ".csv")),
    row.names = FALSE
  )
  cat("\nGrowth factor scores written:", nrow(scores), "rows\n")
  cat(
    "Inflection age (years) -- median:", round(median(scores$inflection_age), 2),
    "| range:", round(range(scores$inflection_age), 2), "\n"
  )

  age_c_range <- range(prep$dat_long$age_c)
  age_grid <- seq(age_c_range[1], age_c_range[2], length.out = 100)
  set.seed(90025)
  samp_idx <- sample(seq_len(prep$ni), min(200, prep$ni))
  traj_df <- map_dfr(samp_idx, function(k) {
    k_i <- exp(log_rates[k])
    alpha_i <- infl_ages_c[k]
    tibble(
      person_idx = k,
      race_grp = scores$race_grp[k],
      age_c = age_grid,
      age = age_grid * prep$age_sd + prep$age_mean,
      eta = 1 + 4 / (1 + exp(-k_i * (age_grid - alpha_i)))
    )
  })

  race_levels <- c("Hispanic", "White", "Black", "Other")
  p_spaghetti <- ggplot(
    traj_df,
    aes(x = age, y = eta, group = person_idx, colour = race_grp)
  ) +
    geom_line(alpha = 0.3, linewidth = 0.5) +
    scale_colour_manual(values = setNames(pal_chains, race_levels)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 2.5))) +
    coord_cartesian(ylim = c(1, 5)) +
    labs(
      title = paste0("Individual sigmoid puberty trajectories by race (SSP) - ", run_label),
      subtitle = paste0("Predicted from growth factors (n = ", length(samp_idx), " sampled)"),
      x = "Age (years)",
      y = "Latent puberty (eta, 1-5 scale)",
      colour = "Race/ethnicity"
    ) +
    theme_minimal(base_size = 13)
  ggsave(
    file.path(trajectories_dir, paste0("trajectories_by_race_", run_tag, ".png")),
    p_spaghetti,
    width = 10,
    height = 6.5,
    dpi = 220
  )

  mu_logk_draws <- as.vector(fitC$draws(variables = "mu_logk", format = "matrix"))
  mu_alpha_draws <- as.vector(fitC$draws(variables = "mu_alpha", format = "matrix"))
  mean_traj <- map_dfr(age_grid, function(a) {
    vals <- 1 + 4 / (1 + exp(-exp(mu_logk_draws) * (a - mu_alpha_draws)))
    tibble(age_c = a, eta_med = median(vals), eta_lo = quantile(vals, 0.05), eta_hi = quantile(vals, 0.95))
  }) %>%
    mutate(age = age_c * prep$age_sd + prep$age_mean)

  p_meantraj <- ggplot(mean_traj, aes(x = age)) +
    geom_ribbon(aes(ymin = eta_lo, ymax = eta_hi), fill = pal_primary_fill, alpha = 0.35) +
    geom_line(aes(y = eta_med), colour = pal_primary, linewidth = 1.8) +
    coord_cartesian(ylim = c(1, 5)) +
    labs(
      title = paste0("Mean sigmoid puberty growth trajectory (SSP) - ", run_label),
      subtitle = "Posterior median + 90% CI; race averaged across groups, WtHR at sample mean",
      x = "Age (years)",
      y = "Latent puberty (eta, 1-5 scale)"
    ) +
    theme_minimal(base_size = 13)
  ggsave(
    file.path(trajectories_dir, paste0("mean_trajectory_", run_tag, ".png")),
    p_meantraj,
    width = 8.5,
    height = 6.5,
    dpi = 220
  )

  a_logk <- fitA$summary(variables = "b_mu_logk")
  a_alpha <- fitA$summary(variables = "b_mu_alpha")
  c_logk <- fitC$summary(variables = "b_mu_logk")
  c_alpha <- fitC$summary(variables = "b_mu_alpha")
  cmp <- bind_rows(
    tibble(covariate = covar_labels, growth_param = "Log-rate (pubertal tempo)", mean = a_logk$mean, q5 = a_logk$q5, q95 = a_logk$q95, stage = "A: Impact only"),
    tibble(covariate = covar_labels, growth_param = "Inflection age (pubertal timing)", mean = a_alpha$mean, q5 = a_alpha$q5, q95 = a_alpha$q95, stage = "A: Impact only"),
    tibble(covariate = covar_labels, growth_param = "Log-rate (pubertal tempo)", mean = c_logk$mean, q5 = c_logk$q5, q95 = c_logk$q95, stage = "C: Impact + DIF"),
    tibble(covariate = covar_labels, growth_param = "Inflection age (pubertal timing)", mean = c_alpha$mean, q5 = c_alpha$q5, q95 = c_alpha$q95, stage = "C: Impact + DIF")
  ) %>%
    mutate(covariate = factor(covariate, levels = rev(covar_labels)))

  p_cmp <- ggplot(cmp, aes(x = mean, y = covariate, colour = stage, shape = stage)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(xmin = q5, xmax = q95), position = position_dodge(width = 0.5), size = 0.5) +
    facet_wrap(~growth_param) +
    scale_colour_manual(values = c("A: Impact only" = pal_two[2], "C: Impact + DIF" = pal_two[1])) +
    scale_shape_manual(values = c("A: Impact only" = pal_shapes_two[2], "C: Impact + DIF" = pal_shapes_two[1])) +
    labs(
      title = paste0("Growth-curve impact, before vs after separating DIF (SSP) - ", run_label),
      subtitle = "Effect-coded deviations from grand mean (race); WtHR is a 1-SD effect; posterior mean + 90% CI",
      x = "Effect on growth-curve parameter",
      y = NULL, colour = NULL, shape = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
  ggsave(
    file.path(impact_comparisons_dir, paste0("impact_comparison_A_vs_C_", run_tag, ".png")),
    p_cmp,
    width = 11,
    height = 6.5,
    dpi = 220
  )

  b_mu_logk_draws <- fitC$draws(variables = "b_mu_logk", format = "matrix")
  b_mu_alpha_draws <- fitC$draws(variables = "b_mu_alpha", format = "matrix")

  any_impact_sig <- FALSE
  for (k in seq_len(kimp)) {
    logk_col <- b_mu_logk_draws[, paste0("b_mu_logk[", k, "]")]
    alpha_col <- b_mu_alpha_draws[, paste0("b_mu_alpha[", k, "]")]
    logk_sig <- unname(quantile(logk_col, 0.05) > 0 || quantile(logk_col, 0.95) < 0)
    alpha_sig <- unname(quantile(alpha_col, 0.05) > 0 || quantile(alpha_col, 0.95) < 0)
    if (!logk_sig && !alpha_sig) {
      cat("\nNo significant impact effect for '", covar_labels[k], "' -- skipping trajectory illustration.\n", sep = "")
      next
    }
    any_impact_sig <- TRUE

    sig_label <- if (logk_sig && alpha_sig) {
      "Log-rate AND inflection-age both significant (90% CI excludes 0)"
    } else if (logk_sig) {
      "Log-rate significant; inflection-age NOT distinguishable from average"
    } else {
      "Inflection-age significant; log-rate NOT distinguishable from average"
    }

    if (k < kimp) {
      level_vals <- setNames(c(0, 1), c("All-race average", covar_labels[k]))
      pal_use <- pal_two
      lty_use <- pal_linetypes_two
      legend_title <- NULL
    } else {
      level_vals <- c("Low (-1 SD)" = -1, "Grand mean" = 0, "High (+1 SD)" = 1)
      pal_use <- pal_three_ref
      lty_use <- pal_linetypes_three
      legend_title <- covar_labels[k]
    }
    traj_by_level <- map_dfr(names(level_vals), function(lvl_name) {
      x_val <- level_vals[[lvl_name]]
      k_i_draws <- exp(mu_logk_draws + x_val * logk_col)
      alpha_i_draws <- mu_alpha_draws + x_val * alpha_col
      map_dfr(age_grid, function(a) {
        vals <- 1 + 4 / (1 + exp(-k_i_draws * (a - alpha_i_draws)))
        tibble(age_c = a, eta_med = median(vals))
      }) %>%
        mutate(level = lvl_name)
    }) %>%
      mutate(age = age_c * prep$age_sd + prep$age_mean, level = factor(level, levels = names(level_vals)))

    p_impact <- ggplot(traj_by_level, aes(x = age, y = eta_med, colour = level, linetype = level)) +
      geom_line(linewidth = 1.8) +
      scale_colour_manual(values = setNames(pal_use, names(level_vals))) +
      scale_linetype_manual(values = setNames(lty_use, names(level_vals))) +
      guides(
        colour = guide_legend(override.aes = list(linewidth = 3)),
        linetype = guide_legend(override.aes = list(linewidth = 3))
      ) +
      coord_cartesian(ylim = c(1, 5)) +
      labs(
        title = paste0("Sigmoid growth trajectory by ", covar_labels[k], " (SSP) - ", run_label),
        subtitle = paste0(sig_label, ". Other covariates at reference."),
        x = "Age (years)",
        y = "Latent puberty (eta, 1-5 scale)",
        colour = legend_title,
        linetype = legend_title
      ) +
      theme_minimal(base_size = 13)
    covar_tag <- gsub("[^A-Za-z0-9]+", "", covar_labels[k])
    ggsave(
      file.path(trajectories_dir, paste0("impact_illustration_", covar_tag, "_", run_tag, ".png")),
      p_impact,
      width = 12,
      height = 6.5,
      dpi = 220
    )
    cat(
      "\nImpact illustration for '", covar_labels[k], "': log-rate significant=", logk_sig,
      ", inflection-age significant=", alpha_sig, "\n",
      sep = ""
    )
  }
  if (!any_impact_sig) {
    cat("\nNo significant impact terms for any covariate -- skipping all impact trajectory illustrations.\n")
  }

  l_dif_summ_all <- fitC$summary(variables = "l_dif")
  n_dif_summ_all <- fitC$summary(variables = "n_dif")
  extract_val <- function(summ, prefix, it, k) {
    val <- summ$mean[summ$variable == paste0(prefix, "[", it, ",", k, "]")]
    if (length(val) == 0) 0 else val
  }

  retained_cells <- which(
    (dif_sel$l_pattern == 1) | (dif_sel$n_pattern == 1),
    arr.ind = TRUE
  )
  if (nrow(retained_cells) == 0) {
    cat("\nNo DIF terms retained for any item/covariate -- skipping all DIF illustration plots.\n")
  }
  for (row_i in seq_len(nrow(retained_cells))) {
    it_idx <- retained_cells[row_i, 1]
    cov_idx <- retained_cells[row_i, 2]
    target_item <- prep$item_names[it_idx]
    target_covar <- dif_covar_names[cov_idx]

    lp_mean <- fitC$summary(variables = "lp")$mean[it_idx]
    np_mean <- fitC$summary(variables = "np")$mean[it_idx]
    tau_summ <- fitC$summary(variables = "tau")
    tau1 <- tau_summ$mean[tau_summ$variable == paste0("tau[", it_idx, ",1]")]
    l_dif_val <- extract_val(l_dif_summ_all, "l_dif", it_idx, cov_idx)
    n_dif_val <- extract_val(n_dif_summ_all, "n_dif", it_idx, cov_idx)

    if (target_covar == "informant") {
      level_vals <- c("Youth-reported" = -1, "Caregiver-reported" = 1)
      pal_use <- pal_two
      lty_use <- pal_linetypes_two
    } else {
      level_vals <- c("Low (-1 SD)" = -1, "Grand mean" = 0, "High (+1 SD)" = 1)
      pal_use <- pal_three_ref
      lty_use <- pal_linetypes_three
    }

    eta_grid_true <- seq(1, 5, by = 0.05)
    curve_df <- map_dfr(
      level_vals,
      function(x_val) {
        nu <- np_mean + n_dif_val * x_val
        lam <- lp_mean * exp(l_dif_val * x_val)
        tibble(eta = eta_grid_true, prob = plogis(nu + lam * (eta_grid_true - 3) - tau1))
      },
      .id = "level"
    )
    curve_df$level <- factor(curve_df$level, levels = names(level_vals))

    p_dif <- ggplot(curve_df, aes(x = eta, y = prob, colour = level, linetype = level)) +
      geom_line(linewidth = 1.1) +
      scale_colour_manual(values = setNames(pal_use, levels(curve_df$level))) +
      scale_linetype_manual(values = setNames(lty_use, levels(curve_df$level))) +
      labs(
        title = paste0("Item characteristic curve for '", target_item, "' by ", target_covar, " (SSP) - ", run_label),
        subtitle = "P(response above lowest category) vs. latent puberty (1-5 scale); other DIF covariates held at reference",
        x = "Latent puberty (eta, 1-5 scale)",
        y = "P(response > lowest category)",
        colour = target_covar,
        linetype = target_covar
      ) +
      theme_minimal(base_size = 13)
    ggsave(
      file.path(
        dif_illustrations_dir,
        paste0("dif_illustration_", target_item, "_", target_covar, "_", run_tag, ".png")
      ),
      p_dif,
      width = 10,
      height = 6.5,
      dpi = 220
    )
    cat(
      "\nDIF illustration for '", target_covar, "': item '", target_item,
      "' | loading DIF:", round(l_dif_val, 3), "| intercept DIF:", round(n_dif_val, 3), "\n"
    )
  }

  cat("\nVisualizations written to:", out_dir, "\n")
}

cat("\nAll outputs written to:", out_dir, "\n")
cat("Done:", sx, "(SSP)\n")
