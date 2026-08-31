## lmnlfa_growth_informant.R
## Next increment beyond the validated growth-only model
## (lmnlfa_growth_only_tanhcor.R): restructures from 8 separate items
## (parent/youth as unrelated items) to 4 shared items with informant
## (parent vs youth) as an explicit DIF covariate -- matching the
## cross-sectional analysis's design and its finding that informant DIF is
## real and pervasive in both sexes.
##
## Growth structure (intercept/slope/correlation via tanh, marker-item
## identification, occasion-specific "wobble" residual, rescaled age) is
## unchanged from the validated growth-only model. Growth-factor impact
## (race, WHtR) and age-varying DIF are deferred to their own later stage.
##
## Data structure: unlike the original 8-item build, which required BOTH
## reporters complete at a wave to keep it, this keeps whichever
## reporter(s) are complete per (person, wave) -- matching the
## cross-sectional model's less restrictive approach, avoiding discarding
## data unnecessarily.
##
## IMPORTANT: this script builds the data and compiles the model, but does
## NOT run the full sampling by default -- set run_fit <- TRUE below first.
## No smoke test was run before this went to HPC (skipped deliberately,
## given how expensive this model class has proven -- growth-only smoke
## tests alone took ~80 min for trivial iteration counts). Watch the real
## run's early log for the known failure signatures instead.
##
## Usage:
##   Rscript lmnlfa_growth_informant.R <sex>
##   sex: female | male
##
## Stan:    scripts/stan/lmnlfa-linear-tanhcor-informant.stan
## Outputs: written to <OUT_DIR>/lmnlfa_growth_informant/

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
# RUN SWITCH -- review before running. Loads from cache if present.
# ---------------------------------------------------------------------------
run_fit <- TRUE

# ---------------------------------------------------------------------------
# ARGUMENTS
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript lmnlfa_growth_informant.R <female|male>")
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
out_dir <- file.path(out_base, "lmnlfa_growth_informant")
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
stan_file <- file.path(
  script_dir,
  "stan",
  "lmnlfa-linear-tanhcor-informant.stan"
)
if (!file.exists(stan_file)) {
  stan_file <- file.path(
    "scripts",
    "stan",
    "lmnlfa-linear-tanhcor-informant.stan"
  )
}
if (!file.exists(stan_file)) {
  stop("Cannot find lmnlfa-linear-tanhcor-informant.stan: ", stan_file)
}

palette_file <- file.path(script_dir, "color_palette.R")
if (!file.exists(palette_file)) palette_file <- file.path("scripts", "color_palette.R")
source(palette_file)

cat("Data dir:  ", data_dir, "\n")
cat("Output dir:", out_dir, "\n")
cat("Stan file: ", stan_file, "\n")

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------
parent_df <- read.csv(file.path(data_dir, paste0(sx, "_parent_long.csv")))
youth_df <- read.csv(file.path(data_dir, paste0(sx, "_youth_long.csv")))

wave_order <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")

# ---------------------------------------------------------------------------
# BUILD LONGITUDINAL DATA: 4 shared items, informant as a DIF covariate,
# keeping whichever reporter(s) are complete per (person, wave)
# ---------------------------------------------------------------------------
build_lmnlfa_data_informant <- function(
  parent_df,
  youth_df,
  sex_label,
  ordinal_items = c("peta", "petb", "petc", "petd"),
  n_subsample = NULL
) {
  cat("\n=== Building longitudinal data (informant DIF) |", sex_label, "===\n")

  clean_reporter <- function(df, reporter_label, informant_val) {
    df %>%
      select(id, wave, age, all_of(ordinal_items)) %>%
      filter(!is.na(age)) %>%
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

  # Optional person-level subsample, e.g. for a diagnostic run under a
  # walltime ceiling: since each person contributes a growth factor pair, a
  # non-centered random-effect pair, and a wobble term per wave (~11
  # parameters/person), cutting the number of people cuts both the
  # dominant parameter count AND the observation count roughly
  # proportionally -- unlike most other levers, which only touch one or
  # the other. Not the final analysis sample; just enough N to check the
  # model converges before committing to the full-sample run.
  all_ids <- sort(unique(dat$id))
  if (!is.null(n_subsample) && n_subsample < length(all_ids)) {
    sub_ids <- sort(sample(all_ids, n_subsample))
    dat <- dat %>% filter(id %in% sub_ids)
    cat(
      "  Subsampled to", n_subsample, "of", length(all_ids), "people (",
      round(100 * n_subsample / length(all_ids), 1), "% )\n"
    )
  }

  ids <- sort(unique(dat$id))
  age_mean <- mean(dat$age, na.rm = TRUE)
  # rescale age to SD units rather than raw years -- see lmnlfa_growth_only_tanhcor.R
  age_sd <- sd(dat$age, na.rm = TRUE)

  dat <- dat %>%
    mutate(
      person_idx = match(id, ids),
      time_idx = as.integer(wave),
      age_c = (age - age_mean) / age_sd
    )

  dat_long <- dat %>%
    select(
      id,
      person_idx,
      time_idx,
      age_c,
      age,
      informant_c,
      reporter,
      all_of(ordinal_items)
    ) %>%
    pivot_longer(
      cols = all_of(ordinal_items),
      names_to = "item",
      values_to = "y_raw"
    ) %>%
    mutate(item_idx = match(item, ordinal_items), y_int = as.integer(y_raw)) %>%
    arrange(person_idx, time_idx, reporter, item_idx)

  cat(
    "  n persons:",
    length(ids),
    "| n items:",
    length(ordinal_items),
    "| n obs:",
    nrow(dat_long),
    "\n"
  )
  cat(
    "  Reporter rows -- parent:",
    sum(dat$reporter == "parent"),
    "| youth:",
    sum(dat$reporter == "youth"),
    "\n"
  )
  cat(
    "  Age range:",
    round(range(dat$age), 1),
    "| mean:",
    round(age_mean, 2),
    "| sd:",
    round(age_sd, 2),
    "\n"
  )

  list(
    dat_long = dat_long,
    item_names = ordinal_items,
    ids = ids,
    age_mean = age_mean,
    age_sd = age_sd,
    ni = length(ids),
    d = 7L,
    p = length(ordinal_items),
    nobs = nrow(dat_long),
    is_binary = rep(0L, length(ordinal_items)),
    k_items = rep(4L, length(ordinal_items)),
    k_max = 4L
  )
}

prep <- build_lmnlfa_data_informant(parent_df, youth_df, sx, n_subsample = 1500)

stan_data <- list(
  nobs = prep$nobs,
  p = prep$p,
  ni = prep$ni,
  d = prep$d,
  person = prep$dat_long$person_idx,
  itm = prep$dat_long$item_idx,
  time = prep$dat_long$time_idx,
  age_c = prep$dat_long$age_c,
  informant_c = prep$dat_long$informant_c,
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

# ---------------------------------------------------------------------------
# COMPILE STAN MODEL (isolated exe -- see lmnlfa_growth_only_tanhcor.R for why)
# ---------------------------------------------------------------------------
bin_dir <- file.path(out_dir, "bin")
dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)
cat("\nCompiling Stan model (local exe)...\n")
stan_model <- cmdstan_model(
  stan_file,
  exe_file = file.path(bin_dir, "lmnlfa-linear-tanhcor-informant-local"),
  force_recompile = TRUE
)
cat("Compiled.\n")

saveRDS(stan_data, file.path(out_dir, paste0("stan_data_", sx, ".rds")))

# ---------------------------------------------------------------------------
# FIT
# ---------------------------------------------------------------------------
rds_fit <- file.path(out_dir, paste0("fit_", sx, ".rds"))

if (file.exists(rds_fit)) {
  cat("\nLoading cached fit:", rds_fit, "\n")
  fit <- readRDS(rds_fit)
} else if (run_fit) {
  cat("\n--- Fitting growth + informant DIF model ---\n")
  t0 <- proc.time()
  fit <- stan_model$sample(
    data = stan_data,
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
  cat("Elapsed:", round((proc.time() - t0)[["elapsed"]] / 60, 1), "min\n")
  fit$save_object(rds_fit)
  diag <- fit$diagnostic_summary(quiet = TRUE)
  cat(
    "Divergences:",
    sum(diag$num_divergent),
    " max treedepth hits:",
    sum(diag$num_max_treedepth),
    "\n"
  )
  print(
    fit$summary(variables = c("mu_slp", "phi_int", "phi_slp", "rho", "eti_sd")),
    digits = 3
  )
} else {
  stop(
    "Not yet run and no cached fit found. Set run_fit <- TRUE at the top ",
    "of this script and re-run when you're ready to sample (expect this ",
    "to be expensive -- see the file header)."
  )
}

# ---------------------------------------------------------------------------
# POST-HOC DIAGNOSTICS + VISUALIZATION
# ---------------------------------------------------------------------------
diag <- fit$diagnostic_summary(quiet = TRUE)
cat(
  "\nDivergences:",
  sum(diag$num_divergent),
  " max treedepth hits:",
  sum(diag$num_max_treedepth),
  "\n"
)

growth_summ <- fit$summary(
  variables = c("mu_slp", "phi_int", "phi_slp", "rho", "eti_sd")
)
cat("\nGrowth parameters:\n")
print(
  growth_summ[, c("variable", "mean", "sd", "q5", "q95", "rhat", "ess_bulk")],
  digits = 3
)
write.csv(
  growth_summ,
  file.path(out_dir, paste0("growth_params_", sx, ".csv")),
  row.names = FALSE
)

# --- (1) Factor scores: fac_gr is saved directly (transformed parameters),
# no reconstruction needed, unlike the cross-sectional model's local eta.
fac_gr_summ <- fit$summary(variables = "fac_gr")
# fac_gr[1,k] = intercept, fac_gr[2,k] = slope, for person k
intercepts <- fac_gr_summ$mean[grepl("^fac_gr\\[1,", fac_gr_summ$variable)]
slopes <- fac_gr_summ$mean[grepl("^fac_gr\\[2,", fac_gr_summ$variable)]

scores <- tibble(
  person_idx = seq_len(prep$ni),
  id = prep$ids,
  intercept = intercepts,
  slope = slopes
)
write.csv(
  scores,
  file.path(out_dir, paste0("growth_factor_scores_", sx, ".csv")),
  row.names = FALSE
)
cat("\nGrowth factor scores written:", nrow(scores), "rows\n")

# per-person predicted trajectory (excludes occasion-specific wobble --
# that's measurement noise, not part of the smooth underlying trajectory)
age_grid <- seq(-3, 3, by = 0.1)
set.seed(90025)
samp_idx <- sample(seq_len(prep$ni), min(200, prep$ni))
traj_df <- map_dfr(samp_idx, function(k) {
  tibble(
    person_idx = k,
    age_c = age_grid,
    age = age_grid * prep$age_sd + prep$age_mean,
    eta = intercepts[k] + slopes[k] * age_grid
  )
})

p_spaghetti <- ggplot(traj_df, aes(x = age, y = eta, group = person_idx)) +
  geom_line(alpha = 0.15, linewidth = 0.3, colour = pal_primary) +
  labs(
    title = paste0("Individual puberty trajectories - ", sx),
    subtitle = paste0(
      "Predicted from growth factors (n = ",
      length(samp_idx),
      " sampled); excludes occasion-specific wobble"
    ),
    x = "Age (years)",
    y = "Latent puberty (eta)"
  ) +
  theme_minimal(base_size = 13)
ggsave(
  file.path(out_dir, paste0("trajectories_sample_", sx, ".png")),
  p_spaghetti,
  width = 8,
  height = 5,
  dpi = 180
)

# population mean trajectory + 90% CI, from posterior draws directly
mu_slp_draws <- fit$draws(variables = "mu_slp", format = "matrix")
mean_traj <- map_dfr(age_grid, function(a) {
  vals <- as.vector(mu_slp_draws) * a # intercept mean fixed at 0 for identification
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
    title = paste0("Mean puberty growth trajectory - ", sx),
    subtitle = "Posterior median + 90% CI",
    x = "Age (years)",
    y = "Latent puberty (eta)"
  ) +
  theme_minimal(base_size = 13)
ggsave(
  file.path(out_dir, paste0("mean_trajectory_", sx, ".png")),
  p_meantraj,
  width = 7,
  height = 5,
  dpi = 180
)

# --- (2) Informant DIF illustration: item characteristic curve for the
# item with the largest informant DIF, parent vs youth
l_dif_summ <- fit$summary(variables = "l_dif_informant")
n_dif_summ <- fit$summary(variables = "n_dif_informant")
target_idx <- which.max(abs(l_dif_summ$mean))
target_item <- prep$item_names[target_idx]

lp_mean <- fit$summary(variables = "lp")$mean[target_idx]
np_mean <- fit$summary(variables = "np")$mean[target_idx]
tau_summ <- fit$summary(variables = "tau")
tau1 <- tau_summ$mean[tau_summ$variable == paste0("tau[", target_idx, ",1]")]
l_dif_val <- l_dif_summ$mean[target_idx]
n_dif_val <- n_dif_summ$mean[target_idx]

eta_grid2 <- seq(-3, 3, by = 0.1)
curve_df <- map_dfr(
  c("Parent" = 1, "Youth" = -1),
  function(inf_val) {
    nu <- np_mean + n_dif_val * inf_val
    lam <- lp_mean * exp(l_dif_val * inf_val)
    tibble(eta = eta_grid2, prob = plogis(nu + lam * eta_grid2 - tau1))
  },
  .id = "informant"
)

p_dif <- ggplot(curve_df, aes(x = eta, y = prob, colour = informant, linetype = informant)) +
  geom_line(linewidth = 1.1) +
  scale_colour_manual(values = c("Parent" = pal_two[1], "Youth" = pal_two[2])) +
  scale_linetype_manual(values = c("Parent" = pal_linetypes_two[1], "Youth" = pal_linetypes_two[2])) +
  labs(
    title = paste0(
      "Item characteristic curve for '",
      target_item,
      "' by informant - ",
      sx
    ),
    subtitle = "P(response above lowest category) vs. latent puberty",
    x = "Latent puberty (eta)",
    y = "P(response > lowest category)",
    colour = "Informant",
    linetype = "Informant"
  ) +
  theme_minimal(base_size = 13)
ggsave(
  file.path(out_dir, paste0("dif_illustration_informant_", sx, ".png")),
  p_dif,
  width = 7.5,
  height = 5,
  dpi = 180
)

cat(
  "\nLargest informant DIF item:",
  target_item,
  "| loading DIF:",
  round(l_dif_val, 3),
  "| intercept DIF:",
  round(n_dif_val, 3),
  "\n"
)

cat("\nAll outputs written to:", out_dir, "\n")
cat("Done:", sx, "\n")
