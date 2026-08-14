## lmnlfa_growth_only_tanhcor.R
## Same isolation test as lmnlfa_growth_only.R (full longitudinal item data,
## random intercept + random slope growth factor, NO DIF, NO covariate
## impact) but using lmnlfa-quad-tanhcor.stan, which fixes two issues found
## via smoke testing:
##   1. cholesky_factor_corr[2] + lkj_corr_cholesky replaced with an
##      unconstrained z_cor mapped through tanh() -- avoids the exact
##      Cholesky-diagonal-hits-zero boundary that caused frequent
##      lkj_corr_cholesky rejections and a chain crash within 25 warmup
##      iterations.
##   2. Item loading lp[1] fixed to 1 (marker item) instead of letting all p
##      loadings float freely against freely-estimated phi_int/phi_slp --
##      the static model's converged fit showed cor(lp, phi0) ~ -0.92, i.e.
##      loadings and factor SD are only identified through their product.
##      With two SDs (phi_int, phi_slp) riding that same ridge plus their
##      correlation, this was the likely root cause of the extreme
##      threshold values and 100% max-treedepth saturation seen even after
##      fix #1.
##
## Usage (local):
##   Rscript lmnlfa_growth_only_tanhcor.R <sex>
##   sex: female | male
##
## Usage (Hoffman2 SGE): see lmnlfa_growth_only_tanhcor_hpc.sh
##
## Stan:    scripts/stan/lmnlfa-quad-tanhcor.stan
## Outputs: written to <OUT_DIR>/lmnlfa_growth_only_tanhcor/

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

# Point to a pre-installed CmdStan if the env var is set (as on Hoffman2 --
# see lmnlfa_hpc.sh)
cmdstan_env <- Sys.getenv("CMDSTAN")
if (nzchar(cmdstan_env) && dir.exists(cmdstan_env)) {
  set_cmdstan_path(cmdstan_env)
}
cat("CmdStan path:", cmdstan_path(), "\n")
options(mc.cores = as.integer(Sys.getenv("NSLOTS", unset = "4")))
set.seed(90025)

# ---------------------------------------------------------------------------
# ARGUMENTS
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript lmnlfa_growth_only_tanhcor.R <female|male>")
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
out_dir <- file.path(out_base, "lmnlfa_growth_only_tanhcor")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# SGE copies the job script, so use SGE_O_WORKDIR (see lmnlfa_hpc.R); fall
# back to the calling frame's file path for local/interactive runs.
script_dir <- Sys.getenv("SGE_O_WORKDIR")
if (!nzchar(script_dir)) {
  script_dir <- tryCatch(
    {
      ofile <- sys.frames()[[1]]$ofile
      if (is.null(ofile)) "scripts" else dirname(ofile)
    },
    error = function(e) "scripts"
  )
}
stan_file <- file.path(script_dir, "stan", "lmnlfa-quad-tanhcor.stan")
if (!file.exists(stan_file)) {
  stan_file <- file.path("scripts", "stan", "lmnlfa-quad-tanhcor.stan")
}
if (!file.exists(stan_file)) {
  stop("Cannot find lmnlfa-quad-tanhcor.stan: ", stan_file)
}

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
# COMPILE STAN MODEL
# ---------------------------------------------------------------------------
# NOTE: scripts/stan/ is shared (sshfs) with the HPC cluster. This is a new
# .stan file the cluster has never compiled, but keep the same isolated-exe
# convention for consistency and to avoid any risk of cross-platform clashes.
bin_dir <- file.path(out_dir, "bin")
dir.create(bin_dir, showWarnings = FALSE, recursive = TRUE)
cat("\nCompiling Stan model (local exe)...\n")
stan_model <- cmdstan_model(
  stan_file,
  exe_file = file.path(bin_dir, "lmnlfa-quad-tanhcor-local"),
  force_recompile = TRUE
)
cat("Compiled.\n")

# ---------------------------------------------------------------------------
# HELPERS (identical to lmnlfa_growth_only.R)
# ---------------------------------------------------------------------------
build_lmnlfa_data <- function(
  parent_df,
  youth_df,
  sex_label,
  ordinal_items = c("peta", "petb", "petc", "petd")
) {
  cat("\n=== Building LMNLFA data |", sex_label, "===\n")

  parent_sel <- parent_df %>%
    select(id, wave, age, any_of(ordinal_items)) %>%
    rename_with(~ paste0(., "_p"), any_of(ordinal_items))

  youth_sel <- youth_df %>%
    select(id, wave, any_of(ordinal_items)) %>%
    rename_with(~ paste0(., "_y"), any_of(ordinal_items))

  all_item_cols <- c(paste0(ordinal_items, "_p"), paste0(ordinal_items, "_y"))

  dat <- inner_join(parent_sel, youth_sel, by = c("id", "wave")) %>%
    filter(!is.na(age)) %>%
    mutate(wave = factor(wave, levels = wave_order)) %>%
    filter(if_all(
      all_of(all_item_cols),
      ~ !is.na(.) & as.integer(.) %in% 1:4
    )) %>%
    arrange(id, wave)

  if (nrow(dat) < 500) {
    stop("Insufficient data for ", sex_label)
  }

  ids <- sort(unique(dat$id))
  age_mean <- mean(dat$age, na.rm = TRUE)
  # Rescale age_c to a smaller SD instead of leaving it in raw years. This is
  # the time metric multiplied into phi_slp/fac_gr[2,] and the DIF terms; per
  # the source paper (met0000685_sm3, p.29) leaving a wide-range covariate
  # unscaled there is a documented cause of non-convergence, since priors are
  # calibrated assuming a compact predictor scale. Their example divides by a
  # fixed constant (~6); we divide by the empirical SD of age instead, which
  # compresses our age_c range from about [-4.2, 4.8] (raw years) to about
  # [-2.1, 2.4] (SD units).
  age_sd <- sd(dat$age, na.rm = TRUE)

  dat <- dat %>%
    mutate(
      person_idx = match(id, ids),
      time_idx = as.integer(wave),
      age_c = (age - age_mean) / age_sd,
      age2_c = age_c^2
    )

  dat_long <- dat %>%
    select(
      id,
      person_idx,
      time_idx,
      age_c,
      age2_c,
      age,
      all_of(all_item_cols)
    ) %>%
    pivot_longer(
      cols = all_of(all_item_cols),
      names_to = "item",
      values_to = "y_raw"
    ) %>%
    filter(!is.na(y_raw)) %>%
    mutate(item_idx = match(item, all_item_cols), y_int = as.integer(y_raw)) %>%
    arrange(person_idx, time_idx, item_idx)

  cat(
    "  n persons:",
    length(ids),
    "| n items:",
    length(all_item_cols),
    "| n obs:",
    nrow(dat_long),
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
    item_names = all_item_cols,
    ids = ids,
    age_mean = age_mean,
    age_sd = age_sd,
    ni = length(ids),
    d = 7L,
    p = length(all_item_cols),
    nobs = nrow(dat_long),
    is_binary = rep(0L, length(all_item_cols)),
    k_items = rep(4L, length(all_item_cols)),
    k_max = 4L
  )
}

make_stan_data <- function(prep, ldf = NULL, baseline_age = TRUE) {
  dat <- prep$dat_long
  if (is.null(ldf)) {
    ldf <- matrix(0L, prep$p, 2)
  }
  mtv <- as.integer(sum(ldf[, 1]))
  mf <- as.integer(sum(ldf[, 2]))

  if (baseline_age) {
    nfpreds <- 1L
    age_by_person <- dat %>%
      group_by(person_idx) %>%
      summarise(mean_age_c = mean(age_c, na.rm = TRUE), .groups = "drop") %>%
      arrange(person_idx)
    xf_person <- matrix(scale(age_by_person$mean_age_c)[, 1], ncol = 1)
    xf <- matrix(xf_person[dat$person_idx, 1], ncol = 1)
  } else {
    nfpreds <- 0L
    xf_person <- matrix(0, prep$ni, 0)
    xf <- matrix(0, prep$nobs, 0)
  }

  list(
    nobs = prep$nobs,
    p = prep$p,
    ni = prep$ni,
    d = prep$d,
    person = dat$person_idx,
    itm = dat$item_idx,
    time = dat$time_idx,
    age_c = dat$age_c,
    age2_c = dat$age2_c,
    y = dat$y_int,
    is_binary = prep$is_binary,
    k_item = prep$k_items,
    k_max = prep$k_max,
    nfpreds = nfpreds,
    ntvpreds = 0L,
    xf_person = xf_person,
    xf = xf,
    xtv = matrix(0, prep$nobs, 0),
    ldf = ldf,
    mtv = mtv,
    mf = mf,
    sigma_l = 1.5,
    sigma_nu = 1.5,
    sigma_cor = 1.0,
    sigma_f = 1.5,
    sigma_di = 0.2
  )
}

extract_growth_params <- function(fit, sex_label) {
  params <- c("mu_slp", "phi_int", "phi_slp", "eti_sd", "rho", "Omega[1,2]")
  draws <- fit$draws(variables = params, format = "df")
  posterior::summarise_draws(
    draws,
    mean,
    sd,
    ~ quantile(.x, c(0.05, 0.25, 0.50, 0.75, 0.95)),
    posterior::default_convergence_measures()
  ) %>%
    mutate(sex = sex_label)
}

save_diagnostics <- function(fit, label, out_dir) {
  key_params <- c("mu_slp", "phi_int", "phi_slp", "eti_sd", "rho", "Omega[1,2]")

  diag <- fit$diagnostic_summary(quiet = TRUE)
  cat("\n--- Convergence:", label, "---\n")
  cat("Divergences per chain:", diag$num_divergent, "\n")
  cat("Max treedepth hits:  ", diag$num_max_treedepth, "\n")
  cat("E-BFMI:              ", round(diag$ebfmi, 3), "\n")

  summ <- fit$summary()
  write.csv(
    summ,
    file.path(out_dir, paste0("convergence_summary_", label, ".csv")),
    row.names = FALSE
  )

  bad <- summ[!is.na(summ$rhat) & summ$rhat > 1.01, ]
  if (nrow(bad) > 0) {
    cat("  Parameters with Rhat > 1.01 (top 10):\n")
    print(
      head(
        bad[order(-bad$rhat), c("variable", "mean", "rhat", "ess_bulk")],
        10
      ),
      digits = 3,
      row.names = FALSE
    )
  } else {
    cat("  All Rhat <= 1.01\n")
  }

  key_summ <- fit$summary(variables = key_params)
  cat("\nKey growth parameters:\n")
  print(
    key_summ[, c("variable", "mean", "sd", "q5", "q95", "rhat", "ess_bulk")],
    digits = 3,
    row.names = FALSE
  )

  draws_long <- tryCatch(
    {
      fit$draws(variables = key_params, format = "df") %>%
        pivot_longer(
          cols = any_of(key_params),
          names_to = "parameter",
          values_to = "value"
        )
    },
    error = function(e) NULL
  )

  if (!is.null(draws_long)) {
    p_trace <- ggplot(
      draws_long,
      aes(x = .iteration, y = value, colour = factor(.chain))
    ) +
      geom_line(alpha = 0.6, linewidth = 0.25) +
      facet_wrap(~parameter, scales = "free_y", ncol = 2) +
      scale_colour_brewer(palette = "Set1") +
      labs(
        title = paste("Trace plots:", label),
        x = "Iteration",
        y = "Value",
        colour = "Chain"
      ) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
    ggsave(
      file.path(out_dir, paste0("trace_", label, ".png")),
      p_trace,
      width = 10,
      height = 8,
      dpi = 150
    )

    p_dens <- ggplot(
      draws_long,
      aes(x = value, fill = factor(.chain))
    ) +
      geom_density(alpha = 0.35) +
      facet_wrap(~parameter, scales = "free", ncol = 2) +
      scale_fill_brewer(palette = "Set1") +
      labs(
        title = paste("Posterior densities:", label),
        x = "Value",
        y = "Density",
        fill = "Chain"
      ) +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
    ggsave(
      file.path(out_dir, paste0("posteriors_", label, ".png")),
      p_dens,
      width = 10,
      height = 8,
      dpi = 150
    )
  }

  invisible(summ)
}

# ---------------------------------------------------------------------------
# DATA PREPARATION
# ---------------------------------------------------------------------------
prep <- build_lmnlfa_data(parent_df, youth_df, sx)

# No DIF (ldf all zero), no covariate impact (baseline_age = FALSE)
ldf_none <- matrix(0L, prep$p, 2)
stan_data <- make_stan_data(prep, ldf = ldf_none, baseline_age = FALSE)

# ---------------------------------------------------------------------------
# FIT
# ---------------------------------------------------------------------------
cat("\n--- Fitting growth-only model, tanh(z_cor) correlation ---\n")
t0 <- proc.time()
fit <- tryCatch(
  {
    stan_model$sample(
      data = stan_data,
      chains = 4,
      parallel_chains = 4,
      iter_warmup = 1000,
      iter_sampling = 1000,
      adapt_delta = 0.95,
      init = 0.1,
      output_dir = out_dir,
      refresh = 100,
      show_messages = TRUE,
      seed = 90025
    )
  },
  error = function(e) {
    cat("\n!!! Sampling failed:", conditionMessage(e), "\n")
    stop(e)
  }
)
cat("Elapsed:", round((proc.time() - t0)[["elapsed"]] / 60, 1), "min\n")

fit$save_object(file.path(out_dir, paste0("fit_growth_only_tanhcor_", sx, ".rds")))

save_diagnostics(fit, paste0("growth_only_tanhcor_", sx), out_dir)

# ---------------------------------------------------------------------------
# GROWTH PARAMETERS
# ---------------------------------------------------------------------------
gp <- extract_growth_params(fit, sx)
cat("\nGrowth parameters [", sx, "]:\n")
print(gp[, c("variable", "mean", "sd", "q5", "q95", "rhat")])
write.csv(
  gp,
  file.path(out_dir, paste0("growth_params_", sx, ".csv")),
  row.names = FALSE
)

cat("\nAll outputs written to:", out_dir, "\n")
cat("Done:", sx, "\n")
