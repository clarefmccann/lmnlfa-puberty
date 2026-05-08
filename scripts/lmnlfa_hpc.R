## lmnlfa_hpc.R
## HPC version of lmnlfa.R — one sex group per job.
##
## Usage:
##   Rscript lmnlfa_hpc.R <sex>
##   sex: female | male
##
## Env vars (set in .sh):
##   DATA_DIR   path to *_long.csv files
##   OUT_DIR    base output directory  (lmnlfa/ subdir created automatically)
##   CMDSTAN    path to CmdStan installation (optional; falls back to default)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(posterior)
})

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("cmdstanr not found. Run: install.packages('cmdstanr', repos=c('https://mc-stan.org/r-packages/', getOption('repos')))")
}
library(cmdstanr)

# Point to a pre-installed CmdStan if the env var is set
cmdstan_env <- Sys.getenv("CMDSTAN")
if (nzchar(cmdstan_env) && dir.exists(cmdstan_env)) {
  set_cmdstan_path(cmdstan_env)
}
cat("CmdStan path:", cmdstan_path(), "\n")

set.seed(90025)

# ---------------------------------------------------------------------------
# ARGUMENTS
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript lmnlfa_hpc.R <female|male>")
sx <- args[1]
if (!sx %in% c("female", "male")) stop("sex must be 'female' or 'male'")
cat("Sex:", sx, "\n")

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
root_path <- Sys.getenv("HOME_DIR")
if (!nzchar(root_path)) root_path <- Sys.getenv("HOME")

data_dir <- Sys.getenv("DATA_DIR")
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  data_dir <- file.path(
    root_path,
    "projects/abcd-projs/abcd-data-release-6.0/cfm/physical-health/puberty"
  )
}
if (!dir.exists(data_dir)) stop("Cannot locate data directory: ", data_dir)

out_base <- Sys.getenv("OUT_DIR")
if (!nzchar(out_base)) {
  out_base <- file.path(
    root_path, "projects/abcd-projs/dissertation/study1/outputs"
  )
}
out_dir <- file.path(out_base, "lmnlfa")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Stan file: SGE copies the job script, so use SGE_O_WORKDIR or a fixed path
script_dir <- Sys.getenv("SGE_O_WORKDIR")
if (!nzchar(script_dir)) script_dir <- "scripts"
stan_file <- file.path(script_dir, "stan", "lmnlfa-quad.stan")
if (!file.exists(stan_file)) {
  stan_file <- file.path("scripts", "stan", "lmnlfa-quad.stan")
}
if (!file.exists(stan_file)) stop("Cannot find lmnlfa-quad.stan: ", stan_file)

cat("Data dir:  ", data_dir, "\n")
cat("Output dir:", out_dir, "\n")
cat("Stan file: ", stan_file, "\n")

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------
parent_df <- read.csv(file.path(data_dir, paste0(sx, "_parent_long.csv")))
youth_df  <- read.csv(file.path(data_dir, paste0(sx, "_youth_long.csv")))

wave_order <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")

# ---------------------------------------------------------------------------
# COMPILE STAN MODEL
# ---------------------------------------------------------------------------
cat("\nCompiling Stan model...\n")
stan_model <- cmdstan_model(stan_file)
cat("Compiled.\n")

# ---------------------------------------------------------------------------
# HELPERS (identical to lmnlfa.R)
# ---------------------------------------------------------------------------
build_lmnlfa_data <- function(parent_df, youth_df, sex_label,
                              ordinal_items = c("peta", "petb", "petc", "petd")) {
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
    filter(if_all(all_of(all_item_cols), ~ !is.na(.) & as.integer(.) %in% 1:4)) %>%
    arrange(id, wave)

  if (nrow(dat) < 500) stop("Insufficient data for ", sex_label)

  ids      <- sort(unique(dat$id))
  age_mean <- mean(dat$age, na.rm = TRUE)

  dat <- dat %>%
    mutate(
      person_idx = match(id, ids),
      time_idx   = as.integer(wave),
      age_c      = age - age_mean,
      age2_c     = age_c^2
    )

  dat_long <- dat %>%
    select(id, person_idx, time_idx, age_c, age2_c, age, all_of(all_item_cols)) %>%
    pivot_longer(cols = all_of(all_item_cols), names_to = "item", values_to = "y_raw") %>%
    filter(!is.na(y_raw)) %>%
    mutate(item_idx = match(item, all_item_cols), y_int = as.integer(y_raw)) %>%
    arrange(person_idx, time_idx, item_idx)

  cat("  n persons:", length(ids),
      "| n items:", length(all_item_cols),
      "| n obs:", nrow(dat_long), "\n")
  cat("  Age range:", round(range(dat$age), 1), "| mean:", round(age_mean, 2), "\n")

  list(
    dat_long   = dat_long,
    item_names = all_item_cols,
    ids        = ids,
    age_mean   = age_mean,
    ni         = length(ids),
    d          = 7L,
    p          = length(all_item_cols),
    nobs       = nrow(dat_long),
    is_binary  = rep(0L, length(all_item_cols)),
    k_items    = rep(4L, length(all_item_cols)),
    k_max      = 4L
  )
}

make_stan_data <- function(prep, ldf = NULL, baseline_age = TRUE) {
  dat <- prep$dat_long
  if (is.null(ldf)) ldf <- matrix(0L, prep$p, 2)
  mtv <- as.integer(sum(ldf[, 1]))
  mf  <- as.integer(sum(ldf[, 2]))

  if (baseline_age) {
    nfpreds <- 1L
    age_by_person <- dat %>%
      group_by(person_idx) %>%
      summarise(mean_age_c = mean(age_c, na.rm = TRUE), .groups = "drop") %>%
      arrange(person_idx)
    xf_person <- matrix(scale(age_by_person$mean_age_c)[, 1], ncol = 1)
    xf        <- matrix(xf_person[dat$person_idx, 1], ncol = 1)
  } else {
    nfpreds   <- 0L
    xf_person <- matrix(0, prep$ni, 0)
    xf        <- matrix(0, prep$nobs, 0)
  }

  list(
    nobs = prep$nobs, p = prep$p, ni = prep$ni, d = prep$d,
    person = dat$person_idx, itm = dat$item_idx, time = dat$time_idx,
    age_c = dat$age_c, age2_c = dat$age2_c, y = dat$y_int,
    is_binary = prep$is_binary, k_item = prep$k_items, k_max = prep$k_max,
    nfpreds = nfpreds, ntvpreds = 0L,
    xf_person = xf_person, xf = xf, xtv = matrix(0, prep$nobs, 0),
    ldf = ldf, mtv = mtv, mf = mf,
    sigma_l = 1.0, sigma_nu = 2.0, sigma_cor = 1.0,
    sigma_f = 1.5, sigma_di = 0.5
  )
}

select_dif <- function(fit1, prep, ci_level = 0.90) {
  p     <- prep$p
  ldf   <- matrix(0L, p, 2)
  alpha <- (1 - ci_level) / 2

  results <- lapply(c("l_diftv", "n_diftv"), function(vname) {
    draws <- tryCatch(fit1$draws(variables = vname, format = "matrix"),
                      error = function(e) NULL)
    if (is.null(draws)) return(rep(FALSE, p))
    apply(draws, 2, function(col) {
      lo <- quantile(col, alpha); hi <- quantile(col, 1 - alpha)
      lo > 0 | hi < 0
    })
  })

  dif_flag     <- Reduce(`|`, results)
  ldf[dif_flag, 1] <- 1L

  cat("\nDIF selection (", ci_level * 100, "% CI):\n", sep = "")
  sel_df <- data.frame(item = prep$item_names,
                       dif_loading   = results[[1]],
                       dif_intercept = results[[2]],
                       selected      = dif_flag)
  print(sel_df, row.names = FALSE)
  list(ldf = ldf, sel_df = sel_df)
}

extract_growth_params <- function(fit2, sex_label) {
  params <- c("mu_slp", "mu_quad", "phi_int", "phi_slp", "eti_sd", "Omega[1,2]")
  draws  <- fit2$draws(variables = params, format = "df")
  posterior::summarise_draws(
    draws, mean, sd,
    ~ quantile(.x, c(0.05, 0.25, 0.50, 0.75, 0.95)),
    posterior::default_convergence_measures()
  ) %>% mutate(sex = sex_label)
}

extract_factor_scores <- function(fit2, prep, ldf_step2, sex_label) {
  draws <- fit2$draws(
    variables = c("mu_slp", "mu_quad", "b_mu", "b_phi",
                  "phi_int", "phi_slp", "L_Omega",
                  "fac_dist", "fac_eti_raw", "eti_sd"),
    format = "df"
  )

  mu_slp  <- mean(draws$mu_slp)
  mu_quad <- mean(draws$mu_quad)
  eti_sd  <- mean(draws$eti_sd)

  fac_dist_mean <- colMeans(as.matrix(draws[, grep("^fac_dist\\[",    names(draws))]))
  fac_eti_mean  <- colMeans(as.matrix(draws[, grep("^fac_eti_raw\\[", names(draws))]))

  ni <- prep$ni
  d  <- prep$d
  fac_dist_mat <- matrix(fac_dist_mean, nrow = 2, ncol = ni)
  fac_eti_mat  <- matrix(fac_eti_mean * eti_sd, nrow = d, ncol = ni)

  b_mu_mean  <- colMeans(as.matrix(draws[, grep("^b_mu\\[",  names(draws)), drop = FALSE]))
  b_phi_mean <- colMeans(as.matrix(draws[, grep("^b_phi\\[", names(draws)), drop = FALSE]))

  age_by_person <- prep$dat_long %>%
    group_by(person_idx) %>%
    summarise(mean_age_c = mean(age_c, na.rm = TRUE), .groups = "drop") %>%
    arrange(person_idx)
  xf_person <- scale(age_by_person$mean_age_c)[, 1]

  L_Omega_mean <- colMeans(as.matrix(draws[, grep("^L_Omega\\[", names(draws))]))
  L_Omega_mat  <- matrix(L_Omega_mean, nrow = 2, ncol = 2)

  phi_eta <- c(mean(draws$phi_int), mean(draws$phi_slp))

  scores <- map_dfr(seq_len(ni), function(k) {
    mu_eta  <- c(0, mu_slp)
    sd_eta  <- phi_eta * exp(b_phi_mean * xf_person[k])
    fac_gr_k <- mu_eta +
      b_mu_mean * xf_person[k] +
      diag(sd_eta) %*% L_Omega_mat %*% fac_dist_mat[, k]

    map_dfr(seq_len(d), function(t) {
      age_obs <- prep$dat_long %>% filter(person_idx == k, time_idx == t) %>% slice(1)
      if (nrow(age_obs) == 0) return(NULL)
      eta_tp <- fac_gr_k[1] +
        fac_gr_k[2] * age_obs$age_c[1] +
        mu_quad    * age_obs$age2_c[1] +
        fac_eti_mat[t, k]
      tibble(person_idx = k, id = prep$ids[k], time_idx = t,
             wave = wave_order[t], age = age_obs$age[1], eta = as.numeric(eta_tp))
    })
  })

  scores %>% mutate(sex = sex_label)
}

# ---------------------------------------------------------------------------
# HELPER: MCMC diagnostics — convergence summary + trace/density plots
# ---------------------------------------------------------------------------
save_diagnostics <- function(fit, label, out_dir) {
  key_params <- c("mu_slp", "mu_quad", "phi_int", "phi_slp", "eti_sd", "Omega[1,2]")

  # Divergence / tree-depth / E-BFMI summary
  diag <- fit$diagnostic_summary(quiet = TRUE)
  cat("\n--- Convergence:", label, "---\n")
  cat("Divergences per chain:", diag$num_divergent, "\n")
  cat("Max treedepth hits:  ", diag$num_max_treedepth, "\n")
  cat("E-BFMI:              ", round(diag$ebfmi, 3), "\n")

  # Full parameter summary → CSV
  summ <- fit$summary()
  write.csv(summ,
            file.path(out_dir, paste0("convergence_summary_", label, ".csv")),
            row.names = FALSE)

  # Worst-Rhat parameters
  bad <- summ[!is.na(summ$rhat) & summ$rhat > 1.01, ]
  if (nrow(bad) > 0) {
    cat("  Parameters with Rhat > 1.01 (top 10):\n")
    print(head(bad[order(-bad$rhat), c("variable", "mean", "rhat", "ess_bulk")], 10),
          digits = 3, row.names = FALSE)
  } else {
    cat("  All Rhat <= 1.01\n")
  }

  # Key parameter summary
  key_summ <- fit$summary(variables = key_params)
  cat("\nKey growth parameters:\n")
  print(key_summ[, c("variable", "mean", "sd", "q5", "q95", "rhat", "ess_bulk")],
        digits = 3, row.names = FALSE)

  # Draws for plotting — use only key params
  draws_long <- tryCatch({
    fit$draws(variables = key_params, format = "df") %>%
      pivot_longer(
        cols      = any_of(key_params),
        names_to  = "parameter",
        values_to = "value"
      )
  }, error = function(e) NULL)

  if (!is.null(draws_long)) {
    # Trace plots
    p_trace <- ggplot(
      draws_long,
      aes(x = .iteration, y = value, colour = factor(.chain))
    ) +
      geom_line(alpha = 0.6, linewidth = 0.25) +
      facet_wrap(~ parameter, scales = "free_y", ncol = 2) +
      scale_colour_brewer(palette = "Set1") +
      labs(title  = paste("Trace plots:", label),
           x = "Iteration", y = "Value", colour = "Chain") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
    ggsave(file.path(out_dir, paste0("trace_", label, ".png")),
           p_trace, width = 10, height = 8, dpi = 150)

    # Posterior density plots (chains overlaid)
    p_dens <- ggplot(
      draws_long,
      aes(x = value, fill = factor(.chain))
    ) +
      geom_density(alpha = 0.35) +
      facet_wrap(~ parameter, scales = "free", ncol = 2) +
      scale_fill_brewer(palette = "Set1") +
      labs(title = paste("Posterior densities:", label),
           x = "Value", y = "Density", fill = "Chain") +
      theme_minimal(base_size = 11) +
      theme(legend.position = "bottom")
    ggsave(file.path(out_dir, paste0("posteriors_", label, ".png")),
           p_dens, width = 10, height = 8, dpi = 150)

    # Scatter: mu_slp vs mu_quad (growth trajectory shape)
    if (all(c("mu_slp", "mu_quad") %in% draws_long$parameter)) {
      scatter_df <- draws_long %>%
        filter(parameter %in% c("mu_slp", "mu_quad")) %>%
        pivot_wider(names_from = parameter, values_from = value)
      p_scatter <- ggplot(scatter_df, aes(x = mu_slp, y = mu_quad,
                                          colour = factor(.chain))) +
        geom_point(alpha = 0.15, size = 0.4) +
        scale_colour_brewer(palette = "Set1") +
        labs(title  = paste("Growth slope vs. quadratic:", label),
             x = "μ_slope", y = "μ_quadratic", colour = "Chain") +
        theme_minimal(base_size = 11) +
        theme(legend.position = "bottom")
      ggsave(file.path(out_dir, paste0("scatter_slp_quad_", label, ".png")),
             p_scatter, width = 6, height = 5, dpi = 150)
    }
  }

  invisible(summ)
}

extract_item_params <- function(fit2, prep, ldf_step2, sex_label,
                                age_grid = seq(-4, 4, by = 0.5)) {
  draws_base  <- fit2$draws(variables = c("lp", "np"), format = "df")
  p           <- prep$p
  items       <- prep$item_names
  lp_mean     <- colMeans(as.matrix(draws_base[, grep("^lp\\[", names(draws_base))]))
  np_mean     <- colMeans(as.matrix(draws_base[, grep("^np\\[", names(draws_base))]))
  ldiftv_mean <- rep(0, p)
  ndiftv_mean <- rep(0, p)
  dif_item_idx <- which(ldf_step2[, 1] == 1)

  if (length(dif_item_idx) > 0) {
    draws_dif   <- fit2$draws(variables = c("l_diftv", "n_diftv"), format = "df")
    ldiftv_cols <- grep("^l_diftv\\[", names(draws_dif))
    ndiftv_cols <- grep("^n_diftv\\[", names(draws_dif))
    if (length(ldiftv_cols) > 0)
      ldiftv_mean[dif_item_idx] <- colMeans(as.matrix(draws_dif[, ldiftv_cols, drop = FALSE]))
    if (length(ndiftv_cols) > 0)
      ndiftv_mean[dif_item_idx] <- colMeans(as.matrix(draws_dif[, ndiftv_cols, drop = FALSE]))
  }

  expand.grid(item_idx = seq_len(p), age_c = age_grid) %>%
    mutate(
      item      = items[item_idx],
      reporter  = ifelse(grepl("_p$", item), "Parent", "Youth"),
      base_item = sub("_(p|y)$", "", item),
      lam       = lp_mean[item_idx] * exp(ldiftv_mean[item_idx] * age_c),
      nu        = np_mean[item_idx] + ndiftv_mean[item_idx] * age_c,
      age       = age_c + prep$age_mean,
      sex       = sex_label
    )
}

# ---------------------------------------------------------------------------
# DATA PREPARATION
# ---------------------------------------------------------------------------
prep <- build_lmnlfa_data(parent_df, youth_df, sx)

# ---------------------------------------------------------------------------
# STEP 1: DIF SCREENING
# ---------------------------------------------------------------------------
rds_fit1 <- file.path(out_dir, paste0("fit1_", sx, ".rds"))

if (file.exists(rds_fit1)) {
  cat("\nLoading cached Step 1 fit:", rds_fit1, "\n")
  fit1 <- readRDS(rds_fit1)
} else {
  cat("\n--- Step 1: DIF screening (all items) ---\n")
  ldf_step1     <- matrix(c(rep(1L, prep$p), rep(0L, prep$p)), ncol = 2)
  stan_data_step1 <- make_stan_data(prep, ldf = ldf_step1)

  t0   <- proc.time()
  fit1 <- stan_model$sample(
    data           = stan_data_step1,
    chains         = 4,
    parallel_chains = 4,
    iter_warmup    = 1000,
    iter_sampling  = 1000,
    adapt_delta    = 0.95,
    init           = "0",
    refresh        = 100,
    show_messages  = TRUE,
    seed           = 90025
  )
  cat("Step 1 elapsed:", round((proc.time() - t0)[["elapsed"]] / 3600, 2), "hr\n")

  diag1 <- fit1$diagnostic_summary(quiet = TRUE)
  cat("Step 1 divergences:", sum(diag1$num_divergent), "\n")
  cat("Step 1 max Rhat:   ", round(max(fit1$summary()$rhat, na.rm = TRUE), 3), "\n")

  fit1$save_object(rds_fit1)
  cat("Step 1 fit saved:", rds_fit1, "\n")
}

save_diagnostics(fit1, paste0("step1_", sx), out_dir)

# DIF selection
dif_result <- select_dif(fit1, prep, ci_level = 0.90)
write.csv(
  dif_result$sel_df %>% mutate(sex = sx),
  file.path(out_dir, paste0("dif_selection_", sx, ".csv")),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# STEP 2: REFIT WITH SELECTED DIF
# ---------------------------------------------------------------------------
rds_fit2 <- file.path(out_dir, paste0("fit2_", sx, ".rds"))

if (file.exists(rds_fit2)) {
  cat("\nLoading cached Step 2 fit:", rds_fit2, "\n")
  fit2 <- readRDS(rds_fit2)
} else {
  cat("\n--- Step 2: Refit with DIF-selected items ---\n")
  ldf_step2     <- dif_result$ldf
  stan_data_step2 <- make_stan_data(prep, ldf = ldf_step2)

  t0   <- proc.time()
  fit2 <- stan_model$sample(
    data           = stan_data_step2,
    chains         = 4,
    parallel_chains = 4,
    iter_warmup    = 1000,
    iter_sampling  = 1000,
    adapt_delta    = 0.95,
    init           = "0",
    refresh        = 100,
    show_messages  = TRUE,
    seed           = 90025
  )
  cat("Step 2 elapsed:", round((proc.time() - t0)[["elapsed"]] / 3600, 2), "hr\n")

  diag2 <- fit2$diagnostic_summary(quiet = TRUE)
  cat("Step 2 divergences:", sum(diag2$num_divergent), "\n")
  cat("Step 2 max Rhat:   ", round(max(fit2$summary()$rhat, na.rm = TRUE), 3), "\n")

  fit2$save_object(rds_fit2)
  cat("Step 2 fit saved:", rds_fit2, "\n")
}

save_diagnostics(fit2, paste0("step2_", sx), out_dir)

ldf_step2 <- dif_result$ldf

# ---------------------------------------------------------------------------
# GROWTH PARAMETERS
# ---------------------------------------------------------------------------
gp <- extract_growth_params(fit2, sx)
cat("\nGrowth parameters [", sx, "]:\n")
print(gp[, c("variable", "mean", "sd", "q5", "q95", "rhat")])
write.csv(gp, file.path(out_dir, paste0("growth_params_", sx, ".csv")),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# FACTOR SCORES
# ---------------------------------------------------------------------------
cat("\nExtracting factor scores...\n")
scores <- tryCatch(
  extract_factor_scores(fit2, prep, ldf_step2, sx),
  error = function(e) { message("Factor scores failed: ", e$message); NULL }
)
if (!is.null(scores)) {
  write.csv(scores,
            file.path(out_dir, paste0("factor_scores_", sx, ".csv")),
            row.names = FALSE)
  cat("Factor scores written:", nrow(scores), "rows\n")

  # Spaghetti plot: sample 200 individuals
  set.seed(90025)
  samp_ids <- sample(unique(scores$person_idx), min(200, length(unique(scores$person_idx))))
  p_scores <- ggplot(
    scores %>% filter(person_idx %in% samp_ids),
    aes(x = age, y = eta, group = id)
  ) +
    geom_line(alpha = 0.15, linewidth = 0.3, colour = "#2166ac") +
    geom_smooth(aes(group = NULL), method = "loess", se = TRUE,
                colour = "black", linewidth = 1.1) +
    labs(title    = paste0("Individual puberty trajectories — ", sx),
         subtitle = "Posterior mean factor scores; n = 200 sampled",
         x = "Age (years)", y = "Latent puberty (η)") +
    theme_minimal(base_size = 13)
  ggsave(file.path(out_dir, paste0("factor_score_trajectories_", sx, ".png")),
         p_scores, width = 8, height = 5, dpi = 180)
}

# ---------------------------------------------------------------------------
# ITEM PARAMETERS + PLOTS
# ---------------------------------------------------------------------------
ip <- extract_item_params(fit2, prep, ldf_step2, sx)
write.csv(ip, file.path(out_dir, paste0("item_params_", sx, ".csv")),
          row.names = FALSE)

p_lam <- ggplot(ip, aes(x = age, y = lam, colour = reporter, linetype = base_item)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(Parent = "#2166ac", Youth = "#d73027")) +
  labs(title = paste0("Item loadings by age — ", sx),
       x = "Age (years)", y = "Loading (λ)",
       colour = "Reporter", linetype = "Item") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, paste0("item_loadings_by_age_", sx, ".png")),
       p_lam, width = 8, height = 5, dpi = 180)

p_nu <- ggplot(ip, aes(x = age, y = nu, colour = reporter, linetype = base_item)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c(Parent = "#2166ac", Youth = "#d73027")) +
  labs(title = paste0("Item intercepts by age — ", sx),
       x = "Age (years)", y = "Intercept (ν)",
       colour = "Reporter", linetype = "Item") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, paste0("item_intercepts_by_age_", sx, ".png")),
       p_nu, width = 8, height = 5, dpi = 180)

# ---------------------------------------------------------------------------
# MEAN GROWTH TRAJECTORY PLOT
# ---------------------------------------------------------------------------
gp_means    <- fit2$draws(variables = c("mu_slp", "mu_quad"), format = "df")
age_grid    <- seq(min(prep$dat_long$age), max(prep$dat_long$age), length.out = 100)
age_c_grid  <- age_grid - prep$age_mean
age2_c_grid <- age_c_grid^2

traj_draws <- map_dfr(seq_len(min(200, nrow(gp_means))), function(i) {
  tibble(age = age_grid,
         eta_mean = gp_means$mu_slp[i] * age_c_grid +
                    gp_means$mu_quad[i] * age2_c_grid,
         draw = i)
})

traj_summ <- traj_draws %>%
  group_by(age) %>%
  summarise(eta_med = median(eta_mean),
            eta_lo  = quantile(eta_mean, 0.05),
            eta_hi  = quantile(eta_mean, 0.95),
            .groups = "drop")

p_traj <- ggplot(traj_summ, aes(x = age)) +
  geom_ribbon(aes(ymin = eta_lo, ymax = eta_hi), fill = "#4393c3", alpha = 0.25) +
  geom_line(aes(y = eta_med), colour = "#2166ac", linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  labs(title = paste0("Mean puberty growth trajectory — ", sx),
       subtitle = "Posterior median + 90% CI",
       x = "Age (years)", y = "Latent puberty (η)") +
  theme_minimal(base_size = 13)
ggsave(file.path(out_dir, paste0("growth_trajectory_", sx, ".png")),
       p_traj, width = 7, height = 5, dpi = 180)

cat("\nAll outputs written to:", out_dir, "\n")
cat("Done:", sx, "\n")
