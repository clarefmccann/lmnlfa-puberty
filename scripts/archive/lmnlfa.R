## lmnlfa.R
## Longitudinal Moderated Nonlinear Factor Analysis (LMNLFA)
## Two-step workflow following Chen & Bauer (2024, Psych Methods).
##
## Measurement model: PDS items from parent AND youth reporters, ordinal 1–4,
##   loading on a single puberty latent factor (justified by high reporter
##   correlation from M2 cross-reporter CFA in 02_psychometrics.R).
##
## Growth model: Quadratic latent growth in puberty as a function of age
##   (centered). Baseline age is a time-invariant predictor of growth factor
##   distribution. Analyses are sex-stratified (females / males run separately).
##
## Two-step DIF workflow:
##   Step 1 — All items have potential time-varying DIF by age.
##             Identify significant DIF via 90% posterior CI excluding 0.
##   Step 2 — Refit with only identified DIF items, weakly informative priors.
##
## Inputs:   female/male parent and youth long CSVs (from 00_data_foundation.R)
## Stan:     scripts/stan/lmnlfa-quad.stan
## Outputs:  dif_selection.csv, growth_params.csv, factor_scores.csv,
##           item_params_by_age_*.png, growth_trajectories_*.png

pacman::p_load(
  dplyr,
  tidyr,
  tibble,
  purrr,
  ggplot2,
  posterior,
  install = TRUE
)
# cmdstanr is loaded separately to allow graceful failure
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop(
    "cmdstanr is required. Install with:\n",
    "  install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))\n",
    "  cmdstanr::install_cmdstan()"
  )
}
library(cmdstanr)

set.seed(90025)

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
    "projects/abcd-projs/abcd-data-release-6.0/cfm/physical-health/puberty"
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
out_dir <- file.path(out_base, "lmnlfa")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

script_dir <- tryCatch(
  {
    ofile <- sys.frames()[[1]]$ofile
    if (is.null(ofile)) "scripts" else dirname(ofile)
  },
  error = function(e) "scripts"
)
stan_file <- file.path(script_dir, "stan", "lmnlfa-quad.stan")
if (!file.exists(stan_file)) {
  stan_file <- file.path("scripts", "stan", "lmnlfa-quad.stan")
}
if (!file.exists(stan_file)) {
  stop("Cannot find lmnlfa-quad.stan: ", stan_file)
}

# ---------------------------------------------------------------------------
# LOAD
# ---------------------------------------------------------------------------
female_parent <- read.csv(file.path(data_dir, "female_parent_long.csv"))
female_youth <- read.csv(file.path(data_dir, "female_youth_long.csv"))
male_parent <- read.csv(file.path(data_dir, "male_parent_long.csv"))
male_youth <- read.csv(file.path(data_dir, "male_youth_long.csv"))

wave_order <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")

# ---------------------------------------------------------------------------
# STAN MODEL COMPILATION
# ---------------------------------------------------------------------------
stan_model <- cmdstan_model(stan_file)

# ---------------------------------------------------------------------------
# HELPER: build long-format data for Stan
# ---------------------------------------------------------------------------
# Returns a list with:
#   dat_long    — tidy long df (one row per person × wave × item)
#   item_names  — ordered item labels (length p)
#   ids         — ordered person IDs (length ni)
#   age_mean    — centering value
#   ni, d, p, nobs, is_binary, k_items, k_max
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

  item_cols <- paste0(ordinal_items, "_p")
  item_cols_y <- paste0(ordinal_items, "_y")
  all_item_cols <- c(item_cols, item_cols_y)

  dat <- inner_join(parent_sel, youth_sel, by = c("id", "wave")) %>%
    filter(!is.na(age)) %>%
    mutate(wave = factor(wave, levels = wave_order)) %>%
    # keep rows where all ordinal items are valid
    filter(if_all(
      all_of(all_item_cols),
      ~ !is.na(.) & as.integer(.) %in% 1:4
    )) %>%
    arrange(id, wave)

  if (nrow(dat) < 500) {
    stop("Insufficient data for ", sex_label)
  }

  # Integer indices
  ids <- sort(unique(dat$id))
  age_mean <- mean(dat$age, na.rm = TRUE)

  dat <- dat %>%
    mutate(
      person_idx = match(id, ids),
      time_idx = as.integer(wave),
      age_c = age - age_mean,
      age2_c = age_c^2
    )

  # Pivot to long (one row per obs)
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
    mutate(
      item_idx = match(item, all_item_cols),
      y_int = as.integer(y_raw)
    ) %>%
    arrange(person_idx, time_idx, item_idx)

  # Item properties — all ordinal 1–4 for the base 4-item PDS set
  is_bin <- rep(0L, length(all_item_cols))
  k_items <- rep(4L, length(all_item_cols))
  k_max <- 4L

  cat(
    "  n persons:",
    length(ids),
    "| n waves: 7 | n items:",
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
    "\n"
  )

  list(
    dat_long = dat_long,
    item_names = all_item_cols,
    ids = ids,
    age_mean = age_mean,
    ni = length(ids),
    d = 7L,
    p = length(all_item_cols),
    nobs = nrow(dat_long),
    is_binary = is_bin,
    k_items = k_items,
    k_max = k_max
  )
}

# ---------------------------------------------------------------------------
# HELPER: assemble Stan data list
# ---------------------------------------------------------------------------
# ldf: p × 2 matrix; col 1 = time-varying DIF (age), col 2 = invariant DIF
# baseline_age: if TRUE, include person-mean age as time-invariant predictor
make_stan_data <- function(prep, ldf = NULL, baseline_age = TRUE) {
  dat <- prep$dat_long

  if (is.null(ldf)) {
    ldf <- matrix(0L, prep$p, 2)
  }

  mtv <- as.integer(sum(ldf[, 1]))
  mf <- as.integer(sum(ldf[, 2]))

  # Time-invariant predictor: baseline age (centered, scaled)
  if (baseline_age) {
    nfpreds <- 1L
    # Person-level mean age (one value per person)
    age_by_person <- dat %>%
      group_by(person_idx) %>%
      summarise(mean_age_c = mean(age_c, na.rm = TRUE), .groups = "drop") %>%
      arrange(person_idx)
    xf_person <- matrix(scale(age_by_person$mean_age_c)[, 1], ncol = 1)
    # obs-level version (same value repeated for all obs of same person)
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
    # prior scales (Chen & Bauer 2024 defaults)
    sigma_l = 1.0,
    sigma_nu = 2.0,
    sigma_cor = 1.0,
    sigma_f = 1.5,
    sigma_di = 0.5
  )
}

# ---------------------------------------------------------------------------
# HELPER: fit Stan model
# ---------------------------------------------------------------------------
fit_stan <- function(
  stan_data,
  model,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  label = ""
) {
  cat("\nFitting Stan model", label, "...\n")
  # init = "0" starts all parameters at 0 in unconstrained space, which maps
  # cholesky_factor_corr to the identity (zero correlation) — avoids the
  # boundary-hit warnings from random initialization
  model$sample(
    data = stan_data,
    chains = chains,
    parallel_chains = min(chains, parallel::detectCores() - 1),
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = 0.95,
    init = "0",
    refresh = 200,
    show_messages = TRUE,
    seed = 90025
  )
}

# ---------------------------------------------------------------------------
# HELPER: DIF selection from Step 1 posterior
# ---------------------------------------------------------------------------
# Returns updated ldf matrix with items flagged where 90% CI excludes 0
select_dif <- function(fit1, prep, ci_level = 0.90) {
  p <- prep$p
  ldf <- matrix(0L, p, 2) # col1 = time-varying, col2 = invariant (not used here)

  # Extract time-varying DIF posteriors for loadings and intercepts
  vars_to_check <- c("l_diftv", "n_diftv")

  # Since Step 1 sets all items with time-varying DIF, mtv = p
  # l_diftv[i] and n_diftv[i] correspond to item i
  alpha <- (1 - ci_level) / 2

  results <- lapply(vars_to_check, function(vname) {
    draws <- tryCatch(
      fit1$draws(variables = vname, format = "matrix"),
      error = function(e) NULL
    )
    if (is.null(draws)) {
      return(rep(FALSE, p))
    }
    # columns are [1]...[p]
    apply(draws, 2, function(col) {
      lo <- quantile(col, alpha)
      hi <- quantile(col, 1 - alpha)
      lo > 0 | hi < 0 # CI excludes 0
    })
  })

  # Flag item if EITHER loading or intercept DIF is significant
  dif_flag <- Reduce(`|`, results)

  ldf[dif_flag, 1] <- 1L

  cat("\nDIF selection (", ci_level * 100, "% CI):\n", sep = "")
  sel_df <- data.frame(
    item = prep$item_names,
    dif_loading = results[[1]],
    dif_intercept = results[[2]],
    selected = dif_flag
  )
  print(sel_df, row.names = FALSE)

  list(ldf = ldf, sel_df = sel_df)
}

# ---------------------------------------------------------------------------
# HELPER: extract growth parameters from Step 2
# ---------------------------------------------------------------------------
extract_growth_params <- function(fit2, sex_label) {
  params <- c("mu_slp", "mu_quad", "phi_int", "phi_slp", "eti_sd", "Omega[1,2]")
  draws <- fit2$draws(variables = params, format = "df")

  summ <- posterior::summarise_draws(
    draws,
    mean,
    sd,
    ~ quantile(.x, c(0.05, 0.25, 0.50, 0.75, 0.95)),
    posterior::default_convergence_measures()
  ) %>%
    mutate(sex = sex_label)

  summ
}

# ---------------------------------------------------------------------------
# HELPER: extract item parameter posteriors and compute implied params by age
# ---------------------------------------------------------------------------
extract_item_params <- function(
  fit2,
  prep,
  ldf_step2,
  sex_label,
  age_grid = seq(-4, 4, by = 0.5)
) {
  draws_base <- fit2$draws(variables = c("lp", "np"), format = "df")

  p <- prep$p
  items <- prep$item_names

  lp_mean <- colMeans(as.matrix(draws_base[, grep(
    "^lp\\[",
    names(draws_base)
  )]))
  np_mean <- colMeans(as.matrix(draws_base[, grep(
    "^np\\[",
    names(draws_base)
  )]))

  # Item-level DIF means, initialized to 0 (non-DIF items stay 0)
  ldiftv_mean <- rep(0, p)
  ndiftv_mean <- rep(0, p)

  # ldf_step2[, 1] == 1 marks which items have time-varying DIF in Step 2.
  # Stan's l_diftv[1..mtv] correspond in order to those flagged items.
  dif_item_idx <- which(ldf_step2[, 1] == 1)

  if (length(dif_item_idx) > 0) {
    draws_dif <- fit2$draws(variables = c("l_diftv", "n_diftv"), format = "df")
    ldiftv_cols <- grep("^l_diftv\\[", names(draws_dif))
    ndiftv_cols <- grep("^n_diftv\\[", names(draws_dif))

    if (length(ldiftv_cols) > 0) {
      ldiftv_mean[dif_item_idx] <- colMeans(
        as.matrix(draws_dif[, ldiftv_cols, drop = FALSE])
      )
    }
    if (length(ndiftv_cols) > 0) {
      ndiftv_mean[dif_item_idx] <- colMeans(
        as.matrix(draws_dif[, ndiftv_cols, drop = FALSE])
      )
    }
  }

  expand.grid(item_idx = seq_len(p), age_c = age_grid) %>%
    mutate(
      item = items[item_idx],
      reporter = ifelse(grepl("_p$", item), "Parent", "Youth"),
      base_item = sub("_(p|y)$", "", item),
      lam = lp_mean[item_idx] * exp(ldiftv_mean[item_idx] * age_c),
      nu = np_mean[item_idx] + ndiftv_mean[item_idx] * age_c,
      age = age_c + prep$age_mean,
      sex = sex_label
    )
}

# ---------------------------------------------------------------------------
# HELPER: extract factor scores (posterior mean of eta per person × wave)
# ---------------------------------------------------------------------------
extract_factor_scores <- function(fit2, prep, ldf_step2, sex_label) {
  # eta_tp = fac_gr[1,p] + fac_gr[2,p]*age_c + fac_gr[3,p]*age2_c + fac_eti[t,p]
  # We reconstruct from saved parameters
  draws <- fit2$draws(
    variables = c(
      "mu_slp",
      "mu_quad",
      "b_mu",
      "b_phi",
      "phi_int",
      "phi_slp",
      "L_Omega",
      "fac_dist",
      "fac_eti_raw",
      "eti_sd"
    ),
    format = "df"
  )

  mu_slp <- mean(draws$mu_slp)
  mu_quad <- mean(draws$mu_quad)
  eti_sd <- mean(draws$eti_sd)

  fac_dist_cols <- grep("^fac_dist\\[", names(draws))
  fac_eti_cols <- grep("^fac_eti_raw\\[", names(draws))

  fac_dist_mean <- colMeans(as.matrix(draws[, fac_dist_cols]))
  fac_eti_mean <- colMeans(as.matrix(draws[, fac_eti_cols]))

  ni <- prep$ni
  d <- prep$d

  # fac_dist is now 2 × ni (row-major in Stan → columns cycle fastest in R)
  fac_dist_mat <- matrix(fac_dist_mean, nrow = 2, ncol = ni)
  fac_eti_mat <- matrix(fac_eti_mean * eti_sd, nrow = d, ncol = ni)

  b_mu_cols <- grep("^b_mu\\[", names(draws))
  b_phi_cols <- grep("^b_phi\\[", names(draws))
  b_mu_mean <- colMeans(as.matrix(draws[, b_mu_cols, drop = FALSE]))
  b_phi_mean <- colMeans(as.matrix(draws[, b_phi_cols, drop = FALSE]))

  age_by_person <- prep$dat_long %>%
    group_by(person_idx) %>%
    summarise(mean_age_c = mean(age_c, na.rm = TRUE), .groups = "drop") %>%
    arrange(person_idx)
  xf_person <- scale(age_by_person$mean_age_c)[, 1]

  # L_Omega is now 2×2 Cholesky
  lomega_cols <- grep("^L_Omega\\[", names(draws))
  L_Omega_mean <- colMeans(as.matrix(draws[, lomega_cols]))
  L_Omega_mat <- matrix(L_Omega_mean, nrow = 2, ncol = 2)

  phi_int <- mean(draws$phi_int)
  phi_slp <- mean(draws$phi_slp)
  phi_eta <- c(phi_int, phi_slp)

  scores <- purrr::map_dfr(seq_len(ni), function(k) {
    mu_eta <- c(0, mu_slp)
    sd_eta <- phi_eta * exp(b_phi_mean * xf_person[k])
    fac_gr_k <- mu_eta +
      b_mu_mean * xf_person[k] +
      diag(sd_eta) %*% L_Omega_mat %*% fac_dist_mat[, k]

    purrr::map_dfr(seq_len(d), function(t) {
      age_obs <- prep$dat_long %>%
        filter(person_idx == k, time_idx == t) %>%
        slice(1)
      if (nrow(age_obs) == 0) {
        return(NULL)
      }

      age_c_val <- age_obs$age_c[1]
      age2_c_val <- age_obs$age2_c[1]
      eta_tp <- fac_gr_k[1] +
        fac_gr_k[2] * age_c_val +
        mu_quad * age2_c_val +
        fac_eti_mat[t, k]

      tibble(
        person_idx = k,
        id = prep$ids[k],
        time_idx = t,
        wave = wave_order[t],
        age = age_obs$age[1],
        eta = as.numeric(eta_tp)
      )
    })
  })

  scores %>% mutate(sex = sex_label)
}

# ---------------------------------------------------------------------------
# MAIN: sex-stratified two-step LMNLFA
# ---------------------------------------------------------------------------
all_dif_sel <- list()
all_growth <- list()
all_scores <- list()

for (sx in c("female", "male")) {
  cat("\n\n", strrep("=", 70), "\n")
  cat("  LMNLFA |", sx, "\n")
  cat(strrep("=", 70), "\n")

  parent_df <- if (sx == "female") female_parent else male_parent
  youth_df <- if (sx == "female") female_youth else male_youth

  # --- Data preparation ---------------------------------------------------
  prep <- build_lmnlfa_data(parent_df, youth_df, sx)

  # --- Step 1: All items have potential age DIF ---------------------------
  cat("\n--- Step 1: DIF screening (all items) ---\n")
  ldf_step1 <- matrix(c(rep(1L, prep$p), rep(0L, prep$p)), ncol = 2)
  stan_data_step1 <- make_stan_data(prep, ldf = ldf_step1, baseline_age = TRUE)

  fit1 <- fit_stan(
    stan_data_step1,
    stan_model,
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    label = paste("Step 1 |", sx)
  )

  # Convergence check
  diag1 <- fit1$diagnostic_summary(quiet = TRUE)
  cat("  Step 1 divergences:", sum(diag1$num_divergent), "\n")
  cat(
    "  Step 1 max Rhat:   ",
    round(max(fit1$summary()$rhat, na.rm = TRUE), 3),
    "\n"
  )

  # --- Diagnostics: growth variance parameters ----------------------------
  cat("\n--- Growth variance diagnostics [Step 1 |", sx, "] ---\n")
  diag_vars <- c("phi_int", "phi_slp", "eti_sd", "mu_slp", "mu_quad")
  diag_summ <- fit1$summary(
    variables = diag_vars,
    mean,
    sd,
    ~ quantile(.x, c(0.05, 0.5, 0.95)),
    posterior::default_convergence_measures()
  )
  print(diag_summ, digits = 3)

  # intercept-slope correlation element of the Cholesky factor
  lomega_summ <- fit1$summary(
    variables = "L_Omega",
    mean,
    sd,
    ~ quantile(.x, c(0.05, 0.5, 0.95))
  )
  cat("\nL_Omega (2x2 Cholesky; [2,1] is the correlation element):\n")
  print(lomega_summ, digits = 3)

  # How many times did phi_slp posterior mass sit below 0.05?
  phi_slp_draws <- fit1$draws("phi_slp", format = "matrix")
  cat(
    sprintf(
      "\n  phi_slp: mean=%.3f, P(phi_slp < 0.05)=%.2f%%\n",
      mean(phi_slp_draws),
      100 * mean(phi_slp_draws < 0.05)
    )
  )
  cat(
    sprintf(
      "  phi_int: mean=%.3f\n",
      mean(fit1$draws("phi_int", format = "matrix"))
    )
  )

  # Worst-Rhat parameters
  s1_summ <- fit1$summary()
  bad_rhat <- s1_summ[!is.na(s1_summ$rhat) & s1_summ$rhat > 1.05, ]
  if (nrow(bad_rhat) > 0) {
    cat("\n  Parameters with Rhat > 1.05:\n")
    print(
      bad_rhat[order(-bad_rhat$rhat), c("variable", "mean", "sd", "rhat")][
        seq_len(min(10, nrow(bad_rhat))),
      ],
      digits = 3
    )
  } else {
    cat("\n  All Rhat <= 1.05\n")
  }
  cat(strrep("-", 50), "\n")

  # --- DIF selection ------------------------------------------------------
  dif_result <- select_dif(fit1, prep, ci_level = 0.90)
  all_dif_sel[[sx]] <- dif_result$sel_df %>% mutate(sex = sx)

  # --- Step 2: Refit with selected DIF only ------------------------------
  cat("\n--- Step 2: Refit with DIF-selected items ---\n")
  ldf_step2 <- dif_result$ldf
  stan_data_step2 <- make_stan_data(prep, ldf = ldf_step2, baseline_age = TRUE)

  fit2 <- fit_stan(
    stan_data_step2,
    stan_model,
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    label = paste("Step 2 |", sx)
  )

  diag2 <- fit2$diagnostic_summary(quiet = TRUE)
  cat("  Step 2 divergences:", sum(diag2$num_divergent), "\n")
  cat(
    "  Step 2 max Rhat:   ",
    round(max(fit2$summary()$rhat, na.rm = TRUE), 3),
    "\n"
  )

  # --- Growth parameters --------------------------------------------------
  gp <- extract_growth_params(fit2, sx)
  cat("\nGrowth parameters [", sx, "]:\n")
  print(gp[, c("variable", "mean", "sd", "q5", "q95", "rhat")])
  all_growth[[sx]] <- gp

  # --- Factor scores -------------------------------------------------------
  cat("\nExtracting factor scores...\n")
  scores <- tryCatch(
    extract_factor_scores(fit2, prep, ldf_step2, sx),
    error = function(e) {
      message("Factor score extraction failed: ", e$message)
      NULL
    }
  )
  if (!is.null(scores)) {
    all_scores[[sx]] <- scores
  }

  # --- Plot: item parameters by age ----------------------------------------
  ip <- extract_item_params(fit2, prep, ldf_step2, sx)
  if (!is.null(ip) && nrow(ip) > 0) {
    p_ip <- ggplot(
      ip,
      aes(x = age, y = lam, colour = reporter, linetype = base_item)
    ) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c(Parent = "#2166ac", Youth = "#d73027")) +
      labs(
        title = paste0("Item loadings by age — ", sx),
        subtitle = "Posterior mean loading (λ) as a function of age",
        x = "Age (years)",
        y = "Loading (λ)",
        colour = "Reporter",
        linetype = "Item"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "bottom",
        plot.subtitle = element_text(colour = "grey45")
      )

    ggsave(
      file.path(out_dir, paste0("item_loadings_by_age_", sx, ".png")),
      p_ip,
      width = 8,
      height = 5,
      dpi = 180
    )

    p_nu <- ggplot(
      ip,
      aes(x = age, y = nu, colour = reporter, linetype = base_item)
    ) +
      geom_line(linewidth = 1) +
      scale_colour_manual(values = c(Parent = "#2166ac", Youth = "#d73027")) +
      labs(
        title = paste0("Item intercepts by age — ", sx),
        subtitle = "Posterior mean intercept (ν) as a function of age",
        x = "Age (years)",
        y = "Intercept (ν)",
        colour = "Reporter",
        linetype = "Item"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        legend.position = "bottom",
        plot.subtitle = element_text(colour = "grey45")
      )

    ggsave(
      file.path(out_dir, paste0("item_intercepts_by_age_", sx, ".png")),
      p_nu,
      width = 8,
      height = 5,
      dpi = 180
    )
  }

  # --- Plot: mean growth trajectory ----------------------------------------
  # Implied mean trajectory: eta = 0 + mu_slp * age_c + mu_quad * age2_c
  gp_means <- fit2$draws(variables = c("mu_slp", "mu_quad"), format = "df")
  age_grid <- seq(
    min(prep$dat_long$age, na.rm = TRUE),
    max(prep$dat_long$age, na.rm = TRUE),
    length.out = 100
  )
  age_c_grid <- age_grid - prep$age_mean
  age2_c_grid <- age_c_grid^2

  # Draw-level trajectories for uncertainty ribbon
  traj_draws <- purrr::map_dfr(seq_len(min(200, nrow(gp_means))), function(i) {
    tibble(
      age = age_grid,
      eta_mean = gp_means$mu_slp[i] *
        age_c_grid +
        gp_means$mu_quad[i] * age2_c_grid,
      draw = i
    )
  })

  traj_summ <- traj_draws %>%
    group_by(age) %>%
    summarise(
      eta_med = median(eta_mean),
      eta_lo = quantile(eta_mean, 0.05),
      eta_hi = quantile(eta_mean, 0.95),
      .groups = "drop"
    )

  p_traj <- ggplot(traj_summ, aes(x = age)) +
    geom_ribbon(
      aes(ymin = eta_lo, ymax = eta_hi),
      fill = "#4393c3",
      alpha = 0.25
    ) +
    geom_line(aes(y = eta_med), colour = "#2166ac", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    labs(
      title = paste0("Mean puberty growth trajectory — ", sx),
      subtitle = "Posterior median + 90% CI; intercept set to mean at mean age",
      x = "Age (years)",
      y = "Latent puberty (η)"
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.subtitle = element_text(colour = "grey45"))

  ggsave(
    file.path(out_dir, paste0("growth_trajectory_", sx, ".png")),
    p_traj,
    width = 7,
    height = 5,
    dpi = 180
  )
  cat("Plots written for", sx, "\n")

  # Save fit objects for diagnostics (optional, large files)
  # fit1$save_object(file.path(out_dir, paste0("fit1_", sx, ".rds")))
  # fit2$save_object(file.path(out_dir, paste0("fit2_", sx, ".rds")))
}

# ---------------------------------------------------------------------------
# OUTPUT: CSVs
# ---------------------------------------------------------------------------
dif_table <- bind_rows(all_dif_sel)
write.csv(dif_table, file.path(out_dir, "dif_selection.csv"), row.names = FALSE)
cat("\nDIF selection table written.\n")

growth_table <- bind_rows(all_growth)
write.csv(
  growth_table,
  file.path(out_dir, "growth_params.csv"),
  row.names = FALSE
)
cat("Growth parameter table written.\n")

if (length(all_scores) > 0) {
  scores_table <- bind_rows(all_scores)
  write.csv(
    scores_table,
    file.path(out_dir, "factor_scores.csv"),
    row.names = FALSE
  )
  cat("Factor scores written (", nrow(scores_table), "rows).\n")
}

cat("\nAll LMNLFA outputs written to:", out_dir, "\n")
