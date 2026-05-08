## 02_psychometrics.R
## Psychometric evaluation of PDS items across waves, reporters, and sex.
## Answers two questions before downstream modelling:
##   (1) Which items are informative? (GRM discrimination + item info)
##   (2) Do we need both parent and youth items, or are they redundant?
##       (cross-reporter bifactor CFA)
## Also runs longitudinal measurement invariance and invariance across
## age tertiles, race/ethnicity, and BMI tertiles. ## to do: add site somehow??
## Requires outputs from 00_data_foundation.R.

# export DATA_DIR="/u/project/silvers/data/ABCD/ABCD-release-6.0/cfm/physical-health/puberty"
# export OUT_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"
# Rscript 02_psychometrics.R

pacman::p_load(
  dplyr,
  tidyr,
  tibble,
  ggplot2,
  psych,
  mirt,
  lavaan,
  semTools,
  install = TRUE
)

set.seed(90025)

# safe print: ensures tibble dispatch so n = Inf works
print_all <- function(x) print(as_tibble(x), n = Inf)

# safe fit-measures row: avoids t() matrix issues with vctrs bind_rows
fit_to_row <- function(fit, fi_names) {
  if (is.null(fit)) {
    return(NULL)
  }
  fi <- tryCatch(lavaan::fitMeasures(fit, fi_names), error = function(e) NULL)
  if (is.null(fi) || length(fi) == 0) {
    return(NULL)
  }
  do.call(data.frame, c(as.list(fi), list(stringsAsFactors = FALSE)))
}

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
pub_root <- data_dir

out_base <- Sys.getenv("OUT_DIR")
if (!nzchar(out_base)) {
  out_base <- file.path(
    root_path,
    "projects/abcd-projs",
    "dissertation/study1/outputs"
  )
}
out_dir <- file.path(out_base, "psychometrics")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# LOAD
# ---------------------------------------------------------------------------
female_parent <- read.csv(file.path(pub_root, "female_parent_long.csv"))
female_youth <- read.csv(file.path(pub_root, "female_youth_long.csv"))
male_parent <- read.csv(file.path(pub_root, "male_parent_long.csv"))
male_youth <- read.csv(file.path(pub_root, "male_youth_long.csv"))

wave_order <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")

# Invariance sections and bifactor use the shared 4-item ordinal set only.
ordinal_items <- c("peta", "petb", "petc", "petd")
female_items <- c("peta", "petb", "petc", "petd", "fpete")
female_binary_items <- "fpete"
male_items <- c("peta", "petb", "petc", "petd", "mpete")
male_binary_items <- character(0)

for (nm in c("female_parent", "female_youth", "male_parent", "male_youth")) {
  d <- get(nm)
  d$wave <- factor(d$wave, levels = wave_order)
  assign(nm, d)
}

# ---------------------------------------------------------------------------
# SECTION 1: PER-WAVE PSYCHOMETRICS
# For each sex × reporter × wave:
#   polychoric correlations, alpha, omega, EFA (1-factor), CFA (WLSMV), GRM
# ---------------------------------------------------------------------------

datasets <- list(
  female_parent = list(
    df = female_parent,
    sex = "female",
    reporter = "parent",
    items = female_items,
    binary = female_binary_items
  ),
  female_youth = list(
    df = female_youth,
    sex = "female",
    reporter = "youth",
    items = female_items,
    binary = female_binary_items
  ),
  male_parent = list(
    df = male_parent,
    sex = "male",
    reporter = "parent",
    items = male_items,
    binary = male_binary_items
  ),
  male_youth = list(
    df = male_youth,
    sex = "male",
    reporter = "youth",
    items = male_items,
    binary = male_binary_items
  )
)

wave_results <- list()
disc_rows <- list() # GRM discrimination params, collected for summary table
iif_plots <- list() # item information function plots

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  df <- ds$df
  sex_lab <- ds$sex
  rep_lab <- ds$reporter
  items_bin <- ds$binary
  items_ord <- setdiff(ds$items, items_bin)

  for (wv in wave_order) {
    key <- paste(ds_name, wv, sep = "_")

    sub_raw <- df %>%
      filter(wave == wv) %>%
      select(any_of(ds$items)) %>%
      mutate(across(everything(), as.integer))

    items_ord_avail <- intersect(items_ord, names(sub_raw))
    items_bin_avail <- intersect(items_bin, names(sub_raw))
    items_avail <- c(items_ord_avail, items_bin_avail)

    sub <- sub_raw %>%
      filter(if_all(all_of(items_ord_avail), ~ !is.na(.) & . %in% 1:4))
    if (length(items_bin_avail) > 0) {
      sub <- sub %>%
        filter(if_all(all_of(items_bin_avail), ~ !is.na(.) & . %in% 1:2))
    }
    sub <- select(sub, all_of(items_avail))

    if (nrow(sub) < 100) {
      message("Skipping ", key, ": n = ", nrow(sub))
      next
    }

    # Drop zero-variance items (e.g., fpete before foundation recode)
    zero_var <- names(sub)[sapply(sub, function(x) {
      length(unique(na.omit(x))) <= 1
    })]
    if (length(zero_var) > 0) {
      message(
        "Wave ",
        wv,
        " [",
        key,
        "]: removing zero-variance items: ",
        paste(zero_var, collapse = ", ")
      )
      sub <- sub[, setdiff(names(sub), zero_var), drop = FALSE]
      items_avail <- setdiff(items_avail, zero_var)
      items_ord_avail <- setdiff(items_ord_avail, zero_var)
      items_bin_avail <- setdiff(items_bin_avail, zero_var)
      if (length(items_avail) == 0 || nrow(sub) < 100) {
        message(
          "Skipping ",
          key,
          ": no valid items after zero-variance removal"
        )
        next
      }
    }

    result <- list(n = nrow(sub))

    # --- polychoric + reliability ---
    pc <- tryCatch(
      psych::polychoric(sub, na.rm = TRUE, correct = 0),
      error = function(e) NULL
    )
    if (!is.null(pc)) {
      result$poly_rho <- pc$rho
      result$alpha <- tryCatch(
        psych::alpha(sub, na.rm = TRUE, check.keys = TRUE)$total,
        error = function(e) NULL
      )
      # omega_h requires >= 2 factors; skip for per-wave (used in bifactor section)
    }

    # --- EFA 1-factor (polychoric) ---
    result$efa <- tryCatch(
      psych::fa(
        sub,
        nfactors = 1,
        fm = "minres",
        rotate = "none",
        cor = "poly"
      ),
      error = function(e) NULL
    )

    # --- CFA 1-factor WLSMV ---
    # WLSMV treats all ordered vars correctly regardless of number of categories.
    cfa_syntax <- paste("pub =~", paste(items_avail, collapse = " + "))
    result$cfa <- tryCatch(
      lavaan::cfa(
        cfa_syntax,
        data = sub,
        ordered = items_avail,
        estimator = "WLSMV",
        missing = "pairwise"
      ),
      error = function(e) NULL
    )

    # --- GRM (ordinal items only; fpete excluded) ---
    grm_fit <- tryCatch(
      {
        mat <- as.matrix(sub[items_ord_avail]) - 1L # mirt expects 0-based
        mirt::mirt(mat, 1, itemtype = "graded", method = "EM", verbose = FALSE)
      },
      error = function(e) {
        message("GRM failed [", key, "]: ", e$message)
        NULL
      }
    )

    if (!is.null(grm_fit)) {
      result$grm <- grm_fit
      coefs <- tryCatch(
        mirt::coef(grm_fit, IRTpars = TRUE, simplify = TRUE)$items,
        error = function(e) NULL
      )
      if (!is.null(coefs)) {
        result$grm_coefs <- coefs
        coef_df <- as.data.frame(coefs)
        # mirt names discrimination 'a' with IRTpars=TRUE (1D); fallback to 'a1'
        a_col <- intersect(c("a", "a1"), names(coef_df))[1]
        if (is.na(a_col)) {
          a_col <- names(coef_df)[1]
        }
        for (col in c("b1", "b2", "b3")) {
          if (!col %in% names(coef_df)) coef_df[[col]] <- NA_real_
        }
        disc_rows[[key]] <- data.frame(
          item = rownames(coef_df),
          a = coef_df[[a_col]],
          b1 = coef_df[["b1"]],
          b2 = coef_df[["b2"]],
          b3 = coef_df[["b3"]],
          wave = wv,
          sex = sex_lab,
          reporter = rep_lab,
          dataset = ds_name,
          stringsAsFactors = FALSE
        )
      }

      # item information function plot
      iif_plots[[key]] <- tryCatch(
        {
          p <- plot(
            grm_fit,
            type = "infotrace",
            main = paste("IIF:", ds_name, wv)
          )
          p
        },
        error = function(e) NULL
      )
    }

    wave_results[[key]] <- result
  }
}

# --- GRM discrimination summary table ---
disc_table <- bind_rows(disc_rows)
cat("\n=== GRM item discrimination (a) by wave × reporter × sex ===\n")
print_all(disc_table)
write.csv(
  disc_table,
  file.path(out_dir, "grm_discrimination.csv"),
  row.names = FALSE
)

# flag items with consistently low discrimination (a < 0.7 at majority of waves)
if (nrow(disc_table) > 0) {
  low_disc <- disc_table %>%
    group_by(item, sex, reporter) %>%
    summarise(
      mean_a = mean(a, na.rm = TRUE),
      pct_low = mean(a < 0.7, na.rm = TRUE),
      n_waves = n(),
      .groups = "drop"
    ) %>%
    arrange(sex, reporter, item)
} else {
  warning(
    "No GRM results — mirt may not be installed correctly. Skipping discrimination summary."
  )
  low_disc <- data.frame()
}

cat(
  "\n=== Item discrimination summary (flag if a < 0.7 at majority of waves) ===\n"
)
print_all(low_disc)
write.csv(
  low_disc,
  file.path(out_dir, "item_discrimination_summary.csv"),
  row.names = FALSE
)

# --- CFA fit summary table ---
fit_indices <- c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")
cfa_fits <- lapply(names(wave_results), function(key) {
  fit <- wave_results[[key]]$cfa
  if (is.null(fit)) {
    return(NULL)
  }
  fi <- tryCatch(lavaan::fitMeasures(fit, fit_indices), error = function(e) {
    NULL
  })
  if (is.null(fi) || length(fi) < length(fit_indices)) {
    return(NULL)
  }
  parts <- strsplit(key, "_")[[1]]
  tibble::tibble(
    dataset = paste(parts[1], parts[2], sep = "_"),
    wave = tail(parts, 1),
    cfi.scaled = unname(fi["cfi.scaled"]),
    tli.scaled = unname(fi["tli.scaled"]),
    rmsea.scaled = unname(fi["rmsea.scaled"]),
    srmr = unname(fi["srmr"])
  )
}) %>%
  dplyr::bind_rows()

cat("\n=== CFA fit indices per wave × reporter × sex ===\n")
print_all(cfa_fits)
write.csv(
  cfa_fits,
  file.path(out_dir, "cfa_fit_per_wave.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# SECTION 1b: POOLED EFA AND RELIABILITY
# Pooled across all waves, per sex × reporter.
# (1) Cronbach's alpha + omega (h and total) with ECV
# (2) Parallel analysis (fa.parallel on polychoric R) to determine factor count
# (3) EFA 1- and 2-factor solutions (oblimin rotation for 2+)
# Outputs: reliability_summary.csv, efa_loadings.csv, parallel_analysis_plot.png
# ---------------------------------------------------------------------------

efa_results <- list()

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  df <- ds$df
  items_bin <- ds$binary
  items_ord <- setdiff(ds$items, items_bin)

  # Pool all waves; apply same validity filters as per-wave loop
  sub_raw <- df %>%
    select(any_of(ds$items)) %>%
    mutate(across(everything(), as.integer))

  sub <- sub_raw %>%
    filter(if_all(all_of(items_ord), ~ !is.na(.) & . %in% 1:4))
  if (length(items_bin) > 0) {
    sub <- sub %>% filter(if_all(all_of(items_bin), ~ !is.na(.) & . %in% 1:2))
  }
  sub <- select(sub, all_of(ds$items))

  cat("\n\n=== EFA / Reliability:", ds_name, "| n =", nrow(sub), "===\n")

  if (nrow(sub) < 200) {
    message("Insufficient data for ", ds_name, " — skipping.")
    next
  }

  # Polychoric correlation matrix (basis for all analyses below)
  pc <- tryCatch(
    psych::polychoric(sub, na.rm = TRUE, correct = 0),
    error = function(e) {
      message("polychoric failed: ", e$message)
      NULL
    }
  )
  if (is.null(pc)) {
    next
  }

  # --- Reliability ---
  alpha_res <- tryCatch(
    psych::alpha(sub, na.rm = TRUE, check.keys = TRUE),
    error = function(e) NULL
  )
  # omega requires >= 2 factors; nfactors=2 is standard for hierarchical omega
  omega_res <- tryCatch(
    psych::omega(pc$rho, nfactors = 2, fm = "minres", plot = FALSE, sl = TRUE),
    error = function(e) {
      message("omega failed: ", e$message)
      NULL
    }
  )

  cat(
    "Alpha  :",
    if (!is.null(alpha_res)) round(alpha_res$total$raw_alpha, 3) else NA,
    "\n"
  )
  if (!is.null(omega_res)) {
    cat("Omega_h:", round(omega_res$omega_h, 3), "  (general factor only)\n")
    cat("Omega_t:", round(omega_res$omega.tot, 3), "  (total reliability)\n")
    cat(
      "ECV    :",
      round(omega_res$ECV, 3),
      "  (explained common variance by g)\n"
    )
  }

  # --- Parallel analysis ---
  cat("\n--- Parallel analysis ---\n")
  pa <- tryCatch(
    psych::fa.parallel(
      pc$rho,
      n.obs = nrow(sub),
      fm = "minres",
      fa = "fa",
      plot = FALSE
    ),
    error = function(e) {
      message("fa.parallel failed: ", e$message)
      NULL
    }
  )
  if (!is.null(pa)) {
    cat("Suggested factors:", pa$nfact, "\n")
    cat("Actual eigenvalues    :", round(pa$fa.values, 3), "\n")
    cat(
      "Simulated eigenvalues :",
      round(pa$fa.sim[seq_along(pa$fa.values)], 3),
      "\n"
    )
  }

  # --- EFA 1-factor ---
  efa1 <- tryCatch(
    psych::fa(sub, nfactors = 1, fm = "minres", rotate = "none", cor = "poly"),
    error = function(e) NULL
  )
  # --- EFA 2-factor (oblimin allows factors to correlate) ---
  efa2 <- tryCatch(
    psych::fa(
      sub,
      nfactors = 2,
      fm = "minres",
      rotate = "oblimin",
      cor = "poly"
    ),
    error = function(e) NULL
  )

  cat("\n--- 1-factor EFA loadings ---\n")
  if (!is.null(efa1)) {
    print(efa1$loadings, cutoff = 0.1)
    cat(
      "RMSEA:",
      round(efa1$RMSEA[1], 3),
      "| TLI:",
      round(efa1$TLI, 3),
      "| RMSR:",
      round(efa1$rms, 3),
      "\n"
    )
  }
  cat("\n--- 2-factor EFA loadings (oblimin) ---\n")
  if (!is.null(efa2)) {
    print(efa2$loadings, cutoff = 0.1)
    cat(
      "RMSEA:",
      round(efa2$RMSEA[1], 3),
      "| TLI:",
      round(efa2$TLI, 3),
      "| RMSR:",
      round(efa2$rms, 3),
      "\n"
    )
    if (!is.null(efa2$Phi)) {
      cat("Factor correlation (phi):", round(efa2$Phi[1, 2], 3), "\n")
    }
  }

  efa_results[[ds_name]] <- list(
    n = nrow(sub),
    alpha = alpha_res,
    omega = omega_res,
    parallel = pa,
    efa1 = efa1,
    efa2 = efa2
  )
}

# --- Reliability summary table ---
rel_rows <- lapply(names(efa_results), function(nm) {
  r <- efa_results[[nm]]
  data.frame(
    dataset = nm,
    n = r$n,
    alpha = if (!is.null(r$alpha)) round(r$alpha$total$raw_alpha, 3) else NA,
    omega_h = if (!is.null(r$omega)) round(r$omega$omega_h, 3) else NA,
    omega_t = if (!is.null(r$omega)) round(r$omega$omega.tot, 3) else NA,
    ecv = if (!is.null(r$omega)) round(r$omega$ECV, 3) else NA,
    n_factors_parallel = if (!is.null(r$parallel)) r$parallel$nfact else NA,
    stringsAsFactors = FALSE
  )
})
rel_table <- bind_rows(rel_rows)
cat("\n=== Reliability summary ===\n")
print(rel_table)
write.csv(
  rel_table,
  file.path(out_dir, "reliability_summary.csv"),
  row.names = FALSE
)

# --- EFA loadings summary table ---
make_loadings_df <- function(efa_obj, n_factors, dataset) {
  if (is.null(efa_obj)) {
    return(NULL)
  }
  L <- as.data.frame(unclass(efa_obj$loadings))
  L$item <- rownames(L)
  L$dataset <- dataset
  L$n_factors <- n_factors
  L$rmsea <- round(efa_obj$RMSEA[1], 3)
  L$tli <- round(efa_obj$TLI, 3)
  L$rmsr <- round(efa_obj$rms, 3)
  L
}

loadings_table <- lapply(names(efa_results), function(nm) {
  r <- efa_results[[nm]]
  bind_rows(
    make_loadings_df(r$efa1, 1, nm),
    make_loadings_df(r$efa2, 2, nm)
  )
}) %>%
  bind_rows()

write.csv(
  loadings_table,
  file.path(out_dir, "efa_loadings.csv"),
  row.names = FALSE
)
cat("\nEFA loadings written to:", file.path(out_dir, "efa_loadings.csv"), "\n")

# --- Parallel analysis scree plot ---
pa_plot_df <- lapply(names(efa_results), function(nm) {
  pa <- efa_results[[nm]]$parallel
  if (is.null(pa)) {
    return(NULL)
  }
  n_ev <- length(pa$fa.values)
  data.frame(
    dataset = nm,
    factor = seq_len(n_ev),
    actual = pa$fa.values,
    simulated = pa$fa.sim[seq_len(n_ev)],
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

if (nrow(pa_plot_df) > 0) {
  pa_long <- pa_plot_df %>%
    tidyr::pivot_longer(
      c(actual, simulated),
      names_to = "type",
      values_to = "eigenvalue"
    ) %>%
    mutate(
      dataset_label = dplyr::recode(
        dataset,
        female_parent = "Female parent",
        female_youth = "Female youth",
        male_parent = "Male parent",
        male_youth = "Male youth"
      )
    )

  p_pa <- ggplot(
    pa_long,
    aes(x = factor(factor), y = eigenvalue, colour = type, group = type)
  ) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    scale_colour_manual(
      values = c(actual = "#2166ac", simulated = "#d73027"),
      labels = c(actual = "Actual", simulated = "Simulated (parallel)")
    ) +
    facet_wrap(~dataset_label, nrow = 1) +
    labs(
      title = "Parallel analysis — factor retention",
      subtitle = "Retain factors where actual eigenvalue > simulated",
      x = "Factor",
      y = "Eigenvalue",
      colour = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(colour = "grey45", size = 11)
    )

  ggsave(
    file.path(out_dir, "parallel_analysis_plot.png"),
    p_pa,
    width = 12,
    height = 4.5,
    dpi = 180
  )
  cat("\nParallel analysis plot written.\n")
}

# ---------------------------------------------------------------------------
# SECTION 2: CROSS-REPORTER BIFACTOR CFA
# One randomly selected wave per participant (preserves independence).
# Per sex. Tests whether a general puberty factor + reporter method factors
# fits better than two separate reporter factors. ## to do: also select one time point randomly for each participant for the above analyses
# ---------------------------------------------------------------------------

build_cr_data <- function(parent_df, youth_df, include_pete = FALSE) {
  base_items <- c("peta", "petb", "petc", "petd")
  # fpete (female) and mpete (male) are sex-specific; any_of() picks whichever exists
  all_items <- if (include_pete) c(base_items, "fpete", "mpete") else base_items

  parent_sel <- parent_df %>%
    select(id, wave, any_of(all_items)) %>%
    rename_with(~ paste0(., "_p"), any_of(all_items))

  youth_sel <- youth_df %>%
    select(id, wave, any_of(all_items)) %>%
    rename_with(~ paste0(., "_y"), any_of(all_items))

  inner_join(parent_sel, youth_sel, by = c("id", "wave")) %>%
    mutate(across(-c(id, wave), as.integer)) %>%
    filter(if_all(-c(id, wave), ~ !is.na(.) & . %in% 1:4)) %>%
    # one randomly selected wave per participant → independent observations
    group_by(id) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(-id, -wave)
}

run_bifactor_section <- function(
  parent_df,
  youth_df,
  sex_label,
  include_pete = FALSE
) {
  dat <- build_cr_data(parent_df, youth_df, include_pete = include_pete)
  cat(
    "\n\n=== Cross-reporter bifactor | sex:",
    sex_label,
    "| n =",
    nrow(dat),
    "===\n"
  )

  if (nrow(dat) < 200) {
    message("Insufficient data for bifactor analysis (", sex_label, ")")
    return(NULL)
  }

  pc <- tryCatch(
    psych::polychoric(dat, na.rm = TRUE, correct = 0)$rho,
    error = function(e) {
      message("polychoric failed")
      NULL
    }
  )
  if (is.null(pc)) {
    return(NULL)
  }

  # Schmid-Leiman / omega to get omega_h
  sl <- tryCatch(
    psych::omega(
      pc,
      nfactors = 2,
      fm = "minres",
      plot = FALSE,
      sl = TRUE,
      key = NULL
    ),
    error = function(e) {
      message("omega/SL failed: ", e$message)
      NULL
    }
  )
  if (!is.null(sl)) {
    cat("Omega_h (general factor):", round(sl$omega_h, 3), "\n")
    cat("Omega_t (total)         :", round(sl$omega.tot, 3), "\n")
    cat("ECV (explained common variance by general):", round(sl$ECV, 3), "\n")
  }

  # Derive item lists from the joined data (automatically includes pete if present)
  cr_items <- names(dat)
  cr_items_p <- grep("_p$", cr_items, value = TRUE)
  cr_items_y <- grep("_y$", cr_items, value = TRUE)

  fi_names <- c(
    "cfi.scaled",
    "tli.scaled",
    "rmsea.scaled",
    "srmr",
    "chisq.scaled",
    "df.scaled"
  )

  safe_cfa <- function(mod, label, theta = FALSE) {
    tryCatch(
      withCallingHandlers(
        lavaan::cfa(
          mod,
          data = dat,
          ordered = cr_items,
          estimator = "WLSMV",
          parameterization = if (theta) "theta" else "delta",
          missing = "pairwise"
        ),
        warning = function(w) {
          message(label, " warning [", sex_label, "]: ", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        message(label, " failed [", sex_label, "]: ", e$message)
        NULL
      }
    )
  }

  # Model 1: single puberty factor (no method effects).
  # Baseline — reporters interchangeable; no systematic reporter variance.
  mod_m1 <- paste("puberty =~", paste(cr_items, collapse = " + "))
  fit_m1 <- safe_cfa(mod_m1, "M1 (single factor)")

  # Model 2: correlated reporter factors, no general trait.
  # reporter_p and reporter_y as separate latent variables, freely correlated.
  # Quantifies latent cross-reporter agreement without imposing a shared construct.
  mod_m2 <- paste0(
    "reporter_p =~ ",
    paste(cr_items_p, collapse = " + "),
    "\n",
    "reporter_y =~ ",
    paste(cr_items_y, collapse = " + ")
  )
  fit_m2 <- safe_cfa(mod_m2, "M2 (correlated reporters)")

  # Model 3: single puberty factor + within-reporter correlated uniqueness.
  # Residuals of same-reporter items are allowed to covary (reporter-specific
  # variance beyond the general factor). CU is within-reporter only — this is
  # identifiable and does not conflict with the general factor.
  cu_within <- c(
    apply(combn(cr_items_p, 2), 2, function(x) paste0(x[1], " ~~ ", x[2])),
    apply(combn(cr_items_y, 2), 2, function(x) paste0(x[1], " ~~ ", x[2]))
  )
  mod_m3 <- paste0(
    "puberty =~ ",
    paste(cr_items, collapse = " + "),
    "\n",
    paste(cu_within, collapse = "\n")
  )
  # theta parameterization frees residual variances — more stable for complex
  # models with many CU parameters (avoids near-singular information matrices)
  fit_m3 <- safe_cfa(
    mod_m3,
    "M3 (single factor + within-reporter CU)",
    theta = TRUE
  )
  fit_m3_diag <- fit_m3 # keep raw fit for diagnostics before nulling
  if (!is.null(fit_m3) && !isTRUE(lavaan::lavInspect(fit_m3, "converged"))) {
    message("M3 did not converge [", sex_label, "] — excluded from fit table")
    fit_m3 <- NULL
  }

  # Model 4: single puberty factor + cross-reporter CU (matched items).
  # The residuals of the *same item* rated by parent and youth are allowed to
  # correlate — item-level convergence beyond the shared trait (CU-MTMM).
  # Match by stripping suffix so peta_p pairs with peta_y regardless of order.
  base_p <- sub("_p$", "", cr_items_p)
  base_y <- sub("_y$", "", cr_items_y)
  matched <- intersect(base_p, base_y)
  cu_cross <- paste0(matched, "_p ~~ ", matched, "_y")
  cat(
    "\nM4 cross-reporter CU pairs [",
    sex_label,
    "]:\n",
    paste(cu_cross, collapse = "\n"),
    "\n"
  )
  mod_m4 <- paste0(
    "puberty =~ ",
    paste(cr_items, collapse = " + "),
    "\n",
    paste(cu_cross, collapse = "\n")
  )
  fit_m4 <- safe_cfa(mod_m4, "M4 (single factor + cross-reporter CU)")
  if (!is.null(fit_m4) && !isTRUE(lavaan::lavInspect(fit_m4, "converged"))) {
    message("M4 did not converge [", sex_label, "] — excluded")
    fit_m4 <- NULL
  }

  # Always show all 4 models; failed models appear as NA rows
  na_row <- as.data.frame(
    setNames(as.list(rep(NA_real_, length(fi_names))), fi_names)
  )
  fit_table <- lapply(
    list(
      m1_single_factor = fit_m1,
      m2_corr_reporters = fit_m2,
      m3_cu_within = fit_m3,
      m4_cu_cross = fit_m4
    ),
    function(fit) if (!is.null(fit)) fit_to_row(fit, fi_names) else na_row
  ) %>%
    bind_rows(.id = "model") %>%
    mutate(sex = sex_label)

  cat("\n--- Model fit comparison ---\n")
  print(fit_table)

  # --- Diagnostics: MIs + Heywood check for all models ---------------------
  diagnose_fit <- function(fit, label) {
    cat("\n>>> Diagnostics:", label, "[", sex_label, "] <<<\n")

    if (is.null(fit)) {
      cat("  Model did not converge or failed — no fit object available.\n")
      return(invisible(NULL))
    }

    converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
    cat("  Converged:", converged, "\n")

    # Heywood check: residual variances < 0 or standardised loading > 1
    pe <- tryCatch(
      lavaan::parameterEstimates(fit, standardized = TRUE),
      error = function(e) NULL
    )
    if (!is.null(pe)) {
      resid_neg <- pe[
        pe$op == "~~" & pe$lhs == pe$rhs & pe$est < 0,
        c("lhs", "est", "std.all")
      ]
      load_high <- pe[
        pe$op == "=~" & !is.na(pe$std.all) & abs(pe$std.all) > 1,
        c("lhs", "rhs", "est", "std.all")
      ]
      if (nrow(resid_neg) > 0) {
        cat("  Heywood: negative residual variances:\n")
        print(resid_neg, row.names = FALSE)
      }
      if (nrow(load_high) > 0) {
        cat("  Heywood: |standardized loading| > 1:\n")
        print(load_high, row.names = FALSE)
      }
      if (nrow(resid_neg) == 0 && nrow(load_high) == 0) {
        cat("  No Heywood cases detected.\n")
      }
    }

    # Modification indices (top 15 by MI value, all ops)
    mi <- tryCatch(
      lavaan::modindices(
        fit,
        sort. = TRUE,
        maximum.number = 15,
        op = c("=~", "~~")
      ),
      error = function(e) {
        cat("  modindices() failed:", e$message, "\n")
        NULL
      }
    )
    if (!is.null(mi) && nrow(mi) > 0) {
      cat("  Top modification indices:\n")
      print(
        mi[, c("lhs", "op", "rhs", "mi", "epc", "sepc.all")],
        row.names = FALSE,
        digits = 3
      )
    }
  }

  diagnose_fit(fit_m1, "M1 (single factor)")
  diagnose_fit(fit_m2, "M2 (correlated reporters)")
  diagnose_fit(fit_m3_diag, "M3 (within-reporter CU)") # use raw fit before NULL
  diagnose_fit(fit_m4, "M4 (cross-reporter CU)")

  # --- M2 characterization -------------------------------------------------
  # Latent reporter correlation: r ≈ 1.0 = interchangeable; r < .90 = distinct
  lat_cor <- if (!is.null(fit_m2)) {
    tryCatch(lavaan::lavInspect(fit_m2, "cor.lv"), error = function(e) NULL)
  } else {
    NULL
  }

  # Latent correlation with SE and 95% CI from parameterEstimates
  lat_cor_row <- NULL
  if (!is.null(fit_m2)) {
    pe <- tryCatch(
      lavaan::parameterEstimates(fit_m2, ci = TRUE, standardized = TRUE),
      error = function(e) NULL
    )
    if (!is.null(pe)) {
      lat_cor_row <- pe[
        pe$op == "~~" & pe$lhs == "reporter_p" & pe$rhs == "reporter_y",
      ]
    }
  }

  if (!is.null(lat_cor)) {
    cat("\nM2 latent reporter correlation [", sex_label, "]:\n")
    print(round(lat_cor, 3))
  }
  if (!is.null(lat_cor_row) && nrow(lat_cor_row) > 0) {
    cat(sprintf(
      "  r = %.3f, SE = %.3f, 95%% CI [%.3f, %.3f]\n",
      lat_cor_row$std.all,
      lat_cor_row$se,
      lat_cor_row$ci.lower,
      lat_cor_row$ci.upper
    ))
  }

  # Standardized loadings for both reporter factors
  m2_loadings <- NULL
  if (!is.null(fit_m2)) {
    pe <- tryCatch(
      lavaan::parameterEstimates(fit_m2, ci = TRUE, standardized = TRUE),
      error = function(e) NULL
    )
    if (!is.null(pe)) {
      m2_loadings <- pe[
        pe$op == "=~",
        c("lhs", "rhs", "std.all", "se", "ci.lower", "ci.upper", "pvalue")
      ] %>%
        dplyr::rename(factor = lhs, item = rhs, std_loading = std.all) %>%
        dplyr::mutate(sex = sex_label)

      cat("\nM2 standardized loadings [", sex_label, "]:\n")
      print(
        m2_loadings[, c(
          "factor",
          "item",
          "std_loading",
          "se",
          "ci.lower",
          "ci.upper"
        )],
        digits = 3,
        row.names = FALSE
      )
    }
  }

  # --- M2 refinement -------------------------------------------------------
  # Add residual correlations with MI > 10 to improve M2 fit.
  # Only ~~ paths between observed items (no cross-loadings, no factor terms).
  fit_m2_refined <- NULL
  mod_m2_refined <- NULL
  m2_ref_loadings <- NULL
  m2_ref_fit_row <- NULL

  if (!is.null(fit_m2)) {
    mi_m2 <- tryCatch(
      lavaan::modindices(fit_m2, sort. = TRUE, op = "~~"),
      error = function(e) NULL
    )
    if (!is.null(mi_m2)) {
      lv_names <- c("reporter_p", "reporter_y")
      mi_sig <- mi_m2[
        mi_m2$mi >= 10 &
          mi_m2$lhs != mi_m2$rhs &
          !(mi_m2$lhs %in% lv_names) &
          !(mi_m2$rhs %in% lv_names),
      ]

      if (nrow(mi_sig) > 0) {
        cat(
          "\nM2 refinement: adding",
          nrow(mi_sig),
          "residual correlations [",
          sex_label,
          "]:\n"
        )
        print(
          mi_sig[, c("lhs", "op", "rhs", "mi", "epc", "sepc.all")],
          row.names = FALSE,
          digits = 3
        )

        extra_cu <- paste0(mi_sig$lhs, " ~~ ", mi_sig$rhs)
        mod_m2_refined <- paste0(
          "reporter_p =~ ",
          paste(cr_items_p, collapse = " + "),
          "\n",
          "reporter_y =~ ",
          paste(cr_items_y, collapse = " + "),
          "\n",
          paste(extra_cu, collapse = "\n")
        )
        fit_m2_refined <- safe_cfa(mod_m2_refined, "M2_refined")

        if (!is.null(fit_m2_refined)) {
          m2_ref_fit_row <- fit_to_row(fit_m2_refined, fi_names) %>%
            mutate(model = "m2_refined", sex = sex_label)
          cat("\nM2_refined fit [", sex_label, "]:\n")
          print(m2_ref_fit_row)

          pe_ref <- tryCatch(
            lavaan::parameterEstimates(
              fit_m2_refined,
              ci = TRUE,
              standardized = TRUE
            ),
            error = function(e) NULL
          )
          if (!is.null(pe_ref)) {
            m2_ref_loadings <- pe_ref[
              pe_ref$op == "=~",
              c("lhs", "rhs", "std.all", "se", "ci.lower", "ci.upper", "pvalue")
            ] %>%
              dplyr::rename(factor = lhs, item = rhs, std_loading = std.all) %>%
              dplyr::mutate(sex = sex_label, model = "m2_refined")
          }
        }
      } else {
        cat(
          "\nM2 fit acceptable — no refinement paths added [",
          sex_label,
          "]\n"
        )
        mod_m2_refined <- paste0(
          "reporter_p =~ ",
          paste(cr_items_p, collapse = " + "),
          "\n",
          "reporter_y =~ ",
          paste(cr_items_y, collapse = " + ")
        )
        fit_m2_refined <- fit_m2
        m2_ref_loadings <- m2_loadings %>% dplyr::mutate(model = "m2_refined")
        m2_ref_fit_row <- fit_to_row(fit_m2, fi_names) %>%
          mutate(model = "m2_refined", sex = sex_label)
      }
    }
  }

  list(
    fit_table = fit_table,
    omega = sl,
    fit_m1 = fit_m1,
    fit_m2 = fit_m2,
    fit_m3 = fit_m3,
    fit_m4 = fit_m4,
    fit_m2_refined = fit_m2_refined,
    mod_m2_refined = mod_m2_refined,
    lat_cor = lat_cor,
    lat_cor_row = lat_cor_row,
    m2_loadings = m2_loadings,
    m2_ref_loadings = m2_ref_loadings,
    m2_ref_fit_row = m2_ref_fit_row
  )
}

bf_female <- run_bifactor_section(
  female_parent,
  female_youth,
  "female",
  include_pete = FALSE
)
bf_male <- run_bifactor_section(
  male_parent,
  male_youth,
  "male",
  include_pete = TRUE
)

bifactor_fits <- bind_rows(
  bf_female$fit_table,
  bf_male$fit_table
)
write.csv(
  bifactor_fits,
  file.path(out_dir, "bifactor_model_comparison.csv"),
  row.names = FALSE
)

# M2 latent correlation summary (both sexes)
lat_cor_summary <- bind_rows(
  bf_female$lat_cor_row,
  bf_male$lat_cor_row
) %>%
  dplyr::mutate(sex = c("female", "male")[seq_len(n())])
write.csv(
  lat_cor_summary,
  file.path(out_dir, "m2_latent_correlation.csv"),
  row.names = FALSE
)

# M2 standardized loadings (both sexes)
m2_loadings_all <- bind_rows(bf_female$m2_loadings, bf_male$m2_loadings)
write.csv(
  m2_loadings_all,
  file.path(out_dir, "m2_standardized_loadings.csv"),
  row.names = FALSE
)

# M2_refined loadings and fit row
m2_ref_loadings_all <- bind_rows(
  bf_female$m2_ref_loadings,
  bf_male$m2_ref_loadings
)
m2_ref_fit_all <- bind_rows(
  bf_female$m2_ref_fit_row,
  bf_male$m2_ref_fit_row
)
if (nrow(m2_ref_loadings_all) > 0) {
  write.csv(
    m2_ref_loadings_all,
    file.path(out_dir, "m2_refined_loadings.csv"),
    row.names = FALSE
  )
}
cat(
  "\nM2 results written to: m2_latent_correlation.csv,",
  "m2_standardized_loadings.csv, m2_refined_loadings.csv\n"
)

# --- Visual: M2_refined standardized loading plot ----------------------------
if (nrow(m2_ref_loadings_all) > 0) {
  load_plot_df <- m2_ref_loadings_all %>%
    mutate(
      reporter = dplyr::recode(
        factor,
        reporter_p = "Parent",
        reporter_y = "Youth"
      ),
      item_label = toupper(sub("_(p|y)$", "", item)),
      sex_label = dplyr::recode(sex, female = "Female", male = "Male"),
      item_label = factor(
        item_label,
        levels = rev(c("PETA", "PETB", "PETC", "PETD", "FPETE", "MPETE"))
      )
    )

  p_loads <- ggplot(
    load_plot_df,
    aes(
      x = std_loading,
      y = item_label,
      colour = reporter,
      xmin = ci.lower,
      xmax = ci.upper
    )
  ) +
    geom_point(size = 3, position = position_dodge(0.5)) +
    geom_errorbarh(height = 0.25, position = position_dodge(0.5)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    scale_colour_manual(
      values = c(Parent = "#2166ac", Youth = "#d73027"),
      name = "Reporter"
    ) +
    facet_wrap(~sex_label) +
    labs(
      title = "M2_refined: standardized factor loadings",
      subtitle = "Point estimate ± 95% CI",
      x = "Standardized loading",
      y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.subtitle = element_text(colour = "grey45", size = 11)
    )

  ggsave(
    file.path(out_dir, "m2_refined_loadings_plot.png"),
    p_loads,
    width = 9,
    height = 5,
    dpi = 180
  )
  cat("M2_refined loading plot written.\n")
}

# --- Visual: model fit comparison (M1, M2, M2_refined) ----------------------
fit_compare_df <- bind_rows(
  bf_female$fit_table %>%
    filter(model %in% c("m1_single_factor", "m2_corr_reporters")),
  m2_ref_fit_all,
  bf_male$fit_table %>%
    filter(model %in% c("m1_single_factor", "m2_corr_reporters"))
) %>%
  filter(!is.na(cfi.scaled)) %>%
  mutate(
    model_label = dplyr::recode(
      model,
      m1_single_factor = "M1: Single factor",
      m2_corr_reporters = "M2: Corr. reporters",
      m2_refined = "M2_refined"
    ),
    model_label = factor(
      model_label,
      levels = c("M1: Single factor", "M2: Corr. reporters", "M2_refined")
    ),
    sex_label = dplyr::recode(sex, female = "Female", male = "Male")
  ) %>%
  tidyr::pivot_longer(
    c(cfi.scaled, rmsea.scaled, srmr),
    names_to = "index",
    values_to = "value"
  ) %>%
  mutate(
    index = dplyr::recode(
      index,
      cfi.scaled = "CFI",
      rmsea.scaled = "RMSEA",
      srmr = "SRMR"
    )
  )

thresholds <- data.frame(
  index = c("CFI", "RMSEA", "SRMR"),
  good = c(0.95, 0.06, 0.08),
  direction = c("above", "below", "below")
)

p_fit <- ggplot(
  fit_compare_df,
  aes(x = model_label, y = value, fill = model_label)
) +
  geom_col(width = 0.6) +
  geom_hline(
    data = thresholds,
    aes(yintercept = good),
    linetype = "dashed",
    colour = "#d73027",
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "M1: Single factor" = "#d1e5f0",
      "M2: Corr. reporters" = "#4393c3",
      "M2_refined" = "#2166ac"
    ),
    guide = "none"
  ) +
  facet_grid(index ~ sex_label, scales = "free_y") +
  labs(
    title = "Model fit comparison: M1, M2, M2_refined",
    subtitle = "Dashed line = conventional threshold",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = 10),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey45", size = 11)
  )

ggsave(
  file.path(out_dir, "model_fit_comparison.png"),
  p_fit,
  width = 9,
  height = 7,
  dpi = 180
)
cat("Model fit comparison plot written.\n")

# ---------------------------------------------------------------------------
# SECTION 2b: CROSS-REPORTER LONGITUDINAL INVARIANCE
# Answers: does the multimethod model hold at each wave, and are the loadings
# stable enough across waves to support a consistent longitudinal model?
# Per-wave fits: M1–M3 at each wave.
# Multigroup invariance (wave as group): M1 (single factor) and M2 (correlated
# reporters) only — M3 too parameter-rich for 7-group invariance testing.
# ---------------------------------------------------------------------------

run_cr_longitudinal_invariance <- function(
  parent_df,
  youth_df,
  sex_label,
  include_pete = FALSE,
  mod_m2_refined = NULL # refined M2 syntax from run_bifactor_section
) {
  cat("\n\n=== Cross-reporter longitudinal invariance |", sex_label, "===\n")

  base_items <- c("peta", "petb", "petc", "petd")
  all_items <- if (include_pete) c(base_items, "fpete", "mpete") else base_items

  parent_sel <- parent_df %>%
    select(id, wave, any_of(all_items)) %>%
    rename_with(~ paste0(., "_p"), any_of(all_items))

  youth_sel <- youth_df %>%
    select(id, wave, any_of(all_items)) %>%
    rename_with(~ paste0(., "_y"), any_of(all_items))

  dat_long <- inner_join(parent_sel, youth_sel, by = c("id", "wave")) %>%
    mutate(
      across(-c(id, wave), as.integer),
      wave = factor(wave, levels = wave_order)
    ) %>%
    filter(if_all(-c(id, wave), ~ !is.na(.) & . %in% 1:4))

  cr_items <- setdiff(names(dat_long), c("id", "wave"))
  cr_items_p <- grep("_p$", cr_items, value = TRUE)
  cr_items_y <- grep("_y$", cr_items, value = TRUE)

  mod_m1 <- paste("puberty =~", paste(cr_items, collapse = " + "))
  mod_m2 <- paste0(
    "reporter_p =~ ",
    paste(cr_items_p, collapse = " + "),
    "\n",
    "reporter_y =~ ",
    paste(cr_items_y, collapse = " + ")
  )
  cu_within <- c(
    apply(combn(cr_items_p, 2), 2, function(x) paste0(x[1], " ~~ ", x[2])),
    apply(combn(cr_items_y, 2), 2, function(x) paste0(x[1], " ~~ ", x[2]))
  )
  mod_m3 <- paste0(
    "puberty =~ ",
    paste(cr_items, collapse = " + "),
    "\n",
    paste(cu_within, collapse = "\n")
  )
  fi_names <- c(
    "cfi.scaled",
    "tli.scaled",
    "rmsea.scaled",
    "srmr",
    "chisq.scaled",
    "df.scaled"
  )

  # --- Per-wave fits --------------------------------------------------------
  cat("\n--- Per-wave fits (M1–M3) ---\n")
  wave_fits <- lapply(wave_order, function(wv) {
    sub <- dat_long %>% filter(wave == wv) %>% select(all_of(cr_items))
    if (nrow(sub) < 100) {
      return(NULL)
    }
    model_list <- list(
      m1_single_factor = list(mod = mod_m1, theta = FALSE),
      m2_corr_reporters = list(mod = mod_m2, theta = FALSE),
      m3_cu_within = list(mod = mod_m3, theta = FALSE)
    )
    lapply(model_list, function(ml) {
      fit <- tryCatch(
        lavaan::cfa(
          ml$mod,
          data = sub,
          ordered = cr_items,
          estimator = "WLSMV",
          parameterization = if (ml$theta) "theta" else "delta",
          missing = "pairwise"
        ),
        error = function(e) NULL
      )
      row <- fit_to_row(fit, fi_names)
      if (!is.null(row)) mutate(row, wave = wv) else NULL
    }) %>%
      bind_rows(.id = "model") %>%
      mutate(sex = sex_label)
  }) %>%
    bind_rows()

  print_all(wave_fits)

  # --- Multigroup invariance across waves -----------------------------------
  # Run for M1 (single factor) and M2 (correlated reporters): the two most
  # stable models for multigroup invariance. M3 carries too many free
  # parameters to test reliably with 7 wave groups.
  dat_grp <- dat_long %>%
    select(-id) %>%
    filter(!is.na(wave)) %>%
    mutate(across(all_of(cr_items), ~ as.ordered(.)))

  fi_inv <- c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")
  inv_levels <- list(
    configural = character(0),
    metric = "loadings",
    scalar = c("loadings", "thresholds")
  )

  run_mg_inv <- function(mod, label, theta = FALSE) {
    fits <- lapply(inv_levels, function(ge) {
      tryCatch(
        lavaan::cfa(
          mod,
          data = dat_grp,
          group = "wave",
          ordered = cr_items,
          estimator = "WLSMV",
          parameterization = if (theta) "theta" else "delta",
          group.equal = ge
        ),
        error = function(e) {
          message(
            sex_label,
            " ",
            label,
            " [",
            paste(ge, collapse = "+"),
            "]: ",
            e$message
          )
          NULL
        }
      )
    })
    tbl <- lapply(fits, function(fit) fit_to_row(fit, fi_inv)) %>%
      bind_rows(.id = "model") %>%
      mutate(sex = sex_label, multimethod_model = label)
    cfg_cfi <- tbl$cfi.scaled[tbl$model == "configural"]
    tbl$delta_cfi <- if (length(cfg_cfi) == 1 && !is.na(cfg_cfi)) {
      tbl$cfi.scaled - cfg_cfi
    } else {
      NA_real_
    }
    tbl
  }

  cat("\n--- Multigroup invariance across waves: M1 (single factor) ---\n")
  inv_m1 <- run_mg_inv(mod_m1, "M1_single_factor")
  print(inv_m1)

  cat(
    "\n--- Multigroup invariance across waves: M2 (correlated reporters) ---\n"
  )
  inv_m2 <- run_mg_inv(mod_m2, "M2_corr_reporters")
  print(inv_m2)

  # M2_refined invariance — if a refined syntax was passed in
  inv_m2_refined <- NULL
  if (!is.null(mod_m2_refined) && !identical(mod_m2_refined, mod_m2)) {
    cat(
      "\n--- Multigroup invariance across waves: M2_refined ---\n"
    )
    inv_m2_refined <- run_mg_inv(mod_m2_refined, "M2_refined")
    print(inv_m2_refined)
  }

  list(
    wave_fits = wave_fits,
    inv_m1 = inv_m1,
    inv_m2 = inv_m2,
    inv_m2_refined = inv_m2_refined
  )
}

# Run with pete items — shows full picture including per-wave failures at
# baseline due to fpete floor effects
cr_long_female <- run_cr_longitudinal_invariance(
  female_parent,
  female_youth,
  "female",
  include_pete = TRUE
)
cr_long_male <- run_cr_longitudinal_invariance(
  male_parent,
  male_youth,
  "male",
  include_pete = TRUE
)

cr_wave_fits <- bind_rows(cr_long_female$wave_fits, cr_long_male$wave_fits)
write.csv(
  cr_wave_fits,
  file.path(out_dir, "cr_per_wave_fits_with_pete.csv"),
  row.names = FALSE
)

# Main invariance run: females exclude fpete ;
# Pass refined M2 syntax so invariance is tested on the final model.
cr_long_female_ord <- run_cr_longitudinal_invariance(
  female_parent,
  female_youth,
  "female",
  include_pete = FALSE,
  mod_m2_refined = bf_female$mod_m2_refined
)
cr_long_male_ord <- run_cr_longitudinal_invariance(
  male_parent,
  male_youth,
  "male",
  include_pete = TRUE,
  mod_m2_refined = bf_male$mod_m2_refined
)

cr_wave_fits_ord <- bind_rows(
  cr_long_female_ord$wave_fits,
  cr_long_male_ord$wave_fits
)
write.csv(
  cr_wave_fits_ord,
  file.path(out_dir, "cr_per_wave_fits_ordinal.csv"),
  row.names = FALSE
)

cr_inv_table <- bind_rows(
  cr_long_female_ord$inv_m1,
  cr_long_female_ord$inv_m2,
  cr_long_female_ord$inv_m2_refined,
  cr_long_male_ord$inv_m1,
  cr_long_male_ord$inv_m2,
  cr_long_male_ord$inv_m2_refined
)
write.csv(
  cr_inv_table,
  file.path(out_dir, "cr_longitudinal_invariance.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# SECTION 3: LONGITUDINAL MEASUREMENT INVARIANCE
# Group = wave (7 annual timepoints). Per sex × reporter.
# Configural → metric → scalar → strict (WLSMV, theta parameterization)
# ---------------------------------------------------------------------------

run_longitudinal_invariance <- function(df, label, items = ordinal_items) {
  cat("\n\n=== Longitudinal invariance:", label, "===\n")

  dat <- df %>%
    select(wave, any_of(items)) %>%
    filter(wave %in% wave_order) %>%
    mutate(
      wave = factor(wave, levels = wave_order),
      across(any_of(items), ~ as.ordered(as.integer(.)))
    )
  items_present <- intersect(items, names(dat))
  dat <- dat %>% filter(if_all(all_of(items_present), ~ !is.na(.)))

  if (nrow(dat) < 500) {
    message("Skipping ", label, ": insufficient data")
    return(NULL)
  }

  syntax <- paste("pub =~", paste(items_present, collapse = " + "))

  models <- list(
    configural = list(group.equal = character(0)),
    metric = list(group.equal = "loadings"),
    scalar = list(group.equal = c("loadings", "thresholds")),
    strict = list(group.equal = c("loadings", "thresholds", "residuals"))
  )

  fits <- lapply(models, function(args) {
    tryCatch(
      lavaan::cfa(
        syntax,
        data = dat,
        group = "wave",
        ordered = items_present,
        estimator = "WLSMV",
        parameterization = "theta",
        group.equal = args$group.equal
      ),
      error = function(e) {
        message(label, " | ", e$message)
        NULL
      }
    )
  })

  fi_names <- c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")
  fit_table <- lapply(fits, function(fit) fit_to_row(fit, fi_names)) %>%
    bind_rows(.id = "model") %>%
    mutate(group = label)

  if (nrow(fit_table) == 0 || !"cfi.scaled" %in% names(fit_table)) {
    return(NULL)
  }

  # delta CFI from configural (threshold: < -.010)
  fit_table <- fit_table %>%
    mutate(delta_cfi = cfi.scaled - cfi.scaled[model == "configural"])

  cat("\n--- Fit indices ---\n")
  print(fit_table)

  # LRT (DIFFTEST for WLSMV)
  non_null <- Filter(Negate(is.null), fits)
  if (length(non_null) >= 2) {
    lrt <- tryCatch(
      do.call(
        lavaan::lavTestLRT,
        c(non_null, list(model.names = names(non_null)))
      ),
      error = function(e) NULL
    )
    if (!is.null(lrt)) {
      cat("\n--- LRT ---\n")
      print(lrt)
    }
  }

  list(fits = fits, fit_table = fit_table)
}

long_inv <- list(
  female_parent = run_longitudinal_invariance(
    female_parent,
    "female_parent",
    items = female_items
  ),
  female_youth = run_longitudinal_invariance(
    female_youth,
    "female_youth",
    items = female_items
  ),
  male_parent = run_longitudinal_invariance(
    male_parent,
    "male_parent",
    items = male_items
  ),
  male_youth = run_longitudinal_invariance(
    male_youth,
    "male_youth",
    items = male_items
  )
)

long_inv_table <- bind_rows(lapply(long_inv, `[[`, "fit_table"))
write.csv(
  long_inv_table,
  file.path(out_dir, "longitudinal_invariance.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# SECTION 4: INVARIANCE ACROSS BACKGROUND COVARIATES
# Groups: age tertile, race/ethnicity, BMI tertile
# Per sex × reporter. Configural vs. scalar only (delta CFI criterion).
# ---------------------------------------------------------------------------

make_tertile <- function(x) {
  cuts <- quantile(x, probs = c(1 / 3, 2 / 3), na.rm = TRUE)
  case_when(
    x <= cuts[1] ~ "low",
    x <= cuts[2] ~ "mid",
    TRUE ~ "high"
  )
}

run_covariate_invariance <- function(
  df,
  covariate,
  group_label,
  ds_label,
  items = ordinal_items
) {
  dat <- df %>%
    filter(!is.na(.data[[covariate]])) %>%
    mutate(
      group_var = .data[[covariate]],
      across(any_of(items), ~ as.ordered(as.integer(.)))
    )
  items_present <- intersect(items, names(dat))
  dat <- dat %>%
    filter(if_all(all_of(items_present), ~ !is.na(.))) %>%
    select(group_var, all_of(items_present))

  if (length(unique(dat$group_var)) < 2) {
    message("Skipping ", ds_label, " | ", group_label, ": < 2 groups")
    return(NULL)
  }

  syntax <- paste("pub =~", paste(items_present, collapse = " + "))

  fit_config <- tryCatch(
    lavaan::cfa(
      syntax,
      data = dat,
      group = "group_var",
      ordered = items_present,
      estimator = "WLSMV",
      parameterization = "theta"
    ),
    error = function(e) NULL
  )
  fit_scalar <- tryCatch(
    lavaan::cfa(
      syntax,
      data = dat,
      group = "group_var",
      ordered = items_present,
      estimator = "WLSMV",
      parameterization = "theta",
      group.equal = c("loadings", "thresholds")
    ),
    error = function(e) NULL
  )

  fi_names <- c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")
  row_config <- fit_to_row(fit_config, fi_names)
  row_scalar <- fit_to_row(fit_scalar, fi_names)

  result <- bind_rows(
    if (!is.null(row_config)) {
      mutate(row_config, model = "configural")
    } else {
      NULL
    },
    if (!is.null(row_scalar)) mutate(row_scalar, model = "scalar") else NULL
  )

  if (nrow(result) == 0 || !"cfi.scaled" %in% names(result)) {
    return(NULL)
  }

  result %>%
    mutate(
      covariate = group_label,
      dataset = ds_label,
      delta_cfi = cfi.scaled - cfi.scaled[model == "configural"]
    )
}

# add tertile columns before testing
add_group_cols <- function(df) {
  df <- df %>%
    mutate(
      age_tertile = make_tertile(age),
      bmi_tertile = make_tertile(bmi),
      race_group = as.character(race)
    )
  if ("site" %in% names(df)) {
    df <- df %>% mutate(site_group = as.character(site))
  } else {
    message(
      "'site' column missing — re-run 00_data_foundation.R to enable site invariance"
    )
    df$site_group <- NA_character_
  }
  df
}

female_parent <- add_group_cols(female_parent)
female_youth <- add_group_cols(female_youth)
male_parent <- add_group_cols(male_parent)
male_youth <- add_group_cols(male_youth)

# Report site availability upfront so it's obvious if it will be skipped
site_available <- !all(is.na(female_parent$site_group))
cat(
  "\nSite invariance:",
  if (site_available) {
    "ENABLED (site column found)"
  } else {
    "SKIPPED — 'site' column missing; re-run 00_data_foundation.R to include it"
  },
  "\n"
)

cov_inv_results <- list()
for (ds_name in names(datasets)) {
  df_aug <- get(ds_name) # already has group cols from add_group_cols above
  ds_items <- datasets[[ds_name]]$items
  for (cov in c("age_tertile", "bmi_tertile", "race_group", "site_group")) {
    res <- run_covariate_invariance(df_aug, cov, cov, ds_name, items = ds_items)
    if (!is.null(res)) cov_inv_results[[paste(ds_name, cov, sep = "_")]] <- res
  }
}

cov_inv_table <- bind_rows(cov_inv_results)
cat("\n\n=== Invariance across background covariates ===\n")
print_all(cov_inv_table)
write.csv(
  cov_inv_table,
  file.path(out_dir, "covariate_invariance.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# SUMMARY PLOT: item discrimination heatmap across waves
# ---------------------------------------------------------------------------
if (nrow(disc_table) > 0) {
  disc_plot_df <- disc_table %>%
    mutate(
      wave = factor(wave, levels = wave_order),
      item = factor(
        item,
        levels = rev(c("peta", "petb", "petc", "petd", "mpete"))
      ),
      reporter_label = dplyr::recode(
        reporter,
        parent = "Parent report",
        youth = "Youth report"
      ),
      sex_label = dplyr::recode(sex, female = "Female", male = "Male"),
      panel = paste0(sex_label, "\n", reporter_label)
    ) %>%
    mutate(
      panel = factor(
        panel,
        levels = c(
          "Female\nParent report",
          "Female\nYouth report",
          "Male\nParent report",
          "Male\nYouth report"
        )
      )
    )

  a_max <- ceiling(max(disc_plot_df$a, na.rm = TRUE))

  p_disc <- ggplot(disc_plot_df, aes(x = wave, y = item, fill = a)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.1f", a)), size = 3.8, colour = "grey15") +
    scale_fill_gradient2(
      low = "#fee090",
      mid = "#1a9850",
      high = "#005a32",
      midpoint = 2,
      limits = c(0, a_max),
      name = "a",
      breaks = c(0, 1, 2, 3, 4),
      labels = c("0", "1", "2", "3", "4")
    ) +
    facet_wrap(~panel, nrow = 1) +
    labs(
      title = "GRM item discrimination across waves",
      subtitle = "a: discrimination parameter  │  higher = better item–trait separation",
      x = "Wave",
      y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.key.width = unit(1.8, "cm"),
      plot.subtitle = element_text(size = 11, colour = "grey45"),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 12, face = "italic")
    )

  ggsave(
    file.path(out_dir, "discrimination_heatmap.png"),
    p_disc,
    width = 13,
    height = 5,
    dpi = 180
  )
}

# ---------------------------------------------------------------------------
# SUMMARY PLOT: covariate invariance heatmap
# ---------------------------------------------------------------------------
if (nrow(cov_inv_table) > 0 && "delta_cfi" %in% names(cov_inv_table)) {
  inv_plot_df <- cov_inv_table %>%
    filter(model == "scalar", !is.na(delta_cfi)) %>%
    mutate(
      fails = delta_cfi < -0.010,
      covariate_label = dplyr::recode(
        covariate,
        age_tertile = "Age tertile",
        bmi_tertile = "BMI tertile",
        race_group = "Race / ethnicity",
        site_group = "Site"
      ),
      dataset_label = dplyr::recode(
        dataset,
        female_parent = "Female\nparent",
        female_youth = "Female\nyouth",
        male_parent = "Male\nparent",
        male_youth = "Male\nyouth"
      ),
      dataset_label = factor(
        dataset_label,
        levels = c(
          "Female\nparent",
          "Female\nyouth",
          "Male\nparent",
          "Male\nyouth"
        )
      )
    )

  p_inv <- ggplot(
    inv_plot_df,
    aes(x = dataset_label, y = covariate_label, fill = delta_cfi)
  ) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(
      aes(
        label = sprintf("%+.3f", delta_cfi),
        colour = fails
      ),
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#d73027",
      mid = "#fee090",
      high = "#1a9850",
      midpoint = -0.010,
      limits = c(min(inv_plot_df$delta_cfi, na.rm = TRUE) - 0.005, 0),
      name = "ΔCFI",
      breaks = c(-0.04, -0.02, -0.01, 0),
      labels = c("−.04", "−.02", "−.01", "0")
    ) +
    scale_colour_manual(
      values = c("FALSE" = "grey25", "TRUE" = "#a50026"),
      guide = "none"
    ) +
    labs(
      title = "Measurement invariance across background covariates",
      subtitle = "ΔCFI (scalar − configural)  │  red bold = |CFI| > .010 (non-invariant)",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(1.8, "cm"),
      plot.subtitle = element_text(size = 11, colour = "grey45"),
      axis.text = element_text(size = 13)
    )

  ggsave(
    file.path(out_dir, "covariate_invariance_heatmap.png"),
    p_inv,
    width = 9,
    height = 5,
    dpi = 180
  )
  cat("\nCovariate invariance heatmap written.\n")
}

# ---------------------------------------------------------------------------
# SUMMARY PLOT: longitudinal invariance heatmap
# ---------------------------------------------------------------------------
if (length(long_inv) > 0) {
  long_inv_plot_df <- long_inv_table %>%
    filter(model != "configural", !is.na(delta_cfi)) %>%
    mutate(
      model = factor(model, levels = c("metric", "scalar", "strict")),
      group_label = dplyr::recode(
        group,
        female_parent = "Female\nparent",
        female_youth = "Female\nyouth",
        male_parent = "Male\nparent",
        male_youth = "Male\nyouth"
      ),
      group_label = factor(
        group_label,
        levels = c(
          "Female\nparent",
          "Female\nyouth",
          "Male\nparent",
          "Male\nyouth"
        )
      ),
      fails = delta_cfi < -0.010
    )

  p_long_inv <- ggplot(
    long_inv_plot_df,
    aes(x = group_label, y = model, fill = delta_cfi)
  ) +
    geom_tile(colour = "white", linewidth = 1) +
    geom_text(
      aes(
        label = sprintf("%+.3f", delta_cfi),
        colour = fails
      ),
      size = 5,
      fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#d73027",
      mid = "#fee090",
      high = "#1a9850",
      midpoint = -0.010,
      limits = c(min(long_inv_plot_df$delta_cfi, na.rm = TRUE) - 0.01, 0),
      name = "ΔCFI vs. configural",
      breaks = c(-0.15, -0.10, -0.05, -0.01, 0),
      labels = c("−.15", "−.10", "−.05", "−.01", "0")
    ) +
    scale_colour_manual(
      values = c("FALSE" = "grey25", "TRUE" = "#a50026"),
      guide = "none"
    ) +
    scale_y_discrete(limits = rev) +
    labs(
      title = "Longitudinal measurement invariance across 7 waves",
      subtitle = "ΔCFI (vs. configural baseline)  │  red bold = |ΔCFI| > .010 (non-invariant)",
      x = NULL,
      y = "Invariance level"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = unit(1.8, "cm"),
      plot.subtitle = element_text(size = 11, colour = "grey45"),
      axis.text = element_text(size = 13)
    )

  ggsave(
    file.path(out_dir, "longitudinal_invariance_heatmap.png"),
    p_long_inv,
    width = 9,
    height = 5,
    dpi = 180
  )
  cat("\nLongitudinal invariance heatmap written.\n")
}

cat("\n\nAll psychometrics output written to:", out_dir, "\n")
