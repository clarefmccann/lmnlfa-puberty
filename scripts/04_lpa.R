## 04_lpa.R
## Latent Profile Analysis on GAMM-derived pubertal timing and tempo.
## Uses tidyLPA (wrapping mclust) for model comparison and profile extraction.
##
## Inputs:  {ds_name}_timing_tempo.csv  (from 03_gamms_hpc.R)
## Outputs: model comparison tables, profile assignments, plots
##
## Runs 6 LPA solutions:
##   female_parent, female_youth, male_parent, male_youth  (sex × reporter)
##   female_averaged, male_averaged                        (reporter-averaged)
##
## Model selection: compare covariance structures × K profiles; BIC-optimal
## solution is saved; comparison plots allow manual override.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tidyLPA)
})

set.seed(90025)

# timing_dev  = timing - cohort mean timing (years; negative = early)  [individual]
# acceleration = d²(PDS)/d(age²) at onset (PDS = 2.5): how fast the growth
#               rate is changing at the moment puberty begins; negative values
#               are meaningful (growth rate already decelerating at onset)  [individual]
# Excluded:
#   obs_tempo / gamm_tempo — OLS slope and GAMM random slope; both conflate
#     overall growth rate with timing and are collinear with timing_dev.
#   peak_velocity — max d(PDS)/d(age); near-identical across individuals once
#     timing is accounted for.
LPA_VARS <- c("timing_dev", "acceleration")
N_PROFILES <- 1:6
# tidyLPA model numbers (mclust parameterisation):
#   1 = EEI  equal variance, zero covariance
#   2 = VVI  varying variance, zero covariance
#   3 = EEE  equal variance + covariance
# Model 6 (VVV) excluded: with only 3 LPA variables the full covariance matrix
# is overparameterized, produces negative BIC, and collapses all variance into
# the covariance structure rather than separating profiles on tempo/velocity.
LPA_MODELS <- c(1, 2, 3)

# ---------------------------------------------------------------------------
# PATHS
# Checks: DATA_DIR / OUT_DIR env vars (set in .sh) → Box local fallback
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
    "projects/abcd-projs",
    "dissertation/study1/outputs"
  )
}
gamm_dir <- file.path(out_base, "gamm")
out_dir <- file.path(out_base, "lpa")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# LOAD timing/tempo from all four datasets
# ---------------------------------------------------------------------------
ds_names <- c("female_parent", "female_youth", "male_parent", "male_youth")

tt_list <- lapply(ds_names, function(ds) {
  path <- file.path(gamm_dir, paste0(ds, "_timing_tempo.csv"))
  if (!file.exists(path)) {
    warning("Not found: ", path)
    return(NULL)
  }
  read.csv(path)
})
names(tt_list) <- ds_names

# ---------------------------------------------------------------------------
# BUILD AVERAGED DATASETS (parent + youth averaged per person per sex)
# ---------------------------------------------------------------------------
build_averaged <- function(sex) {
  p <- tt_list[[paste0(sex, "_parent")]]
  y <- tt_list[[paste0(sex, "_youth")]]

  if (is.null(p) && is.null(y)) {
    return(NULL)
  }
  if (is.null(p)) {
    return(y %>% mutate(reporter = "youth_only"))
  }
  if (is.null(y)) {
    return(p %>% mutate(reporter = "parent_only"))
  }

  joined <- full_join(
    p %>%
      select(id, all_of(LPA_VARS)) %>%
      rename_with(~ paste0(., "_p"), all_of(LPA_VARS)),
    y %>%
      select(id, all_of(LPA_VARS)) %>%
      rename_with(~ paste0(., "_y"), all_of(LPA_VARS)),
    by = "id"
  )

  # Average each LPA variable across reporters without hardcoding column names
  for (.v in LPA_VARS) {
    joined[[.v]] <- rowMeans(
      cbind(joined[[paste0(.v, "_p")]], joined[[paste0(.v, "_y")]]),
      na.rm = TRUE
    )
  }

  joined %>%
    mutate(dataset = paste0(sex, "_averaged"), reporter = "averaged") %>%
    select(id, dataset, reporter, all_of(LPA_VARS))
}

tt_list[["female_averaged"]] <- build_averaged("female")
tt_list[["male_averaged"]] <- build_averaged("male")

all_labels <- c(ds_names, "female_averaged", "male_averaged")

# ---------------------------------------------------------------------------
# Extract class assignments from a tidyLPA fit without relying on get_data()
# storing the original data (some tidyLPA versions don't attach it).
# ---------------------------------------------------------------------------
extract_assignments <- function(fit) {
  # Try the official route first
  result <- tryCatch(
    get_data(fit) %>% select(Class, starts_with("CPROB")),
    error = function(e) NULL
  )
  if (!is.null(result)) {
    return(result)
  }

  # Fall back: pull the mclust object from the tidyProfile/tidyLPA structure
  mc <- NULL
  for (.acc in list(
    function(f) f$fit,
    function(f) f[[1]]$fit,
    function(f) f$model_object,
    function(f) f[[1]]$model_object
  )) {
    mc <- tryCatch(.acc(fit), error = function(e) NULL)
    if (!is.null(mc) && inherits(mc, "Mclust")) {
      break
    }
    mc <- NULL
  }
  if (is.null(mc)) {
    stop("Cannot extract class assignments from tidyLPA fit.")
  }

  probs <- as.data.frame(mc$z)
  names(probs) <- paste0("CPROB", seq_len(ncol(mc$z)))
  bind_cols(data.frame(Class = as.integer(mc$classification)), probs)
}

# ---------------------------------------------------------------------------
# LPA HELPER
# ---------------------------------------------------------------------------
run_lpa <- function(df, label) {
  df_cc <- df %>% filter(complete.cases(.[LPA_VARS]))
  cat("\n", strrep("=", 60), "\n")
  cat("LPA:", label, "| complete cases:", nrow(df_cc), "\n")
  cat(strrep("=", 60), "\n")

  if (nrow(df_cc) < 50) {
    warning("Too few complete cases for ", label, " — skipping.")
    return(invisible(NULL))
  }

  # Exclude participants with negative OLS slope — Tanner stages do not reverse;
  # negative values are measurement artifacts (sparse data, caregiver switches).
  # NOTE: "acceleration" (d²PDS/dt² at onset) is intentionally excluded from
  # this filter — negative values there are biologically meaningful (growth rate
  # already decelerating at onset) and should not be dropped.
  slope_var <- grep(
    "^obs_tempo$|^ols_slope$|^tempo$",
    LPA_VARS,
    value = TRUE,
    ignore.case = TRUE
  )
  if (length(slope_var) > 0) {
    n_neg <- sum(df_cc[[slope_var[1]]] < 0, na.rm = TRUE)
    if (n_neg > 0) {
      cat(
        "Excluding",
        n_neg,
        "participants with negative",
        slope_var[1],
        "(OLS < 0)\n"
      )
      df_cc <- df_cc %>% filter(.data[[slope_var[1]]] >= 0)
    }
  }

  if (nrow(df_cc) < 50) {
    warning("Too few cases after slope exclusion for ", label, " — skipping.")
    return(invisible(NULL))
  }

  # Winsorize at ±3 SD before z-scoring.
  # obs_tempo has extreme right outliers (max ~10 SD above median); without
  # winsorizing, scale() maps the bulk of the distribution to z ≈ 0 and LPA
  # cannot distinguish profiles on tempo.
  winsorise <- function(x, k = 3) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    pmax(pmin(x, m + k * s), m - k * s)
  }
  df_cc <- df_cc %>%
    mutate(across(all_of(LPA_VARS), winsorise))

  # Standardise within this sample (LPA on z-scores; raw means saved separately)
  df_z <- df_cc %>%
    mutate(across(all_of(LPA_VARS), ~ as.numeric(scale(.))))

  # --- Model comparison ---------------------------------------------------
  # Estimate each model/k individually so failed mclust fits don't crash the
  # whole comparison (compare_solutions() chokes on NULL fits with empty names).
  cat("Comparing models (this may take a few minutes)...\n")
  successful_fits <- list()
  for (.m in LPA_MODELS) {
    for (.k in N_PROFILES) {
      .key <- paste0("m", .m, "_k", .k)
      .fit <- tryCatch(
        suppressWarnings(
          df_z %>%
            select(all_of(LPA_VARS)) %>%
            estimate_profiles(n_profiles = .k, models = .m)
        ),
        error = function(e) NULL
      )
      if (!is.null(.fit)) {
        .stats <- tryCatch(get_fit(.fit), error = function(e) NULL)
        if (!is.null(.stats) && nrow(.stats) > 0) {
          # Column names vary by tidyLPA/tibble version; match by pattern.
          .bic_col <- grep("^BIC", names(.stats), value = TRUE)[1]
          .ent_col <- grep(
            "(?i)^entropy",
            names(.stats),
            value = TRUE,
            perl = TRUE
          )[1]
          .bic_val <- if (!is.na(.bic_col)) {
            as.numeric(.stats[[.bic_col]][1])
          } else {
            NA_real_
          }
          .ent_val <- if (!is.na(.ent_col)) {
            as.numeric(.stats[[.ent_col]][1])
          } else {
            NA_real_
          }
          successful_fits[[.key]] <- list(
            fit = .fit,
            model = .m,
            k = .k,
            BIC = .bic_val,
            Entropy = .ent_val
          )
        }
      }
    }
  }

  if (length(successful_fits) == 0) {
    warning("No models converged for ", label, " — skipping.")
    return(invisible(NULL))
  }

  fits_df <- bind_rows(lapply(successful_fits, function(x) {
    data.frame(
      Model = x$model,
      Classes = x$k,
      BIC = x$BIC,
      Entropy = x$Entropy,
      stringsAsFactors = FALSE
    )
  }))

  cat("\nFit statistics:\n")
  print(fits_df)

  write.csv(
    fits_df,
    file.path(out_dir, paste0(label, "_model_comparison.csv")),
    row.names = FALSE
  )

  # BIC plot (standard BIC: lower = better)
  p_bic <- ggplot(
    fits_df,
    aes(x = Classes, y = BIC, colour = factor(Model), group = factor(Model))
  ) +
    geom_line() +
    geom_point(size = 2) +
    scale_x_continuous(breaks = N_PROFILES) +
    labs(
      title = paste("BIC by profiles and model —", label),
      x = "Number of profiles",
      y = "BIC (lower = better)",
      colour = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  # Entropy plot
  p_ent <- ggplot(
    fits_df %>% filter(Classes > 1),
    aes(x = Classes, y = Entropy, colour = factor(Model), group = factor(Model))
  ) +
    geom_line() +
    geom_point(size = 2) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
    scale_x_continuous(breaks = N_PROFILES) +
    labs(
      title = paste("Entropy by profiles and model —", label),
      x = "Number of profiles",
      y = "Entropy (>0.80 = good separation)",
      colour = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(
    file.path(out_dir, paste0(label, "_bic_plot.png")),
    p_bic,
    width = 8,
    height = 5,
    dpi = 150
  )
  ggsave(
    file.path(out_dir, paste0(label, "_entropy_plot.png")),
    p_ent,
    width = 8,
    height = 5,
    dpi = 150
  )

  # --- Best solution (BIC-optimal) ----------------------------------------
  # get_fit() returns standard BIC (lower = better); exclude K=1 since a
  # single-component solution is not a meaningful LPA result.
  valid_fits <- Filter(function(x) !is.na(x$BIC) && x$k > 1, successful_fits)
  if (length(valid_fits) == 0) {
    warning("No valid multi-profile BIC values for ", label, " — skipping.")
    return(invisible(NULL))
  }
  best_key <- names(which.min(sapply(valid_fits, `[[`, "BIC")))
  best_k <- valid_fits[[best_key]]$k
  best_model <- valid_fits[[best_key]]$model
  fit_best <- valid_fits[[best_key]]$fit
  cat("\nBIC-optimal solution: K =", best_k, "| model =", best_model, "\n")

  assigns <- extract_assignments(fit_best) %>%
    bind_cols(df_cc %>% select(id, all_of(LPA_VARS)))

  write.csv(
    assigns,
    file.path(out_dir, paste0(label, "_k", best_k, "_assignments.csv")),
    row.names = FALSE
  )

  # --- Profile means plot (z-score scale, faceted) ------------------------
  # Raw-scale plot is misleading: timing (~8-15 yr) and obs_tempo (~0.2-0.4
  # PDS/yr) are on incompatible scales; tempo differences are invisible on a
  # shared y-axis. Faceted z-score plot gives each variable its own axis.
  profile_means_raw <- assigns %>%
    group_by(Class) %>%
    summarise(
      n = n(),
      across(
        all_of(LPA_VARS),
        list(mean = mean, se = ~ sd(.) / sqrt(n())),
        .names = "{.col}__{.fn}"
      ),
      .groups = "drop"
    )

  # Also compute profile means on the z-scored data (what LPA actually saw)
  assigns_z <- assigns %>%
    mutate(across(all_of(LPA_VARS), ~ as.numeric(scale(.))))
  profile_means_z <- assigns_z %>%
    group_by(Class) %>%
    summarise(
      n = n(),
      across(
        all_of(LPA_VARS),
        list(mean = mean, se = ~ sd(.) / sqrt(n())),
        .names = "{.col}__{.fn}"
      ),
      .groups = "drop"
    )

  write.csv(
    profile_means_raw,
    file.path(out_dir, paste0(label, "_k", best_k, "_profile_means.csv")),
    row.names = FALSE
  )

  var_labels <- c(
    "Timing deviation (yr from cohort mean)",
    "Acceleration (d²PDS/dt² at onset)"
  )

  # Pivot z-score means for faceted plot
  pm_long_z <- profile_means_z %>%
    pivot_longer(
      cols = ends_with("__mean") | ends_with("__se"),
      names_to = c("variable", ".value"),
      names_sep = "__"
    ) %>%
    mutate(
      variable = factor(variable, levels = LPA_VARS, labels = var_labels),
      Class = factor(Class)
    )

  p_prof <- ggplot(
    pm_long_z,
    aes(x = Class, y = mean, colour = Class, group = Class)
  ) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = mean - 2 * se, ymax = mean + 2 * se),
      width = 0.2,
      linewidth = 0.7
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    facet_wrap(~variable, scales = "free_y") +
    labs(
      title = paste0(
        "Pubertal profiles — ",
        label,
        " (K = ",
        best_k,
        ", model ",
        best_model,
        ")"
      ),
      x = "Profile",
      y = "Mean (z-score ± 2 SE)",
      colour = "Profile",
      caption = paste0(
        "n per profile: ",
        paste(profile_means_raw$n, collapse = " / ")
      )
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")

  ggsave(
    file.path(out_dir, paste0(label, "_k", best_k, "_profile_plot.png")),
    p_prof,
    width = 8,
    height = 5,
    dpi = 150
  )

  # --- Also save K = 2,3,4 of best model type for inspection --------------
  for (k in setdiff(2:4, best_k)) {
    .key_k <- paste0("m", best_model, "_k", k)
    fit_k <- if (!is.null(successful_fits[[.key_k]])) {
      successful_fits[[.key_k]]$fit
    } else {
      tryCatch(
        suppressWarnings(
          df_z %>%
            select(all_of(LPA_VARS)) %>%
            estimate_profiles(n_profiles = k, models = best_model)
        ),
        error = function(e) NULL
      )
    }
    if (is.null(fit_k)) {
      next
    }

    asgn_k <- extract_assignments(fit_k) %>%
      bind_cols(df_cc %>% select(id, all_of(LPA_VARS)))

    write.csv(
      asgn_k,
      file.path(out_dir, paste0(label, "_k", k, "_assignments.csv")),
      row.names = FALSE
    )
  }

  cat("Outputs written for", label, "\n")
  invisible(list(
    comparison = fits_df,
    best_fit = fit_best,
    best_k = best_k,
    best_model = best_model
  ))
}

# ---------------------------------------------------------------------------
# RUN ALL SIX LPAs
# ---------------------------------------------------------------------------
results <- list()
for (label in all_labels) {
  df <- tt_list[[label]]
  if (is.null(df)) {
    next
  }
  results[[label]] <- run_lpa(df, label)
}

# ---------------------------------------------------------------------------
# SUMMARY TABLE across all solutions
# ---------------------------------------------------------------------------
summary_rows <- lapply(names(results), function(label) {
  r <- results[[label]]
  if (is.null(r)) {
    return(NULL)
  }
  data.frame(
    label = label,
    best_K = r$best_k,
    best_model = r$best_model,
    stringsAsFactors = FALSE
  )
})
summary_tbl <- bind_rows(Filter(Negate(is.null), summary_rows))
cat("\n=== BIC-optimal solutions across all runs ===\n")
print(summary_tbl)
write.csv(
  summary_tbl,
  file.path(out_dir, "lpa_best_solutions_summary.csv"),
  row.names = FALSE
)

cat("\nDone. All outputs written to:", out_dir, "\n")
