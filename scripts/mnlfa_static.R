## mnlfa_static.R
## Stage 1+2 convergence check: single-wave (baseline), cross-sectional MNLFA
## with mean + variance impact of age on the latent factor. No DIF, no
## growth structure. Run this before lmnlfa_hpc.R to confirm the basic
## measurement + impact model converges; add complexity from here in stages
## (random-intercept-only longitudinal -> full growth -> DIF).
##
## Usage:
##   Rscript mnlfa_static.R <sex>
##   sex: female | male
##
## Stan:    scripts/stan/mnlfa-static.stan
## Outputs: written to <OUT_DIR>/lmnlfa_static/

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
  stop("Usage: Rscript mnlfa_static.R <female|male>")
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
out_dir <- file.path(out_base, "lmnlfa_static")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

script_dir <- tryCatch(
  {
    ofile <- sys.frames()[[1]]$ofile
    if (is.null(ofile)) "scripts" else dirname(ofile)
  },
  error = function(e) "scripts"
)
stan_file <- file.path(script_dir, "stan", "mnlfa-static.stan")
if (!file.exists(stan_file)) {
  stan_file <- file.path("scripts", "stan", "mnlfa-static.stan")
}
if (!file.exists(stan_file)) {
  stop("Cannot find mnlfa-static.stan: ", stan_file)
}

cat("Data dir:  ", data_dir, "\n")
cat("Output dir:", out_dir, "\n")
cat("Stan file: ", stan_file, "\n")

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------
parent_df <- read.csv(file.path(data_dir, paste0(sx, "_parent_long.csv")))
youth_df <- read.csv(file.path(data_dir, paste0(sx, "_youth_long.csv")))

# ---------------------------------------------------------------------------
# COMPILE STAN MODEL
# ---------------------------------------------------------------------------
cat("\nCompiling Stan model...\n")
stan_model <- cmdstan_model(stan_file)
cat("Compiled.\n")

# ---------------------------------------------------------------------------
# BUILD CROSS-SECTIONAL (BASELINE-ONLY) DATA
# ---------------------------------------------------------------------------
build_static_data <- function(
  parent_df,
  youth_df,
  sex_label,
  ordinal_items = c("peta", "petb", "petc", "petd")
) {
  cat("\n=== Building static MNLFA data |", sex_label, "===\n")

  parent_sel <- parent_df %>%
    filter(wave == "bl") %>%
    select(id, age, any_of(ordinal_items)) %>%
    rename_with(~ paste0(., "_p"), any_of(ordinal_items))

  youth_sel <- youth_df %>%
    filter(wave == "bl") %>%
    select(id, any_of(ordinal_items)) %>%
    rename_with(~ paste0(., "_y"), any_of(ordinal_items))

  all_item_cols <- c(paste0(ordinal_items, "_p"), paste0(ordinal_items, "_y"))

  dat <- inner_join(parent_sel, youth_sel, by = "id") %>%
    filter(!is.na(age)) %>%
    filter(if_all(
      all_of(all_item_cols),
      ~ !is.na(.) & as.integer(.) %in% 1:4
    ))

  if (nrow(dat) < 200) {
    stop("Insufficient baseline data for ", sex_label)
  }

  ids <- sort(unique(dat$id))
  age_mean <- mean(dat$age, na.rm = TRUE)

  dat <- dat %>%
    mutate(
      person_idx = match(id, ids),
      age_c = age - age_mean
    ) %>%
    arrange(person_idx)

  age_by_person <- dat %>%
    distinct(person_idx, age_c) %>%
    arrange(person_idx)

  dat_long <- dat %>%
    select(id, person_idx, all_of(all_item_cols)) %>%
    pivot_longer(
      cols = all_of(all_item_cols),
      names_to = "item",
      values_to = "y_raw"
    ) %>%
    mutate(item_idx = match(item, all_item_cols), y_int = as.integer(y_raw)) %>%
    arrange(person_idx, item_idx)

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
    "\n"
  )

  list(
    dat_long = dat_long,
    item_names = all_item_cols,
    ids = ids,
    age_mean = age_mean,
    age_c_person = age_by_person$age_c,
    ni = length(ids),
    p = length(all_item_cols),
    nobs = nrow(dat_long),
    k_items = rep(4L, length(all_item_cols)),
    k_max = 4L
  )
}

prep <- build_static_data(parent_df, youth_df, sx)

stan_data <- list(
  nobs = prep$nobs,
  p = prep$p,
  ni = prep$ni,
  person = prep$dat_long$person_idx,
  itm = prep$dat_long$item_idx,
  y = prep$dat_long$y_int,
  k_item = prep$k_items,
  k_max = prep$k_max,
  age_c = prep$age_c_person,
  sigma_l = 1.5,
  sigma_nu = 1.5,
  sigma_f = 1.5
)

# ---------------------------------------------------------------------------
# FIT
# ---------------------------------------------------------------------------
cat("\n--- Fitting static MNLFA (impact only, no DIF) ---\n")
t0 <- proc.time()
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
cat("Elapsed:", round((proc.time() - t0)[["elapsed"]] / 60, 1), "min\n")

fit$save_object(file.path(out_dir, paste0("fit_static_", sx, ".rds")))

# ---------------------------------------------------------------------------
# DIAGNOSTICS
# ---------------------------------------------------------------------------
diag <- fit$diagnostic_summary(quiet = TRUE)
cat("\nDivergences per chain:", diag$num_divergent, "\n")
cat("Max treedepth hits:  ", diag$num_max_treedepth, "\n")
cat("E-BFMI:              ", round(diag$ebfmi, 3), "\n")

summ <- fit$summary()
write.csv(
  summ,
  file.path(out_dir, paste0("convergence_summary_static_", sx, ".csv")),
  row.names = FALSE
)

bad <- summ[!is.na(summ$rhat) & summ$rhat > 1.01, ]
if (nrow(bad) > 0) {
  cat("\nParameters with Rhat > 1.01 (top 10):\n")
  print(
    head(bad[order(-bad$rhat), c("variable", "mean", "rhat", "ess_bulk")], 10),
    digits = 3,
    row.names = FALSE
  )
} else {
  cat("\nAll Rhat <= 1.01\n")
}

key_summ <- fit$summary(variables = c("b_mu", "b_phi", "phi0"))
cat("\nImpact parameters:\n")
print(key_summ[, c("variable", "mean", "sd", "q5", "q95", "rhat", "ess_bulk")], digits = 3, row.names = FALSE)

item_summ <- fit$summary(variables = c("lp", "np"))
item_summ$item <- prep$item_names[
  as.integer(gsub(".*\\[|\\]", "", item_summ$variable))
]
write.csv(
  item_summ,
  file.path(out_dir, paste0("item_params_static_", sx, ".csv")),
  row.names = FALSE
)
cat("\nItem parameters:\n")
print(item_summ[, c("item", "variable", "mean", "rhat", "ess_bulk")], digits = 3, row.names = FALSE)

# ---------------------------------------------------------------------------
# PLOT: posterior mean latent factor vs. age
# ---------------------------------------------------------------------------
eta_summ <- fit$summary(variables = "eta")
eta_df <- tibble(
  person_idx = seq_len(prep$ni),
  id = prep$ids,
  age_c = prep$age_c_person,
  age = prep$age_c_person + prep$age_mean,
  eta = eta_summ$mean
)
write.csv(
  eta_df,
  file.path(out_dir, paste0("factor_scores_static_", sx, ".csv")),
  row.names = FALSE
)

p_eta <- ggplot(eta_df, aes(x = age, y = eta)) +
  geom_point(alpha = 0.15, size = 0.6, colour = "#2166ac") +
  geom_smooth(method = "loess", se = TRUE, colour = "black") +
  labs(
    title = paste0("Static MNLFA factor scores vs. age — ", sx),
    subtitle = "Baseline wave only; posterior means",
    x = "Age (years)",
    y = "Latent puberty (η)"
  ) +
  theme_minimal(base_size = 13)
ggsave(
  file.path(out_dir, paste0("factor_scores_static_", sx, ".png")),
  p_eta,
  width = 7,
  height = 5,
  dpi = 180
)

cat("\nAll outputs written to:", out_dir, "\n")
cat("Done:", sx, "\n")
