## lmnlfa_sigmoid_premeeting_diagnostics.R
## Pre-meeting diagnostics for the staged sigmoidal L-MNLFA of pubertal
## development (lmnlfa_growth_sigmoid_staged.R). This script does NOT run
## any new Stan sampling -- it only post-processes the cached fitA/fitB/fitC
## .rds objects (fits/fitA_impact_<sex>_n1500.rds, etc.) and computes
## descriptives directly from the raw parent/youth CSVs.
##
## It does NOT modify lmnlfa_growth_sigmoid_staged.R or any Stan file.
## build_lmnlfa_data_staged() and select_dif() below are copy-pasted
## VERBATIM from lmnlfa_growth_sigmoid_staged.R (not source()'d, since that
## script has top-level side effects -- CLI parsing, Stan compilation/
## sampling calls -- that can't safely run as a library import). If the
## main script's data-prep or DIF-selection logic ever changes, these two
## copies need to be updated by hand to match.
##
## Usage:
##   Rscript lmnlfa_sigmoid_premeeting_diagnostics.R <sex>
##   sex: female | male
##
## Outputs: written to <OUT_DIR>/lmnlfa_growth_sigmoid_staged/premeeting/

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(posterior)
})

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("cmdstanr not found.")
}
library(cmdstanr)

options(mc.cores = as.integer(Sys.getenv("NSLOTS", unset = "4")))

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
# ARGUMENTS
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript lmnlfa_sigmoid_premeeting_diagnostics.R <female|male>")
}
sx <- args[1]
if (!sx %in% c("female", "male")) {
  stop("sex must be 'female' or 'male'")
}
n_subsample <- 1500
run_tag <- paste0(sx, "_n", n_subsample)
run_label <- paste0(sx, " (n=", n_subsample, " subsample)")
cat("Sex:", sx, "\n")

# ---------------------------------------------------------------------------
# PATHS -- sshfs-mounted path preferred locally (HOME_DIR/HOME resolve to a
# stale Box mirror on this machine); DATA_DIR/OUT_DIR env vars (set on HPC)
# always take precedence when present.
# ---------------------------------------------------------------------------
sshfs_data_default <- "/private/tmp/sshfs/projects/abcd-projs/dissertation/lmnlfa-puberty/data"
sshfs_out_default <- "/private/tmp/sshfs/projects/abcd-projs/dissertation/lmnlfa-puberty/outputs"

data_dir <- Sys.getenv("DATA_DIR")
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  if (dir.exists(sshfs_data_default)) {
    data_dir <- sshfs_data_default
  } else {
    root_path <- Sys.getenv("HOME_DIR")
    if (!nzchar(root_path)) root_path <- Sys.getenv("HOME")
    data_dir <- file.path(root_path, "projects/abcd-projs/dissertation/lmnlfa-puberty/data")
  }
}
if (!dir.exists(data_dir)) stop("Cannot locate data directory: ", data_dir)

out_base <- Sys.getenv("OUT_DIR")
if (!nzchar(out_base) || !dir.exists(out_base)) {
  if (dir.exists(sshfs_out_default)) {
    out_base <- sshfs_out_default
  } else {
    root_path <- Sys.getenv("HOME_DIR")
    if (!nzchar(root_path)) root_path <- Sys.getenv("HOME")
    out_base <- file.path(root_path, "projects/abcd-projs/dissertation/lmnlfa-puberty/outputs")
  }
}

model_dir <- file.path(out_base, "lmnlfa_growth_sigmoid_staged")
fits_dir <- file.path(model_dir, "fits")
if (!dir.exists(fits_dir)) stop("Cannot locate fits directory: ", fits_dir)

out_dir <- file.path(model_dir, "premeeting")
audit_dir <- file.path(out_dir, "data-audit")
ceiling_dir <- file.path(out_dir, "ceiling-floor")
conv_dir <- file.path(out_dir, "convergence")
sens_dir <- file.path(out_dir, "dif-sensitivity")
impact_dir <- file.path(out_dir, "impact-comparison")
agree_dir <- file.path(out_dir, "agreement")
for (d in c(out_dir, audit_dir, ceiling_dir, conv_dir, sens_dir, impact_dir, agree_dir)) {
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
if (!file.exists(palette_file)) palette_file <- file.path("scripts", "color_palette.R")
source(palette_file)

cat("Data dir:  ", data_dir, "\n")
cat("Model dir: ", model_dir, "\n")
cat("Output dir:", out_dir, "\n")

written_files <- character(0)
track <- function(path) {
  written_files <<- c(written_files, path)
  path
}

# ---------------------------------------------------------------------------
# LOAD RAW DATA
# ---------------------------------------------------------------------------
parent_df <- read.csv(file.path(data_dir, paste0(sx, "_parent_long.csv")))
youth_df <- read.csv(file.path(data_dir, paste0(sx, "_youth_long.csv")))

wave_order <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")
ordinal_items <- c("peta", "petb", "petc", "petd")
whtr_min <- 0.25
whtr_max <- 0.75

# ---------------------------------------------------------------------------
# COPIED VERBATIM from lmnlfa_growth_sigmoid_staged.R -- see header note
# above. DO NOT edit the logic here without also checking the main script.
# ---------------------------------------------------------------------------
build_lmnlfa_data_staged <- function(
  parent_df,
  youth_df,
  sex_label,
  ordinal_items = c("peta", "petb", "petc", "petd"),
  n_subsample = NULL
) {
  cat(
    "\n=== Building longitudinal data (staged impact + DIF, sigmoid) |",
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

# same seed, same call order as the main script -- reproduces the identical
# n=1500 subsample used by the cached fits
set.seed(90025)
prep <- build_lmnlfa_data_staged(parent_df, youth_df, sx, n_subsample = n_subsample)

covar_labels <- c("Hispanic", "White", "Black", "WtHR (baseline)")
dif_covar_names <- c("age", "race_c1", "race_c2", "race_c3", "whtr", "informant")
covar_display <- setNames(
  c("Age", "Hispanic", "White", "Black", "WtHR", "Informant"),
  dif_covar_names
)
item_display <- c(
  peta = "Growth spurt",
  petb = "Body hair",
  petc = "Skin changes",
  petd = if (sx == "male") "Voice changes" else "Breast growth"
)[prep$item_names]

ximp <- prep$person_covars %>%
  select(race_c1, race_c2, race_c3, whtr_baseline_c) %>%
  as.matrix()
xdif <- prep$dat_long %>%
  select(age_c, race_c1, race_c2, race_c3, whtr_c, informant_c) %>%
  as.matrix()
kimp <- ncol(ximp)
kdif <- ncol(xdif)

# ---------------------------------------------------------------------------
# COPIED VERBATIM from lmnlfa_growth_sigmoid_staged.R
# ---------------------------------------------------------------------------
select_dif <- function(fitB, p, kdif, covar_names, fdr_q = 0.05, mag_floor = 0.05) {
  l_draws <- fitB$draws(variables = "l_dif", format = "matrix")
  n_draws <- fitB$draws(variables = "n_dif", format = "matrix")

  tail_prob <- function(draws, varname, p, kdif) {
    m <- matrix(NA_real_, p, kdif)
    for (i in 1:p) {
      for (k in 1:kdif) {
        col <- draws[, paste0(varname, "[", i, ",", k, "]")]
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

  l_adj <- matrix(p.adjust(as.vector(l_p), method = "BH"), nrow = p, ncol = kdif)
  l_pattern <- ((l_adj < fdr_q) & (abs(l_mean) > mag_floor)) * 1

  remaining <- l_pattern == 0
  n_adj_remaining <- p.adjust(n_p[remaining], method = "BH")
  n_pattern <- l_pattern
  n_pattern[remaining] <- ((n_adj_remaining < fdr_q) & (abs(n_mean[remaining]) > mag_floor)) * 1

  rownames(l_pattern) <- rownames(n_pattern) <- prep$item_names
  colnames(l_pattern) <- colnames(n_pattern) <- covar_names

  list(l_pattern = l_pattern, n_pattern = n_pattern, l_mean = l_mean, n_mean = n_mean)
}

# ---------------------------------------------------------------------------
# LOAD CACHED FITS (skip gracefully if missing)
# ---------------------------------------------------------------------------
load_fit <- function(letter, tag) {
  f <- file.path(fits_dir, paste0(tag, "_", run_tag, ".rds"))
  if (!file.exists(f)) {
    cat("Stage", letter, "fit not found at", f, "-- skipping.\n")
    return(NULL)
  }
  cat("Loading Stage", letter, "fit:", f, "\n")
  readRDS(f)
}
fitA <- load_fit("A", "fitA_impact")
fitB <- load_fit("B", "fitB_difscreen")
fitC <- load_fit("C", "fitC_final")

summary_lines <- character(0)
add_summary <- function(...) {
  summary_lines <<- c(summary_lines, paste0(...))
}

cat("\n############################################################\n")
cat("SECTION 1: DATA AUDIT\n")
cat("############################################################\n")

audit_reporter <- function(df, label) {
  steps <- list()
  d <- df
  steps[["1_raw_rows"]] <- d
  d <- d %>% filter(!is.na(age))
  steps[["2_drop_missing_age"]] <- d
  d <- d %>% filter(!is.na(race))
  steps[["3_drop_missing_race"]] <- d
  d <- d %>% filter(!is.na(whtr))
  steps[["4_drop_missing_whtr"]] <- d
  d <- d %>% filter(whtr >= whtr_min, whtr <= whtr_max)
  steps[["5_drop_whtr_out_of_range"]] <- d
  d <- d %>% filter(if_all(all_of(ordinal_items), ~ !is.na(.) & as.integer(.) %in% 1:4))
  steps[["6_drop_item_missing_or_oob"]] <- d
  d <- d %>%
    mutate(race_grp = case_when(
      race == 1 ~ "Hispanic", race == 2 ~ "White", race == 3 ~ "Black",
      race %in% c(7, 11, 12, 13) ~ "Other", TRUE ~ NA_character_
    )) %>%
    filter(!is.na(race_grp))
  steps[["7_drop_race_not_in_four_groups"]] <- d

  rows_remaining <- sapply(steps, nrow)
  persons_remaining <- sapply(steps, function(x) length(unique(x$id)))
  tibble(
    informant = label,
    step = names(steps),
    rows_remaining = rows_remaining,
    persons_remaining = persons_remaining
  ) %>%
    mutate(rows_dropped = lag(rows_remaining, default = first(rows_remaining)) - rows_remaining) %>%
    relocate(rows_dropped, .after = step)
}

waterfall_parent <- audit_reporter(parent_df, "parent")
waterfall_youth <- audit_reporter(youth_df, "youth")
waterfall <- bind_rows(waterfall_parent, waterfall_youth)
write.csv(waterfall, track(file.path(audit_dir, paste0("waterfall_", run_tag, ".csv"))), row.names = FALSE)
cat("Waterfall (parent + youth):\n")
print(waterfall)

# persons retained per wave, by informant (post full cleaning, step 7, not
# subsampled -- reflects the full data-availability picture)
persons_per_wave <- bind_rows(
  waterfall_parent[waterfall_parent$step == "7_drop_race_not_in_four_groups", ] %>% select(informant),
  .id = NULL
)
wave_counts <- function(df, label) {
  d <- df %>%
    filter(!is.na(age), !is.na(race), !is.na(whtr), whtr >= whtr_min, whtr <= whtr_max) %>%
    filter(if_all(all_of(ordinal_items), ~ !is.na(.) & as.integer(.) %in% 1:4)) %>%
    mutate(race_grp = case_when(
      race == 1 ~ "Hispanic", race == 2 ~ "White", race == 3 ~ "Black",
      race %in% c(7, 11, 12, 13) ~ "Other", TRUE ~ NA_character_
    )) %>%
    filter(!is.na(race_grp))
  d %>%
    mutate(wave = factor(wave, levels = wave_order)) %>%
    group_by(wave) %>%
    summarize(n_persons = n_distinct(id), .groups = "drop") %>%
    mutate(informant = label)
}
persons_per_wave <- bind_rows(wave_counts(parent_df, "parent"), wave_counts(youth_df, "youth"))
write.csv(persons_per_wave, track(file.path(audit_dir, paste0("persons_per_wave_", run_tag, ".csv"))), row.names = FALSE)
cat("\nPersons retained per wave, by informant:\n")
print(persons_per_wave)

# baseline-occasion wave origin, in the actual n=1500 ANALYSIS sample (prep)
baseline_wave_tab <- prep$dat_long %>%
  distinct(person_idx, time_idx) %>%
  group_by(person_idx) %>%
  summarize(min_time = min(time_idx), .groups = "drop") %>%
  mutate(wave = wave_order[min_time]) %>%
  count(wave, name = "n_persons") %>%
  mutate(is_baseline_wave = wave == "bl")
write.csv(baseline_wave_tab, track(file.path(audit_dir, paste0("baseline_wave_origin_", run_tag, ".csv"))), row.names = FALSE)
n_baseline_not_bl <- sum(baseline_wave_tab$n_persons[!baseline_wave_tab$is_baseline_wave])
pct_baseline_not_bl <- round(100 * n_baseline_not_bl / sum(baseline_wave_tab$n_persons), 1)
cat("\nBaseline-occasion wave origin (analysis sample):\n")
print(baseline_wave_tab)
cat("Persons whose 'baseline' WtHR comes from a wave other than bl:", n_baseline_not_bl,
  "(", pct_baseline_not_bl, "% )\n")

# race group counts in the retained analysis sample
race_counts <- prep$person_raw %>% count(race_grp, name = "n_persons")
write.csv(race_counts, track(file.path(audit_dir, paste0("race_group_counts_", run_tag, ".csv"))), row.names = FALSE)
cat("\nRace group counts (analysis sample):\n")
print(race_counts)

add_summary(
  "SECTION 1 -- DATA AUDIT. Parent data starts at ", waterfall_parent$rows_remaining[1],
  " rows and retains ", waterfall_parent$rows_remaining[nrow(waterfall_parent)],
  " (", round(100 * waterfall_parent$rows_remaining[nrow(waterfall_parent)] / waterfall_parent$rows_remaining[1], 1),
  "%) after all filters; youth data starts at ", waterfall_youth$rows_remaining[1], " rows and retains ",
  waterfall_youth$rows_remaining[nrow(waterfall_youth)], " (",
  round(100 * waterfall_youth$rows_remaining[nrow(waterfall_youth)] / waterfall_youth$rows_remaining[1], 1), "%)."
)
add_summary(
  "In the n=", n_subsample, " analysis sample, ", n_baseline_not_bl, " people (", pct_baseline_not_bl,
  "%) have their 'baseline' WtHR value drawn from a wave other than the true baseline (bl), because bl itself ",
  "was missing/filtered for those people -- baseline-WtHR-based impact estimates for them reflect a later occasion."
)
add_summary(
  "Race group counts in the analysis sample: ",
  paste(paste0(race_counts$race_grp, "=", race_counts$n_persons), collapse = ", "), "."
)

cat("\n############################################################\n")
cat("SECTION 2: CEILING / FLOOR SUPPORT\n")
cat("############################################################\n")

age_band_breaks <- c(-Inf, 10.9999, 12.9999, 14.9999, Inf)
age_band_labels <- c("9-10", "11-12", "13-14", "15+")

item_resp <- prep$dat_long %>%
  mutate(
    age = age_c * prep$age_sd + prep$age_mean,
    age_band = cut(age, breaks = age_band_breaks, labels = age_band_labels),
    informant = ifelse(informant_c == 1, "parent", "youth"),
    item = item_display[item]
  )
item_levels <- unname(item_display[prep$item_names])

# Facial hair (mpete) is decided to be added as a genuine 5th item for
# males on the next full model rerun (it's already ordinal 1-4, same as
# peta-petd, and other scripts in this project already treat it that way --
# see the fifth-item exploration note below). It isn't in prep$dat_long
# (prep matches the CACHED Stage A/B/C fits, which don't include it), so
# it's pulled from the raw parent/youth files here and matched to the exact
# same modeled person-occasions by id x wave x informant, purely so it can
# be shown on the same ceiling/floor scale as the other four ahead of that
# rerun -- this does NOT change prep, the cached fits, or any Stage A-C output.
if (sx == "male") {
  mpete_lookup <- bind_rows(
    parent_df %>% select(id, wave, mpete) %>% mutate(informant_c = 1),
    youth_df %>% select(id, wave, mpete) %>% mutate(informant_c = -1)
  ) %>%
    mutate(time_idx = match(wave, wave_order)) %>%
    filter(!is.na(mpete), as.integer(mpete) %in% 1:4)

  occasion_lookup <- prep$dat_long %>% distinct(id, person_idx, time_idx, informant_c, age_c)
  mpete_resp <- occasion_lookup %>%
    inner_join(mpete_lookup, by = c("id", "time_idx", "informant_c")) %>%
    mutate(
      age = age_c * prep$age_sd + prep$age_mean,
      age_band = cut(age, breaks = age_band_breaks, labels = age_band_labels),
      informant = ifelse(informant_c == 1, "parent", "youth"),
      item = "Facial hair",
      y_int = as.integer(mpete)
    ) %>%
    select(person_idx, time_idx, age_c, informant_c, item, age, age_band, informant, y_int)
  cat("Facial hair (mpete) matched to", nrow(mpete_resp), "of", nrow(occasion_lookup),
    "modeled person-occasions (", round(100 * nrow(mpete_resp) / nrow(occasion_lookup), 1),
    "% ) -- added to the ceiling/floor diagnostics below as a 5th item, ahead of being added to the full model.\n")
  item_resp <- bind_rows(item_resp, mpete_resp)
  item_levels <- c(item_levels, "Facial hair")
}
item_resp <- item_resp %>% mutate(item = factor(item, levels = item_levels))

ceiling_floor_tab <- item_resp %>%
  group_by(item, age_band, informant) %>%
  summarize(
    n = n(),
    prop_cat1 = mean(y_int == 1),
    prop_cat4 = mean(y_int == 4),
    .groups = "drop"
  )
write.csv(ceiling_floor_tab, track(file.path(ceiling_dir, paste0("ceiling_floor_by_item_", run_tag, ".csv"))), row.names = FALSE)
cat("Ceiling/floor support by item x age band x informant:\n")
print(ceiling_floor_tab, n = 100)

cf_plot_df <- ceiling_floor_tab %>%
  pivot_longer(c(prop_cat1, prop_cat4), names_to = "category", values_to = "proportion") %>%
  mutate(category = recode(category, prop_cat1 = "Category 1 (floor)", prop_cat4 = "Category 4 (ceiling)"))
p_cf <- ggplot(cf_plot_df, aes(x = age_band, y = proportion, fill = informant)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_grid(category ~ item) +
  scale_fill_manual(values = setNames(pal_two, c("parent", "youth"))) +
  labs(
    title = paste0("Floor/ceiling response proportions by item, age band, and informant - ", run_label),
    subtitle = "Proportion of item responses at the lowest (1) or highest (4) category",
    x = "Age band", y = "Proportion", fill = "Informant"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
ggsave(track(file.path(ceiling_dir, paste0("ceiling_floor_by_item_", run_tag, ".png"))), p_cf, width = if (sx == "male") 13 else 11, height = 7, dpi = 220)

# person-occasions where ALL 4 items are at ceiling / floor
occ_all <- prep$dat_long %>%
  mutate(
    age = age_c * prep$age_sd + prep$age_mean,
    age_band = cut(age, breaks = age_band_breaks, labels = age_band_labels),
    informant = ifelse(informant_c == 1, "parent", "youth")
  ) %>%
  group_by(person_idx, time_idx, age_band, informant) %>%
  summarize(all_ceiling = all(y_int == 4), all_floor = all(y_int == 1), .groups = "drop")

occ_all_tab <- occ_all %>%
  group_by(age_band, informant) %>%
  summarize(
    n_occasions = n(),
    prop_all_ceiling = mean(all_ceiling),
    prop_all_floor = mean(all_floor),
    .groups = "drop"
  )
write.csv(occ_all_tab, track(file.path(ceiling_dir, paste0("all_items_ceiling_floor_", run_tag, ".csv"))), row.names = FALSE)
cat("\nPerson-occasions with ALL items at ceiling/floor, by age band x informant:\n")
print(occ_all_tab)

# inflection age per person, from Stage C if available else Stage A
infl_source <- NULL
infl_ages_years <- NULL
if (!is.null(fitC)) {
  infl_source <- "Stage C"
  fac_gr_summ <- fitC$summary(variables = "fac_gr")
  infl_ages_c <- fac_gr_summ$mean[grepl("^fac_gr\\[2,", fac_gr_summ$variable)]
  infl_ages_years <- infl_ages_c * prep$age_sd + prep$age_mean
} else if (!is.null(fitA)) {
  infl_source <- "Stage A"
  fac_gr_summ <- fitA$summary(variables = "fac_gr")
  infl_ages_c <- fac_gr_summ$mean[grepl("^fac_gr\\[2,", fac_gr_summ$variable)]
  infl_ages_years <- infl_ages_c * prep$age_sd + prep$age_mean
}

if (!is.null(infl_ages_years)) {
  obs_age_range <- range(prep$dat_long$age_c * prep$age_sd + prep$age_mean)
  infl_df <- tibble(person_idx = seq_along(infl_ages_years), inflection_age = infl_ages_years)
  infl_summary_tab <- tibble(
    source = infl_source,
    n = length(infl_ages_years),
    mean = mean(infl_ages_years),
    median = median(infl_ages_years),
    sd = sd(infl_ages_years),
    q5 = quantile(infl_ages_years, 0.05),
    q95 = quantile(infl_ages_years, 0.95),
    min = min(infl_ages_years),
    max = max(infl_ages_years),
    obs_age_min = obs_age_range[1],
    obs_age_max = obs_age_range[2],
    share_above_obs_max = mean(infl_ages_years > obs_age_range[2]),
    share_below_obs_min = mean(infl_ages_years < obs_age_range[1])
  )
  write.csv(infl_summary_tab, track(file.path(ceiling_dir, paste0("inflection_age_summary_", run_tag, ".csv"))), row.names = FALSE)
  write.csv(infl_df, track(file.path(ceiling_dir, paste0("inflection_age_per_person_", run_tag, ".csv"))), row.names = FALSE)
  cat("\nInflection age distribution (", infl_source, "):\n", sep = "")
  print(infl_summary_tab)

  p_infl <- ggplot(infl_df, aes(x = inflection_age)) +
    geom_histogram(bins = 40, fill = pal_primary, colour = "white") +
    annotate("rect", xmin = obs_age_range[1], xmax = obs_age_range[2], ymin = -Inf, ymax = Inf,
      alpha = 0.12, fill = sunset_anchors["honey"]) +
    geom_vline(xintercept = obs_age_range, linetype = "dashed", colour = sunset_anchors["brick"]) +
    labs(
      title = paste0("Posterior-mean inflection age per person - ", run_label, " (", infl_source, ")"),
      subtitle = paste0(
        "Shaded band = observed age range [", round(obs_age_range[1], 1), ", ", round(obs_age_range[2], 1), "]; ",
        round(100 * infl_summary_tab$share_above_obs_max, 1), "% above, ",
        round(100 * infl_summary_tab$share_below_obs_min, 1), "% below"
      ),
      x = "Inflection age (years)", y = "Number of people"
    ) +
    theme_minimal(base_size = 13)
  ggsave(track(file.path(ceiling_dir, paste0("inflection_age_histogram_", run_tag, ".png"))), p_infl, width = 9, height = 6, dpi = 220)
} else {
  cat("\nNeither Stage A nor Stage C fit available -- skipping inflection-age diagnostics.\n")
}

# population mean curve at ximp = 0 -- predicted eta at youngest/oldest
# observed age (uses Stage C's population parameters if available, else A)
pop_curve_tab <- NULL
pop_fit <- if (!is.null(fitC)) fitC else fitA
pop_fit_label <- if (!is.null(fitC)) "Stage C" else if (!is.null(fitA)) "Stage A" else NA
if (!is.null(pop_fit)) {
  mu_logk_draws <- as.vector(pop_fit$draws(variables = "mu_logk", format = "matrix"))
  mu_alpha_draws <- as.vector(pop_fit$draws(variables = "mu_alpha", format = "matrix"))
  age_c_range <- range(prep$dat_long$age_c)
  eta_at <- function(a) {
    vals <- 1 + 4 / (1 + exp(-exp(mu_logk_draws) * (a - mu_alpha_draws)))
    c(mean = mean(vals), q5 = unname(quantile(vals, 0.05)), q95 = unname(quantile(vals, 0.95)))
  }
  young <- eta_at(age_c_range[1])
  old <- eta_at(age_c_range[2])
  pop_curve_tab <- tibble(
    source = pop_fit_label,
    age_years = c(age_c_range[1] * prep$age_sd + prep$age_mean, age_c_range[2] * prep$age_sd + prep$age_mean),
    which = c("youngest observed", "oldest observed"),
    eta_mean = c(young["mean"], old["mean"]),
    eta_q5 = c(young["q5"], old["q5"]),
    eta_q95 = c(young["q95"], old["q95"])
  )
  write.csv(pop_curve_tab, track(file.path(ceiling_dir, paste0("population_mean_curve_range_", run_tag, ".csv"))), row.names = FALSE)
  cat("\nPopulation-mean curve at ximp = 0, at observed age extremes (", pop_fit_label, "):\n", sep = "")
  print(pop_curve_tab)
}

if (!is.null(pop_curve_tab)) {
  range_traversed <- pop_curve_tab$eta_mean[2] - pop_curve_tab$eta_mean[1]
  add_summary(
    "SECTION 2 -- CEILING/FLOOR SUPPORT. The population-mean curve (", pop_fit_label, ") moves from eta = ",
    round(pop_curve_tab$eta_mean[1], 2), " at the youngest observed age (", round(pop_curve_tab$age_years[1], 1),
    " yrs) to eta = ", round(pop_curve_tab$eta_mean[2], 2), " at the oldest observed age (",
    round(pop_curve_tab$age_years[2], 1), " yrs) -- traversing ", round(range_traversed, 2),
    " of the full 1-5 range (", round(100 * range_traversed / 4, 1), "%)."
  )
}
if (!is.null(infl_ages_years)) {
  add_summary(
    "Individual inflection ages (", infl_source, ") have a median of ", round(median(infl_ages_years), 2),
    " years (range ", round(min(infl_ages_years), 1), "-", round(max(infl_ages_years), 1), "); ",
    round(100 * infl_summary_tab$share_above_obs_max, 1), "% of people are predicted to reach their midpoint ",
    "AFTER the oldest age observed in the sample, and ", round(100 * infl_summary_tab$share_below_obs_min, 1),
    "% BEFORE the youngest age observed -- both indicate extrapolation beyond directly observed data for those people."
  )
}
max_ceiling <- occ_all_tab$prop_all_ceiling[which.max(occ_all_tab$prop_all_ceiling)]
max_ceiling_band <- occ_all_tab$age_band[which.max(occ_all_tab$prop_all_ceiling)]
add_summary(
  "The highest observed rate of ALL FOUR items simultaneously at the top category is ", round(100 * max_ceiling, 1),
  "%, in the ", max_ceiling_band, " / ", occ_all_tab$informant[which.max(occ_all_tab$prop_all_ceiling)],
  " cell -- i.e. even in the oldest/most-developed cell, a hard ceiling in the RAW items is not the norm, ",
  "consistent with (but not proof of) the model's fixed ceiling being reasonable rather than truncating real variation."
)

cat("\n############################################################\n")
cat("SECTION 3: CONVERGENCE (DIF-vs-GROWTH CONFOUND)\n")
cat("############################################################\n")

relabel_var <- function(vars, item_names, dif_names, imp_names) {
  vapply(vars, function(v) {
    m <- regmatches(v, regexec("^([a-zA-Z_]+)\\[(\\d+)(?:,(\\d+))?\\]$", v))[[1]]
    if (length(m) == 0) return(v)
    base <- m[2]
    i1 <- as.integer(m[3])
    i2 <- if (is.na(m[4]) || m[4] == "") NA_integer_ else as.integer(m[4])
    if (base %in% c("b_mu_logk", "b_mu_alpha")) return(paste0(base, "_", imp_names[i1]))
    if (base %in% c("lp", "np")) return(paste0(base, "_", item_names[i1]))
    if (base %in% c("l_dif", "n_dif")) return(paste0(base, "_", item_names[i1], "_x_", dif_names[i2]))
    v
  }, character(1), USE.NAMES = FALSE)
}

conv_table <- function(fit, stage_name, has_impact = FALSE, has_dif = FALSE) {
  vars <- c("mu_logk", "mu_alpha", "phi_logk", "phi_alpha", "rho", "eti_sd")
  if (has_impact) vars <- c(vars, "b_mu_logk", "b_mu_alpha")
  vars <- c(vars, "lp", "np")
  if (has_dif) vars <- c(vars, "l_dif", "n_dif")
  summ <- fit$summary(variables = vars)
  summ$label <- relabel_var(summ$variable, prep$item_names, dif_covar_names, covar_labels)
  summ$stage <- stage_name
  summ
}

conv_tabs <- list()
if (!is.null(fitA)) conv_tabs[["A"]] <- conv_table(fitA, "A: Impact only", has_impact = TRUE, has_dif = FALSE)
if (!is.null(fitB)) conv_tabs[["B"]] <- conv_table(fitB, "B: DIF screening", has_impact = FALSE, has_dif = TRUE)
if (!is.null(fitC)) conv_tabs[["C"]] <- conv_table(fitC, "C: Final", has_impact = TRUE, has_dif = TRUE)

conv_all <- bind_rows(conv_tabs) %>%
  select(stage, variable, label, mean, sd, q5, q95, rhat, ess_bulk, ess_tail)
write.csv(conv_all, track(file.path(conv_dir, paste0("convergence_all_params_", run_tag, ".csv"))), row.names = FALSE)
cat("Convergence table written (", nrow(conv_all), " rows across ", length(conv_tabs), " stage(s) ).\n", sep = "")

diag_tab <- purrr::map_dfr(names(conv_tabs), function(nm) {
  fit <- list(A = fitA, B = fitB, C = fitC)[[nm]]
  ds <- fit$diagnostic_summary(quiet = TRUE)
  tibble(stage = nm, n_divergent = sum(ds$num_divergent), n_max_treedepth = sum(ds$num_max_treedepth))
})
write.csv(diag_tab, track(file.path(conv_dir, paste0("divergences_treedepth_", run_tag, ".csv"))), row.names = FALSE)
cat("\nDivergences / max-treedepth by stage:\n")
print(diag_tab)

flagged <- conv_all %>% filter(rhat > 1.01 | ess_bulk < 400)
write.csv(flagged, track(file.path(conv_dir, paste0("flagged_params_", run_tag, ".csv"))), row.names = FALSE)
cat("\nFlagged parameters (rhat > 1.01 or ess_bulk < 400):", nrow(flagged), "of", nrow(conv_all), "\n")
if (nrow(flagged) > 0) print(flagged)

growth_param_prefixes <- c("mu_logk", "mu_alpha", "phi_logk", "phi_alpha", "rho", "eti_sd", "b_mu_logk_", "b_mu_alpha_")
flagged_growth <- flagged %>% filter(variable %in% c("mu_logk", "mu_alpha", "phi_logk", "phi_alpha", "rho", "eti_sd") |
  grepl("^b_mu_logk\\[|^b_mu_alpha\\[", variable))
flagged_age_dif <- flagged %>% filter(grepl("_x_age$", label))
cat("Flagged growth-curve parameters:", nrow(flagged_growth), "| Flagged age-DIF terms:", nrow(flagged_age_dif), "\n")

# Stage B: posterior correlations between mu_alpha/mu_logk and age l_dif/n_dif
dif_growth_cor <- NULL
if (!is.null(fitB)) {
  draws_mat <- fitB$draws(variables = c("mu_logk", "mu_alpha", "l_dif", "n_dif"), format = "matrix")
  age_idx <- which(dif_covar_names == "age")
  age_cols <- c(
    paste0("l_dif[", seq_len(prep$p), ",", age_idx, "]"),
    paste0("n_dif[", seq_len(prep$p), ",", age_idx, "]")
  )
  sub_mat <- draws_mat[, c("mu_logk", "mu_alpha", age_cols)]
  colnames(sub_mat) <- c(
    "mu_logk", "mu_alpha",
    paste0("l_dif_age_", item_display[prep$item_names]),
    paste0("n_dif_age_", item_display[prep$item_names])
  )
  cor_mat <- cor(sub_mat)
  cor_df <- as.data.frame(cor_mat) %>% rownames_to_column("param")
  write.csv(cor_df, track(file.path(conv_dir, paste0("stageB_growth_age_dif_correlations_", run_tag, ".csv"))), row.names = FALSE)

  cor_long <- as.data.frame(cor_mat) %>%
    rownames_to_column("param1") %>%
    pivot_longer(-param1, names_to = "param2", values_to = "correlation")
  p_cor <- ggplot(cor_long, aes(x = param2, y = param1, fill = correlation)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
    sunset_diverging(limits = c(-1, 1)) +
    labs(
      title = paste0("Stage B posterior correlations: growth curve vs age DIF - ", run_label),
      subtitle = "High magnitude correlation = direct evidence of the age-DIF vs growth-curve tradeoff",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(track(file.path(conv_dir, paste0("stageB_growth_age_dif_cor_heatmap_", run_tag, ".png"))), p_cor, width = 9, height = 7.5, dpi = 220)

  age_dif_labels <- setdiff(colnames(cor_mat), c("mu_logk", "mu_alpha"))
  dif_growth_cor <- cor_mat[c("mu_logk", "mu_alpha"), age_dif_labels, drop = FALSE]
  max_abs_cor <- max(abs(dif_growth_cor))
  cat("\nMax magnitude correlation between growth params and age-DIF terms:", round(max_abs_cor, 3), "\n")
}

add_summary(
  "SECTION 3 -- CONVERGENCE. Across all available stages (", paste(names(conv_tabs), collapse = ", "),
  "), ", nrow(flagged), " of ", nrow(conv_all), " reported parameters are flagged (rhat > 1.01 or ess_bulk < 400); ",
  nrow(flagged_growth), " of those are growth-curve parameters and ", nrow(flagged_age_dif), " are age-DIF terms. ",
  "Total divergences across stages: ", sum(diag_tab$n_divergent), "; max-treedepth hits: ", sum(diag_tab$n_max_treedepth), "."
)
if (!is.null(dif_growth_cor)) {
  add_summary(
    "The largest-magnitude posterior correlation between a growth-curve mean parameter (log-rate or inflection age) ",
    "and an age loading/intercept DIF term in Stage B is ", round(max_abs_cor, 3),
    " -- ", if (max_abs_cor > 0.4) "large enough to be a real identifiability concern for separating age-related growth from age-related measurement bias." else "modest, suggesting age-DIF and the growth curve's mean structure are reasonably separable in this sample.", ""
  )
}

cat("\n############################################################\n")
cat("SECTION 4: DIF SCREENING SENSITIVITY\n")
cat("############################################################\n")

sens_summary <- NULL
sens_cells <- NULL
raw_dif_ci <- NULL
if (!is.null(fitB)) {
  fdr_grid <- c(0.01, 0.05, 0.10)
  mag_grid <- c(0, 0.05, 0.10, 0.20)
  grid_results <- list()
  cell_results <- list()
  for (q in fdr_grid) {
    for (m in mag_grid) {
      sel <- select_dif(fitB, prep$p, kdif, dif_covar_names, fdr_q = q, mag_floor = m)
      grid_results[[length(grid_results) + 1]] <- tibble(
        fdr_q = q, mag_floor = m,
        n_loading_retained = sum(sel$l_pattern),
        n_intercept_retained = sum(sel$n_pattern)
      )
      cell_results[[length(cell_results) + 1]] <- as.data.frame(sel$l_pattern) %>%
        rownames_to_column("item") %>%
        pivot_longer(-item, names_to = "covariate", values_to = "loading_retained") %>%
        left_join(
          as.data.frame(sel$n_pattern) %>%
            rownames_to_column("item") %>%
            pivot_longer(-item, names_to = "covariate", values_to = "intercept_retained"),
          by = c("item", "covariate")
        ) %>%
        mutate(fdr_q = q, mag_floor = m)
    }
  }
  sens_summary <- bind_rows(grid_results)
  sens_cells <- bind_rows(cell_results)
  write.csv(sens_summary, track(file.path(sens_dir, paste0("dif_sensitivity_grid_summary_", run_tag, ".csv"))), row.names = FALSE)
  write.csv(sens_cells, track(file.path(sens_dir, paste0("dif_sensitivity_grid_cells_", run_tag, ".csv"))), row.names = FALSE)
  cat("DIF selection sensitivity grid:\n")
  print(sens_summary)

  sens_plot_df <- sens_cells %>%
    mutate(
      status = case_when(
        loading_retained == 1 & intercept_retained == 1 ~ "Both",
        loading_retained == 1 ~ "Loading only",
        intercept_retained == 1 ~ "Intercept only",
        TRUE ~ "Not retained"
      ),
      item_label = factor(item_display[item], levels = rev(unname(item_display))),
      covar_label = factor(covar_display[covariate], levels = unname(covar_display)),
      fdr_label = paste0("FDR q=", fdr_q),
      mag_label = paste0("mag floor=", mag_floor)
    )
  p_sens <- ggplot(sens_plot_df, aes(x = covar_label, y = item_label, fill = status)) +
    geom_tile(colour = "white") +
    facet_grid(fdr_label ~ mag_label) +
    scale_fill_manual(values = c(
      "Both" = unname(sunset_anchors["bordeaux"]),
      "Loading only" = unname(sunset_anchors["slate"]),
      "Intercept only" = unname(sunset_anchors["honey"]),
      "Not retained" = unname(sunset_anchors["cream"])
    )) +
    labs(
      title = paste0("DIF selection sensitivity to threshold choices - ", run_label),
      subtitle = "Stage B posterior, re-screened over a grid of BH-FDR q and magnitude-floor thresholds",
      x = NULL, y = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.text.y = element_text(size = 7),
      legend.position = "bottom",
      strip.text = element_text(size = 8)
    )
  ggsave(track(file.path(sens_dir, paste0("dif_sensitivity_heatmap_", run_tag, ".png"))), p_sens, width = 13, height = 9, dpi = 220)

  l_summ <- fitB$summary(variables = "l_dif")
  n_summ <- fitB$summary(variables = "n_dif")
  l_summ$label <- relabel_var(l_summ$variable, prep$item_names, dif_covar_names, covar_labels)
  n_summ$label <- relabel_var(n_summ$variable, prep$item_names, dif_covar_names, covar_labels)
  raw_dif_ci <- bind_rows(
    l_summ %>% mutate(parameter = "loading") %>% select(parameter, label, mean, q5, q95),
    n_summ %>% mutate(parameter = "intercept") %>% select(parameter, label, mean, q5, q95)
  )
  write.csv(raw_dif_ci, track(file.path(sens_dir, paste0("stageB_raw_dif_posteriors_", run_tag, ".csv"))), row.names = FALSE)
  cat("\nRaw Stage B DIF posterior means + 90% CIs written:", nrow(raw_dif_ci), "rows\n")
}

if (!is.null(sens_summary)) {
  default_row <- sens_summary %>% filter(fdr_q == 0.05, mag_floor == 0.05)
  loading_range <- range(sens_summary$n_loading_retained)
  add_summary(
    "SECTION 4 -- DIF SCREENING SENSITIVITY. Across the 12-cell grid (FDR q in .01/.05/.10, magnitude floor in 0/.05/.10/.20), ",
    "the number of retained loading-DIF terms ranges from ", loading_range[1], " to ", loading_range[2],
    " (default q=.05/floor=.05 setting retains ", default_row$n_loading_retained[1], " loading and ",
    default_row$n_intercept_retained[1], " intercept terms), out of a maximum possible ", prep$p * kdif, " cells each."
  )
}

cat("\n############################################################\n")
cat("SECTION 5: STAGE A vs STAGE C IMPACT\n")
cat("############################################################\n")

impact_cmp_covar <- NULL
impact_cmp_scalar <- NULL
if (!is.null(fitA) && !is.null(fitC)) {
  a_logk_d <- fitA$draws(variables = "b_mu_logk", format = "matrix")
  a_alpha_d <- fitA$draws(variables = "b_mu_alpha", format = "matrix")
  c_logk_d <- fitC$draws(variables = "b_mu_logk", format = "matrix")
  c_alpha_d <- fitC$draws(variables = "b_mu_alpha", format = "matrix")

  build_cmp <- function(a_summ, c_summ, a_draws, param_name) {
    a_sd <- apply(a_draws, 2, sd)
    tibble(
      growth_param = param_name,
      covariate = covar_labels,
      A_mean = a_summ$mean, A_q5 = a_summ$q5, A_q95 = a_summ$q95,
      C_mean = c_summ$mean, C_q5 = c_summ$q5, C_q95 = c_summ$q95,
      diff_mean = c_summ$mean - a_summ$mean,
      diff_as_frac_A_sd = (c_summ$mean - a_summ$mean) / a_sd,
      A_excludes_zero = (a_summ$q5 > 0) | (a_summ$q95 < 0),
      C_excludes_zero = (c_summ$q5 > 0) | (c_summ$q95 < 0)
    ) %>%
      mutate(ci_status_changed = A_excludes_zero != C_excludes_zero)
  }
  impact_cmp_covar <- bind_rows(
    build_cmp(fitA$summary(variables = "b_mu_logk"), fitC$summary(variables = "b_mu_logk"), a_logk_d, "Log-rate (tempo)"),
    build_cmp(fitA$summary(variables = "b_mu_alpha"), fitC$summary(variables = "b_mu_alpha"), a_alpha_d, "Inflection age (timing)")
  )
  write.csv(impact_cmp_covar, track(file.path(impact_dir, paste0("stageA_vs_C_impact_covariates_", run_tag, ".csv"))), row.names = FALSE)
  cat("Stage A vs C impact (by covariate):\n")
  print(impact_cmp_covar)

  scalar_vars <- c("mu_logk", "mu_alpha", "phi_logk", "phi_alpha", "rho", "eti_sd")
  a_scal <- fitA$summary(variables = scalar_vars)
  c_scal <- fitC$summary(variables = scalar_vars)
  a_scal_draws <- fitA$draws(variables = scalar_vars, format = "matrix")
  a_scal_sd <- apply(a_scal_draws, 2, sd)
  impact_cmp_scalar <- tibble(
    parameter = scalar_vars,
    A_mean = a_scal$mean, A_q5 = a_scal$q5, A_q95 = a_scal$q95,
    C_mean = c_scal$mean, C_q5 = c_scal$q5, C_q95 = c_scal$q95,
    diff_mean = c_scal$mean - a_scal$mean,
    diff_as_frac_A_sd = (c_scal$mean - a_scal$mean) / a_scal_sd[scalar_vars],
    A_excludes_zero = (a_scal$q5 > 0) | (a_scal$q95 < 0),
    C_excludes_zero = (c_scal$q5 > 0) | (c_scal$q95 < 0)
  ) %>%
    mutate(ci_status_changed = A_excludes_zero != C_excludes_zero)
  write.csv(impact_cmp_scalar, track(file.path(impact_dir, paste0("stageA_vs_C_impact_scalars_", run_tag, ".csv"))), row.names = FALSE)
  cat("\nStage A vs C impact (scalar growth params):\n")
  print(impact_cmp_scalar)

  p_impact_cmp <- ggplot(
    impact_cmp_covar %>% mutate(covariate = factor(covariate, levels = rev(covar_labels))),
    aes(y = covariate)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_pointrange(aes(x = A_mean, xmin = A_q5, xmax = A_q95, colour = "A: Impact only"),
      position = position_nudge(y = 0.15), size = 0.5) +
    geom_pointrange(aes(x = C_mean, xmin = C_q5, xmax = C_q95, colour = "C: Impact + DIF"),
      position = position_nudge(y = -0.15), size = 0.5) +
    facet_wrap(~growth_param) +
    scale_colour_manual(values = c("A: Impact only" = pal_two[2], "C: Impact + DIF" = pal_two[1])) +
    labs(
      title = paste0("Stage A vs Stage C impact estimates - ", run_label),
      x = "Effect on growth-curve parameter", y = NULL, colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
  ggsave(track(file.path(impact_dir, paste0("stageA_vs_C_impact_plot_", run_tag, ".png"))), p_impact_cmp, width = 10, height = 6, dpi = 220)
}

if (!is.null(impact_cmp_covar)) {
  max_frac <- impact_cmp_covar$diff_as_frac_A_sd[which.max(abs(impact_cmp_covar$diff_as_frac_A_sd))]
  n_changed <- sum(impact_cmp_covar$ci_status_changed) + sum(impact_cmp_scalar$ci_status_changed)
  add_summary(
    "SECTION 5 -- STAGE A vs C IMPACT. The largest Stage A-to-C shift in an impact estimate, expressed as a fraction ",
    "of Stage A's own posterior SD, is ", round(max_frac, 2), " SDs (",
    impact_cmp_covar$covariate[which.max(abs(impact_cmp_covar$diff_as_frac_A_sd))], " / ",
    impact_cmp_covar$growth_param[which.max(abs(impact_cmp_covar$diff_as_frac_A_sd))], "); ",
    n_changed, " of ", nrow(impact_cmp_covar) + nrow(impact_cmp_scalar),
    " parameters changed whether their 90% CI excludes zero between stages."
  )
} else {
  add_summary("SECTION 5 -- STAGE A vs C IMPACT. Not computed (Stage A and/or Stage C fit unavailable).")
}

cat("\n############################################################\n")
cat("SECTION 6: PARENT-YOUTH AGREEMENT\n")
cat("############################################################\n")

clean_full <- function(df) {
  df %>%
    select(id, wave, age, race, whtr, all_of(ordinal_items)) %>%
    filter(!is.na(age), !is.na(race), !is.na(whtr), whtr >= whtr_min, whtr <= whtr_max) %>%
    filter(if_all(all_of(ordinal_items), ~ !is.na(.) & as.integer(.) %in% 1:4)) %>%
    mutate(race_grp = case_when(
      race == 1 ~ "Hispanic", race == 2 ~ "White", race == 3 ~ "Black",
      race %in% c(7, 11, 12, 13) ~ "Other", TRUE ~ NA_character_
    )) %>%
    filter(!is.na(race_grp)) %>%
    mutate(wave = factor(wave, levels = wave_order))
}
parent_clean <- clean_full(parent_df)
youth_clean <- clean_full(youth_df)

paired <- inner_join(parent_clean, youth_clean, by = c("id", "wave"), suffix = c("_p", "_y")) %>%
  mutate(age_band = cut(age_p, breaks = age_band_breaks, labels = age_band_labels))
cat("Paired parent-youth occasions:", nrow(paired), "\n")

agree_tab <- purrr::map_dfr(ordinal_items, function(it) {
  col_p <- paste0(it, "_p")
  col_y <- paste0(it, "_y")
  paired %>%
    group_by(age_band) %>%
    summarize(
      n = n(),
      exact_agreement = mean(.data[[col_p]] == .data[[col_y]]),
      within_one = mean(abs(.data[[col_p]] - .data[[col_y]]) <= 1),
      spearman_r = suppressWarnings(cor(.data[[col_p]], .data[[col_y]], method = "spearman")),
      polychoric_r = tryCatch(
        suppressWarnings(polycor::polychor(.data[[col_p]], .data[[col_y]])),
        error = function(e) NA_real_
      ),
      mean_signed_diff = mean(.data[[col_p]] - .data[[col_y]]),
      .groups = "drop"
    ) %>%
    mutate(item = item_display[it], .before = 1)
})
write.csv(agree_tab, track(file.path(agree_dir, paste0("parent_youth_agreement_", run_tag, ".csv"))), row.names = FALSE)
cat("Parent-youth agreement by item x age band:\n")
print(agree_tab, n = 50)

p_agree_diff <- ggplot(agree_tab, aes(x = age_band, y = mean_signed_diff, colour = item, group = item)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = setNames(pal_chains, unname(item_display))) +
  labs(
    title = paste0("Mean signed parent-youth difference by age band - ", run_label),
    subtitle = "Parent minus youth response category; positive = parent rates higher",
    x = "Age band", y = "Mean signed difference (parent - youth)", colour = "Item"
  ) +
  theme_minimal(base_size = 13)
ggsave(track(file.path(agree_dir, paste0("parent_youth_mean_diff_", run_tag, ".png"))), p_agree_diff, width = 9, height = 6, dpi = 220)

overall_exact <- weighted.mean(agree_tab$exact_agreement, agree_tab$n)
overall_within1 <- weighted.mean(agree_tab$within_one, agree_tab$n)
overall_spear <- mean(agree_tab$spearman_r, na.rm = TRUE)
add_summary(
  "SECTION 6 -- PARENT-YOUTH AGREEMENT. Across ", nrow(paired), " paired parent-youth occasions, average exact ",
  "item agreement is ", round(100 * overall_exact, 1), "% and within-one-category agreement is ",
  round(100 * overall_within1, 1), "%; mean Spearman correlation across items/age bands is ", round(overall_spear, 2),
  ". Mean signed differences (parent minus youth) are plotted by age band in parent_youth_mean_diff_", run_tag, ".png; ",
  "a value consistently near zero across bands is consistent with parent and youth reports tracking the same latent quantity, ",
  "while a systematic trend across age bands would suggest reporter divergence changes over development."
)

cat("\n############################################################\n")
cat("SECTION 7: MENARCHE (fpete) -- CONVERGENCE DIAGNOSIS (female only)\n")
cat("Facial hair (mpete, male) is decided -- it's now folded into Section 2's\n")
cat("ceiling/floor diagnostics as a genuine 5th item instead of a separate\n")
cat("exploration. This section is menarche-only: menarche previously caused\n")
cat("Stage A/B/C convergence problems when included, so these descriptives\n")
cat("exist to diagnose why, not to ask whether to include it. Descriptive\n")
cat("only -- no model changes made here.\n")
cat("############################################################\n")

if (sx == "female") {
item5_dir <- file.path(out_dir, "item5-exploration")
dir.create(item5_dir, showWarnings = FALSE, recursive = TRUE)

item5_name <- "fpete"
item5_label <- "Menarche (fpete)"
item5_type <- "binary"

# valid-value set determined empirically from the raw data, not assumed --
# fpete's on-disk coding doesn't match 00_data_foundation.R's own comment
# ("shift to 1/2"); observed values are {2,3}, so whatever the intended
# labels, we treat the two codes as "not yet" (lower) vs "reached" (higher)
# and confirm that direction empirically below (proportion at the higher
# code should rise with age).
item5_valid_vals <- sort(unique(na.omit(c(parent_df[[item5_name]], youth_df[[item5_name]]))))
cat("Item 5 = '", item5_name, "' (", item5_type, "), observed valid values: ",
  paste(item5_valid_vals, collapse = ", "), "\n", sep = "")

# completeness of item 5 vs. the four items already in the model, among
# rows that already pass the age/race/whtr filters (i.e. before requiring
# any PDS item to be non-missing) -- an apples-to-apples completeness check
completeness_calc <- function(df, reporter_label) {
  d <- df %>% filter(!is.na(age), !is.na(race), !is.na(whtr), whtr >= whtr_min, whtr <= whtr_max)
  n <- nrow(d)
  pct <- sapply(c(ordinal_items, item5_name), function(it) {
    if (it == item5_name) {
      round(100 * mean(!is.na(d[[it]]) & d[[it]] %in% item5_valid_vals), 1)
    } else {
      round(100 * mean(!is.na(d[[it]]) & as.integer(d[[it]]) %in% 1:4), 1)
    }
  })
  tibble(reporter = reporter_label, item = c(ordinal_items, item5_name), pct_valid = pct, n_base = n)
}
completeness_tab <- bind_rows(completeness_calc(parent_df, "parent"), completeness_calc(youth_df, "youth"))
write.csv(completeness_tab, track(file.path(item5_dir, paste0("item5_completeness_vs_other_items_", run_tag, ".csv"))), row.names = FALSE)
cat("\nCompleteness (% valid, among age/race/whtr-filtered rows), item 5 vs the four modeled items:\n")
print(completeness_tab)

# cost of REQUIRING item 5 too: among person-occasions that already pass
# the current 4-item filter (i.e. the universe build_lmnlfa_data_staged()
# samples from), what fraction also have item 5 valid?
clean_with_item5 <- function(df, reporter_label) {
  df %>%
    select(id, wave, age, race, whtr, all_of(ordinal_items), all_of(item5_name)) %>%
    filter(!is.na(age), !is.na(race), !is.na(whtr), whtr >= whtr_min, whtr <= whtr_max) %>%
    filter(if_all(all_of(ordinal_items), ~ !is.na(.) & as.integer(.) %in% 1:4)) %>%
    mutate(race_grp = case_when(
      race == 1 ~ "Hispanic", race == 2 ~ "White", race == 3 ~ "Black",
      race %in% c(7, 11, 12, 13) ~ "Other", TRUE ~ NA_character_
    )) %>%
    filter(!is.na(race_grp)) %>%
    mutate(
      wave = factor(wave, levels = wave_order),
      reporter = reporter_label,
      age_band = cut(age, breaks = age_band_breaks, labels = age_band_labels),
      item5_valid = !is.na(.data[[item5_name]]) & .data[[item5_name]] %in% item5_valid_vals
    )
}
item5_all <- bind_rows(clean_with_item5(parent_df, "parent"), clean_with_item5(youth_df, "youth"))

item5_missing_tab <- item5_all %>%
  group_by(reporter) %>%
  summarize(
    n_4item_retained = n(),
    n_item5_also_valid = sum(item5_valid),
    pct_item5_valid = round(100 * mean(item5_valid), 1),
    .groups = "drop"
  )
write.csv(item5_missing_tab, track(file.path(item5_dir, paste0("item5_missingness_cost_", run_tag, ".csv"))), row.names = FALSE)
cat("\nItem 5 validity among person-occasions that already pass the CURRENT 4-item filter (i.e. cost of also requiring item 5):\n")
print(item5_missing_tab)

# response rate by age band x informant -- confirms coding direction
# (higher code should become more common with age) and shows whether
# menarche behaves like a normal PDS item or transitions far more abruptly
# than the other four (either could independently explain a convergence
# failure: a positive-constrained loading (lp >= 0) fighting a NEGATIVE raw
# relationship, or a near-step-function item whose implied loading/threshold
# combination sits far outside what lp ~ normal(1, sigma_l) and tau ~
# normal(0, 1.5) were calibrated for)
#
# determine which of the two codes empirically means "reached/developed" by
# checking which one becomes MORE common with age, rather than assuming the
# higher numeric code does -- 00_data_foundation.R's own comment says fpete
# should be 1=no/2=yes after its +1 recode, but the values actually on disk
# are {2,3}, so the numeric-code-to-meaning mapping can't be trusted without
# checking
rate_by_max <- item5_all %>%
  filter(item5_valid, reporter == "parent") %>%
  group_by(age_band) %>%
  summarize(prop = mean(.data[[item5_name]] == max(item5_valid_vals)), .groups = "drop") %>%
  arrange(age_band)
d <- diff(rate_by_max$prop)
direction <- if (all(d >= -0.01)) {
  "increasing"
} else if (all(d <= 0.01)) {
  "decreasing"
} else {
  "non-monotonic"
}
developed_code <- if (direction == "decreasing") min(item5_valid_vals) else max(item5_valid_vals)
coding_matches_assumption <- direction == "increasing"

item5_rate_tab <- item5_all %>%
  filter(item5_valid) %>%
  group_by(age_band, reporter) %>%
  summarize(n = n(), prop_developed = mean(.data[[item5_name]] == developed_code), .groups = "drop")
write.csv(item5_rate_tab, track(file.path(item5_dir, paste0("item5_rate_by_age_", run_tag, ".csv"))), row.names = FALSE)
cat("\n", item5_label, " -- coding direction check: raw trend is '", direction,
  "' as age increases, so code ", developed_code, " (of observed values ",
  paste(item5_valid_vals, collapse = "/"), ") is being treated as 'reached/developed' below.\n", sep = "")
cat(item5_label, " -- proportion at the 'developed' code, by age band x informant:\n", sep = "")
print(item5_rate_tab)

p_item5 <- ggplot(item5_rate_tab, aes(x = age_band, y = prop_developed, colour = reporter, group = reporter)) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = setNames(pal_two, c("parent", "youth"))) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = paste0(item5_label, " -- proportion reached, by age band - ", run_label),
    subtitle = paste0(
      "Being diagnosed for Stage A/B/C convergence failure; code ", developed_code, " treated as 'reached' (raw trend: ", direction, " with age)"
    ),
    x = "Age band", y = "Proportion reached", colour = "Informant"
  ) +
  theme_minimal(base_size = 13)
ggsave(track(file.path(item5_dir, paste0("item5_rate_by_age_", run_tag, ".png"))), p_item5, width = 8, height = 6, dpi = 220)

# steepness check: menarche's biggest single-age-band jump in "developed"
# proportion, vs. the same statistic for the four modeled items' ceiling
# (category-4) proportion -- quantifies whether menarche transitions far
# more abruptly than anything the model currently has to represent
item5_max_jump <- max(abs(diff(
  item5_rate_tab %>% filter(reporter == "parent") %>% arrange(age_band) %>% pull(prop_developed)
)))
modeled_items_max_jump <- item_resp %>%
  filter(informant == "parent", item != "Facial hair") %>%
  group_by(item, age_band) %>%
  summarize(prop_cat4 = mean(y_int == 4), .groups = "drop") %>%
  arrange(item, age_band) %>%
  group_by(item) %>%
  summarize(max_jump = max(abs(diff(prop_cat4))), .groups = "drop")
cat("\nSteepness check -- largest single-age-band-transition jump in 'developed' proportion:\n")
cat("  ", item5_label, ": ", round(item5_max_jump, 3), "\n", sep = "")
print(modeled_items_max_jump)

cost_pct <- 100 - mean(item5_missing_tab$pct_item5_valid)
item5_direction_tab <- tibble(
  item5_name = item5_name, item5_type = item5_type, direction = direction,
  developed_code = developed_code, coding_matches_assumption = coding_matches_assumption,
  missingness_cost_pct = round(cost_pct, 1)
)
write.csv(item5_direction_tab, track(file.path(item5_dir, paste0("item5_direction_summary_", run_tag, ".csv"))), row.names = FALSE)
write.csv(modeled_items_max_jump, track(file.path(item5_dir, paste0("item5_steepness_comparison_", run_tag, ".csv"))), row.names = FALSE)

steepness_ratio <- item5_max_jump / mean(modeled_items_max_jump$max_jump)
add_summary(
  "SECTION 7 -- MENARCHE (fpete) CONVERGENCE DIAGNOSIS. Menarche previously caused Stage A/B/C convergence problems when included; ",
  "these descriptives point to two candidate causes, in order of how likely each is to be the actual bug. (1) CODING DIRECTION: ",
  if (coding_matches_assumption) {
    "fpete's raw coding matches the usual PDS convention (higher code = more developed), so this is NOT the cause here."
  } else {
    paste0(
      "fpete's raw response rate is ", direction, " with age, meaning the numeric coding is the OPPOSITE of the other PDS items' ",
      "'higher code = more developed' convention -- fed in as-is, the model would be forced to fit a POSITIVE-constrained loading ",
      "(lp >= 0) to an item with a NEGATIVE true relationship to the latent trait, which is a plausible direct cause of non-convergence ",
      "on its own. Fix: recode so the 'developed' response maps to 1 (bernoulli_logit) before adding it back."
    )
  },
  " (2) STEEPNESS: menarche's largest single-age-band jump in 'developed' proportion is ", round(item5_max_jump, 2),
  " (", round(100 * item5_max_jump, 0), " percentage points in a 2-year band), vs. an average of ", round(mean(modeled_items_max_jump$max_jump), 2),
  " across the four modeled items' ceiling proportions -- roughly ", round(steepness_ratio, 1), "x steeper. Even with the coding fixed, ",
  "this near-step-function transition may need a loading prior wider than lp ~ normal(1, sigma_l) (calibrated for the smoother four items) ",
  "to be sampled without difficulty. Completeness is not a concern either way: fpete is valid in ",
  paste(completeness_tab$pct_valid[completeness_tab$item == item5_name], collapse = "%/"),
  "% of parent/youth rows already passing the age/race/WtHR filters (vs. ",
  paste(round(rowMeans(matrix(completeness_tab$pct_valid[completeness_tab$item != item5_name], nrow = 2)), 1), collapse = "%/"),
  "% average for the four modeled items), and requiring it would only additionally drop about ", round(cost_pct, 1),
  "% of currently-retained person-occasions. No Stan-file changes are needed to re-add it (is_binary/k_item are already per-item vectors) -- ",
  "only the coding fix, a possibly-wider loading prior, and a sex-conditional item list in build_lmnlfa_data_staged(). Not implemented or ",
  "fit here, since it requires new Stan sampling."
)
} # end if (sx == "female")

cat("\n############################################################\n")
cat("SECTION 8: SUMMARY TEXT\n")
cat("############################################################\n")

summary_path <- file.path(out_dir, paste0("premeeting_summary_", sx, ".txt"))
header <- c(
  paste0("Pre-meeting diagnostics summary -- ", run_label),
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  paste0("Stages available: ", paste(c("A", "B", "C")[c(!is.null(fitA), !is.null(fitB), !is.null(fitC))], collapse = ", ")),
  ""
)
writeLines(c(header, summary_lines), track(summary_path))
cat("Summary written to:", summary_path, "\n")
cat(paste(c(header, summary_lines), collapse = "\n"), "\n")

cat("\n############################################################\n")
cat("DONE:", sx, "\n")
cat("Files written (", length(written_files), "):\n", sep = "")
for (f in written_files) cat(" ", f, "\n")
cat("############################################################\n")
