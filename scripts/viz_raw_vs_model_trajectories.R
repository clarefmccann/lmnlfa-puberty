## viz_raw_vs_model_trajectories.R
## Companion figures to lmnlfa_growth_sigmoid_staged.R: RAW observed
## pubertal-item trajectories (no measurement model at all) alongside the
## model-implied latent trajectory, to visually show what the staged
## sigmoidal L-MNLFA is correcting for.
##
## Three figure types, for both sexes:
##   1. "raw_trajectories_<sex>_<reporter>.png" -- individual raw
##      trajectories (mean of the 4 pubertal items, 1-4 scale) by race,
##      one file per sex x reporter (4 total), sampling 200 of the FULL
##      filtered pool (~5,500-6,500 people/sex). Direct raw-scale
##      companion to the model's own trajectories_by_race_<run_tag>.png.
##   1b. "raw_trajectories_modeled_sample_<sex>_<reporter>.png" -- the
##      same plot, but showing EVERY person (not a 200-person draw) in
##      the actual n=1500-per-sex sample the staged Stan models were fit
##      on -- "what does my modeled sample's raw data actually look like."
##   2. "raw_vs_model_<sex>.png" -- THE headline comparison, one file per
##      sex (2 total). Left panel: raw population mean trajectory,
##      caregiver vs youth report, smoothed with 90% CI ribbons -- this is
##      where reporter (informant) measurement bias shows up as two
##      diverging curves. Right panel: the model's single, informant-
##      independent latent trajectory (from whichever staged fit is most
##      advanced on disk for that sex: Stage C > B > A, clearly labeled).
##      Because this model only ever lets race/WtHR shift the growth
##      curve's mean (never informant -- informant enters only as an
##      item-level DIF covariate), a single coherent latent trajectory
##      where the two raw reporter curves diverge IS the visual evidence
##      that reporter bias is being absorbed into item parameters rather
##      than distorting the estimated trajectory itself.
##
## Sample: matches lmnlfa_growth_sigmoid_staged.R's own filtering (age/
## race/WtHR non-missing, WtHR in [.25,.75], all 4 items scored 1-4) so
## the "raw" picture reflects the same people the model actually saw.
## Figure 1 uses the FULL filtered sample (no need to subsample just to
## draw a scatter/smooth). Figure 2's raw panel is restricted to
## whichever subsample the overlaid model fit was actually trained on
## (reconstructed via the identical seed + build order used in
## lmnlfa_growth_sigmoid_staged.R), so the two panels are a genuine
## before/after on the SAME people -- overlaying a full-sample raw curve
## against a subsample-trained model curve would not be apples-to-apples.
##
## Usage:  Rscript viz_raw_vs_model_trajectories.R
## (no arguments -- always does both sexes x both reporters)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(patchwork)
})

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  stop("cmdstanr not found -- needed to read cached model .rds fits.")
}
library(cmdstanr)

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

# DATA_DIR (explicit override, e.g. set by an HPC job script) wins if set.
# Otherwise prefer the known sshfs project path over HOME_DIR/HOME-derived
# guesses -- HOME_DIR is sometimes aliased locally to an unrelated synced
# mirror that does NOT receive freshly-written HPC model outputs, which
# would silently (and wrongly) make this script think no model has been
# fit yet. Only fall back to HOME_DIR/HOME construction if the sshfs path
# isn't present at all (e.g. running directly on the cluster).
data_dir <- Sys.getenv("DATA_DIR")
sshfs_data_dir <- "/private/tmp/sshfs/projects/abcd-projs/dissertation/study1/outputs"
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  if (dir.exists(sshfs_data_dir)) {
    data_dir <- sshfs_data_dir
  } else {
    root_path <- Sys.getenv("HOME_DIR")
    if (!nzchar(root_path)) root_path <- Sys.getenv("HOME")
    data_dir <- file.path(root_path, "projects/abcd-projs/dissertation/study1/outputs")
  }
}
if (!dir.exists(data_dir)) stop("Cannot locate data directory: ", data_dir)

model_out_dir <- file.path(data_dir, "lmnlfa_growth_sigmoid_staged")
model_fits_dir <- file.path(model_out_dir, "fits") # cached .rds fits live here, not flat in model_out_dir

out_dir <- file.path(data_dir, "viz_raw_vs_model")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Data dir: ", data_dir, "\n")
cat("Model dir:", model_out_dir, "\n")
cat("Out dir:  ", out_dir, "\n")

ordinal_items <- c("peta", "petb", "petc", "petd")
wave_order <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")
whtr_min <- 0.25
whtr_max <- 0.75

# ---------------------------------------------------------------------------
# BUILD FILTERED RAW SAMPLE -- same filtering rule as
# lmnlfa_growth_sigmoid_staged.R's build_lmnlfa_data_staged(), but keeping
# reporter-separated raw item scores instead of pivoting for Stan. Reusing
# the exact same filter (not just "similar") matters: it's what makes
# figure 2's raw panel a fair "before" picture of what the model's
# "after" panel is actually correcting.
# ---------------------------------------------------------------------------
build_raw_sample <- function(parent_df, youth_df, n_subsample = NULL) {
  clean_reporter <- function(df, reporter_label, informant_val) {
    df %>%
      select(id, wave, age, race, whtr, all_of(ordinal_items)) %>%
      filter(!is.na(age), !is.na(race), !is.na(whtr), whtr >= whtr_min, whtr <= whtr_max) %>%
      filter(if_all(all_of(ordinal_items), ~ !is.na(.) & as.integer(.) %in% 1:4)) %>%
      mutate(
        wave = factor(wave, levels = wave_order),
        reporter = reporter_label,
        informant_c = informant_val,
        raw_score = rowMeans(across(all_of(ordinal_items), as.numeric))
      )
  }

  dat <- bind_rows(
    clean_reporter(parent_df, "Caregiver-reported", 1),
    clean_reporter(youth_df, "Youth-reported", -1)
  ) %>%
    mutate(race_grp = case_when(
      race == 1 ~ "Hispanic", race == 2 ~ "White", race == 3 ~ "Black",
      race %in% c(7, 11, 12, 13) ~ "Other", TRUE ~ NA_character_
    )) %>%
    filter(!is.na(race_grp)) %>%
    mutate(race_grp = factor(race_grp, levels = c("Hispanic", "White", "Black", "Other")))

  all_ids <- sort(unique(dat$id))
  age_mean_full <- mean(dat$age, na.rm = TRUE)
  age_sd_full <- sd(dat$age, na.rm = TRUE)

  if (!is.null(n_subsample) && n_subsample < length(all_ids)) {
    sub_ids <- sort(sample(all_ids, n_subsample))
    dat <- dat %>% filter(id %in% sub_ids)
  }

  list(dat = dat %>% arrange(id, wave, reporter), age_mean = age_mean_full, age_sd = age_sd_full)
}

# ---------------------------------------------------------------------------
# LOCATE + LOAD THE MOST-ADVANCED CACHED SIGMOID FIT FOR A SEX. Returns
# NULL (with a message) if nothing is cached yet for that sex.
# ---------------------------------------------------------------------------
find_best_fit <- function(sex) {
  stages <- list(
    c("fitC_final_", "Stage C (final, DIF-adjusted)"),
    c("fitB_difscreen_", "Stage B (DIF screening; impact fixed at Stage A)"),
    c("fitA_impact_", "Stage A (impact only; informant DIF not yet estimated)")
  )
  short_labels <- c("Stage C (final)", "Stage B (DIF screening)", "Stage A (impact only)")
  for (i in seq_along(stages)) {
    s <- stages[[i]]
    hits <- list.files(model_fits_dir, pattern = paste0("^", s[1], sex, ".*\\.rds$"), full.names = TRUE)
    if (length(hits) > 0) {
      f <- hits[1]
      n_match <- regmatches(basename(f), regexpr("(?<=_n)[0-9]+", basename(f), perl = TRUE))
      n_sub <- if (length(n_match) > 0) as.integer(n_match) else NULL
      cat("  Using", s[2], "for", sex, "->", basename(f), "\n")
      return(list(fit = readRDS(f), stage_label = s[2], stage_label_short = short_labels[i], n_subsample = n_sub))
    }
  }
  cat("  No cached sigmoid fit found for", sex, "-- model panel will be skipped.\n")
  NULL
}

# ---------------------------------------------------------------------------
# FIGURE 1: individual raw trajectories by race, one per sex x reporter
# ---------------------------------------------------------------------------
race_levels <- c("Hispanic", "White", "Black", "Other")

make_raw_spaghetti <- function(dat, sex, reporter_label, run_label, n_display = 200, file_tag = NULL, subtitle_note = NULL) {
  sub <- dat %>% filter(reporter == reporter_label)
  ids_here <- unique(sub$id)
  set.seed(90025)
  n_show <- min(n_display, length(ids_here))
  samp_ids <- if (n_show >= length(ids_here)) ids_here else sample(ids_here, n_show)
  plot_df <- sub %>% filter(id %in% samp_ids)

  subtitle_text <- if (!is.null(subtitle_note)) {
    subtitle_note
  } else {
    paste0(
      "Mean of the 4 raw pubertal items (no measurement model applied); n = ",
      length(samp_ids), " sampled of ", length(ids_here), " people"
    )
  }

  p <- ggplot(plot_df, aes(x = age, y = raw_score, group = id, colour = race_grp)) +
    geom_line(alpha = 0.2, linewidth = 0.3) +
    scale_colour_manual(values = setNames(pal_chains, race_levels)) +
    coord_cartesian(ylim = c(1, 4)) +
    labs(
      title = paste0("Raw ", tolower(reporter_label), " puberty trajectories by race - ", run_label),
      subtitle = subtitle_text,
      x = "Age (years)",
      y = "Raw mean item score (1-4 scale)",
      colour = "Race/ethnicity"
    ) +
    theme_minimal(base_size = 13)

  reporter_tag <- ifelse(grepl("Caregiver", reporter_label), "caregiver", "youth")
  name_tag <- if (!is.null(file_tag)) paste0(file_tag, "_", sex, "_", reporter_tag) else paste0(sex, "_", reporter_tag)
  ggsave(
    file.path(out_dir, paste0("raw_trajectories_", name_tag, ".png")),
    p, width = 10, height = 6.5, dpi = 220
  )
  invisible(p)
}

# ---------------------------------------------------------------------------
# FIGURE 2: raw (caregiver vs youth) vs model-implied trajectory, per sex
# ---------------------------------------------------------------------------
make_raw_vs_model <- function(raw_dat_model_sample, best_fit, sex, run_label) {
  age_range_raw <- range(raw_dat_model_sample$age, na.rm = TRUE)

  p_raw <- ggplot(raw_dat_model_sample, aes(x = age, y = raw_score, colour = reporter, linetype = reporter, fill = reporter)) +
    geom_smooth(method = "loess", se = TRUE, level = 0.90, span = 0.75, linewidth = 1.1) +
    scale_colour_manual(values = setNames(pal_two, c("Caregiver-reported", "Youth-reported"))) +
    scale_linetype_manual(values = setNames(pal_linetypes_two, c("Caregiver-reported", "Youth-reported"))) +
    scale_fill_manual(values = setNames(pal_two, c("Caregiver-reported", "Youth-reported"))) +
    coord_cartesian(xlim = age_range_raw, ylim = c(1, 4)) +
    labs(
      title = "Raw observed",
      subtitle = "Loess-smoothed mean item score + 90% CI, by reporter",
      x = "Age (years)",
      y = "Raw mean item score (1-4 scale)",
      colour = NULL, linetype = NULL, fill = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")

  if (!is.null(best_fit)) {
    # restrict the model curve to the SAME age range actually observed in
    # this sample -- extrapolating out to +/-3 SD would draw the sigmoid
    # well past any data support and mislead the raw-vs-model comparison
    age_c_range <- (age_range_raw - best_fit$age_mean) / best_fit$age_sd
    age_grid <- seq(age_c_range[1], age_c_range[2], length.out = 100)
    mu_logk_draws <- as.vector(best_fit$fit$draws(variables = "mu_logk", format = "matrix"))
    mu_alpha_draws <- as.vector(best_fit$fit$draws(variables = "mu_alpha", format = "matrix"))
    mean_traj <- map_dfr(age_grid, function(a) {
      vals <- 1 + 4 / (1 + exp(-exp(mu_logk_draws) * (a - mu_alpha_draws)))
      tibble(age_c = a, eta_med = median(vals), eta_lo = quantile(vals, 0.05), eta_hi = quantile(vals, 0.95))
    }) %>%
      mutate(age = age_c * best_fit$age_sd + best_fit$age_mean)

    p_model <- ggplot(mean_traj, aes(x = age)) +
      geom_ribbon(aes(ymin = eta_lo, ymax = eta_hi), fill = pal_primary_fill, alpha = 0.35) +
      geom_line(aes(y = eta_med), colour = pal_primary, linewidth = 1.2) +
      coord_cartesian(xlim = age_range_raw, ylim = c(1, 5)) +
      labs(
        title = "Model-implied (informant-independent)",
        subtitle = paste0(best_fit$stage_label_short, "; median + 90% CI"),
        x = "Age (years)",
        y = "Latent puberty (eta, 1-5 scale)"
      ) +
      theme_minimal(base_size = 13)
  } else {
    p_model <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No cached model fit yet", size = 5, colour = "grey40") +
      theme_void() +
      labs(title = "Model-implied (informant-independent)", subtitle = "Run lmnlfa_growth_sigmoid_staged.R first")
  }

  combined <- p_raw + p_model +
    plot_annotation(
      title = paste0("Raw reporter divergence vs. model-corrected trajectory - ", run_label),
      subtitle = "If the model is absorbing reporter measurement bias, the two raw curves at left should diverge more than the single curve at right"
    )

  ggsave(file.path(out_dir, paste0("raw_vs_model_", sex, ".png")), combined, width = 14, height = 7, dpi = 220)
  invisible(combined)
}

# ---------------------------------------------------------------------------
# RUN FOR BOTH SEXES
# ---------------------------------------------------------------------------
for (sx in c("female", "male")) {
  cat("\n=== ", sx, " ===\n")
  parent_df <- read.csv(file.path(data_dir, paste0(sx, "_parent_long.csv")))
  youth_df <- read.csv(file.path(data_dir, paste0(sx, "_youth_long.csv")))

  # Figure 1: full filtered sample, no subsampling needed for a descriptive plot
  full_build <- build_raw_sample(parent_df, youth_df, n_subsample = NULL)
  cat("  Full filtered sample:", length(unique(full_build$dat$id)), "people |", nrow(full_build$dat), "person-occasion-reporter rows\n")

  for (rep_label in c("Caregiver-reported", "Youth-reported")) {
    make_raw_spaghetti(full_build$dat, sx, rep_label, run_label = sx)
  }

  # Figure 2: best cached model fit + matching subsample
  best_fit <- find_best_fit(sx)
  n_sub <- if (!is.null(best_fit)) best_fit$n_subsample else 1500

  set.seed(90025) # matches lmnlfa_growth_sigmoid_staged.R's seed + build order exactly
  model_sample_build <- build_raw_sample(parent_df, youth_df, n_subsample = n_sub)
  if (!is.null(best_fit)) {
    best_fit$age_mean <- model_sample_build$age_mean
    best_fit$age_sd <- model_sample_build$age_sd
  }
  n_sub_label <- if (is.null(n_sub)) "full" else n_sub
  cat("  Model-matched sample:", length(unique(model_sample_build$dat$id)), "people (n_subsample =", n_sub_label, ")\n")

  make_raw_vs_model(model_sample_build$dat, best_fit, sx, run_label = sx)

  # raw trajectories for the ACTUAL n=1500-per-sex sample the staged Stan
  # models were fit on (not a 200-person draw from the full ~5,500-6,500
  # person pool above) -- every person in the modeled sample is drawn, to
  # show the real messiness of the raw data the model is working with
  for (rep_label in c("Caregiver-reported", "Youth-reported")) {
    reporter_ids <- unique(model_sample_build$dat$id[model_sample_build$dat$reporter == rep_label])
    make_raw_spaghetti(
      model_sample_build$dat,
      sx,
      rep_label,
      run_label = paste0(sx, " (n=", n_sub_label, " staged-model sample)"),
      n_display = Inf,
      file_tag = "modeled_sample",
      subtitle_note = paste0(
        "Mean of the 4 raw pubertal items (no measurement model applied); ALL ",
        length(reporter_ids), " people in the staged-model sample"
      )
    )
  }
}

cat("\nAll figures written to:", out_dir, "\n")
