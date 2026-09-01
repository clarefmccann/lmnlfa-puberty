## fmri_puberty_association.R
## Associates N-back betas (fmri_nback_foundation.R) with the validated
## cross-sectional MNLFA pubertal factor scores (mnlfa_crosssectional_
## staged.R, Stage C: factor_scores_<sex>.csv), for every contrast in
## fmri_contrasts.R. Purely descriptive (scatter + correlation) -- no
## modeling.
##
## IMPORTANT CAVEAT: the cross-sectional factor score is ONE value per
## person, from a single randomly-selected wave (see build_crosssectional_
## data()'s slice_sample(n = 1) -- one wave is sampled per person
## specifically so the cross-sectional model doesn't treat repeated
## observations of the same person as independent). factor_scores_<sex>.csv
## has no wave identifier saying which occasion that score came from, so
## the person-level join below attaches the same eta_mean to every fMRI
## wave that person has. Read what follows as "how does this person's
## overall (single-occasion) pubertal level relate to their fMRI
## activation" for that version.
##
## A wave-matched version is also produced: factor_scores_<sex>.csv's
## age/whtr are the exact raw values from whichever wave was sampled, so
## the wave can be recovered by matching that age back to the puberty long
## data's per-wave ages -- this lets beta and puberty score be compared at
## the exact same occasion. The puberty-score side of this (loading
## factor scores, recovering the wave) is contrast-independent and done
## once; only the fMRI side is repeated per contrast.
##
## Usage: Rscript fmri_puberty_association.R
## Output: outputs/fmri_puberty_association_<contrast>/*.png, *.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

root_path <- Sys.getenv("HOME_DIR")
if (!nzchar(root_path)) root_path <- Sys.getenv("HOME")

data_dir <- Sys.getenv("DATA_DIR")
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  data_dir <- file.path(root_path, "projects/abcd-projs/dissertation/study1/outputs")
}
out_base <- Sys.getenv("OUT_DIR")
if (!nzchar(out_base)) out_base <- data_dir

script_dir <- tryCatch(
  {
    ofile <- sys.frames()[[1]]$ofile
    if (is.null(ofile)) "scripts" else dirname(ofile)
  },
  error = function(e) "scripts"
)
palette_file <- file.path(script_dir, "color_palette.R")
if (!file.exists(palette_file)) palette_file <- file.path("scripts", "color_palette.R")
source(palette_file)
contrasts_file <- file.path(script_dir, "fmri_contrasts.R")
if (!file.exists(contrasts_file)) contrasts_file <- file.path("scripts", "fmri_contrasts.R")
source(contrasts_file)

atlas_labeller <- as_labeller(c(aseg = "Subcortical (aseg)", dsk = "Cortical (Desikan)"))
puberty_caveat <- "Cross-sectional puberty score (one value per person, single random wave) vs. this wave's beta"
# re-key the shared wave palette/shapes (bl, fu1, ...) to the display labels
# ("Baseline", "Year 1", ...) used by wave_label
wave_palette_display <- setNames(pal_waves, wave_display_labels[names(pal_waves)])
wave_shapes_display <- setNames(pal_shapes_waves, wave_display_labels[names(pal_shapes_waves)])

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok]))
}

# ---------------------------------------------------------------------------
# PUBERTY FACTOR SCORES + WAVE RECOVERY -- contrast-independent, done once
# ---------------------------------------------------------------------------
eta_path_f <- file.path(data_dir, "mnlfa_crosssectional_staged", "factor_scores_female.csv")
eta_path_m <- file.path(data_dir, "mnlfa_crosssectional_staged", "factor_scores_male.csv")
if (!file.exists(eta_path_f) || !file.exists(eta_path_m)) {
  stop(
    "Cross-sectional factor scores not found -- run mnlfa_crosssectional_staged.R ",
    "(Stage C) for both sexes first. Expected:\n  ", eta_path_f, "\n  ", eta_path_m
  )
}
eta_all <- bind_rows(
  read.csv(eta_path_f) %>% select(id, age, whtr, eta_mean, eta_sd),
  read.csv(eta_path_m) %>% select(id, age, whtr, eta_mean, eta_sd)
)
if (any(duplicated(eta_all$id))) {
  stop("Duplicate ids across factor_scores_female.csv / factor_scores_male.csv -- investigate before proceeding.")
}

puberty_path <- file.path(data_dir, "all_long.csv")
wave_match <- NULL
if (file.exists(puberty_path)) {
  puberty_waves <- read.csv(puberty_path) %>%
    distinct(id, wave, age)

  wm <- eta_all %>%
    select(id, age_eta = age, eta_mean, eta_sd) %>%
    left_join(puberty_waves, by = "id") %>%
    mutate(age_diff = abs(age - age_eta)) %>%
    filter(!is.na(age_diff)) %>%
    group_by(id) %>%
    slice_min(age_diff, n = 1, with_ties = FALSE) %>%
    ungroup()

  n_unresolved <- length(unique(eta_all$id)) - nrow(wm)
  n_ambiguous <- sum(wm$age_diff > 0.01)
  cat(
    "Wave recovery: resolved for", nrow(wm), "of",
    length(unique(eta_all$id)), "people (", n_unresolved,
    "had no puberty-data row to match against;", n_ambiguous,
    "had age_diff > 0.01 yr and are treated as unresolved)\n"
  )
  cat("Recovered wave distribution (all people, not just those with fMRI):\n")
  print(table(wm$wave, useNA = "ifany"))

  wave_match <- wm %>%
    filter(age_diff <= 0.01) %>%
    select(id, matched_wave = wave, eta_mean, eta_sd)
} else {
  cat("Skipping wave recovery -- all_long.csv not found in", data_dir, "\n")
}

# ---------------------------------------------------------------------------
# PER-CONTRAST ASSOCIATION
# ---------------------------------------------------------------------------
run_association_for_contrast <- function(cst) {
  clabel <- fmri_contrasts[[cst]]
  cat("\n==================== Contrast:", cst, "(", clabel, ") ====================\n")

  out_dir <- file.path(out_base, paste0("fmri_puberty_association_", cst))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  fmri_path <- file.path(data_dir, paste0("all_fmri_nback_", cst, "_long.csv"))
  lookup_path <- file.path(data_dir, paste0("fmri_nback_", cst, "_region_lookup.csv"))
  if (!file.exists(fmri_path) || !file.exists(lookup_path)) {
    cat("Skipping", cst, "-- missing", fmri_path, "or", lookup_path, "\n")
    return(invisible(NULL))
  }

  fmri <- read.csv(fmri_path)
  region_lookup <- read.csv(lookup_path)

  region_cols <- region_lookup$region_key
  region_cols <- region_cols[region_cols %in% names(fmri)]
  if (length(region_cols) == 0) {
    cat("Skipping", cst, "-- no region columns from the lookup table found in", fmri_path, "\n")
    return(invisible(NULL))
  }

  fmri <- fmri %>%
    mutate(
      wave = factor(wave, levels = wave_codes),
      wave_label = factor(wave_display_labels[as.character(wave)], levels = wave_display_labels),
      sex_label = case_when(
        sex == 1 ~ "Male",
        sex == 2 ~ "Female",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(wave), !is.na(sex_label))

  merged <- fmri %>% inner_join(eta_all %>% select(id, eta_mean, eta_sd), by = "id")
  cat(
    "fMRI person-waves:", nrow(fmri),
    "| with a matched cross-sectional puberty score:", nrow(merged),
    "(", round(100 * nrow(merged) / nrow(fmri), 1), "% )\n"
  )
  cat("Matched rows by wave:\n")
  print(table(merged$wave_label, useNA = "ifany"))

  write.csv(
    merged,
    file.path(out_dir, "fmri_puberty_merged.csv"),
    row.names = FALSE
  )

  matched_merged <- NULL
  matched_beta_long <- NULL
  if (!is.null(wave_match)) {
    matched_merged <- fmri %>%
      inner_join(wave_match, by = c("id", "wave" = "matched_wave"))

    cat(
      "People with a resolved wave AND fMRI data at that exact wave:",
      nrow(matched_merged), "of", nrow(wave_match),
      "people with a resolved wave (",
      round(100 * nrow(matched_merged) / nrow(wave_match), 1),
      "% -- the gap is people whose cross-sectional wave wasn't an imaging",
      "wave, or didn't pass fMRI QC that wave)\n"
    )
    cat("Wave-matched rows by wave:\n")
    print(table(matched_merged$wave_label, useNA = "ifany"))

    write.csv(
      matched_merged,
      file.path(out_dir, "fmri_puberty_wave_matched.csv"),
      row.names = FALSE
    )

    matched_beta_long <- matched_merged %>%
      select(id, wave, wave_label, sex_label, eta_mean, all_of(region_cols)) %>%
      pivot_longer(all_of(region_cols), names_to = "region_key", values_to = "beta") %>%
      filter(!is.na(beta)) %>%
      left_join(region_lookup, by = "region_key")
  }

  # long format: one row per person-wave-region, with eta_mean carried along
  beta_long <- merged %>%
    select(id, wave, wave_label, sex_label, eta_mean, all_of(region_cols)) %>%
    pivot_longer(all_of(region_cols), names_to = "region_key", values_to = "beta") %>%
    filter(!is.na(beta)) %>%
    left_join(region_lookup, by = "region_key")

  # -------------------------------------------------------------------------
  # (1) Correlation table: Pearson r(beta, eta_mean) by region x wave x sex,
  # and pooled across wave/sex -- a lightweight numeric summary of association
  # -------------------------------------------------------------------------
  region_wave_sex_cor <- beta_long %>%
    group_by(region_key, atlas, hemi, region_label, wave_label, sex_label) %>%
    summarise(r = safe_cor(beta, eta_mean), n = n(), .groups = "drop") %>%
    arrange(desc(abs(r)))
  write.csv(
    region_wave_sex_cor,
    file.path(out_dir, "region_wave_sex_puberty_correlation.csv"),
    row.names = FALSE
  )

  region_pooled_cor <- beta_long %>%
    group_by(region_key, atlas, hemi, region_label) %>%
    summarise(r = safe_cor(beta, eta_mean), n = n(), .groups = "drop") %>%
    arrange(desc(abs(r)))
  write.csv(
    region_pooled_cor,
    file.path(out_dir, "region_puberty_correlation_pooled.csv"),
    row.names = FALSE
  )

  cat("Top 10 |r(beta, puberty factor score)| regions (pooled across wave/sex):\n")
  print(as.data.frame(head(region_pooled_cor, 10)), digits = 3)

  # -------------------------------------------------------------------------
  # (2) Global mean beta (averaged across all regions) vs. puberty factor
  # score, faceted by wave, colour+shape+linetype by sex
  # -------------------------------------------------------------------------
  person_wave_mean <- beta_long %>%
    group_by(id, wave_label, sex_label, eta_mean) %>%
    summarise(mean_beta = mean(beta), .groups = "drop")

  p_global <- ggplot(
    person_wave_mean,
    aes(x = eta_mean, y = mean_beta, colour = sex_label, shape = sex_label)
  ) +
    geom_point(alpha = 0.25, size = 1) +
    geom_smooth(aes(linetype = sex_label), method = "loess", se = TRUE, linewidth = 1) +
    facet_wrap(~wave_label) +
    scale_colour_manual(values = c(Female = pal_two[1], Male = pal_two[2])) +
    scale_linetype_manual(values = c(Female = pal_linetypes_two[1], Male = pal_linetypes_two[2])) +
    scale_shape_manual(values = c(Female = pal_shapes_two[1], Male = pal_shapes_two[2])) +
    labs(
      title = paste0("Mean N-back ", clabel, " beta vs. pubertal factor score"),
      subtitle = puberty_caveat,
      x = "Puberty factor score (cross-sectional eta)",
      y = paste0("Mean beta (averaged across all ", length(region_cols), " regions)"),
      colour = "Sex",
      shape = "Sex",
      linetype = "Sex"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
  ggsave(
    file.path(out_dir, "mean_beta_vs_puberty_score.png"),
    p_global,
    width = 9,
    height = 6,
    dpi = 180
  )

  # -------------------------------------------------------------------------
  # (3) Region-specific scatter: amygdala, insula, rostral middle frontal --
  # left/right hemispheres averaged per person-wave-region here, since the
  # question is "does this region track puberty," not a hemisphere contrast
  # -------------------------------------------------------------------------
  region_scatter <- function(label, filename) {
    d <- beta_long %>%
      filter(region_label == label) %>%
      group_by(id, wave_label, sex_label, eta_mean) %>%
      summarise(beta = mean(beta), .groups = "drop")

    if (nrow(d) == 0) {
      cat("Skipping", filename, "-- no rows for", label, "\n")
      return(invisible(NULL))
    }

    p <- ggplot(d, aes(x = eta_mean, y = beta, colour = sex_label, shape = sex_label)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_point(alpha = 0.3, size = 1) +
      geom_smooth(aes(linetype = sex_label), method = "loess", se = TRUE, linewidth = 1) +
      facet_wrap(~wave_label) +
      scale_colour_manual(values = c(Female = pal_two[1], Male = pal_two[2])) +
      scale_linetype_manual(values = c(Female = pal_linetypes_two[1], Male = pal_linetypes_two[2])) +
      scale_shape_manual(values = c(Female = pal_shapes_two[1], Male = pal_shapes_two[2])) +
      labs(
        title = paste0("N-back ", clabel, " beta (", label, ") vs. pubertal factor score"),
        subtitle = paste0(puberty_caveat, "; left/right hemispheres averaged"),
        x = "Puberty factor score (cross-sectional eta)",
        y = "Beta",
        colour = "Sex",
        shape = "Sex",
        linetype = "Sex"
      ) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
    ggsave(file.path(out_dir, filename), p, width = 9, height = 6, dpi = 180)
  }

  region_scatter("amygdala", "amygdala_vs_puberty_score.png")
  region_scatter("insula", "insula_vs_puberty_score.png")
  region_scatter("rostral middle frontal", "rostral_middle_frontal_vs_puberty_score.png")

  # -------------------------------------------------------------------------
  # (4) Wave-matched region-specific scatter: beta and puberty score from the
  # exact same occasion, one point per person -- points coloured/shaped by
  # which wave that occasion was, since wave now varies BETWEEN people
  # rather than within, and a single loess per sex panel (not split further
  # by wave -- per-wave-per-sex cells get thin).
  # -------------------------------------------------------------------------
  region_scatter_matched <- function(label, filename) {
    if (is.null(matched_merged)) {
      cat("Skipping", filename, "-- no wave-matched data available\n")
      return(invisible(NULL))
    }
    d <- matched_beta_long %>%
      filter(region_label == label) %>%
      group_by(id, wave_label, sex_label, eta_mean) %>%
      summarise(beta = mean(beta), .groups = "drop")

    if (nrow(d) == 0) {
      cat("Skipping", filename, "-- no rows for", label, "\n")
      return(invisible(NULL))
    }

    p <- ggplot(d, aes(x = eta_mean, y = beta)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_point(aes(colour = wave_label, shape = wave_label), alpha = 0.6, size = 1.6) +
      geom_smooth(method = "loess", se = TRUE, colour = "black", linewidth = 0.8) +
      facet_wrap(~sex_label) +
      scale_colour_manual(values = wave_palette_display, breaks = as.character(wave_display_labels)) +
      scale_shape_manual(values = wave_shapes_display, breaks = as.character(wave_display_labels)) +
      labs(
        title = paste0("N-back ", clabel, " beta (", label, ") vs. pubertal factor score - same wave"),
        subtitle = paste0(
          "Beta and puberty score both from the exact wave the cross-sectional model used; ",
          "n = ", nrow(d), " people; left/right hemispheres averaged"
        ),
        x = "Puberty factor score (cross-sectional eta)",
        y = "Beta",
        colour = "Wave",
        shape = "Wave"
      ) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")
    ggsave(file.path(out_dir, filename), p, width = 9, height = 5.5, dpi = 180)
  }

  region_scatter_matched("amygdala", "amygdala_vs_puberty_score_wave_matched.png")
  region_scatter_matched("insula", "insula_vs_puberty_score_wave_matched.png")
  region_scatter_matched("rostral middle frontal", "rostral_middle_frontal_vs_puberty_score_wave_matched.png")

  cat("Outputs written to:", out_dir, "\n")
}

for (cst in names(fmri_contrasts)) {
  run_association_for_contrast(cst)
}
