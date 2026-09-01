## fmri_nback_descriptives.R
## Descriptive summaries and plots of the cleaned N-Back beta data
## (fmri_nback_foundation.R), for every contrast in fmri_contrasts.R: how
## betas change by wave, differ by sex, and vary across regions. Purely
## descriptive -- no modeling.
##
## Usage: Rscript fmri_nback_descriptives.R
## Output: outputs/fmri_nback_descriptives_<contrast>/*.png, *.csv

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
# re-key the shared wave palette/shapes (bl, fu1, ...) to the display labels
# ("Baseline", "Year 1", ...) used by wave_label
wave_palette_display <- setNames(pal_waves, wave_display_labels[names(pal_waves)])

# Boxplot y-limits, computed from the plotted data itself: with beta's true
# range dominated by rare outliers, the default full-range y-axis squashes
# the boxes themselves down to a sliver. coord_cartesian() only changes the
# visible window -- boxplot statistics (and outlier points beyond the
# window) are still computed from the complete data, nothing is dropped.
tight_ylim <- function(x, probs = c(0.02, 0.98), pad = 0.15) {
  r <- quantile(x, probs, na.rm = TRUE)
  p <- diff(r) * pad
  c(r[1] - p, r[2] + p)
}

run_descriptives_for_contrast <- function(cst) {
  clabel <- fmri_contrasts[[cst]]
  cat("\n==================== Contrast:", cst, "(", clabel, ") ====================\n")

  out_dir <- file.path(out_base, paste0("fmri_nback_descriptives_", cst))
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

  cat("Rows:", nrow(fmri), "| people:", length(unique(fmri$id)), "\n")
  cat("Rows by wave:\n")
  print(table(fmri$wave_label, useNA = "ifany"))
  cat("\nRows by sex:\n")
  print(table(fmri$sex_label, useNA = "ifany"))

  # long format: one row per person-wave-region
  beta_long <- fmri %>%
    select(id, wave, wave_label, sex_label, all_of(region_cols)) %>%
    pivot_longer(all_of(region_cols), names_to = "region_key", values_to = "beta") %>%
    filter(!is.na(beta)) %>%
    left_join(region_lookup, by = "region_key")

  # ---------------------------------------------------------------------------
  # (1) SUMMARY TABLE: mean/sd/n beta by region x wave x sex
  # ---------------------------------------------------------------------------
  region_wave_sex_summ <- beta_long %>%
    group_by(region_key, atlas, hemi, region_label, wave_label, sex_label) %>%
    summarise(
      mean_beta = mean(beta),
      sd_beta = sd(beta),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(atlas, region_label, hemi, wave_label, sex_label)
  write.csv(
    region_wave_sex_summ,
    file.path(out_dir, "region_wave_sex_summary.csv"),
    row.names = FALSE
  )

  # ---------------------------------------------------------------------------
  # (2) Distribution of beta by wave, faceted by atlas -- all regions pooled.
  # Boxes sit at distinct x positions (one per wave), so position alone
  # already disambiguates them; colour here is a visual/consistency aid, not
  # the only channel carrying wave identity.
  # ---------------------------------------------------------------------------
  p_wave_box <- ggplot(beta_long, aes(x = wave_label, y = beta, colour = wave_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_boxplot(outlier.alpha = 0.05, outlier.size = 0.3, linewidth = 0.4) +
    facet_wrap(~atlas, labeller = atlas_labeller) +
    scale_colour_manual(values = wave_palette_display, guide = "none") +
    coord_cartesian(ylim = tight_ylim(beta_long$beta)) +
    labs(
      title = paste0("N-back ", clabel, " beta by wave"),
      subtitle = "Distribution pooled across all regions and people",
      x = NULL,
      y = "Beta"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(
    file.path(out_dir, "beta_by_wave_boxplot.png"),
    p_wave_box,
    width = 9,
    height = 5,
    dpi = 180
  )

  # ---------------------------------------------------------------------------
  # (3) Distribution of beta by wave x sex, faceted by atlas -- overlapping/
  # adjacent dodged groups within each wave cluster, so sex also gets a
  # redundant linetype (outline dash pattern) on top of colour+dodge position.
  # ---------------------------------------------------------------------------
  p_wave_sex_box <- ggplot(
    beta_long,
    aes(x = wave_label, y = beta, colour = sex_label, linetype = sex_label)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_boxplot(outlier.alpha = 0.03, outlier.size = 0.2, linewidth = 0.5) +
    facet_wrap(~atlas, labeller = atlas_labeller, ncol = 1) +
    scale_colour_manual(values = c(Female = pal_two[1], Male = pal_two[2])) +
    scale_linetype_manual(values = c(Female = pal_linetypes_two[1], Male = pal_linetypes_two[2])) +
    coord_cartesian(ylim = tight_ylim(beta_long$beta)) +
    labs(
      title = paste0("N-back ", clabel, " beta by wave and sex"),
      subtitle = "Distribution pooled across all regions",
      x = NULL,
      y = "Beta",
      colour = "Sex",
      linetype = "Sex"
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
  ggsave(
    file.path(out_dir, "beta_by_wave_sex_boxplot.png"),
    p_wave_sex_box,
    width = 8,
    height = 8,
    dpi = 180
  )

  # ---------------------------------------------------------------------------
  # (3b) Pubertal status: join in the puberty pipeline's own pds_categ (ABCD's
  # 1-5 categorical stage: prepubertal/early/mid/late/postpubertal). One
  # status per person-wave is needed, but puberty data has a row per
  # reporter (parent/youth); prefer youth self-report, falling back to
  # parent report when youth is missing, rather than inventing a new
  # composite -- this reuses ABCD's own already-computed category as-is.
  # ---------------------------------------------------------------------------
  puberty_path <- file.path(data_dir, "all_long.csv")
  beta_long_stage <- NULL
  if (file.exists(puberty_path)) {
    puberty_status <- read.csv(puberty_path) %>%
      select(id, wave, reporter, pds_categ) %>%
      filter(!is.na(pds_categ)) %>%
      mutate(reporter_rank = if_else(reporter == "youth", 1, 2)) %>%
      arrange(id, wave, reporter_rank) %>%
      distinct(id, wave, .keep_all = TRUE) %>%
      mutate(
        wave = as.character(wave),
        pds_categ = as.integer(pds_categ),
        stage_label = factor(
          pds_categ,
          levels = 1:5,
          labels = c("Prepubertal", "Early", "Mid", "Late", "Postpubertal")
        )
      ) %>%
      select(id, wave, stage_label)

    beta_long_stage <- beta_long %>%
      mutate(wave = as.character(wave)) %>%
      inner_join(puberty_status, by = c("id", "wave"))

    cat(
      "\nRows with a matched pubertal status:", nrow(beta_long_stage),
      "of", nrow(beta_long), "\n"
    )

    pal_stages <- setNames(
      sunset_ramp(5),
      c("Prepubertal", "Early", "Mid", "Late", "Postpubertal")
    )

    p_stage_box <- ggplot(
      beta_long_stage,
      aes(x = stage_label, y = beta, colour = stage_label)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_boxplot(outlier.alpha = 0.05, outlier.size = 0.3, linewidth = 0.4) +
      facet_wrap(~atlas, labeller = atlas_labeller) +
      scale_colour_manual(values = pal_stages, guide = "none") +
      coord_cartesian(ylim = tight_ylim(beta_long_stage$beta)) +
      labs(
        title = paste0("N-back ", clabel, " beta by pubertal status"),
        subtitle = "Distribution pooled across all regions and people; pds_categ, youth-preferred",
        x = NULL,
        y = "Beta"
      ) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1))
    ggsave(
      file.path(out_dir, "beta_by_pubertal_status_boxplot.png"),
      p_stage_box,
      width = 9,
      height = 5,
      dpi = 180
    )

    p_stage_sex_box <- ggplot(
      beta_long_stage,
      aes(x = stage_label, y = beta, colour = sex_label, linetype = sex_label)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_boxplot(outlier.alpha = 0.03, outlier.size = 0.2, linewidth = 0.5) +
      facet_wrap(~atlas, labeller = atlas_labeller, ncol = 1) +
      scale_colour_manual(values = c(Female = pal_two[1], Male = pal_two[2])) +
      scale_linetype_manual(values = c(Female = pal_linetypes_two[1], Male = pal_linetypes_two[2])) +
      coord_cartesian(ylim = tight_ylim(beta_long_stage$beta)) +
      labs(
        title = paste0("N-back ", clabel, " beta by pubertal status and sex"),
        subtitle = "Distribution pooled across all regions; pds_categ, youth-preferred",
        x = NULL,
        y = "Beta",
        colour = "Sex",
        linetype = "Sex"
      ) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "bottom")
    ggsave(
      file.path(out_dir, "beta_by_pubertal_status_sex_boxplot.png"),
      p_stage_sex_box,
      width = 8,
      height = 8,
      dpi = 180
    )
  } else {
    cat("\nSkipping pubertal-status plots -- all_long.csv not found in", data_dir, "\n")
  }

  # ---------------------------------------------------------------------------
  # (4) Global mean beta trajectory across waves, by sex (with 90% CI ribbon).
  # One value per person-wave = mean beta across all regions, then
  # mean/CI across people -- the two sex lines/ribbons overlap in the same
  # space, so this is exactly the case that needs colour+linetype redundancy.
  # ---------------------------------------------------------------------------
  person_wave_mean <- beta_long %>%
    group_by(id, wave_label, sex_label) %>%
    summarise(mean_beta = mean(beta), .groups = "drop")

  traj_summ <- person_wave_mean %>%
    group_by(wave_label, sex_label) %>%
    summarise(
      est = mean(mean_beta),
      se = sd(mean_beta) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(lo = est - 1.645 * se, hi = est + 1.645 * se)

  p_traj <- ggplot(
    traj_summ,
    aes(x = wave_label, y = est, colour = sex_label, linetype = sex_label, group = sex_label)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_ribbon(
      aes(ymin = lo, ymax = hi, fill = sex_label),
      alpha = 0.2,
      colour = NA
    ) +
    geom_line(linewidth = 1.1) +
    geom_point(aes(shape = sex_label), size = 2.2) +
    scale_colour_manual(values = c(Female = pal_two[1], Male = pal_two[2])) +
    scale_fill_manual(values = c(Female = pal_two[1], Male = pal_two[2])) +
    scale_linetype_manual(values = c(Female = pal_linetypes_two[1], Male = pal_linetypes_two[2])) +
    scale_shape_manual(values = c(Female = pal_shapes_two[1], Male = pal_shapes_two[2])) +
    labs(
      title = paste0("Mean N-back ", clabel, " beta across waves, by sex"),
      subtitle = paste0(
        "Averaged across all ", length(region_cols), " regions and all people; mean +/- 90% CI"
      ),
      x = NULL,
      y = "Mean beta",
      colour = "Sex",
      fill = "Sex",
      linetype = "Sex",
      shape = "Sex"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")
  ggsave(
    file.path(out_dir, "mean_beta_trajectory_by_sex.png"),
    p_traj,
    width = 7.5,
    height = 5,
    dpi = 180
  )

  # ---------------------------------------------------------------------------
  # (5) Region x wave heatmap, faceted by atlas x sex -- mean beta per cell.
  # Continuous fill (not a categorical grouping), so the accessibility fix
  # here is a genuinely diverging, colorblind-legible palette rather than a
  # redundant shape/linetype channel.
  # ---------------------------------------------------------------------------
  region_order <- region_lookup %>%
    arrange(atlas, region_label, hemi) %>%
    mutate(region_tag = paste0(region_label, ifelse(is.na(hemi), "", paste0(" (", hemi, ")")))) %>%
    pull(region_tag) %>%
    unique()

  heat_df <- region_wave_sex_summ %>%
    mutate(
      region_tag = paste0(region_label, ifelse(is.na(hemi), "", paste0(" (", hemi, ")"))),
      region_tag = factor(region_tag, levels = rev(region_order))
    )

  for (a in c("aseg", "dsk")) {
    hd <- heat_df %>% filter(atlas == a)
    n_regions <- length(unique(hd$region_tag))
    p_heat <- ggplot(hd, aes(x = wave_label, y = region_tag, fill = mean_beta)) +
      geom_tile() +
      facet_wrap(~sex_label) +
      sunset_diverging(name = "Mean\nbeta") +
      labs(
        title = paste0(
          "Mean N-back ", clabel, " beta by region and wave - ",
          ifelse(a == "aseg", "subcortical (aseg)", "cortical (Desikan)")
        ),
        x = NULL,
        y = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))
    ggsave(
      file.path(out_dir, paste0("region_wave_heatmap_", a, ".png")),
      p_heat,
      width = 8,
      height = max(4, 0.16 * n_regions),
      dpi = 180
    )
  }

  # ---------------------------------------------------------------------------
  # (6) Region-pair comparison boxplots: amygdala vs. rostral middle frontal,
  # amygdala vs. insula. Both regions are bilateral, so hemisphere is kept as
  # a facet (left/right) rather than averaged away.
  # ---------------------------------------------------------------------------
  region_pair_boxplot <- function(
    label_a,
    label_b,
    filename,
    data = beta_long,
    x_col = "wave_label",
    subtitle = "Distribution across all people, by wave"
  ) {
    d <- data %>% filter(region_label %in% c(label_a, label_b)) %>%
      mutate(region_label = factor(region_label, levels = c(label_a, label_b)))

    if (nrow(d) == 0) {
      cat("\nSkipping", filename, "-- no rows for", label_a, "/", label_b, "\n")
      return(invisible(NULL))
    }

    pal_pair <- setNames(pal_two, c(label_a, label_b))
    linetype_pair <- setNames(pal_linetypes_two, c(label_a, label_b))

    p <- ggplot(
      d,
      aes(x = .data[[x_col]], y = beta, colour = region_label, linetype = region_label)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_boxplot(outlier.alpha = 0.05, outlier.size = 0.3, linewidth = 0.5) +
      facet_wrap(~hemi, labeller = as_labeller(c(left = "Left hemisphere", right = "Right hemisphere"))) +
      scale_colour_manual(values = pal_pair) +
      scale_linetype_manual(values = linetype_pair) +
      coord_cartesian(ylim = tight_ylim(d$beta)) +
      labs(
        title = paste0("N-back ", clabel, " beta: ", label_a, " vs. ", label_b),
        subtitle = subtitle,
        x = NULL,
        y = "Beta",
        colour = "Region",
        linetype = "Region"
      ) +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
    ggsave(file.path(out_dir, filename), p, width = 9, height = 5, dpi = 180)
  }

  region_pair_boxplot("amygdala", "rostral middle frontal", "amygdala_vs_rostral_middle_frontal_boxplot.png")
  region_pair_boxplot("amygdala", "insula", "amygdala_vs_insula_boxplot.png")

  if (!is.null(beta_long_stage)) {
    stage_subtitle <- "Distribution across all people, by pubertal status (pds_categ, youth-preferred)"
    region_pair_boxplot(
      "amygdala", "rostral middle frontal",
      "amygdala_vs_rostral_middle_frontal_by_pubertal_status_boxplot.png",
      data = beta_long_stage, x_col = "stage_label", subtitle = stage_subtitle
    )
    region_pair_boxplot(
      "amygdala", "insula",
      "amygdala_vs_insula_by_pubertal_status_boxplot.png",
      data = beta_long_stage, x_col = "stage_label", subtitle = stage_subtitle
    )
  } else {
    cat("\nSkipping pubertal-status region-pair plots -- no pubertal status match available\n")
  }

  cat("\nOutputs written to:", out_dir, "\n")
}

for (cst in names(fmri_contrasts)) {
  run_descriptives_for_contrast(cst)
}
