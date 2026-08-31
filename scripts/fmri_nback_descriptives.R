## fmri_nback_descriptives.R
## Descriptive summaries and plots of the cleaned N-Back face-vs-place beta
## data (fmri_nback_foundation.R): how betas change by wave, differ by sex,
## and vary across regions. Purely descriptive -- no modeling.
##
## Usage: Rscript fmri_nback_descriptives.R
## Output: outputs/fmri_nback_descriptives/*.png, *.csv

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
out_dir <- file.path(out_base, "fmri_nback_descriptives")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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

# ---------------------------------------------------------------------------
# LOAD
# ---------------------------------------------------------------------------
fmri <- read.csv(file.path(data_dir, "all_fmri_nback_long.csv"))
region_lookup <- read.csv(file.path(data_dir, "fmri_nback_region_lookup.csv"))

region_cols <- region_lookup$region_key
region_cols <- region_cols[region_cols %in% names(fmri)]
if (length(region_cols) == 0) {
  stop("No region columns from the lookup table found in all_fmri_nback_long.csv")
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

atlas_labeller <- as_labeller(c(aseg = "Subcortical (aseg)", dsk = "Cortical (Desikan)"))

# re-key the shared wave palette/shapes (bl, fu1, ...) to the display labels
# ("Baseline", "Year 1", ...) used by wave_label
wave_palette_display <- setNames(pal_waves, wave_display_labels[names(pal_waves)])

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
  labs(
    title = "N-back face-vs-place beta by wave",
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
  labs(
    title = "N-back face-vs-place beta by wave and sex",
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
    title = "Mean N-back face-vs-place beta across waves, by sex",
    subtitle = "Averaged across all 98 regions and all people; mean +/- 90% CI",
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
        "Mean N-back face-vs-place beta by region and wave - ",
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

cat("\nOutputs written to:", out_dir, "\n")
