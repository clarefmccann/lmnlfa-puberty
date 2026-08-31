## viz_sample_timeline.R
## Descriptive visualization of data availability across the whole analytic
## sample (both sexes, all 11,860 participants in all_long.csv) -- NOT the
## HPC-run subsample. Shows, per participant, which waves they have any
## valid pubertal-item data for (parent OR youth report), plotted against
## their actual age at each occasion.
##
## "Valid" at a wave = at least one reporter (parent or youth) has all four
## PDS items (peta-petd) non-missing and in 1-4 -- matches the completeness
## criterion used to build the analytic samples in lmnlfa_growth_informant.R
## and friends, so this plot reflects the data actually available to those
## models rather than raw row presence (which includes partially-missing
## reporter rows).
##
## Usage: Rscript viz_sample_timeline.R
## Output: outputs/sample_timeline/timeline_<all|female|male>.png + .csv

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
out_dir <- file.path(out_base, "sample_timeline")
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

df <- read.csv(file.path(data_dir, "all_long.csv"))

items <- c("peta", "petb", "petc", "petd")
wave_order <- wave_codes
wave_labels <- wave_display_labels
# re-key the shared wave palette/shapes (bl, fu1, ...) to the display labels
# ("Baseline", "Year 1", ...) used by wave_label below
wave_palette <- setNames(pal_waves, wave_display_labels[names(pal_waves)])
wave_shapes <- setNames(pal_shapes_waves, wave_display_labels[names(pal_shapes_waves)])

# one row per (id, wave, reporter) that clears the same completeness bar
# used when building the analytic samples
valid_rows <- df %>%
  filter(!is.na(age)) %>%
  filter(if_all(all_of(items), ~ !is.na(.) & as.integer(.) %in% 1:4))

# collapse reporters: a wave counts as "has data" if EITHER reporter is valid
wave_avail <- valid_rows %>%
  distinct(id, wave, sex, .keep_all = TRUE) %>%
  group_by(id, wave) %>%
  summarise(age = mean(age), sex = first(sex), .groups = "drop") %>%
  mutate(
    wave = factor(wave, levels = wave_order),
    sex_label = case_match(sex, 1 ~ "Female", 2 ~ "Male", .default = NA_character_)
  ) %>%
  filter(!is.na(wave))

# baseline age = age at the participant's earliest valid occasion (not
# necessarily the "bl" wave itself, if that one happens to be missing)
baseline_age <- wave_avail %>%
  group_by(id) %>%
  summarise(baseline_age = min(age), n_waves = n(), sex_label = first(sex_label))

plot_df <- wave_avail %>%
  left_join(baseline_age %>% select(id, baseline_age, n_waves), by = "id") %>%
  group_by(sex_label) %>%
  mutate(id_rank = dense_rank(baseline_age)) %>%
  ungroup() %>%
  mutate(
    wave_label = factor(wave_labels[as.character(wave)], levels = wave_labels)
  ) %>%
  arrange(id, wave) %>%
  group_by(id) %>%
  mutate(age_prev = lag(age)) %>%
  ungroup()

cat("n participants with any valid data:", nrow(baseline_age), "\n")
cat("n waves per person -- summary:\n")
print(summary(baseline_age$n_waves))
cat("\nby sex:\n")
print(table(baseline_age$sex_label, useNA = "always"))

write.csv(
  baseline_age,
  file.path(out_dir, "participant_wave_counts.csv"),
  row.names = FALSE
)

make_timeline_plot <- function(pdat, title_suffix) {
  ggplot(pdat, aes(x = age, y = id_rank)) +
    geom_segment(
      data = pdat %>% filter(!is.na(age_prev)),
      aes(x = age_prev, xend = age, y = id_rank, yend = id_rank, colour = wave_label),
      alpha = 0.08,
      linewidth = 0.2
    ) +
    geom_point(
      aes(colour = wave_label, shape = wave_label),
      alpha = 0.5,
      size = 0.7,
      stroke = 0.4
    ) +
    scale_colour_manual(
      values = wave_palette,
      breaks = as.character(wave_labels),
      name = "Wave"
    ) +
    scale_shape_manual(
      values = wave_shapes,
      breaks = as.character(wave_labels),
      name = "Wave"
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
    labs(
      title = paste0("Data availability across the sample", title_suffix),
      subtitle = paste0(
        "Rows = participants, sorted by age at first valid occasion\n",
        "Solid shapes = MRI waves (baseline, yr 2/4/6); open shapes = annual-only (yr 1/3/5)"
      ),
      x = "Age (years)",
      y = "Participants, ranked by baseline age"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}

p_all <- make_timeline_plot(plot_df, " (both sexes)") +
  facet_wrap(~sex_label, ncol = 2)
ggsave(
  file.path(out_dir, "timeline_all.png"),
  p_all,
  width = 12,
  height = 7,
  dpi = 180
)

for (sx_lab in c("Female", "Male")) {
  p_sx <- make_timeline_plot(
    plot_df %>% filter(sex_label == sx_lab),
    paste0(" - ", sx_lab)
  )
  ggsave(
    file.path(out_dir, paste0("timeline_", tolower(sx_lab), ".png")),
    p_sx,
    width = 8,
    height = 7,
    dpi = 180
  )
}

cat("\nOutputs written to:", out_dir, "\n")
