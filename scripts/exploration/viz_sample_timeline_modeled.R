## viz_sample_timeline_modeled.R
## Same wave-availability timeline as viz_sample_timeline.R (per-participant
## rows, ranked by age at first valid occasion, points/segments colored by
## wave) -- but restricted to the EXACT n=1500-per-sex sample the staged
## sigmoid Stan models were actually fit on, not the full ~11,860-person
## ABCD analytic sample.
##
## Sample reconstruction matches lmnlfa_growth_sigmoid_staged.R's own
## build_lmnlfa_data_staged() filtering exactly (age/race/WHtR non-missing,
## WHtR in [.25,.75], all 4 PDS items non-missing for a given reporter row)
## plus the identical seed (90025) + subsample size (1500), so this is
## genuinely the same people, not merely a similarly-sized draw.
##
## "Valid" at a wave (for the purposes of THIS plot) collapses reporters --
## a wave counts as available if EITHER parent OR youth report is valid --
## matching viz_sample_timeline.R's own definition, applied within the
## modeled person-set.
##
## Usage: Rscript viz_sample_timeline_modeled.R
## Output: outputs/sample_timeline/timeline_modeled_<all|female|male>.png

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

data_dir <- Sys.getenv("DATA_DIR")
sshfs_data_dir <- "/private/tmp/sshfs/projects/abcd-projs/dissertation/lmnlfa-puberty/data"
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  if (dir.exists(sshfs_data_dir)) {
    data_dir <- sshfs_data_dir
  } else {
    root_path <- Sys.getenv("HOME_DIR")
    if (!nzchar(root_path)) root_path <- Sys.getenv("HOME")
    data_dir <- file.path(root_path, "projects/abcd-projs/dissertation/lmnlfa-puberty/data")
  }
}
if (!dir.exists(data_dir)) stop("Cannot locate data directory: ", data_dir)

out_base <- Sys.getenv("OUT_DIR")
sshfs_out_base <- "/private/tmp/sshfs/projects/abcd-projs/dissertation/lmnlfa-puberty/outputs"
if (!nzchar(out_base) || !dir.exists(out_base)) {
  if (dir.exists(sshfs_out_base)) {
    out_base <- sshfs_out_base
  } else {
    root_path <- Sys.getenv("HOME_DIR")
    if (!nzchar(root_path)) root_path <- Sys.getenv("HOME")
    out_base <- file.path(root_path, "projects/abcd-projs/dissertation/lmnlfa-puberty/outputs")
  }
}
out_dir <- file.path(out_base, "exploration", "sample_timeline")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

script_dir <- Sys.getenv("SGE_O_WORKDIR")
if (!nzchar(script_dir)) {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", cmd_args[grep("^--file=", cmd_args)])
  script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg[1])) else "scripts"
}
palette_file <- file.path(script_dir, "color_palette.R")
if (!file.exists(palette_file)) palette_file <- file.path("scripts", "color_palette.R")
source(palette_file)

ordinal_items <- c("peta", "petb", "petc", "petd")
wave_order <- wave_codes
wave_labels <- wave_display_labels
wave_palette <- setNames(pal_waves, wave_display_labels[names(pal_waves)])
wave_shapes <- setNames(pal_shapes_waves, wave_display_labels[names(pal_shapes_waves)])

whtr_min <- 0.25
whtr_max <- 0.75
n_subsample <- 1500

# ---------------------------------------------------------------------------
# reconstruct the EXACT modeled sample -- same filter, same seed/order as
# lmnlfa_growth_sigmoid_staged.R's build_lmnlfa_data_staged()
# ---------------------------------------------------------------------------
build_modeled_sample <- function(sx) {
  parent_df <- read.csv(file.path(data_dir, paste0(sx, "_parent_long.csv")))
  youth_df <- read.csv(file.path(data_dir, paste0(sx, "_youth_long.csv")))

  clean_reporter <- function(df) {
    df %>%
      select(id, wave, age, race, whtr, all_of(ordinal_items)) %>%
      filter(!is.na(age), !is.na(race), !is.na(whtr), whtr >= whtr_min, whtr <= whtr_max) %>%
      filter(if_all(all_of(ordinal_items), ~ !is.na(.) & as.integer(.) %in% 1:4))
  }

  dat <- bind_rows(clean_reporter(parent_df), clean_reporter(youth_df))
  all_ids <- sort(unique(dat$id))

  set.seed(90025) # matches lmnlfa_growth_sigmoid_staged.R's seed + build order exactly
  sub_ids <- sort(sample(all_ids, min(n_subsample, length(all_ids))))

  dat %>%
    filter(id %in% sub_ids) %>%
    mutate(wave = factor(wave, levels = wave_order)) %>%
    filter(!is.na(wave)) %>%
    distinct(id, wave, .keep_all = TRUE) %>% # collapse reporters: either valid -> wave available
    select(id, wave, age) %>%
    mutate(sex_label = tools::toTitleCase(sx))
}

wave_avail <- bind_rows(lapply(c("female", "male"), build_modeled_sample))

cat("n participants in modeled samples (post wave-collapse):\n")
print(table(wave_avail$sex_label[!duplicated(paste(wave_avail$sex_label, wave_avail$id))]))

baseline_age <- wave_avail %>%
  group_by(id) %>%
  summarise(baseline_age = min(age), n_waves = n(), sex_label = first(sex_label))

cat("\nn waves per person -- summary:\n")
print(summary(baseline_age$n_waves))

write.csv(
  baseline_age,
  file.path(out_dir, "participant_wave_counts_modeled.csv"),
  row.names = FALSE
)

plot_df <- wave_avail %>%
  left_join(baseline_age %>% select(id, baseline_age, n_waves), by = "id") %>%
  group_by(sex_label) %>%
  mutate(id_rank = dense_rank(baseline_age)) %>%
  ungroup() %>%
  mutate(wave_label = factor(wave_labels[as.character(wave)], levels = wave_labels)) %>%
  arrange(id, wave) %>%
  group_by(id) %>%
  mutate(age_prev = lag(age)) %>%
  ungroup()

make_timeline_plot <- function(pdat, title_suffix) {
  ggplot(pdat, aes(x = age, y = id_rank)) +
    geom_segment(
      data = pdat %>% filter(!is.na(age_prev)),
      aes(x = age_prev, xend = age, y = id_rank, yend = id_rank, colour = wave_label),
      alpha = 0.1,
      linewidth = 0.35
    ) +
    geom_point(
      aes(colour = wave_label, shape = wave_label),
      alpha = 0.5,
      size = 0.9,
      stroke = 0.5
    ) +
    scale_colour_manual(values = wave_palette, breaks = as.character(wave_labels), name = "Wave") +
    scale_shape_manual(values = wave_shapes, breaks = as.character(wave_labels), name = "Wave") +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3, linewidth = 2))) +
    labs(
      title = paste0("Data availability, staged-model sample (n=", n_subsample, "/sex)", title_suffix),
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
ggsave(file.path(out_dir, "timeline_modeled_all.png"), p_all, width = 12, height = 7, dpi = 180)

for (sx_lab in c("Female", "Male")) {
  p_sx <- make_timeline_plot(plot_df %>% filter(sex_label == sx_lab), paste0(" - ", sx_lab))
  ggsave(
    file.path(out_dir, paste0("timeline_modeled_", tolower(sx_lab), ".png")),
    p_sx,
    width = 8,
    height = 7,
    dpi = 180
  )
}

cat("\nOutputs written to:", out_dir, "\n")
