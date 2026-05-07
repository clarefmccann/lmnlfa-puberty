## 01_regression_characterization.R
## Characterizes participants who show item-level regression (score decrease
## across consecutive annual waves) in parent and youth PDS reports.
## Accounts for caregiver switching between waves.
## Requires outputs from 00_data_foundation.R.

library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(broom.mixed)

set.seed(90025)

# export DATA_DIR="/u/project/silvers/data/ABCD/ABCD-release-6.0/cfm/physical-health/puberty"
# export OUT_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"
# Rscript 01_regression_characteristics.R

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
out_dir <- file.path(out_base, "regression")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# LOAD
# ---------------------------------------------------------------------------
female_parent <- read.csv(file.path(pub_root, "female_parent_long.csv"))
female_youth <- read.csv(file.path(pub_root, "female_youth_long.csv"))
male_parent <- read.csv(file.path(pub_root, "male_parent_long.csv"))
male_youth <- read.csv(file.path(pub_root, "male_youth_long.csv"))

wave_order <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")
for (df_name in c(
  "female_parent",
  "female_youth",
  "male_parent",
  "male_youth"
)) {
  d <- get(df_name)
  d$wave <- factor(d$wave, levels = wave_order)
  assign(df_name, d)
}

# ---------------------------------------------------------------------------
# SECTION 0: RAW TRAJECTORY CHARACTERIZATION
# Descriptives + plots for individual items, PDS composite, categorical stage.
# ---------------------------------------------------------------------------
char_dir <- file.path(out_base, "characterization")
dir.create(char_dir, showWarnings = FALSE, recursive = TRUE)

# mpete is ordinal 1-4 (genital development), same scale as peta-petd.
# Only fpete (menarche) is binary (1 = no, 2 = yes after foundation recode).
ord_f <- c("peta", "petb", "petc", "petd")
ord_m <- c("peta", "petb", "petc", "petd", "mpete")

group_labels <- c(
  female_parent = "Female · Parent",
  female_youth = "Female · Youth",
  male_parent = "Male · Parent",
  male_youth = "Male · Youth"
)

datasets_char <- list(
  female_parent = list(df = female_parent, ord = ord_f, bin = "fpete"),
  female_youth = list(df = female_youth, ord = ord_f, bin = "fpete"),
  male_parent = list(df = male_parent, ord = ord_m, bin = NULL),
  male_youth = list(df = male_youth, ord = ord_m, bin = NULL)
)

as_wave_fct <- function(x) factor(x, levels = wave_order)

# ---- 0a. Ordinal item descriptives by wave ----------------------------------
describe_ordinal <- function(df, items, label) {
  df %>%
    select(wave, all_of(items)) %>%
    pivot_longer(all_of(items), names_to = "item", values_to = "value") %>%
    filter(!is.na(wave)) %>%
    group_by(wave, item) %>%
    summarise(
      n_valid = sum(!is.na(value)),
      n_total = n(),
      pct_miss = round(100 * mean(is.na(value)), 1),
      mean = round(mean(value, na.rm = TRUE), 3),
      sd = round(sd(value, na.rm = TRUE), 3),
      median = median(value, na.rm = TRUE),
      p25 = quantile(value, 0.25, na.rm = TRUE),
      p75 = quantile(value, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(group = label)
}

item_desc <- bind_rows(lapply(names(datasets_char), function(nm) {
  d <- datasets_char[[nm]]
  describe_ordinal(d$df, d$ord, nm)
}))

write.csv(
  item_desc,
  file.path(char_dir, "item_descriptives_by_wave.csv"),
  row.names = FALSE
)
cat("\n=== Ordinal item descriptives by wave ===\n")
print(item_desc, n = Inf)

# ---- 0b. Binary item (fpete: menarche) proportion-yes by wave ---------------
# Males have no binary item; only female datasets have fpete.
describe_binary <- function(df, bin_item, label) {
  if (is.null(bin_item) || !bin_item %in% names(df)) {
    return(NULL)
  }
  df %>%
    select(wave, value = all_of(bin_item)) %>%
    filter(!is.na(wave)) %>%
    group_by(wave) %>%
    summarise(
      n_valid = sum(!is.na(value)),
      n_total = n(),
      pct_miss = round(100 * mean(is.na(value)), 1),
      pct_yes = round(100 * mean(value == 2, na.rm = TRUE), 1),
      .groups = "drop"
    ) %>%
    mutate(item = bin_item, group = label)
}

binary_desc <- bind_rows(lapply(names(datasets_char), function(nm) {
  d <- datasets_char[[nm]]
  describe_binary(d$df, d$bin, nm)
}))

if (nrow(binary_desc) > 0) {
  write.csv(
    binary_desc,
    file.path(char_dir, "binary_item_descriptives_by_wave.csv"),
    row.names = FALSE
  )
  cat("\n=== Menarche (fpete) % yes by wave ===\n")
  print(binary_desc, n = Inf)
}

# ---- 0c. PDS composite descriptives by wave ---------------------------------
describe_composite <- function(df, label) {
  if (!"pds_comp" %in% names(df)) {
    return(NULL)
  }
  df %>%
    filter(!is.na(wave), !is.na(pds_comp)) %>%
    group_by(wave) %>%
    summarise(
      n = n(),
      mean = round(mean(pds_comp, na.rm = TRUE), 3),
      sd = round(sd(pds_comp, na.rm = TRUE), 3),
      median = round(median(pds_comp, na.rm = TRUE), 2),
      p25 = round(quantile(pds_comp, 0.25, na.rm = TRUE), 2),
      p75 = round(quantile(pds_comp, 0.75, na.rm = TRUE), 2),
      .groups = "drop"
    ) %>%
    mutate(group = label)
}

comp_desc <- bind_rows(lapply(names(datasets_char), function(nm) {
  describe_composite(datasets_char[[nm]]$df, nm)
}))

write.csv(
  comp_desc,
  file.path(char_dir, "composite_descriptives_by_wave.csv"),
  row.names = FALSE
)
cat("\n=== PDS composite descriptives by wave ===\n")
print(comp_desc, n = Inf)

# ---- 0d. Categorical stage frequencies by wave ------------------------------
describe_categ <- function(df, label) {
  if (!"pds_categ" %in% names(df)) {
    return(NULL)
  }
  df %>%
    filter(!is.na(wave), !is.na(pds_categ)) %>%
    mutate(
      stage = factor(
        as.integer(pds_categ),
        levels = 1:5,
        labels = c("prepubertal", "early", "mid", "late", "postpubertal")
      )
    ) %>%
    count(wave, stage, .drop = FALSE) %>%
    group_by(wave) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup() %>%
    mutate(group = label)
}

categ_desc <- bind_rows(lapply(names(datasets_char), function(nm) {
  describe_categ(datasets_char[[nm]]$df, nm)
}))

if (nrow(categ_desc) > 0) {
  write.csv(
    categ_desc,
    file.path(char_dir, "categ_stage_by_wave.csv"),
    row.names = FALSE
  )
  cat("\n=== Categorical stage distribution by wave ===\n")
  print(categ_desc, n = Inf)
}

# ---- 0e. PLOTS --------------------------------------------------------------

# 1. Ordinal item mean trajectories — ribbon = ±1 SD
p_item_traj <- ggplot(
  item_desc %>%
    mutate(
      wave = as_wave_fct(wave),
      group_label = factor(
        dplyr::recode(group, !!!group_labels),
        levels = group_labels
      )
    ),
  aes(x = wave, y = mean, colour = item, group = item)
) +
  geom_ribbon(
    aes(ymin = mean - sd, ymax = mean + sd, fill = item),
    alpha = 0.10,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~group_label, nrow = 1) +
  scale_y_continuous(limits = c(1, 4), breaks = 1:4) +
  labs(
    title = "PDS ordinal item mean trajectories across waves",
    subtitle = "Ribbon = ±1 SD",
    x = "Wave",
    y = "Mean score (1–4)",
    colour = "Item",
    fill = "Item"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 9)
  )

ggsave(
  file.path(char_dir, "item_mean_trajectories.png"),
  p_item_traj,
  width = 14,
  height = 5.5,
  dpi = 150
)

# 2. Menarche proportion-yes trajectory (females only)
if (nrow(binary_desc) > 0) {
  p_binary <- ggplot(
    binary_desc %>%
      mutate(
        wave = as_wave_fct(wave),
        group_label = factor(
          dplyr::recode(group, !!!group_labels),
          levels = group_labels
        )
      ),
    aes(x = wave, y = pct_yes, colour = group_label, group = group_label)
  ) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 2.5) +
    scale_y_continuous(limits = c(0, 100), labels = function(x) {
      paste0(x, "%")
    }) +
    labs(
      title = "% reporting menarche (fpete = yes) by wave",
      x = "Wave",
      y = "% yes",
      colour = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom")

  ggsave(
    file.path(char_dir, "menarche_trajectory.png"),
    p_binary,
    width = 8,
    height = 5,
    dpi = 150
  )
}

# 3. PDS composite mean trajectory — ribbon = ±1 SE
p_comp_traj <- ggplot(
  comp_desc %>%
    mutate(
      wave = as_wave_fct(wave),
      group_label = factor(
        dplyr::recode(group, !!!group_labels),
        levels = group_labels
      )
    ),
  aes(x = wave, y = mean, colour = group_label, group = group_label)
) +
  geom_ribbon(
    aes(
      ymin = mean - sd / sqrt(n),
      ymax = mean + sd / sqrt(n),
      fill = group_label
    ),
    alpha = 0.15,
    colour = NA
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_y_continuous(limits = c(1, 4)) +
  labs(
    title = "PDS composite mean trajectory across waves",
    subtitle = "Ribbon = ±1 SE",
    x = "Wave",
    y = "Mean PDS composite (1–4)",
    colour = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(
  file.path(char_dir, "composite_trajectory.png"),
  p_comp_traj,
  width = 9,
  height = 5,
  dpi = 150
)

# 4. Categorical stage distribution — stacked bar
if (nrow(categ_desc) > 0) {
  stage_cols <- c(
    prepubertal = "#d0e9f7",
    early = "#74b9d4",
    mid = "#2980b9",
    late = "#1a5276",
    postpubertal = "#0b2545"
  )

  p_categ <- ggplot(
    categ_desc %>%
      filter(!is.na(stage)) %>%
      mutate(
        wave = as_wave_fct(wave),
        group_label = factor(
          dplyr::recode(group, !!!group_labels),
          levels = group_labels
        ),
        stage = factor(stage, levels = names(stage_cols))
      ),
    aes(x = wave, y = pct, fill = stage)
  ) +
    geom_bar(stat = "identity", width = 0.85) +
    scale_fill_manual(values = stage_cols) +
    facet_wrap(~group_label, nrow = 1) +
    labs(
      title = "ABCD categorical pubertal stage distribution across waves",
      x = "Wave",
      y = "% of participants",
      fill = "Stage"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold", size = 11),
      axis.text.x = element_text(size = 9)
    )

  ggsave(
    file.path(char_dir, "categ_stage_distribution.png"),
    p_categ,
    width = 14,
    height = 5.5,
    dpi = 150
  )
}

cat("\nCharacterization outputs written to:", char_dir, "\n")

# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

# Items to check for regression (ordinal, 1-4 scale)
ordinal_items_f <- c("peta", "petb", "petc", "petd")
ordinal_items_m <- c("peta", "petb", "petc", "petd", "mpete")

# Only fpete (menarche) is binary — mpete is ordinal 1-4, handled with peta-petd
binary_item_f <- "fpete"

flag_regression <- function(df, ordinal_items, binary_item = NULL) {
  df <- df %>%
    arrange(id, wave) %>%
    group_by(id) %>%
    mutate(
      across(
        all_of(ordinal_items),
        ~ . < lag(.),
        .names = "reg_{col}"
      )
    )

  if (!is.null(binary_item) && binary_item %in% names(df)) {
    # menarche: flag any decrease (should be impossible biologically)
    df <- df %>%
      mutate(
        !!paste0("reg_", binary_item) := get(binary_item) <
          lag(get(binary_item))
      )
  }

  reg_cols <- grep("^reg_", names(df), value = TRUE)

  df <- df %>%
    mutate(
      # at least one item regressed this wave transition
      any_regression = if_any(all_of(reg_cols), ~ . == TRUE & !is.na(.)),
      # count of items that regressed
      n_items_regressed = rowSums(
        across(all_of(reg_cols), ~ . == TRUE & !is.na(.)),
        na.rm = TRUE
      ),
      # flag if regression co-occurred with a caregiver switch
      regression_with_cg_switch = any_regression & caregiver_switched
    ) %>%
    ungroup()

  df
}

# ---------------------------------------------------------------------------
# FLAG REGRESSIONS
# ---------------------------------------------------------------------------
female_parent <- flag_regression(female_parent, ordinal_items_f, binary_item_f)
female_youth <- flag_regression(female_youth, ordinal_items_f, binary_item_f)
male_parent <- flag_regression(male_parent, ordinal_items_m, NULL)
male_youth <- flag_regression(male_youth, ordinal_items_m, NULL)

# ---------------------------------------------------------------------------
# SUMMARY TABLES
# ---------------------------------------------------------------------------

summarise_regression <- function(df, label) {
  reg_cols <- grep("^reg_", names(df), value = TRUE)

  # per-item regression rates across all wave transitions
  item_rates <- df %>%
    filter(!is.na(wave)) %>%
    summarise(
      across(
        all_of(reg_cols),
        ~ mean(. == TRUE, na.rm = TRUE),
        .names = "rate_{col}"
      )
    ) %>%
    pivot_longer(
      everything(),
      names_to = "item",
      values_to = "regression_rate"
    ) %>%
    mutate(
      item = sub("rate_reg_", "", item),
      group = label
    )

  # per-wave regression rates
  wave_rates <- df %>%
    group_by(wave) %>%
    summarise(
      n_transitions = sum(!is.na(any_regression)),
      pct_any_regression = mean(any_regression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(group = label)

  # caregiver switch × regression cross-tab
  cg_cross <- df %>%
    filter(!is.na(caregiver_switched)) %>%
    group_by(caregiver_switched) %>%
    summarise(
      pct_any_regression = mean(any_regression, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(group = label)

  list(item_rates = item_rates, wave_rates = wave_rates, cg_cross = cg_cross)
}

res_fp <- summarise_regression(female_parent, "female_parent")
res_fy <- summarise_regression(female_youth, "female_youth")
res_mp <- summarise_regression(male_parent, "male_parent")
res_my <- summarise_regression(male_youth, "male_youth")

item_rate_table <- bind_rows(
  res_fp$item_rates,
  res_fy$item_rates,
  res_mp$item_rates,
  res_my$item_rates
)

wave_rate_table <- bind_rows(
  res_fp$wave_rates,
  res_fy$wave_rates,
  res_mp$wave_rates,
  res_my$wave_rates
)

cg_cross_table <- bind_rows(
  res_fp$cg_cross,
  res_fy$cg_cross,
  res_mp$cg_cross,
  res_my$cg_cross
)

cat("=== Item-level regression rates ===\n")
print(item_rate_table, n = Inf)
cat("\n=== Wave-level regression rates ===\n")
print(wave_rate_table, n = Inf)
cat("\n=== Regression rate by caregiver switch ===\n")
print(cg_cross_table, n = Inf)

# ---------------------------------------------------------------------------
# PARTICIPANT-LEVEL REGRESSION SUMMARY
# ---------------------------------------------------------------------------
# how many participants regress at least once (ever-regressen), and on how
# many items / waves?

participant_summary <- function(df, label) {
  df %>%
    group_by(id) %>%
    summarise(
      ever_regressed = any(any_regression, na.rm = TRUE),
      total_regressions = sum(n_items_regressed, na.rm = TRUE),
      n_waves_w_regression = sum(any_regression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(group = label)
}

part_summary <- bind_rows(
  participant_summary(female_parent, "female_parent"),
  participant_summary(female_youth, "female_youth"),
  participant_summary(male_parent, "male_parent"),
  participant_summary(male_youth, "male_youth")
)

cat("\n=== Participant-level regression prevalence ===\n")
part_summary %>%
  group_by(group) %>%
  summarise(
    n_participants = n(),
    pct_ever_regressed = mean(ever_regressed, na.rm = TRUE),
    mean_total_reg = mean(total_regressions, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# ---------------------------------------------------------------------------
# PREDICTORS OF REGRESSION
# Mixed logistic: any_regression ~ age + bmi_z + race + caregiver_switched
# + (1 | id), run separately per sex × reporter
# ---------------------------------------------------------------------------

run_regression_model <- function(df, label) {
  model_data <- df %>%
    filter(!is.na(any_regression), !is.na(age), !is.na(bmi_z), !is.na(race)) %>%
    mutate(
      any_regression = as.integer(any_regression),
      caregiver_switched = as.integer(coalesce(caregiver_switched, FALSE)),
      race = factor(race)
    )

  if (nrow(model_data) < 100) {
    message("Skipping ", label, ": insufficient data after NA removal.")
    return(NULL)
  }

  fit <- tryCatch(
    glmer(
      any_regression ~ age + bmi_z + race + caregiver_switched + (1 | id),
      data = model_data,
      family = binomial,
      control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    ),
    error = function(e) {
      message(label, " model failed: ", e$message)
      NULL
    }
  )

  if (is.null(fit)) {
    return(NULL)
  }

  tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    mutate(group = label)
}

cat("\n=== Mixed logistic models: predictors of item-level regression ===\n")
reg_models <- bind_rows(
  run_regression_model(female_parent, "female_parent"),
  run_regression_model(female_youth, "female_youth"),
  run_regression_model(male_parent, "male_parent"),
  run_regression_model(male_youth, "male_youth")
)
print(reg_models, n = Inf)

# ---------------------------------------------------------------------------
# CAREGIVER SWITCH CHARACTERIZATION
# ---------------------------------------------------------------------------
cg_summary <- bind_rows(
  female_parent,
  female_youth,
  male_parent,
  male_youth
) %>%
  distinct(id, sex, reporter, n_cg_switches, ever_switched_cg) %>%
  group_by(sex, reporter) %>%
  summarise(
    n = n(),
    pct_ever_switched = mean(ever_switched_cg, na.rm = TRUE),
    mean_switches = mean(n_cg_switches, na.rm = TRUE),
    pct_switched_once = mean(n_cg_switches == 1, na.rm = TRUE),
    pct_switched_2plus = mean(n_cg_switches >= 2, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Caregiver switching summary ===\n")
print(cg_summary)

# ---------------------------------------------------------------------------
# PLOTS
# ---------------------------------------------------------------------------

# item regression rates by group
p_item <- ggplot(
  item_rate_table,
  aes(x = item, y = regression_rate, fill = group)
) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Item-level regression rates across wave transitions",
    x = "PDS item",
    y = "% wave transitions showing regression",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(
  file.path(out_dir, "item_regression_rates.png"),
  p_item,
  width = 10,
  height = 5,
  dpi = 150
)

# wave-level regression rates
p_wave <- ggplot(
  wave_rate_table,
  aes(x = wave, y = pct_any_regression, colour = group, group = group)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Regression rate by wave transition",
    x = "Wave (transition from prior wave)",
    y = "% with any item regression",
    colour = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(
  file.path(out_dir, "wave_regression_rates.png"),
  p_wave,
  width = 9,
  height = 5,
  dpi = 150
)

# caregiver switch × regression
p_cg <- ggplot(
  cg_cross_table,
  aes(
    x = factor(caregiver_switched, labels = c("No switch", "Switch")),
    y = pct_any_regression,
    fill = group
  )
) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Regression rate by caregiver switch",
    x = NULL,
    y = "% wave transitions with regression",
    fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(
  file.path(out_dir, "caregiver_switch_regression.png"),
  p_cg,
  width = 8,
  height = 5,
  dpi = 150
)

# ---------------------------------------------------------------------------
# SAVE TABLES
# ---------------------------------------------------------------------------
write.csv(
  item_rate_table,
  file.path(out_dir, "reg_item_rates.csv"),
  row.names = FALSE
)
write.csv(
  wave_rate_table,
  file.path(out_dir, "reg_wave_rates.csv"),
  row.names = FALSE
)
write.csv(
  cg_cross_table,
  file.path(out_dir, "reg_cg_cross.csv"),
  row.names = FALSE
)
write.csv(
  part_summary,
  file.path(out_dir, "reg_participant_summary.csv"),
  row.names = FALSE
)
write.csv(
  reg_models,
  file.path(out_dir, "reg_predictors_models.csv"),
  row.names = FALSE
)
write.csv(
  cg_summary,
  file.path(out_dir, "caregiver_switch_summary.csv"),
  row.names = FALSE
)

cat("\nPlots and tables written to:", out_dir, "\n")
