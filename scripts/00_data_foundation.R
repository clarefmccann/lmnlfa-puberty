## 00_data_foundation.R
## Loads raw ABCD PDS data (annual waves only), attaches covariates
## (age, sex, BMI, race/ethnicity, caregiver), and writes clean
## long-format datasets for each sex × reporter combination.

# export DATA_DIR="/u/project/silvers/data/ABCD/ABCD-release-6.0/cfm/physical-health/puberty"
# export OUT_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"
# Rscript 01_data_foundation.R

library(dplyr)
library(tidyr)
library(NBDCtools)

set.seed(90025)

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
root_path <- Sys.getenv("HOME_DIR")
if (!nzchar(root_path)) {
  root_path <- Sys.getenv("HOME")
}

data_root <- file.path(
  root_path,
  "projects/abcd-projs/abcd-data-release-6.0/nbdc-tools-data"
)
if (!dir.exists(data_root)) {
  stop("Cannot locate nbdc-tools-data: ", data_root)
}

out_root <- file.path(
  root_path,
  "projects/abcd-projs/abcd-data-release-6.0/cfm/physical-health/puberty"
)

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------
vars <- c(
  "ab_g_dyn__visit_age",
  "ab_g_stc__cohort_sex",
  "ab_g_stc__cohort_ethnrace__mhisp",
  "ab_g_dyn__visit__day1_inform",
  "ab_g_dyn__design_site",
  "ph_y_anthr__height_mean", # inches
  "ph_y_anthr__weight_mean", # lbs
  # PDS items — parent report
  "ph_p_pds_001",
  "ph_p_pds_002",
  "ph_p_pds_003",
  "ph_p_pds__f_001",
  "ph_p_pds__f_002",
  "ph_p_pds__m_001",
  "ph_p_pds__m_002",
  # PDS items — youth report
  "ph_y_pds_001",
  "ph_y_pds_002",
  "ph_y_pds_003",
  "ph_y_pds__f_001",
  "ph_y_pds__f_002",
  "ph_y_pds__m_001",
  "ph_y_pds__m_002",
  # ABCD-calculated PDS composites
  "ph_p_pds__f_mean",
  "ph_p_pds__m_mean",
  "ph_y_pds__f_mean",
  "ph_y_pds__m_mean",
  # ABCD-calculated PDS categorical stage
  "ph_p_pds__f_categ",
  "ph_p_pds__m_categ",
  "ph_y_pds__f_categ",
  "ph_y_pds__m_categ"
)

# ---------------------------------------------------------------------------
# CHARACTERIZE "DON'T KNOW" (999) AND "REFUSED" (777) RESPONSES
# ---------------------------------------------------------------------------
pds_item_vars <- c(
  "ph_p_pds_001",
  "ph_p_pds_002",
  "ph_p_pds_003",
  "ph_p_pds__f_001",
  "ph_p_pds__f_002",
  "ph_p_pds__m_001",
  "ph_p_pds__m_002",
  "ph_y_pds_001",
  "ph_y_pds_002",
  "ph_y_pds_003",
  "ph_y_pds__f_001",
  "ph_y_pds__f_002",
  "ph_y_pds__m_001",
  "ph_y_pds__m_002"
)

raw_uncleaned <- create_dataset(
  dir_data = data_root,
  study = "abcd",
  vars = c("ab_g_stc__cohort_sex", pds_item_vars),
  value_to_na = FALSE,
  bind_shadow = FALSE
)

wave_labels <- c(
  "ses-00A" = "bl",
  "ses-01A" = "fu1",
  "ses-02A" = "fu2",
  "ses-03A" = "fu3",
  "ses-04A" = "fu4",
  "ses-05A" = "fu5",
  "ses-06A" = "fu6"
)

dk_long <- raw_uncleaned %>%
  filter(session_id %in% names(wave_labels)) %>%
  mutate(wave = recode(session_id, !!!wave_labels)) %>%
  rename(sex = ab_g_stc__cohort_sex) %>%
  select(participant_id, wave, sex, all_of(pds_item_vars)) %>%
  mutate(across(all_of(pds_item_vars), ~ as.numeric(as.character(.)))) %>%
  tidyr::pivot_longer(
    all_of(pds_item_vars),
    names_to = "item_raw",
    values_to = "value"
  ) %>%
  mutate(
    reporter = if_else(grepl("^ph_p", item_raw), "parent", "youth"),
    response_type = case_when(
      value == 777 ~ "refused",
      value == 999 ~ "dont_know",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(response_type))

dk_summary <- dk_long %>%
  count(item_raw, reporter, sex, wave, response_type) %>%
  arrange(item_raw, sex, reporter, wave, response_type)

dk_overall <- dk_long %>%
  count(item_raw, reporter, response_type) %>%
  tidyr::pivot_wider(
    names_from = response_type,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(item_raw, reporter)

out_root_dk <- file.path(
  root_path,
  "projects/abcd-projs/dissertation/study1/outputs/data_quality"
)
dir.create(out_root_dk, showWarnings = FALSE, recursive = TRUE)
write.csv(
  dk_summary,
  file.path(out_root_dk, "dk_refused_by_wave.csv"),
  row.names = FALSE
)
write.csv(
  dk_overall,
  file.path(out_root_dk, "dk_refused_overall.csv"),
  row.names = FALSE
)

cat("\n=== Don't-know / refused response summary ===\n")
print(dk_overall)
rm(raw_uncleaned, dk_long)

raw <- create_dataset(
  dir_data = data_root,
  study = "abcd",
  vars = vars,
  value_to_na = TRUE,
  bind_shadow = FALSE
)

# ---------------------------------------------------------------------------
# ANNUAL WAVES ONLY + WAVE LABEL
# ---------------------------------------------------------------------------

data <- raw %>%
  filter(session_id %in% names(wave_labels)) %>%
  mutate(
    wave = recode(session_id, !!!wave_labels),
    wave = factor(wave, levels = wave_labels)
  ) %>%
  rename(
    id = participant_id,
    age = ab_g_dyn__visit_age,
    sex = ab_g_stc__cohort_sex,
    race = ab_g_stc__cohort_ethnrace__mhisp,
    caregiver = ab_g_dyn__visit__day1_inform,
    site = ab_g_dyn__design_site,
    height_in = ph_y_anthr__height_mean,
    weight_lb = ph_y_anthr__weight_mean,
    peta_p = ph_p_pds_001,
    petb_p = ph_p_pds_002,
    petc_p = ph_p_pds_003,
    petdf_p = ph_p_pds__f_001,
    fpete_p = ph_p_pds__f_002,
    petdm_p = ph_p_pds__m_001,
    mpete_p = ph_p_pds__m_002,
    peta_y = ph_y_pds_001,
    petb_y = ph_y_pds_002,
    petc_y = ph_y_pds_003,
    petdf_y = ph_y_pds__f_001,
    fpete_y = ph_y_pds__f_002,
    petdm_y = ph_y_pds__m_001,
    mpete_y = ph_y_pds__m_002,
    abcd_comp_fp = ph_p_pds__f_mean,
    abcd_comp_mp = ph_p_pds__m_mean,
    abcd_comp_fy = ph_y_pds__f_mean,
    abcd_comp_my = ph_y_pds__m_mean,
    pds_categ_fp = ph_p_pds__f_categ,
    pds_categ_mp = ph_p_pds__m_categ,
    pds_categ_fy = ph_y_pds__f_categ,
    pds_categ_my = ph_y_pds__m_categ
  )

# ---------------------------------------------------------------------------
# DERIVED VARIABLES
# ---------------------------------------------------------------------------

# BMI
data <- data %>%
  mutate(bmi = 703 * weight_lb / height_in^2)

# age-standardized BMI z-score within sex
data <- data %>%
  group_by(sex, wave) %>%
  mutate(bmi_z = as.numeric(scale(bmi))) %>%
  ungroup()

# caregiver switch flags (relative to each participant's prior annual wave)
data <- data %>%
  arrange(id, wave) %>%
  group_by(id) %>%
  mutate(
    caregiver_switched = caregiver != lag(caregiver) &
      !is.na(caregiver) &
      !is.na(lag(caregiver)),
    n_cg_switches = sum(caregiver_switched, na.rm = TRUE),
    ever_switched_cg = n_cg_switches > 0
  ) %>%
  ungroup()

# ---------------------------------------------------------------------------
# RECODE INVALID VALUES TO NA
# (999 = don't know, 777 = refused in ABCD)
# ---------------------------------------------------------------------------
pds_items <- c(
  "peta_p",
  "petb_p",
  "petc_p",
  "petdf_p",
  "fpete_p",
  "petdm_p",
  "mpete_p",
  "peta_y",
  "petb_y",
  "petc_y",
  "petdf_y",
  "fpete_y",
  "petdm_y",
  "mpete_y"
)

data <- data %>%
  mutate(across(
    all_of(pds_items),
    ~ if_else(. %in% c(777, 999), NA_real_, as.numeric(.))
  ))

# fpete raw coding is 0 = no menarche, 1 = yes menarche.
# Shift to 1/2 so it sits on the same range as the ordinal items (1–4)
# and passes the %in% 1:2 validity filter in downstream scripts.
data <- data %>%
  mutate(
    fpete_p = if_else(!is.na(fpete_p), fpete_p + 1L, NA_real_),
    fpete_y = if_else(!is.na(fpete_y), fpete_y + 1L, NA_real_)
  )

# ---------------------------------------------------------------------------
# PDS COMPOSITES
# ---------------------------------------------------------------------------
use_abcd_comp <- function(abcd_col, fallback_cols, data) {
  if (all(is.na(data[[abcd_col]]))) {
    cat("  ABCD composite", abcd_col, "unavailable — using manual mean\n")
    rowMeans(data[fallback_cols], na.rm = FALSE)
  } else {
    data[[abcd_col]]
  }
}

data <- data %>%
  mutate(
    pds_comp_f_p = use_abcd_comp(
      "abcd_comp_fp",
      c("peta_p", "petb_p", "petc_p", "petdf_p"),
      .
    ),
    pds_comp_f_y = use_abcd_comp(
      "abcd_comp_fy",
      c("peta_y", "petb_y", "petc_y", "petdf_y"),
      .
    ),
    pds_comp_m_p = use_abcd_comp(
      "abcd_comp_mp",
      c("peta_p", "petb_p", "petc_p", "petdm_p"),
      .
    ),
    pds_comp_m_y = use_abcd_comp(
      "abcd_comp_my",
      c("peta_y", "petb_y", "petc_y", "petdm_y"),
      .
    )
  )

# ---------------------------------------------------------------------------
# CHARACTERIZE ABCD CATEGORICAL PUBERTAL STAGE
# Frequency of each stage (1–5) by sex × reporter × wave.
# Categories: 1=prepubertal, 2=early puberty, 3=mid puberty,
#             4=late puberty, 5=postpubertal
# ---------------------------------------------------------------------------
categ_long <- data %>%
  select(
    id,
    wave,
    sex,
    parent_female = pds_categ_fp,
    parent_male = pds_categ_mp,
    youth_female = pds_categ_fy,
    youth_male = pds_categ_my
  ) %>%
  tidyr::pivot_longer(
    c(parent_female, parent_male, youth_female, youth_male),
    names_to = "reporter_sex",
    values_to = "pds_categ"
  ) %>%
  tidyr::separate(reporter_sex, into = c("reporter", "item_sex"), sep = "_") %>%
  filter(!is.na(pds_categ)) %>%
  mutate(
    pds_categ = as.integer(pds_categ),
    stage_label = factor(
      pds_categ,
      levels = 1:5,
      labels = c("prepubertal", "early", "mid", "late", "postpubertal")
    )
  )

categ_freq <- categ_long %>%
  count(sex, item_sex, reporter, wave, stage_label) %>%
  group_by(sex, item_sex, reporter, wave) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(sex, item_sex, reporter, wave, stage_label)

categ_overall <- categ_long %>%
  count(sex, item_sex, reporter, stage_label) %>%
  group_by(sex, item_sex, reporter) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  ungroup() %>%
  arrange(sex, item_sex, reporter, stage_label)

write.csv(
  categ_freq,
  file.path(out_root_dk, "pds_categ_by_wave.csv"),
  row.names = FALSE
)
write.csv(
  categ_overall,
  file.path(out_root_dk, "pds_categ_overall.csv"),
  row.names = FALSE
)

cat("\n=== ABCD categorical pubertal stage (pooled across waves) ===\n")
print(categ_overall, n = Inf)

# ---------------------------------------------------------------------------
# SEX × REPORTER LONG DATASETS
# ---------------------------------------------------------------------------
covariate_cols <- c(
  "id",
  "wave",
  "age",
  "sex",
  "race",
  "site",
  "bmi",
  "bmi_z",
  "caregiver",
  "caregiver_switched",
  "n_cg_switches",
  "ever_switched_cg"
)

female_parent <- data %>%
  filter(sex == 2) %>%
  select(
    all_of(covariate_cols),
    peta = peta_p,
    petb = petb_p,
    petc = petc_p,
    petd = petdf_p,
    fpete = fpete_p,
    pds_comp = pds_comp_f_p,
    pds_categ = pds_categ_fp
  ) %>%
  mutate(reporter = "parent")

female_youth <- data %>%
  filter(sex == 2) %>%
  select(
    all_of(covariate_cols),
    peta = peta_y,
    petb = petb_y,
    petc = petc_y,
    petd = petdf_y,
    fpete = fpete_y,
    pds_comp = pds_comp_f_y,
    pds_categ = pds_categ_fy
  ) %>%
  mutate(reporter = "youth")

male_parent <- data %>%
  filter(sex == 1) %>%
  select(
    all_of(covariate_cols),
    peta = peta_p,
    petb = petb_p,
    petc = petc_p,
    petd = petdm_p,
    mpete = mpete_p,
    pds_comp = pds_comp_m_p,
    pds_categ = pds_categ_mp
  ) %>%
  mutate(reporter = "parent")

male_youth <- data %>%
  filter(sex == 1) %>%
  select(
    all_of(covariate_cols),
    peta = peta_y,
    petb = petb_y,
    petc = petc_y,
    petd = petdm_y,
    mpete = mpete_y,
    pds_comp = pds_comp_m_y,
    pds_categ = pds_categ_my
  ) %>%
  mutate(reporter = "youth")

all_long <- bind_rows(female_parent, female_youth, male_parent, male_youth)

# ---------------------------------------------------------------------------
# WRITE OUTPUTS
# ---------------------------------------------------------------------------
write.csv(
  female_parent,
  file.path(out_root, "female_parent_long.csv"),
  row.names = FALSE
)
write.csv(
  female_youth,
  file.path(out_root, "female_youth_long.csv"),
  row.names = FALSE
)
write.csv(
  male_parent,
  file.path(out_root, "male_parent_long.csv"),
  row.names = FALSE
)
write.csv(
  male_youth,
  file.path(out_root, "male_youth_long.csv"),
  row.names = FALSE
)
write.csv(all_long, file.path(out_root, "all_long.csv"), row.names = FALSE)

cat("Written to:", out_root, "\n")
cat(
  "Rows per dataset — female parent:",
  nrow(female_parent),
  "| female youth:",
  nrow(female_youth),
  "| male parent:",
  nrow(male_parent),
  "| male youth:",
  nrow(male_youth),
  "\n"
)
