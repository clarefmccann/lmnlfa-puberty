## fmri_nback_foundation.R
## Loads raw ABCD task fMRI data (Emotional N-Back task) for every contrast
## in fmri_contrasts.R and writes clean long-format datasets in the same
## id/wave/covariate shape as the puberty data (00_data_foundation.R), so
## the two can be joined or modeled with the same downstream tooling.
##
## Inclusion: mr_y_qc__incl__tfmri__nback_indicator == 1 (ABCD's own QC flag
## for whether this person-wave's n-back data is recommended for use --
## task-level, so the same flag applies to every contrast below). Rows that
## are 0 or NA (imaging wasn't collected that wave -- MRI only happens
## every other year, unlike the annual PDS assessments) are dropped, since
## there's no usable beta data for them either way.
##
## Regions: for each contrast, every mr_y_tfmri__nback__<contrast>__
## {aseg,dsk}__*_beta variable is pulled directly from the ABCD data
## dictionary at run time (not hardcoded), so this stays correct if a
## future NBDCtoolsData release adds/renames regions. aseg = subcortical
## (30 vars: 13 paired left/right structures + 4 unpaired midline ones like
## brain stem/CSF). dsk = cortical, Desikan-Killiany atlas (68 vars: 34
## regions x 2 hemispheres). Column names are shortened from the full ABCD
## variable name to "<atlas>_<region>_<hemi>_beta" (e.g. "aseg_hc_lh_beta"),
## with a companion lookup CSV per contrast mapping each short key back to
## its full anatomical label for plotting.
##
## All contrasts share one create_dataset() call (covariates + QC + every
## contrast's region vars, loaded once) and are then split into per-
## contrast output files.
##
## Outputs per contrast <c> (see fmri_contrasts.R for the tag list):
##   all_fmri_nback_<c>_long.csv, female_fmri_nback_<c>_long.csv,
##   male_fmri_nback_<c>_long.csv, fmri_nback_<c>_region_lookup.csv

# export DATA_DIR="/u/project/silvers/data/ABCD/ABCD-release-6.0/cfm/physical-health/puberty"
# export OUT_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"
# Rscript fmri_nback_foundation.R

library(dplyr)
library(tidyr)
library(NBDCtools)

set.seed(90025)

# ---------------------------------------------------------------------------
# PATHS (identical convention to 00_data_foundation.R)
# ---------------------------------------------------------------------------
root_path <- Sys.getenv("HOME_DIR")
if (!nzchar(root_path)) {
  root_path <- Sys.getenv("HOME")
}

data_root <- file.path(
  root_path,
  "projects/abcd-projs/abcd-data-release-7.0/phenotype"
)
if (!dir.exists(data_root)) {
  stop("Cannot locate nbdc-tools-data: ", data_root)
}

out_root <- file.path(
  root_path,
  "projects/abcd-projs/dissertation/study1/outputs"
)
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

script_dir <- Sys.getenv("SGE_O_WORKDIR")
if (!nzchar(script_dir)) {
  script_dir <- tryCatch(
    {
      ofile <- sys.frames()[[1]]$ofile
      if (is.null(ofile)) "scripts" else dirname(ofile)
    },
    error = function(e) "scripts"
  )
}
contrasts_file <- file.path(script_dir, "fmri_contrasts.R")
if (!file.exists(contrasts_file)) contrasts_file <- file.path("scripts", "fmri_contrasts.R")
source(contrasts_file)

# ---------------------------------------------------------------------------
# RESOLVE REGION VARIABLE NAMES FROM THE DATA DICTIONARY, PER CONTRAST
# ---------------------------------------------------------------------------
dd <- get_dd_abcd()

qc_var <- "mr_y_qc__incl__tfmri__nback_indicator"
if (!qc_var %in% dd$name) {
  stop("QC variable not found in data dictionary: ", qc_var)
}

resolve_contrast_vars <- function(contrast) {
  beta_vars <- dd %>%
    filter(grepl(paste0("^mr_y_tfmri__nback__", contrast, "__(aseg|dsk)__.*_beta$"), name)) %>%
    pull(name) %>%
    sort()

  if (length(beta_vars) == 0) {
    stop("No ", contrast, " aseg/dsk beta variables found -- check the data dictionary.")
  }

  # short column key: strip the shared prefix, collapse remaining "__" to "_"
  region_key <- beta_vars %>%
    sub(paste0("^mr_y_tfmri__nback__", contrast, "__"), "", .) %>%
    gsub("__", "_", .)

  if (any(duplicated(region_key))) {
    stop("Region key collision for contrast ", contrast, " -- inspect beta_vars manually.")
  }

  # region lookup: key -> atlas, hemisphere, anatomical label (parsed from
  # the dictionary's free-text label, e.g. "... in Subcortical ROI:
  # hippocampus (left hemisphere)" -> "hippocampus")
  region_lookup <- tibble(var = beta_vars, region_key = region_key) %>%
    left_join(dd %>% select(name, label), by = c("var" = "name")) %>%
    mutate(
      atlas = if_else(grepl("__aseg__", var), "aseg", "dsk"),
      hemi = case_when(
        grepl("__lh_beta$", var) ~ "left",
        grepl("__rh_beta$", var) ~ "right",
        TRUE ~ NA_character_
      ),
      region_label = sub(".*ROI: ", "", label),
      region_label = sub(" \\((left|right) hemisphere\\)$", "", region_label)
    ) %>%
    select(region_key, atlas, hemi, region_label)

  list(beta_vars = beta_vars, region_key = region_key, region_lookup = region_lookup)
}

contrast_info <- lapply(names(fmri_contrasts), resolve_contrast_vars)
names(contrast_info) <- names(fmri_contrasts)

for (cst in names(fmri_contrasts)) {
  ci <- contrast_info[[cst]]
  n_aseg <- sum(grepl("^aseg_", ci$region_key))
  n_dsk <- sum(grepl("^dsk_", ci$region_key))
  cat(
    sprintf(
      "%-10s (%s): %d beta variables ( %d aseg, %d dsk )\n",
      cst, fmri_contrasts[[cst]], length(ci$beta_vars), n_aseg, n_dsk
    )
  )
}

all_beta_vars <- unique(unlist(lapply(contrast_info, `[[`, "beta_vars")))

# ---------------------------------------------------------------------------
# LOAD DATA -- one call across all contrasts
# ---------------------------------------------------------------------------
covariate_vars <- c(
  "ab_g_dyn__visit_age",
  "ab_g_stc__cohort_sex",
  "ab_g_stc__cohort_ethnrace__mhisp",
  "ab_g_dyn__design_site"
)

cat("\nLoading", length(all_beta_vars), "beta variables across", length(fmri_contrasts), "contrasts...\n")
raw <- create_dataset(
  dir_data = data_root,
  study = "abcd",
  vars = c(covariate_vars, qc_var, all_beta_vars),
  value_to_na = TRUE,
  bind_shadow = FALSE
)

# ---------------------------------------------------------------------------
# ANNUAL WAVE LABELS (same mapping as 00_data_foundation.R; only bl/fu2/
# fu4/fu6 will actually have imaging data, since MRI is collected every
# other year, but session_id filtering + the QC/completeness filters below
# handle that naturally without hardcoding which waves have scans)
# ---------------------------------------------------------------------------
wave_labels <- c(
  "ses-00A" = "bl",
  "ses-01A" = "fu1",
  "ses-02A" = "fu2",
  "ses-03A" = "fu3",
  "ses-04A" = "fu4",
  "ses-05A" = "fu5",
  "ses-06A" = "fu6"
)

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
    site = ab_g_dyn__design_site,
    qc_nback_incl = !!qc_var
  )

n_total <- nrow(data)
n_qc_pass <- sum(data$qc_nback_incl == 1, na.rm = TRUE)
n_qc_fail <- sum(data$qc_nback_incl == 0, na.rm = TRUE)
n_qc_na <- sum(is.na(data$qc_nback_incl))
cat(
  "\nQC inclusion (mr_y_qc__incl__tfmri__nback_indicator) across all",
  "person-waves with a session in bl-fu6 (shared across every contrast):\n"
)
cat("  Total rows:      ", n_total, "\n")
cat("  Pass (1, kept):  ", n_qc_pass, "\n")
cat("  Fail (0):        ", n_qc_fail, "\n")
cat("  Missing (NA):    ", n_qc_na, "\n")

data <- data %>% filter(qc_nback_incl == 1)
cat("\nRows after QC filter (shared across every contrast):", nrow(data), "\n")
cat("Rows by wave:\n")
print(table(data$wave, useNA = "ifany"))

covariate_cols <- c("id", "wave", "age", "sex", "race", "site")

# ---------------------------------------------------------------------------
# PER-CONTRAST OUTPUTS
# ---------------------------------------------------------------------------
for (cst in names(fmri_contrasts)) {
  ci <- contrast_info[[cst]]

  all_fmri_long <- data %>%
    select(all_of(covariate_cols), all_of(ci$beta_vars)) %>%
    rename_with(~ci$region_key, all_of(ci$beta_vars))

  female_fmri <- all_fmri_long %>% filter(sex == 2)
  male_fmri <- all_fmri_long %>% filter(sex == 1)

  write.csv(
    all_fmri_long,
    file.path(out_root, paste0("all_fmri_nback_", cst, "_long.csv")),
    row.names = FALSE
  )
  write.csv(
    female_fmri,
    file.path(out_root, paste0("female_fmri_nback_", cst, "_long.csv")),
    row.names = FALSE
  )
  write.csv(
    male_fmri,
    file.path(out_root, paste0("male_fmri_nback_", cst, "_long.csv")),
    row.names = FALSE
  )
  write.csv(
    ci$region_lookup,
    file.path(out_root, paste0("fmri_nback_", cst, "_region_lookup.csv")),
    row.names = FALSE
  )

  cat(
    sprintf(
      "\n[%s] written -- all: %d | female: %d | male: %d | regions: %d\n",
      cst, nrow(all_fmri_long), nrow(female_fmri), nrow(male_fmri), length(ci$region_key)
    )
  )
}

cat("\nAll contrasts written to:", out_root, "\n")
