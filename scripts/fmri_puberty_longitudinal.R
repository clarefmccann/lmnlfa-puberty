## fmri_puberty_longitudinal.R
## Models the LONGITUDINAL relationship between N-back betas
## (fmri_nback_foundation.R) and pubertal development, for every contrast
## in fmri_contrasts.R. Unlike fmri_puberty_association.R (which uses the
## cross-sectional MNLFA factor score -- one value per person, from a
## single random wave), this uses pds_comp: ABCD's continuous PDS
## composite (1-4 scale), computed at EVERY annual wave. Joined to fMRI by
## (id, wave), this gives genuine repeated-measures data -- multiple
## (pds_comp, beta) pairs per person across development.
##
## pds_comp has a row per reporter (parent/youth); one value per
## person-wave is needed, so youth self-report is preferred, falling back
## to parent report when youth is missing -- same convention used for
## pds_categ in fmri_nback_descriptives.R.
##
## Model: a GAMM (mgcv), matching the convention already used elsewhere in
## this project for pubertal trajectories (03_gamms_hpc.R). Male and female
## are fit as COMPLETELY SEPARATE models (not a shared model with a sex
## interaction/by-factor term) -- each sex's data is subset first, then:
##   beta ~ s(pds_comp, k = 5) + s(id_fac, bs = "re") + s(id_fac, pds_comp, bs = "re")
## s(pds_comp):              population-level, possibly nonlinear curve
##                           (not forced to be a straight line)
## s(id_fac, bs = "re"):     random intercept per person
## s(id_fac, pds_comp, bs="re"): random slope per person (individual tempo)
## mgcv estimates smoothness via REML penalization -- no external nonlinear
## optimizer needed (unlike lme4/nloptr), so this has no extra HPC install
## dependency beyond mgcv itself.
## Fit for the global mean beta (averaged across all regions) and for
## amygdala/insula/rostral middle frontal specifically (main plots), plus
## a simpler sweep across every region (coefficient table only: population
## smooth + random intercept, no random slope -- keeps ~98 x 2 sexes
## per-contrast model fits fast and robust). Every model below is fit
## per-sex from the ground up, so there is no sex covariate anywhere.
##
## Plots are longitudinal spaghetti plots: each dot is one person-wave
## observation (coloured/shaped by wave), thin lines connect one
## participant's own observations in chronological (wave) order -- via
## geom_path() on data pre-sorted by (id, wave), not geom_line(), since
## pds_comp isn't guaranteed to increase monotonically wave-to-wave and
## geom_line() would otherwise connect points in x-sorted order and draw
## the wrong trajectory. The population GAMM curve (+ 95% CI ribbon,
## random effects excluded) is overlaid per sex panel.
##
## Usage: Rscript fmri_puberty_longitudinal.R
## Output: outputs/fmri_puberty_longitudinal_<contrast>/*.png, *.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(parallel)
})

# ~1,000 independent GAM fits total (per-region x per-sex sweep, dominates;
# plus a handful of richer key-region/global fits) -- embarrassingly
# parallel, so fan out across whatever slots the job actually has
# (Sys.getenv("NSLOTS") is set by SGE's -pe shared N; defaults to 1 for a
# local/interactive run). parallel::mclapply forks, so this only works on
# Unix (fine for Hoffman2 and macOS/Linux local testing, not Windows).
n_cores <- max(1L, as.integer(Sys.getenv("NSLOTS", unset = "1")))
cat("Parallel model fitting across", n_cores, "core(s) (NSLOTS)\n")

root_path <- Sys.getenv("HOME_DIR")
if (!nzchar(root_path)) root_path <- Sys.getenv("HOME")

data_dir <- Sys.getenv("DATA_DIR")
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  data_dir <- file.path(root_path, "projects/abcd-projs/dissertation/study1/outputs")
}
out_base <- Sys.getenv("OUT_DIR")
if (!nzchar(out_base)) out_base <- data_dir

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
contrasts_file <- file.path(script_dir, "fmri_contrasts.R")
if (!file.exists(contrasts_file)) contrasts_file <- file.path("scripts", "fmri_contrasts.R")
source(contrasts_file)

wave_palette_display <- setNames(pal_waves, wave_display_labels[names(pal_waves)])
wave_shapes_display <- setNames(pal_shapes_waves, wave_display_labels[names(pal_shapes_waves)])
key_regions <- c("amygdala", "insula", "rostral middle frontal")

# ---------------------------------------------------------------------------
# PDS COMPOSITE (pds_comp), ONE VALUE PER PERSON-WAVE -- contrast-independent
# ---------------------------------------------------------------------------
puberty_path <- file.path(data_dir, "all_long.csv")
if (!file.exists(puberty_path)) {
  stop("Cannot find all_long.csv in ", data_dir, " -- run 00_data_foundation.R first.")
}
pds_comp_by_wave <- read.csv(puberty_path) %>%
  select(id, wave, reporter, pds_comp) %>%
  filter(!is.na(pds_comp)) %>%
  mutate(reporter_rank = if_else(reporter == "youth", 1, 2)) %>%
  arrange(id, wave, reporter_rank) %>%
  distinct(id, wave, .keep_all = TRUE) %>%
  mutate(wave = as.character(wave)) %>%
  select(id, wave, pds_comp)

cat(
  "PDS composite resolved for", nrow(pds_comp_by_wave), "person-waves",
  "(youth-preferred, parent fallback)\n"
)

# -- sweep model: population smooth + random intercept only, fit on
# already-single-sex data (grouping includes sex_label below, so each call
# only ever sees one sex). Returns the pds_comp smooth term's edf/F/p
# (edf > 1 indicates a nonlinear population curve; edf = 1 is a straight
# line).
fit_sweep_gam <- function(data, label) {
  data <- data %>% mutate(id_fac = factor(id))
  m <- tryCatch(
    suppressWarnings(gam(
      beta ~ s(pds_comp, k = 5) + s(id_fac, bs = "re"),
      data = data, method = "REML"
    )),
    error = function(e) NULL
  )
  if (is.null(m)) {
    return(tibble(
      edf = NA_real_, ref_df = NA_real_, stat = NA_real_, p = NA_real_,
      region = label, n_obs = nrow(data), n_people = length(unique(data$id))
    ))
  }
  st <- as.data.frame(summary(m)$s.table)
  tibble(
    edf = st["s(pds_comp)", "edf"],
    ref_df = st["s(pds_comp)", "Ref.df"],
    stat = st["s(pds_comp)", "F"],
    p = st["s(pds_comp)", "p-value"],
    region = label,
    n_obs = nrow(data),
    n_people = length(unique(data$id))
  )
}

run_longitudinal_for_contrast <- function(cst) {
  clabel <- fmri_contrasts[[cst]]
  cat("\n==================== Contrast:", cst, "(", clabel, ") ====================\n")

  out_dir <- file.path(out_base, paste0("fmri_puberty_longitudinal_", cst))
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
    cat("Skipping", cst, "-- no region columns found in", fmri_path, "\n")
    return(invisible(NULL))
  }

  fmri <- fmri %>%
    mutate(
      wave = as.character(wave),
      wave_label = factor(wave_display_labels[wave], levels = wave_display_labels),
      sex_label = factor(case_when(
        sex == 1 ~ "Male",
        sex == 2 ~ "Female",
        TRUE ~ NA_character_
      ), levels = c("Female", "Male"))
    ) %>%
    filter(!is.na(wave_label), !is.na(sex_label))

  merged <- fmri %>% inner_join(pds_comp_by_wave, by = c("id", "wave"))
  cat(
    "fMRI person-waves:", nrow(fmri),
    "| with a matched pds_comp:", nrow(merged),
    "(", round(100 * nrow(merged) / nrow(fmri), 1), "% )\n"
  )
  n_multi_wave <- merged %>% count(id) %>% filter(n > 1) %>% nrow()
  cat(
    "People with >1 imaging wave in this contrast (true repeated measures):",
    n_multi_wave, "of", length(unique(merged$id)), "\n"
  )

  beta_long <- merged %>%
    select(id, wave_label, sex_label, pds_comp, all_of(region_cols)) %>%
    pivot_longer(all_of(region_cols), names_to = "region_key", values_to = "beta") %>%
    filter(!is.na(beta)) %>%
    left_join(region_lookup, by = "region_key")

  # -------------------------------------------------------------------------
  # (1) Coefficient sweep across every region -- the bulk of this script's
  # compute (~196 independent fits per contrast: every region x both
  # sexes), so fanned out across n_cores via mclapply rather than a serial
  # group_modify(). group_split()/group_keys() are guaranteed to return
  # groups in the same order, so the keys can be safely bound back onto
  # the parallel results afterward.
  # -------------------------------------------------------------------------
  sweep_grouped <- beta_long %>% group_by(region_key, atlas, hemi, region_label, sex_label)
  sweep_keys <- group_keys(sweep_grouped)
  sweep_groups <- group_split(sweep_grouped)

  sweep_results <- mclapply(
    sweep_groups,
    function(g) fit_sweep_gam(g, g$region_key[1]),
    mc.cores = n_cores
  )
  region_coefs <- bind_cols(sweep_keys, bind_rows(sweep_results) %>% select(-region)) %>%
    arrange(p)
  write.csv(
    region_coefs,
    file.path(out_dir, "region_pds_comp_coefficients.csv"),
    row.names = FALSE
  )
  cat("Top 10 region x sex by pds_comp smooth p-value (population-curve model):\n")
  print(as.data.frame(head(region_coefs, 10)), digits = 3)

  # -------------------------------------------------------------------------
  # (2)/(3) Global mean beta + key regions: a completely separate GAMM per
  # sex (fit on that sex's own subset, no sex covariate at all -- not one
  # shared model with a by-sex smooth), longitudinal spaghetti plot with
  # each sex's own population curve overlaid in its own panel.
  # -------------------------------------------------------------------------
  fit_one_sex <- function(d_sx, region_label, sx) {
    m <- tryCatch(
      suppressWarnings(gam(
        beta ~ s(pds_comp, k = 5) + s(id_fac, bs = "re") + s(id_fac, pds_comp, bs = "re"),
        data = d_sx, method = "REML"
      )),
      error = function(e) NULL
    )
    if (is.null(m)) {
      return(list(coef_tbl = NULL, pred = NULL))
    }

    st <- as.data.frame(summary(m)$s.table)
    pt <- as.data.frame(summary(m)$p.table)
    coef_tbl <- bind_rows(
      tibble(
        term = rownames(st), edf = st[, "edf"], ref_df = st[, "Ref.df"],
        stat = st[, "F"], p = st[, "p-value"], estimate = NA_real_, se = NA_real_
      ),
      tibble(
        term = rownames(pt), edf = NA_real_, ref_df = NA_real_,
        stat = pt[, "t value"], p = pt[, "Pr(>|t|)"],
        estimate = pt[, "Estimate"], se = pt[, "Std. Error"]
      )
    ) %>%
      mutate(
        region = region_label, sex_label = sx,
        n_obs = nrow(d_sx), n_people = length(unique(d_sx$id))
      )

    re_terms <- vapply(m$smooth, function(s) s$label, character(1))
    re_terms <- re_terms[grepl("^s\\(id_fac", re_terms)]

    pred <- expand.grid(pds_comp = seq(min(d_sx$pds_comp), max(d_sx$pds_comp), length.out = 50))
    pred$id_fac <- d_sx$id_fac[1]
    pr <- predict(m, newdata = pred, type = "response", se.fit = TRUE, exclude = re_terms)
    pred$fit <- pr$fit
    pred$lo <- pr$fit - 1.96 * pr$se.fit
    pred$hi <- pr$fit + 1.96 * pr$se.fit
    pred$sex_label <- sx

    list(coef_tbl = coef_tbl, pred = pred)
  }

  # Fitting only (no ggplot/ggsave) -- this is the part that's actually
  # expensive and worth parallelizing. Keeping graphics calls out of the
  # forked worker matters beyond just mgcv's own cost: ggsave()/graphics
  # devices are not reliably fork-safe on every platform (macOS's
  # Objective-C/CoreText font machinery in particular crashes forked
  # children if touched mid-fork), so plotting is done afterward in the
  # main process by render_plot() instead.
  fit_region_model <- function(d, region_label) {
    if (nrow(d) < 20) {
      return(list(coef_tbl = NULL, pred_grid = NULL, d = NULL, region_label = region_label))
    }
    d <- d %>% mutate(id_fac = factor(id))

    fits <- lapply(intersect(c("Female", "Male"), unique(d$sex_label)), function(sx) {
      fit_one_sex(d %>% filter(sex_label == sx), region_label, sx)
    })

    coef_tbl <- bind_rows(lapply(fits, `[[`, "coef_tbl"))
    if (nrow(coef_tbl) == 0) coef_tbl <- NULL
    pred_grid <- bind_rows(lapply(fits, `[[`, "pred"))
    if (nrow(pred_grid) == 0) pred_grid <- NULL

    list(coef_tbl = coef_tbl, pred_grid = pred_grid, d = d, region_label = region_label)
  }

  render_plot <- function(fit_result, title, filename) {
    if (is.null(fit_result$d)) {
      cat("Skipping", filename, "-- too few rows\n")
      return(invisible(NULL))
    }
    d <- fit_result$d
    pred_grid <- fit_result$pred_grid
    d_sorted <- d %>% arrange(id, wave_label)

    p <- ggplot(d_sorted, aes(x = pds_comp, y = beta)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_path(aes(group = id), alpha = 0.15, linewidth = 0.3, colour = "grey40") +
      geom_point(aes(colour = wave_label, shape = wave_label), alpha = 0.6, size = 1.4) +
      scale_colour_manual(values = wave_palette_display, breaks = as.character(wave_display_labels)) +
      scale_shape_manual(values = wave_shapes_display, breaks = as.character(wave_display_labels)) +
      facet_wrap(~sex_label) +
      labs(
        title = title,
        subtitle = paste0(
          "Longitudinal GAMM, fit separately per sex (beta ~ s(pds_comp) + s(id, re) + s(id, pds_comp, re)); ",
          "dots = observations, lines = one participant; ",
          "n = ", nrow(d), " person-waves, ", length(unique(d$id)), " people"
        ),
        x = "PDS composite (pds_comp)",
        y = "Beta",
        colour = "Wave",
        shape = "Wave"
      ) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom")

    if (!is.null(pred_grid)) {
      p <- p +
        geom_ribbon(
          data = pred_grid, aes(x = pds_comp, ymin = lo, ymax = hi),
          inherit.aes = FALSE, fill = pal_primary_fill, alpha = 0.3
        ) +
        geom_line(
          data = pred_grid, aes(x = pds_comp, y = fit),
          inherit.aes = FALSE, colour = pal_primary, linewidth = 1.1
        )
    }

    ggsave(file.path(out_dir, filename), p, width = 9, height = 5.5, dpi = 180)
    invisible(NULL)
  }

  person_wave_mean <- beta_long %>%
    group_by(id, wave_label, sex_label, pds_comp) %>%
    summarise(beta = mean(beta), .groups = "drop")

  # Global mean + the 3 key regions: only 4 independent fits, so run the
  # (expensive, graphics-free) fitting step across up to 4 cores -- never
  # request more workers than jobs -- then render+save every plot
  # sequentially afterward in the main process.
  key_jobs <- c(
    list(list(
      data = person_wave_mean,
      title = paste0("Mean N-back ", clabel, " beta vs. PDS composite, longitudinal"),
      filename = "mean_beta_vs_pds_comp.png",
      region_label = "Global mean"
    )),
    lapply(key_regions, function(rl) {
      list(
        data = beta_long %>%
          filter(region_label == rl) %>%
          group_by(id, wave_label, sex_label, pds_comp) %>%
          summarise(beta = mean(beta), .groups = "drop"),
        title = paste0("N-back ", clabel, " beta (", rl, ") vs. PDS composite, longitudinal"),
        filename = paste0(gsub(" ", "_", rl), "_vs_pds_comp.png"),
        region_label = rl
      )
    })
  )

  key_fits <- mclapply(
    key_jobs,
    function(j) fit_region_model(j$data, j$region_label),
    mc.cores = min(n_cores, length(key_jobs))
  )

  for (i in seq_along(key_jobs)) {
    render_plot(key_fits[[i]], key_jobs[[i]]$title, key_jobs[[i]]$filename)
  }

  key_coef_tbl <- bind_rows(lapply(key_fits, `[[`, "coef_tbl"))
  write.csv(
    key_coef_tbl,
    file.path(out_dir, "key_region_pds_comp_model_coefficients.csv"),
    row.names = FALSE
  )

  cat("Outputs written to:", out_dir, "\n")
}

for (cst in names(fmri_contrasts)) {
  run_longitudinal_for_contrast(cst)
}
