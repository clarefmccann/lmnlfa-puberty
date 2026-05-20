## 05_gmm_hpc.R
## Growth Mixture Models for pubertal trajectories (Phase 1: PDS composite)
##
## Identifies latent classes of developmental trajectories using lcmm::hlme().
## Phase 1 (this script): PDS composite, nwg=FALSE (equal within-class variance),
##   quadratic trajectory in age, random intercept + slope within class.
##
## Phase 2 (future): item-level trajectories via multlcmm(), and/or LMNLFA
##   factor scores as outcome after those are computed.
##
## Usage:
##   Rscript 05_gmm_hpc.R <dataset_name> [max_k] [n_starts]
##   dataset_name : female_parent | female_youth | male_parent | male_youth
##   max_k        : max latent classes to fit (default 6)
##   n_starts     : random starts per K>1 model via gridsearch (default 10)
##
## Outputs (all in OUT_DIR/gmm_trajectories/):
##   {ds}_gmm_k{k}.rds             — saved hlme model per K
##   {ds}_gmm_fit_indices.csv       — BIC, AIC, loglik, entropy, AvePP per K
##   {ds}_gmm_class_assignments.csv — id, class, posterior probs
##   {ds}_gmm_trajectories_k{k}.png — trajectory profiles for each K
##   {ds}_gmm_model_selection.png   — BIC + entropy plot across K

suppressPackageStartupMessages({
  library(lcmm)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(90025)

# ---------------------------------------------------------------------------
# ARGUMENTS
# ---------------------------------------------------------------------------
args    <- commandArgs(trailingOnly = TRUE)
ds_name <- args[1]
max_k   <- if (length(args) >= 2) as.integer(args[2]) else 6L
n_starts <- if (length(args) >= 3) as.integer(args[3]) else 10L

valid_ds <- c("female_parent", "female_youth", "male_parent", "male_youth")
if (!ds_name %in% valid_ds) {
  stop("dataset_name must be one of: ", paste(valid_ds, collapse = ", "))
}
cat("Dataset:", ds_name, " max_k:", max_k, " n_starts:", n_starts, "\n")

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
root_path <- Sys.getenv("HOME_DIR")
if (!nzchar(root_path)) root_path <- Sys.getenv("HOME")

data_dir <- Sys.getenv("DATA_DIR")
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  data_dir <- file.path(
    root_path,
    "projects/abcd-projs/abcd-data-release-6.0/cfm/physical-health/puberty"
  )
}
if (!dir.exists(data_dir)) stop("Cannot locate data directory: ", data_dir)

out_base <- Sys.getenv("OUT_DIR")
if (!nzchar(out_base)) {
  out_base <- file.path(
    root_path, "projects/abcd-projs/dissertation/study1/outputs"
  )
}
out_dir <- file.path(out_base, "gmm_trajectories")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat("Output dir:", out_dir, "\n")

# ---------------------------------------------------------------------------
# LOAD DATA
# ---------------------------------------------------------------------------
csv_path <- file.path(data_dir, paste0(ds_name, "_long.csv"))
if (!file.exists(csv_path)) stop("File not found: ", csv_path)

df <- read.csv(csv_path) %>%
  filter(!is.na(pds_comp), !is.na(age)) %>%
  arrange(id, age)

cat("Loaded:", nrow(df), "rows,", n_distinct(df$id), "participants\n")
cat("Age range:", round(range(df$age), 1), "\n")
cat("PDS range:", round(range(df$pds_comp), 2), "\n")

# hlme needs a numeric subject id
df <- df %>% mutate(id_num = as.integer(factor(id)))

# ---------------------------------------------------------------------------
# HELPERS
# ---------------------------------------------------------------------------

# Relative entropy (0–1, higher = better class separation)
calc_entropy <- function(model) {
  if (model$ng == 1) return(NA_real_)
  pp <- as.matrix(model$pprob[, -(1:2)])  # drop id + assigned class columns
  pp[pp <= 0] <- .Machine$double.eps
  n <- nrow(pp)
  k <- ncol(pp)
  1 + sum(pp * log(pp)) / (n * log(k))
}

# Average posterior probability for assigned class (measure of certainty)
calc_avepp <- function(model) {
  if (model$ng == 1) return(NA_real_)
  pp  <- model$pprob
  cls <- pp$class
  probs <- as.matrix(pp[, -(1:2)])
  mean(sapply(seq_len(nrow(pp)), function(i) probs[i, cls[i]]))
}

# Class sizes as percentages
class_pct <- function(model) {
  if (model$ng == 1) return("100%")
  pp  <- model$pprob
  tab <- table(pp$class)
  pct <- round(100 * tab / sum(tab), 1)
  paste(paste0("C", names(pct), ":", pct, "%"), collapse = " | ")
}

# Predicted class trajectories on a fine age grid
predict_trajectories <- function(model, age_min, age_max, n = 150) {
  newdata <- data.frame(age = seq(age_min, age_max, length.out = n))
  pred    <- predictY(model, newdata = newdata, var.time = "age")
  traj    <- as.data.frame(pred$pred)
  traj$age <- newdata$age
  k       <- model$ng
  if (k == 1) {
    traj <- traj %>% rename(Ypred_class1 = Ypred)
  }
  traj %>%
    pivot_longer(
      cols      = starts_with("Ypred_class"),
      names_to  = "class",
      values_to = "pds_pred"
    ) %>%
    mutate(class = sub("Ypred_class", "Class ", class))
}

# Class sizes from posterior modal assignment for legend labelling
class_labels <- function(model) {
  if (model$ng == 1) return(c("Class 1" = "Class 1 (100%)"))
  pp  <- model$pprob
  tab <- table(pp$class)
  pct <- round(100 * tab / sum(tab), 1)
  setNames(
    paste0("Class ", names(tab), " (", pct, "%)"),
    paste0("Class ", names(tab))
  )
}

# ---------------------------------------------------------------------------
# FIT K=1 REFERENCE MODEL
# ---------------------------------------------------------------------------
rds_k1 <- file.path(out_dir, paste0(ds_name, "_gmm_k1.rds"))

if (file.exists(rds_k1)) {
  cat("\nLoading K=1 from cache:", rds_k1, "\n")
  m1 <- readRDS(rds_k1)
} else {
  cat("\nFitting K=1 reference model...\n")
  t0 <- proc.time()
  m1 <- hlme(
    pds_comp ~ age + I(age^2),
    random  = ~ age,
    subject = "id_num",
    ng      = 1,
    data    = df
  )
  cat("K=1 done in", round((proc.time() - t0)[["elapsed"]] / 60, 1), "min\n")
  cat("  loglik:", round(m1$loglik, 2), " BIC:", round(m1$BIC, 2), "\n")
  saveRDS(m1, rds_k1)
}

# ---------------------------------------------------------------------------
# FIT K=2..max_k WITH GRIDSEARCH RANDOM STARTS
# ---------------------------------------------------------------------------
models      <- vector("list", max_k)
models[[1]] <- m1

for (k in 2:max_k) {
  rds_k <- file.path(out_dir, paste0(ds_name, "_gmm_k", k, ".rds"))

  if (file.exists(rds_k)) {
    cat("\nLoading K =", k, "from cache:", rds_k, "\n")
    models[[k]] <- readRDS(rds_k)
  } else {
    cat("\nFitting K =", k, "(", n_starts, "random starts via gridsearch)...\n")
    t0 <- proc.time()

    models[[k]] <- tryCatch(
      gridsearch(
        hlme(
          pds_comp ~ age + I(age^2),
          mixture = ~ age + I(age^2),
          random  = ~ age,
          subject = "id_num",
          ng      = k,
          nwg     = FALSE,
          data    = df
        ),
        rep     = n_starts,
        maxiter = 30,
        minit   = m1
      ),
      error = function(e) {
        message("K=", k, " failed: ", e$message)
        NULL
      }
    )

    if (!is.null(models[[k]])) {
      elapsed <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
      cat("  K =", k, "done in", elapsed, "min\n")
      cat("  loglik:", round(models[[k]]$loglik, 2),
          " BIC:", round(models[[k]]$BIC, 2),
          " entropy:", round(calc_entropy(models[[k]]), 3), "\n")
      cat("  Class sizes:", class_pct(models[[k]]), "\n")
      saveRDS(models[[k]], rds_k)
    }
  }
}

# Remove NULL (failed fits)
models <- Filter(Negate(is.null), models)

# ---------------------------------------------------------------------------
# FIT INDICES TABLE
# ---------------------------------------------------------------------------
fit_table <- bind_rows(lapply(models, function(m) {
  tibble(
    dataset  = ds_name,
    K        = m$ng,
    loglik   = m$loglik,
    AIC      = m$AIC,
    BIC      = m$BIC,
    entropy  = calc_entropy(m),
    avepp    = calc_avepp(m),
    class_pct = class_pct(m)
  )
}))

cat("\n=== Model fit indices ===\n")
print(fit_table[, c("K", "loglik", "BIC", "entropy", "avepp", "class_pct")],
      digits = 3, row.names = FALSE)

write.csv(fit_table,
          file.path(out_dir, paste0(ds_name, "_gmm_fit_indices.csv")),
          row.names = FALSE)
cat("Fit indices saved.\n")

# ---------------------------------------------------------------------------
# MODEL SELECTION PLOT  (BIC + entropy across K)
# ---------------------------------------------------------------------------
fit_long <- fit_table %>%
  select(K, BIC, entropy) %>%
  pivot_longer(c(BIC, entropy), names_to = "metric", values_to = "value") %>%
  mutate(
    metric = recode(metric, BIC = "BIC (lower = better)",
                            entropy = "Entropy (higher = better)")
  )

p_sel <- ggplot(fit_long, aes(x = K, y = value)) +
  geom_line(linewidth = 1, colour = "#2166ac") +
  geom_point(size = 3, colour = "#2166ac") +
  facet_wrap(~ metric, scales = "free_y") +
  scale_x_continuous(breaks = seq_len(max_k)) +
  labs(
    title    = paste("GMM model selection:", ds_name),
    subtitle = paste0(n_starts, " random starts per K; quadratic trajectory, nwg=FALSE"),
    x = "Number of classes (K)", y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, paste0(ds_name, "_gmm_model_selection.png")),
       p_sel, width = 9, height = 4, dpi = 180)
cat("Model selection plot saved.\n")

# ---------------------------------------------------------------------------
# TRAJECTORY PLOTS FOR EACH FITTED K
# ---------------------------------------------------------------------------
age_min <- min(df$age, na.rm = TRUE)
age_max <- max(df$age, na.rm = TRUE)

pal <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628")

for (m in models) {
  k     <- m$ng
  traj  <- predict_trajectories(m, age_min, age_max)
  lbls  <- class_labels(m)

  p_traj <- ggplot(traj, aes(x = age, y = pds_pred, colour = class)) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual(values = setNames(pal[seq_len(k)], paste0("Class ", seq_len(k))),
                        labels = lbls) +
    scale_x_continuous(breaks = seq(floor(age_min), ceiling(age_max))) +
    labs(
      title    = paste0("GMM K=", k, " trajectory profiles — ", ds_name),
      subtitle = paste0("Posterior modal class assignment; ", class_pct(m)),
      x = "Age (years)", y = "PDS composite",
      colour = "Class"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())

  ggsave(
    file.path(out_dir, paste0(ds_name, "_gmm_trajectories_k", k, ".png")),
    p_traj, width = 8, height = 5, dpi = 180
  )
}
cat("Trajectory plots saved.\n")

# ---------------------------------------------------------------------------
# CLASS ASSIGNMENTS for ALL K (save the full posterior prob table)
# ---------------------------------------------------------------------------
all_assignments <- bind_rows(lapply(models, function(m) {
  pp <- m$pprob
  # pprob columns: id_num, class, prob_class1, prob_class2, ...
  prob_cols <- grep("^prob", names(pp), value = TRUE)
  pp %>%
    select(id_num, class, all_of(prob_cols)) %>%
    mutate(K = m$ng, dataset = ds_name)
}))

# Merge back the original string id
id_lookup <- df %>% distinct(id, id_num)
all_assignments <- all_assignments %>%
  left_join(id_lookup, by = "id_num") %>%
  select(dataset, id, id_num, K, class, everything())

write.csv(all_assignments,
          file.path(out_dir, paste0(ds_name, "_gmm_class_assignments.csv")),
          row.names = FALSE)
cat("Class assignments saved.\n")

cat("\nDone. All GMM outputs written to:", out_dir, "\n")
