#!/bin/bash
#$ -cwd
#$ -N fmriLongit
#$ -l h_rt=12:00:00
#$ -l h_data=8G
#$ -pe shared 4
#$ -j y
#$ -o logs/fmriLongit_$JOB_ID.log
#$ -m bea
#$ -M clarefmccann@g.ucla.edu

# ---------------------------------------------------------------------------
# Hoffman2 SGE job: fmri_puberty_longitudinal.R
#
# GAMM (mgcv) association between N-back betas and pds_comp, fit separately
# per sex, looped across every contrast in fmri_contrasts.R and every
# region (~98 regions x 2 sexes per contrast for the sweep, plus 4 richer
# key-region/global models x 2 sexes per contrast) -- around 1,000 GAM fits
# total. That's real compute, not a quick interactive script, hence this
# job wrapper instead of running it on a login node.
#
# mgcv ships as part of base R (one of R's own "recommended" packages), so
# unlike lme4/nloptr this needs no extra package install on Hoffman2.
#
# No task array / arguments: the R script loops over every contrast and
# both sexes internally in one run.
#
# 12h / 8G x 4 slots is a first guess -- this job has not yet been
# benchmarked at full ABCD scale (only tested locally on a small synthetic
# dataset). The R script parallelizes across NSLOTS via parallel::mclapply
# (fork-based, Unix only -- fine on Hoffman2) for the ~196 independent
# region x sex fits per contrast, so raising -pe shared here directly buys
# wall-clock speed, not just memory headroom. Fitting and plotting are
# deliberately kept separate (all GAM fitting happens inside the forked
# workers; every ggplot/ggsave call runs afterward in the main process),
# since graphics-device code is not reliably fork-safe on every platform.
# Watch the early log: if it's clearly going to blow past 12h, kill it and
# increase h_rt rather than let it hit a walltime kill, and note the
# actual per-contrast timing here for next time.
#
# Before first submission:
#   1. Verify DATA_DIR/OUT_DIR paths below.
#   2. Create logs/ directory:  mkdir -p logs
#   3. Confirm all_fmri_nback_<contrast>_long.csv / fmri_nback_<contrast>_
#      region_lookup.csv exist for every contrast (fmri_nback_foundation.R)
#      and all_long.csv exists (00_data_foundation.R).
#   4. Submit:  qsub fmri_puberty_longitudinal_hpc.sh
# ---------------------------------------------------------------------------

. /u/local/Modules/default/init/bash

echo "============================================"
echo "Host: $HOSTNAME   Cores: $NSLOTS"
echo "Start: $(date)"
echo "============================================"

module load R/4.2.2
module load gcc/10.2.0

export DATA_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"
export OUT_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"

SCRIPT_DIR="$SGE_O_WORKDIR"
mkdir -p "${SCRIPT_DIR}/logs"

export R_MAX_VSIZE=100Gb

Rscript "${SCRIPT_DIR}/fmri_puberty_longitudinal.R"

echo "============================================"
echo "End: $(date)"
echo "============================================"
