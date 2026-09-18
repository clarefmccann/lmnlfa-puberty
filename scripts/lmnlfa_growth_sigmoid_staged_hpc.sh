#!/bin/bash
#$ -cwd
#$ -N lmnlfaSigStg
#$ -t 1-2:1
#$ -l h_rt=24:00:00
#$ -l h_data=32G
#$ -pe shared 4
#$ -j y
#$ -o logs/lmnlfaSigStg_$JOB_ID.$TASK_ID.log
#$ -m bea
#$ -M clarefmccann@g.ucla.edu

# ---------------------------------------------------------------------------
# Hoffman2 SGE job: staged SIGMOIDAL longitudinal MNLFA
# (lmnlfa_growth_sigmoid_staged.R). Same staging pattern as
# lmnlfa_growth_staged_hpc.sh, applied to the bounded logistic growth curve
# (fixed floor=1/ceiling=5, person-specific rate + inflection age) instead
# of the linear one -- a cached stage (its RDS already on disk) loads
# instantly regardless of its flag, so re-submitting after flipping the
# next stage's flag picks up right where you left off.
#
# THIS IS A NEW, UNBENCHMARKED MODEL. It has never been run end-to-end --
# only compiled and dry-run on the data-build step. Treat Stage A's first
# run as a genuine diagnostic pass, not a known quantity: watch divergences,
# max-treedepth hits, and Rhat/ESS on mu_logk/phi_logk/mu_alpha/phi_alpha
# closely before trusting it enough to move on to Stage B.
#
# COST: parameter count and per-observation work are comparable to
# lmnlfa_growth_staged.R's linear model (same p x kdif DIF term count in
# Stage B, same xdif matrix-vector product per observation) -- the only
# difference is the per-observation growth-curve formula itself (a couple
# of extra transcendental ops: exp() for k_i, inv_logit() for the sigmoid),
# which is unlikely to change runtime much. Still, treat the linear model's
# runtime history as a rough guide only, not a guarantee.
#
# 24h is the hard ceiling on the standard queue here regardless of what's
# requested (h_rt is a fixed queue attribute, not just a job setting), so
# there's no larger walltime to fall back on if a stage doesn't fit --
# STRONGLY prefer the current n_subsample <- 1500 default in
# lmnlfa_growth_sigmoid_staged.R (the value that worked cleanly for the
# linear staged model) for this first diagnostic pass through all three
# stages before attempting any stage at full sample.
#
# Before first submission:
#   1. Confirm run_stage_A <- TRUE (B, C FALSE) and n_subsample <- 1500 in
#      lmnlfa_growth_sigmoid_staged.R for the first diagnostic pass.
#   2. Verify DATA_DIR/OUT_DIR paths below.
#   3. Create logs/ directory:  mkdir -p logs
#   4. Submit from scripts/:  cd scripts && qsub lmnlfa_growth_sigmoid_staged_hpc.sh
#   5. Once Stage A completes and looks healthy, flip run_stage_A <- FALSE,
#      run_stage_B <- TRUE, resubmit -- Stage A's cached fit loads
#      instantly. Repeat for Stage C.
#
# Task 1 -> female only for now (SEXES array below already supports male;
# change "-t 1-1" to "-t 1-2" once you're ready to add it -- but see the
# note above about run_stage_A/B/C being script-level, not per-sex, flags:
# don't add male until female's own staged run is fully done, or a male
# submission sharing the same flags could pick up the wrong stage for
# whichever sex is still in flight).
# ---------------------------------------------------------------------------

. /u/local/Modules/default/init/bash

SEXES=(female male)
SX=${SEXES[$((SGE_TASK_ID - 1))]}

echo "============================================"
echo "Task $SGE_TASK_ID  ->  sex: $SX"
echo "Host: $HOSTNAME   Cores: $NSLOTS"
echo "Start: $(date)"
echo "============================================"

module load R/4.2.2
module load gcc/10.2.0

export DATA_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/data"
export OUT_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"

export CMDSTAN="${HOME}/.cmdstan/cmdstan-2.38.0"
export TMPDIR="/u/scratch/c/clarefmc/tmp"
mkdir -p "$TMPDIR"
SCRIPT_DIR="$SGE_O_WORKDIR"

mkdir -p "${OUT_DIR}/lmnlfa_growth_sigmoid_staged"
mkdir -p "${SCRIPT_DIR}/logs"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_MAX_VSIZE=100G

Rscript "${SCRIPT_DIR}/lmnlfa_growth_sigmoid_staged.R" "$SX"

echo "============================================"
echo "End: $(date)"
echo "============================================"
