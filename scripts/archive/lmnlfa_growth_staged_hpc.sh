#!/bin/bash
#$ -cwd
#$ -N lmnlfaStaged
#$ -t 1-1
#$ -l h_rt=24:00:00
#$ -l h_data=32G
#$ -pe shared 4
#$ -j y
#$ -o logs/lmnlfaStaged_$JOB_ID.$TASK_ID.log
#$ -m bea
#$ -M clarefmccann@g.ucla.edu

# ---------------------------------------------------------------------------
# Hoffman2 SGE job: staged longitudinal MNLFA (lmnlfa_growth_staged.R).
#
# Runs whichever stage(s) are currently flagged TRUE (run_stage_A/B/C) at
# the top of lmnlfa_growth_staged.R -- a cached stage (its RDS already on
# disk) loads instantly regardless of its flag, so re-submitting after
# flipping the next stage's flag picks up right where you left off, same
# pattern as mnlfa_crosssectional_staged_hpc.sh.
#
# COST: this is a bigger model than lmnlfa_growth_informant.R. Stage A has
# a comparable parameter count to the informant-only model (impact terms
# instead of DIF terms). Stage B is heavier -- up to p x kdif x 2 = 48 DIF
# parameters (vs. 8 for informant-only) plus an xdif matrix-vector product
# computed for every one of ~187,000 (female, full sample) observations.
# Real runtime for this exact design is not yet benchmarked.
#
# 24h is the hard ceiling on the standard queue here regardless of what's
# requested (h_rt is a fixed queue attribute, not just a job setting), so
# there's no larger walltime to fall back on if a stage doesn't fit --
# STRONGLY prefer setting n_subsample in lmnlfa_growth_staged.R (e.g. 1500,
# the value that worked cleanly for the informant-only model) for the
# first diagnostic pass through ALL THREE stages before attempting any
# stage at full sample. Two prior full-sample submissions of the simpler
# informant-only model hit exactly this wall before subsampling fixed it.
#
# Before first submission:
#   1. Set run_stage_A <- TRUE (and B, C FALSE) in lmnlfa_growth_staged.R,
#      and set n_subsample for a first diagnostic pass.
#   2. Verify DATA_DIR/OUT_DIR paths below.
#   3. Create logs/ directory:  mkdir -p logs
#   4. Submit:  qsub lmnlfa_growth_staged_hpc.sh
#   5. Once Stage A completes, flip run_stage_A <- FALSE, run_stage_B <-
#      TRUE, resubmit -- Stage A's cached fit loads instantly. Repeat for
#      Stage C.
#
# Task 1 -> female only for now. Change "-t 1-1" to "-t 1-2" once you're
# ready to add male (SEXES array below already supports it).
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

export DATA_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"
export OUT_DIR="/u/home/c/clarefmc/projects/abcd-projs/dissertation/study1/outputs"

export CMDSTAN="${HOME}/.cmdstan/cmdstan-2.38.0"
export TMPDIR="/u/scratch/c/clarefmc/tmp"
mkdir -p "$TMPDIR"
SCRIPT_DIR="$SGE_O_WORKDIR"

mkdir -p "${OUT_DIR}/lmnlfa_growth_staged"
mkdir -p "${SCRIPT_DIR}/logs"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_MAX_VSIZE=100G

Rscript "${SCRIPT_DIR}/lmnlfa_growth_staged.R" "$SX"

echo "============================================"
echo "End: $(date)"
echo "============================================"
