#!/bin/bash
#$ -cwd
#$ -N growthonly
#$ -t 1-1
#$ -l h_rt=12:00:00
#$ -l h_data=32G
#$ -pe shared 4
#$ -j y
#$ -o logs/growthonly_$JOB_ID.$TASK_ID.log
#$ -m bea
#$ -M clarefmccann@g.ucla.edu

# ---------------------------------------------------------------------------
# Hoffman2 SGE job: growth-only convergence isolation test (Stage 4/5).
#
# Full longitudinal item data, random intercept + random slope + their
# correlation -- NO DIF, NO covariate impact. Uses lmnlfa-quad-tanhcor.stan
# (tanh-parameterized correlation + marker-item loading fix), meant to
# confirm the growth structure converges cleanly before returning to the
# full DIF/impact model (lmnlfa_hpc.R / lmnlfa-quad.stan).
#
# Task 1 -> female only for now (matches what's been validated locally).
# Change "-t 1-1" to "-t 1-2" (and SEXES below already supports it) to add
# male once this stage looks good.
#
# 12h is a starting guess -- this is a smaller model than the full one (no
# DIF, no impact), but local smoke tests showed 100% max-treedepth
# saturation before the marker-item fix, so watch early log output and
# adjust h_rt if needed.
#
# Before first submission on Hoffman2 (skip if already done for lmnlfa_hpc.sh):
#   1. Install cmdstanr (once, on login node):
#        module load R/4.2.2 gcc/10.2.0
#        Rscript -e "install.packages('cmdstanr', repos=c('https://mc-stan.org/r-packages/', getOption('repos')))"
#        Rscript -e "cmdstanr::install_cmdstan(dir='~/.cmdstan')"
#   2. Set CMDSTAN below to the exact versioned path that install_cmdstan created.
#   3. Verify DATA_DIR and OUT_DIR paths.
#   4. Create logs/ directory:  mkdir -p logs
#   5. Submit:  qsub lmnlfa_growth_only_tanhcor_hpc.sh
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

export DATA_DIR="/u/home/c/clarefmc/projects/dissertation/study1/outputs"
export OUT_DIR="/u/project/silvers/data/ABCD/cfm-dissertation-output/study1/outputs"

# Path to the CmdStan installation created by cmdstanr::install_cmdstan()
export CMDSTAN="${HOME}/.cmdstan/cmdstan-2.38.0"
export TMPDIR="/u/scratch/${USER}/tmp"
mkdir -p "$TMPDIR"
SCRIPT_DIR="$SGE_O_WORKDIR"

mkdir -p "${OUT_DIR}/lmnlfa_growth_only_tanhcor"
mkdir -p "${SCRIPT_DIR}/logs"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_MAX_VSIZE=100G

Rscript "${SCRIPT_DIR}/lmnlfa_growth_only_tanhcor.R" "$SX"

echo "============================================"
echo "End: $(date)"
echo "============================================"
