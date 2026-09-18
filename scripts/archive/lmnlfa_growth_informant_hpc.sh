#!/bin/bash
#$ -cwd
#$ -N lmnlfaInf
#$ -t 1-1
#$ -l h_rt=24:00:00
#$ -l h_data=32G
#$ -pe shared 4
#$ -j y
#$ -o logs/lmnlfaInf_$JOB_ID.$TASK_ID.log
#$ -m bea
#$ -M clarefmccann@g.ucla.edu

# ---------------------------------------------------------------------------
# Hoffman2 SGE job: growth model + informant DIF (lmnlfa_growth_informant.R).
#
# First real run of this restructured model (4 shared items, informant as a
# DIF covariate, growth structure carried over from the validated
# growth-only model). No smoke test was run beforehand -- data build and
# Stan compile were verified locally, but the actual sampling behavior at
# full scale is untested. Watch the early log for the failure signatures
# from earlier debugging: divergences, 100% max-treedepth, or extreme
# ordered_logistic cut-point values in the first few hundred warmup
# iterations. Clean early behavior there is a good sign; a repeat of any of
# those symptoms means stopping and diagnosing rather than waiting it out.
#
# 24h is a generous margin given real uncertainty about runtime -- this
# dataset (244K obs, 5,656 people) is somewhat larger than the growth-only
# model's, and this exact configuration has never been benchmarked. Watch
# early log pace and adjust future submissions if it's clearly over- or
# under-provisioned.
#
# run_fit is FALSE by default at the top of lmnlfa_growth_informant.R --
# flip it to TRUE before submitting, or this job will build the data,
# compile, and then stop with an informative error instead of sampling.
#
# Task 1 -> female only for now. Change "-t 1-1" to "-t 1-2" once this
# looks good (SEXES array below already supports it).
#
# Before first submission (skip if already done for your other *_hpc.sh jobs):
#   1. Install cmdstanr (once, on login node):
#        module load R/4.2.2 gcc/10.2.0
#        Rscript -e "install.packages('cmdstanr', repos=c('https://mc-stan.org/r-packages/', getOption('repos')))"
#        Rscript -e "cmdstanr::install_cmdstan(dir='~/.cmdstan')"
#   2. Set CMDSTAN below to the exact versioned path that install_cmdstan created.
#   3. Verify DATA_DIR and OUT_DIR paths.
#   4. Create logs/ directory:  mkdir -p logs
#   5. Set run_fit <- TRUE in lmnlfa_growth_informant.R.
#   6. Submit:  qsub lmnlfa_growth_informant_hpc.sh
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

# Path to the CmdStan installation created by cmdstanr::install_cmdstan()
export CMDSTAN="${HOME}/.cmdstan/cmdstan-2.38.0"
export TMPDIR="/u/scratch/c/clarefmc/tmp"
mkdir -p "$TMPDIR"
SCRIPT_DIR="$SGE_O_WORKDIR"

mkdir -p "${OUT_DIR}/lmnlfa_growth_informant"
mkdir -p "${SCRIPT_DIR}/logs"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_MAX_VSIZE=100G

Rscript "${SCRIPT_DIR}/lmnlfa_growth_informant.R" "$SX"

echo "============================================"
echo "End: $(date)"
echo "============================================"
