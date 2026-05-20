#!/bin/bash
#$ -cwd
#$ -N lmnlfa
#$ -t 1-2
#$ -l h_rt=20:00:00
#$ -l h_data=32G
#$ -pe shared 4
#$ -j y
#$ -o logs/lmnlfa_$JOB_ID.$TASK_ID.log
#$ -m bea
#$ -M clarefmccann@g.ucla.edu

# ---------------------------------------------------------------------------
# Hoffman2 SGE array job: Longitudinal MNLFA for puberty
#
# Task 1 → female   Task 2 → male
# 4 chains run in parallel (one core each, hence -pe shared 4).
# 36h wall-clock is conservative; expect ~20h per sex group.
#
# Before first submission on Hoffman2:
#   1. Install cmdstanr (once, on login node):
#        module load R/4.2.2 gcc/10.2.0
#        Rscript -e "install.packages('cmdstanr', repos=c('https://mc-stan.org/r-packages/', getOption('repos')))"
#        Rscript -e "cmdstanr::install_cmdstan(dir='~/.cmdstan')"
#   2. Set CMDSTAN below to the exact versioned path that install_cmdstan created.
#   3. Verify DATA_DIR and OUT_DIR paths.
#   4. Create logs/ directory:  mkdir -p logs
#   5. Submit:  qsub lmnlfa_hpc.sh
# ---------------------------------------------------------------------------

. /u/local/Modules/default/init/bash

SEXES=(female male)
SX=${SEXES[$((SGE_TASK_ID - 1))]}

echo "============================================"
echo "Task $SGE_TASK_ID  →  sex: $SX"
echo "Host: $HOSTNAME   Cores: $NSLOTS"
echo "Start: $(date)"
echo "============================================"

module load R/4.2.2
module load gcc/10.2.0

export DATA_DIR="/u/project/silvers/data/ABCD/ABCD-release-6.0/cfm/physical-health/puberty"
export OUT_DIR="/u/project/silvers/data/ABCD/cfm-dissertation-output/study1/outputs"

# Path to the CmdStan installation created by cmdstanr::install_cmdstan()
# Update the version number if you install a newer CmdStan.
export CMDSTAN="${HOME}/.cmdstan/cmdstan-2.38.0"
export TMPDIR="/u/scratch/${USER}/tmp"
mkdir -p "$TMPDIR"
SCRIPT_DIR="$SGE_O_WORKDIR"

mkdir -p "${OUT_DIR}/lmnlfa"
mkdir -p "${SCRIPT_DIR}/logs"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_MAX_VSIZE=100G

Rscript "${SCRIPT_DIR}/lmnlfa_hpc.R" "$SX"

echo "============================================"
echo "End: $(date)"
echo "============================================"
