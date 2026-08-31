#!/bin/bash
#$ -cwd
#$ -N mnlfaXsec
#$ -t 2-2
#$ -l h_rt=15:00:00
#$ -l h_data=32G
#$ -pe shared 4
#$ -j y
#$ -o logs/mnlfaXsec_$JOB_ID.$TASK_ID.log
#$ -m bea
#$ -M clarefmccann@g.ucla.edu

# ---------------------------------------------------------------------------
# Hoffman2 SGE job: staged cross-sectional MNLFA (mnlfa_crosssectional_staged.R).
#
# Runs whichever stages are currently flagged TRUE (run_stage_A/B/C) at the
# top of mnlfa_crosssectional_staged.R -- a cached stage (its RDS already on
# disk) loads instantly regardless of its flag, so re-submitting after
# flipping the next stage's flag picks up right where you left off rather
# than re-fitting anything already done.
#
# Task 1 -> female only for now. Change "-t 1-1" to "-t 1-2" once you're
# ready to add male (SEXES array below already supports it).
#
# 15h is generous headroom given Stage A alone took ~3.75h; Stage B has more
# free parameters (48 DIF terms vs. 10 impact terms) so may run comparably
# or somewhat longer. Adjust if a stage's early log suggests otherwise.
#
# Before first submission (skip if already done for your other *_hpc.sh jobs):
#   1. Install cmdstanr (once, on login node):
#        module load R/4.2.2 gcc/10.2.0
#        Rscript -e "install.packages('cmdstanr', repos=c('https://mc-stan.org/r-packages/', getOption('repos')))"
#        Rscript -e "cmdstanr::install_cmdstan(dir='~/.cmdstan')"
#   2. Set CMDSTAN below to the exact versioned path that install_cmdstan created.
#   3. Verify DATA_DIR and OUT_DIR paths.
#   4. Create logs/ directory:  mkdir -p logs
#   5. Review run_stage_A/B/C flags in mnlfa_crosssectional_staged.R.
#   6. Submit:  qsub mnlfa_crosssectional_staged_hpc.sh
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

mkdir -p "${OUT_DIR}/mnlfa_crosssectional_staged"
mkdir -p "${SCRIPT_DIR}/logs"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_MAX_VSIZE=100G

Rscript "${SCRIPT_DIR}/mnlfa_crosssectional_staged.R" "$SX"

echo "============================================"
echo "End: $(date)"
echo "============================================"
