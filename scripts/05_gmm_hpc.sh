#!/bin/bash
#$ -cwd
#$ -N gmm_traj
#$ -t 1-4
#$ -l h_rt=6:00:00
#$ -l h_data=16G
#$ -pe shared 1
#$ -j y
#$ -o logs/gmm_traj_$JOB_ID.$TASK_ID.log
#$ -m bea
#$ -M clarefmccann@g.ucla.edu

# ---------------------------------------------------------------------------
# Hoffman2 SGE array job: Growth Mixture Models (PDS composite)
#
# Task 1 → female_parent   Task 3 → male_parent
# Task 2 → female_youth    Task 4 → male_youth
#
# lcmm::hlme() is single-threaded so we request 1 slot.
# 6h wall-clock: K=1..6 × 10 random starts each ≈ 1–3h per dataset.
# ---------------------------------------------------------------------------

. /u/local/Modules/default/init/bash

DATASETS=(female_parent female_youth male_parent male_youth)
DS=${DATASETS[$((SGE_TASK_ID - 1))]}

echo "============================================"
echo "Task $SGE_TASK_ID  →  dataset: $DS"
echo "Host: $HOSTNAME"
echo "Start: $(date)"
echo "============================================"

module load R/4.2.2
module load gcc/10.2.0

export DATA_DIR="/u/project/silvers/data/ABCD/ABCD-release-6.0/cfm/physical-health/puberty"
export OUT_DIR="/u/project/silvers/data/ABCD/cfm-dissertation-output/study1/outputs"

SCRIPT_DIR="$SGE_O_WORKDIR"

mkdir -p "${OUT_DIR}/gmm_trajectories"
mkdir -p "${SCRIPT_DIR}/logs"

export R_MAX_VSIZE=32G

# Args: dataset_name  max_k  n_starts
Rscript "${SCRIPT_DIR}/05_gmm_hpc.R" "$DS" 6 10

echo "============================================"
echo "End: $(date)"
echo "============================================"
