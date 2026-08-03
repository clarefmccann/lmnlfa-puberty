#!/bin/bash
#$ -cwd
#$ -N gamm_fs
#$ -t 1-4
#$ -l h_rt=8:00:00
#$ -l h_data=32G
#$ -pe shared 4
#$ -j y
#$ -o logs/gamm_fs_$JOB_ID.$TASK_ID.log
#$ -m bea
#$ -M clarefmccann@g.ucla.edu
    
# ---------------------------------------------------------------------------
# Hoffman2 SGE array job: factor-smooth GAMMs for PDS composite
# Submits 4 tasks (one per sex × reporter dataset), each using 4 cores.
#
# Before submitting:
#   1. Set DATA_DIR below to wherever the *_long.csv files live on Hoffman2
#   2. Verify the R module name with:  module avail R
#   3. Create the logs/ directory:    mkdir -p logs
#   4. Submit:  qsub 03_gamms_hpc.sh
# ---------------------------------------------------------------------------

# Initialize the module system (required in SGE jobs on Hoffman2)
. /u/local/Modules/default/init/bash

# Map array task ID → dataset name
DATASETS=(male_parent male_youth female_parent female_youth)
DS=${DATASETS[$((SGE_TASK_ID - 1))]}

echo "============================================"
echo "Task $SGE_TASK_ID  →  dataset: $DS"
echo "Nodes: $HOSTNAME   Cores: $NSLOTS"
echo "Start: $(date)"
echo "============================================"

# Load R — check available versions with: module avail R
module load R/4.2.2
# Load GCC to provide a newer libstdc++ (GLIBCXX_3.4.21+) required by tidyr
module load gcc/10.2.0

# Set data directory (edit this to match where your puberty CSVs live on H2)
export DATA_DIR="/u/project/silvers/data/ABCD/cfm-dissertation-output/study1/outputs"
export OUT_DIR="/u/project/silvers/data/ABCD/cfm-dissertation-output/study1/outputs"

# SGE copies the script to a temp location, so dirname "$0" doesn't work.
# $SGE_O_WORKDIR is the directory qsub was called from (set by -cwd).
SCRIPT_DIR="$SGE_O_WORKDIR"

mkdir -p "${OUT_DIR}/gamm"

export R_MAX_VSIZE=100Gb
Rscript "${SCRIPT_DIR}/03_gamms_hpc.R" "$DS" "$NSLOTS"

echo "============================================"
echo "End: $(date)"
echo "============================================"
