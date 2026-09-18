#!/bin/bash
#$ -cwd
#$ -N lmnlfaSigSSP
#$ -t 1-2:1
#$ -l h_rt=24:00:00
#$ -l h_data=32G
#$ -pe shared 4
#$ -j y
#$ -o logs/lmnlfaSigSSP_$JOB_ID.$TASK_ID.log
#$ -m bea
#$ -M clarefmccann@g.ucla.edu

# ---------------------------------------------------------------------------
# Hoffman2 SGE job: spike-and-slab DIF-selection alternative
# (lmnlfa_growth_sigmoid_staged_ssp.R), for comparison against the
# BH-FDR-based lmnlfa_growth_sigmoid_staged.R on the SAME growth curve and
# SAME Stage A impact estimates.
#
# PREREQUISITE: lmnlfa_growth_sigmoid_staged.R must have ALREADY completed
# Stage A (run_stage_A <- TRUE) for this sex -- this script LOADS that
# cached fit rather than refitting it (Stage A doesn't involve DIF at all,
# so it's identical either way). If fits/fitA_impact_<sex>_n1500.rds
# doesn't exist yet under lmnlfa_growth_sigmoid_staged/, this job will fail
# immediately with a clear message telling you to run that first.
#
# THIS IS ALSO A NEW, UNBENCHMARKED MODEL/METHOD COMBINATION. The
# spike-and-slab prior (Laplace x Beta(.5,.5) mixture) has a more
# irregular, near-discontinuous posterior geometry than the BH-FDR model's
# plain diffuse priors, and is more prone to divergences/poor mixing --
# adapt_delta is already bumped to 0.97 (vs 0.95) in the R script for this
# reason, but watch Stage B's diagnostics closely regardless: divergences,
# max-treedepth hits, and whether r_incl's posterior distribution actually
# separates into "near 0.5" (not retained) vs "near 1" (retained) values
# rather than sitting uniformly near 0.5 for everything (which would mean
# the lasso penalty is too tight, or too loose, for this data's effect
# sizes -- see the SSP hyperparameter comments at the top of the R script).
#
# Before first submission:
#   1. Confirm lmnlfa_growth_sigmoid_staged.R's Stage A has already
#      completed for this sex (check for fits/fitA_impact_<sex>_n1500.rds
#      under lmnlfa_growth_sigmoid_staged/).
#   2. Confirm run_stage_B <- TRUE (C FALSE) in
#      lmnlfa_growth_sigmoid_staged_ssp.R for the first diagnostic pass,
#      and review ssp_scale/phi_shape/phi_rate/ssp_threshold.
#   3. Verify DATA_DIR/OUT_DIR paths below (must match the BH-FDR job's,
#      since Stage A is read from there).
#   4. Create logs/ directory:  mkdir -p logs
#   5. Submit from scripts/:  cd scripts && qsub lmnlfa_growth_sigmoid_staged_ssp_hpc.sh
#   6. Once Stage B looks healthy, flip run_stage_B <- FALSE,
#      run_stage_C <- TRUE, resubmit -- Stage B's cached fit loads
#      instantly. Stage C reuses the SAME final.stan as the BH-FDR
#      pipeline, just with the SSP-derived DIF pattern.
#
# Both tasks (female + male) enabled by default, since each sex's Stage A
# is independently already cached (unlike the BH-FDR script's original
# female-only start) -- if that's not true for both sexes yet, change
# "-t 1-2" to "-t 1-1" (female) or "-t 2-2" (male) as needed.
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

mkdir -p "${OUT_DIR}/lmnlfa_growth_sigmoid_staged_ssp"
mkdir -p "${SCRIPT_DIR}/logs"

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_MAX_VSIZE=100G

Rscript "${SCRIPT_DIR}/lmnlfa_growth_sigmoid_staged_ssp.R" "$SX"

echo "============================================"
echo "End: $(date)"
echo "============================================"
