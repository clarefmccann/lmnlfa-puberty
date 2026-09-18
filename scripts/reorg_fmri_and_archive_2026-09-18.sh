#!/bin/bash
## reorg_fmri_and_archive_2026-09-18.sh
## One-time reorg: (1) move fMRI scripts out of study1 into a new sibling
## dissertation/study3/scripts/ directory, (2) archive the linear/informant
## lmnlfa lineage (superseded by the sigmoid model) into scripts/archive/.
##
## SCOPE, DELIBERATELY LIMITED:
##  - fMRI OUTPUT DATA (~502MB under outputs/*fmri*) is NOT moved here --
##    scripts only, per your call. Move it separately whenever you're ready
##    (same-filesystem mv is a server-side rename, not a real copy, so it's
##    cheap whenever you do it).
##  - The study1 -> lmnlfa-puberty directory rename is NOT part of this
##    script -- that's a separate, higher-blast-radius pass (it touches
##    ~27 active scripts/job files that hardcode "dissertation/study1" in
##    absolute paths) once this lower-risk pass is reviewed and committed.
##  - mnlfa_crosssectional_staged.R (+ its 3 stan files) and the
##    01_regression_characterization.R / 02_psychometrics.R /
##    03_gamms_hpc.R / 05_gmm_hpc.R / 05_slides.qmd pipeline are NOT
##    touched -- fmri_puberty_association.R reads mnlfa_crosssectional_staged's
##    OUTPUT factor-score CSVs (not the script itself), and 05_slides.qmd
##    reads output from 01/02, so none of these were confirmed safe to move.
##
## GIT NOTE: dissertation/study3 is OUTSIDE the study1 git repo (study1/
## IS the repo root -- study3 is a new sibling directory, not part of this
## repo). Because of that, the fMRI files use plain `mv` (not `git mv`,
## which only works for renames WITHIN one repo) -- they're staged in
## study1's git index as deletions by the `git add -A` at the end. Whether
## study3 becomes its own git repo is entirely your call, separate from
## this script. The archive moves (linear/informant lineage), by contrast,
## stay INSIDE study1, so those use `git mv` to preserve rename history.
##
## SAFETY: DRY_RUN=1 by default -- prints every command without running it.
## Review the printed plan, then re-run with DRY_RUN=0 to actually execute.
## Nothing is committed by this script -- review `git status` / `git diff
## --stat` yourself and commit when you're satisfied.
##
## Usage:
##   cd ~/projects/abcd-projs/dissertation/study1   # repo root
##   bash scripts/reorg_fmri_and_archive_2026-09-18.sh          # dry run
##   DRY_RUN=0 bash scripts/reorg_fmri_and_archive_2026-09-18.sh  # for real

set -euo pipefail

DRY_RUN="${DRY_RUN:-1}"

STUDY1_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DISSERTATION_ROOT="$(dirname "$STUDY1_ROOT")"
STUDY3_ROOT="$DISSERTATION_ROOT/study3"

echo "STUDY1_ROOT  = $STUDY1_ROOT"
echo "STUDY3_ROOT  = $STUDY3_ROOT (will be created)"
echo "DRY_RUN      = $DRY_RUN"
echo

if [ ! -d "$STUDY1_ROOT/.git" ]; then
  echo "ERROR: $STUDY1_ROOT does not look like the study1 git repo root (.git not found). Aborting." >&2
  exit 1
fi

cd "$STUDY1_ROOT"

run() {
  echo "+ $*"
  if [ "$DRY_RUN" != "1" ]; then
    "$@"
  fi
}

echo "=== STEP 1: create study3 (sibling to study1, outside this git repo) ==="
run mkdir -p "$STUDY3_ROOT/scripts"

echo
echo "=== STEP 2: move fMRI scripts to study3 (plain mv -- crosses repo boundary) ==="
FMRI_FILES=(
  "scripts/fmri_contrasts.R"
  "scripts/fmri_nback_descriptives.R"
  "scripts/fmri_nback_foundation.R"
  "scripts/fmri_puberty_association.R"
  "scripts/fmri_puberty_longitudinal.R"
  "scripts/fmri_puberty_longitudinal_hpc.sh"
)
for f in "${FMRI_FILES[@]}"; do
  if [ -f "$f" ]; then
    run mv "$f" "$STUDY3_ROOT/scripts/$(basename "$f")"
  else
    echo "  (skip, not found: $f)"
  fi
done

echo
echo "=== STEP 3: archive linear/informant lmnlfa lineage (git mv -- stays inside study1) ==="
ARCHIVE_FILES=(
  "scripts/lmnlfa_growth_staged.R"
  "scripts/lmnlfa_growth_staged_hpc.sh"
  "scripts/lmnlfa_growth_informant.R"
  "scripts/lmnlfa_growth_informant_hpc.sh"
  "scripts/stan/lmnlfa-growth-difscreen.stan"
  "scripts/stan/lmnlfa-growth-final.stan"
  "scripts/stan/lmnlfa-growth-impact.stan"
  "scripts/stan/lmnlfa-linear-tanhcor-informant.stan"
  "scripts/stan/lmnlfa-quad.stan"
  "scripts/stan/lmnlfa-quad-tanhcor.stan"
)
for f in "${ARCHIVE_FILES[@]}"; do
  if [ -f "$f" ]; then
    run git mv "$f" "scripts/archive/$(basename "$f")"
  else
    echo "  (skip, not found: $f)"
  fi
done

echo
echo "=== STEP 4: stage the fMRI deletions in study1's git index ==="
run git add -A -- scripts/fmri_contrasts.R scripts/fmri_nback_descriptives.R \
  scripts/fmri_nback_foundation.R scripts/fmri_puberty_association.R \
  scripts/fmri_puberty_longitudinal.R scripts/fmri_puberty_longitudinal_hpc.sh

echo
echo "=== DONE (dry run: $DRY_RUN) ==="
echo "Review before committing:"
echo "  git status"
echo "  git diff --cached --stat"
echo
echo "MANUAL FOLLOW-UP NEEDED (not automated -- a real code change, not a file move):"
echo "  $STUDY3_ROOT/scripts/fmri_puberty_association.R reads mnlfa_crosssectional_staged's"
echo "  OUTPUT factor-score CSVs (outputs/mnlfa_crosssectional_staged/factor_scores_*.csv)."
echo "  Those stay in study1's outputs/, not study3's. After this move, that script's"
echo "  data_dir (or eta_path_f/eta_path_m) needs to point back at study1's outputs dir"
echo "  explicitly, not wherever study3's own DATA_DIR/OUT_DIR resolves to."
echo
echo "NOT done by this script (separate, later passes):"
echo "  - Moving the ~502MB of fMRI output data under outputs/*fmri*"
echo "  - Renaming study1/ -> lmnlfa-puberty/ and updating ~27 hardcoded paths"
