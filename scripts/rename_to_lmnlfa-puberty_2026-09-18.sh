#!/bin/bash
## rename_to_lmnlfa-puberty_2026-09-18.sh
## Renames dissertation/study1 -> dissertation/lmnlfa-puberty (the directory
## that IS the git repo root), then updates every active script/job file
## that hardcodes "dissertation/study1" in an absolute path so the pipeline
## keeps working afterward.
##
## GIT NOTE: this is safe for git history. git tracks paths RELATIVE TO
## the repo root (wherever .git lives) -- it does not care what the
## containing checkout directory is named. Renaming the folder is a plain
## filesystem mv, .git moves with it intact, and `git log`/`git blame`/
## history are completely unaffected. This script does NOT touch git at
## all beyond that implicit consequence -- nothing is staged or committed.
##
## THE REMOTE (GitHub/GitLab/etc.) IS A SEPARATE, MANUAL STEP, NOT COVERED
## HERE -- renaming the remote repo itself is done on the hosting
## platform's own settings page:
##   GitHub: repo -> Settings -> repository name -> Rename
##   GitLab: repo -> Settings -> General -> Advanced -> Rename/transfer
## Most platforms auto-redirect the old URL for a while, but update your
## local remote afterward regardless:
##   git remote set-url origin <new-url>
##   git remote -v   # confirm
## Do this whenever you're ready -- independent of, and can happen before
## or after, this script.
##
## SCOPE: only ACTIVE files are touched (main scripts/, scripts/exploration/,
## and the 4 cross-referencing files in study3/scripts/). scripts/archive/
## and scripts/logs/ are left untouched -- they're dead/historical and
## don't need to keep resolving paths correctly. Also NOT touched:
## scripts/exploration/03_gamms_hpc.sh and 05_gmm_hpc.sh's OUT_DIR (lab
## group storage at /u/project/silvers/..., flagged previously as unclear
## and out of scope again here), and their DATA_DIR's pre-existing
## "abcd-projs"-missing typo -- same reasoning as the previous reorg pass.
##
## SAFETY: DRY_RUN=1 by default. Review the printed plan, then re-run with
## DRY_RUN=0 to actually execute. Nothing is committed -- review `git
## status` / `git diff --stat` yourself afterward and commit when ready.
##
## Usage:
##   cd ~/projects/abcd-projs/dissertation   # the PARENT of study1
##   bash study1/scripts/rename_to_lmnlfa-puberty_2026-09-18.sh          # dry run
##   DRY_RUN=0 bash study1/scripts/rename_to_lmnlfa-puberty_2026-09-18.sh  # for real

set -euo pipefail

DRY_RUN="${DRY_RUN:-1}"

OLD_NAME="study1"
NEW_NAME="lmnlfa-puberty"

# resolve paths before any rename happens
SCRIPT_PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # .../dissertation/study1
DISSERTATION_ROOT="$(dirname "$SCRIPT_PATH")"                    # .../dissertation
OLD_ROOT="$DISSERTATION_ROOT/$OLD_NAME"
NEW_ROOT="$DISSERTATION_ROOT/$NEW_NAME"
STUDY3_ROOT="$DISSERTATION_ROOT/study3"

echo "OLD_ROOT = $OLD_ROOT"
echo "NEW_ROOT = $NEW_ROOT"
echo "DRY_RUN  = $DRY_RUN"
echo

if [ ! -d "$OLD_ROOT/.git" ]; then
  echo "ERROR: $OLD_ROOT does not look like the study1 git repo root (.git not found). Aborting." >&2
  exit 1
fi
if [ -e "$NEW_ROOT" ]; then
  echo "ERROR: $NEW_ROOT already exists. Aborting to avoid clobbering it." >&2
  exit 1
fi

run() {
  echo "+ $*"
  if [ "$DRY_RUN" != "1" ]; then
    "$@"
  fi
}

echo "=== STEP 1: rename the repo directory itself ==="
run mv "$OLD_ROOT" "$NEW_ROOT"

# In a real run, STEP 1 already renamed the directory by the time we get
# here, so files live under NEW_ROOT. In a dry run, STEP 1 was only
# printed (not executed), so files still live under OLD_ROOT -- check
# existence against whichever one is actually real right now, so the dry
# run's "not found" skips are accurate instead of universally false.
CHECK_ROOT="$NEW_ROOT"
if [ "$DRY_RUN" = "1" ]; then
  CHECK_ROOT="$OLD_ROOT"
fi

echo
echo "=== STEP 2: update hardcoded paths inside the renamed repo ==="
# files identified by: grep -rl "dissertation/study1" <repo>, excluding
# archive/, logs/, .git/, and this script + the prior reorg script
STUDY1_FILES=(
  "scripts/exploration/00_data_foundation.R"
  "scripts/exploration/01_regression_characterization.R"
  "scripts/exploration/02_psychometrics.R"
  "scripts/exploration/03_gamms_hpc.R"
  "scripts/exploration/03_gamms_hpc.sh"
  "scripts/exploration/05_gmm_hpc.R"
  "scripts/exploration/05_gmm_hpc.sh"
  "scripts/exploration/05_slides.qmd"
  "scripts/exploration/viz_sample_timeline.R"
  "scripts/exploration/viz_sample_timeline_modeled.R"
  "scripts/lmnlfa_growth_sigmoid_report.Rmd"
  "scripts/lmnlfa_growth_sigmoid_staged.R"
  "scripts/lmnlfa_growth_sigmoid_staged_hpc.sh"
  "scripts/lmnlfa_growth_sigmoid_staged_ssp.R"
  "scripts/lmnlfa_growth_sigmoid_staged_ssp_hpc.sh"
  "scripts/lmnlfa_sigmoid_premeeting_diagnostics.R"
  "scripts/lmnlfa_sigmoid_premeeting_report.Rmd"
  "scripts/mnlfa_crosssectional_staged.R"
  "scripts/mnlfa_crosssectional_staged_hpc.sh"
  "scripts/viz_raw_vs_model_trajectories.R"
)
for f in "${STUDY1_FILES[@]}"; do
  check_path="$CHECK_ROOT/$f"
  real_path="$NEW_ROOT/$f"
  if [ -f "$check_path" ]; then
    run sed -i "s#dissertation/${OLD_NAME}#dissertation/${NEW_NAME}#g" "$real_path"
  else
    echo "  (skip, not found: $f)"
  fi
done

echo
echo "=== STEP 3: update the 4 cross-referencing files in study3 (outside the renamed repo) ==="
STUDY3_FILES=(
  "scripts/fmri_nback_descriptives.R"
  "scripts/fmri_puberty_longitudinal_hpc.sh"
  "scripts/fmri_puberty_association.R"
  "scripts/fmri_puberty_longitudinal.R"
)
for f in "${STUDY3_FILES[@]}"; do
  path="$STUDY3_ROOT/$f"
  if [ -f "$path" ]; then
    run sed -i "s#dissertation/${OLD_NAME}#dissertation/${NEW_NAME}#g" "$path"
  else
    echo "  (skip, not found: $f)"
  fi
done

echo
echo "=== DONE (dry run: $DRY_RUN) ==="
echo "Review before committing (from inside the renamed repo):"
echo "  cd $NEW_ROOT"
echo "  git status"
echo "  git diff --stat"
echo
echo "NOT done by this script (separate, manual, your call on timing):"
echo "  - Renaming the remote (GitHub/GitLab settings page), then:"
echo "      git remote set-url origin <new-url>"
echo "  - scripts/exploration/03_gamms_hpc.sh / 05_gmm_hpc.sh's DATA_DIR (pre-existing"
echo "    typo) and OUT_DIR (lab group storage) -- left exactly as before, still unclear"
echo "  - Deleting/archiving this script and reorg_fmri_and_archive_2026-09-18.sh"
echo "    now that both have served their purpose"
