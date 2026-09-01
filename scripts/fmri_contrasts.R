## fmri_contrasts.R
## Shared list of N-back Emotional task contrasts processed by the
## fmri_nback_foundation.R / fmri_nback_descriptives.R /
## fmri_puberty_association.R pipeline. Source this near the top of each so
## the three stages never drift out of sync on which contrasts exist.
##
## Names are the short tags used in ABCD variable names
## (mr_y_tfmri__nback__<tag>__{aseg,dsk}__*_beta) and in this pipeline's own
## output filenames; values are human-readable labels for plot titles.
##
## Two tags were corrected from initially-requested "emot"/"emotvntf" --
## the data dictionary uses "emo"/"emovntf" (confirmed via get_dd_abcd()).

fmri_contrasts <- c(
  emo     = "Emotion",
  emovntf = "Emotion vs. neutral face",
  ngfvntf = "Negative face vs. neutral face",
  psfvntf = "Positive face vs. neutral face",
  fvplc   = "Face vs. place"
)
