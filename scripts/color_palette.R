## color_palette.R
## Shared "sunset" palette (slate -> cream -> honey -> brick -> bordeaux) for
## every plot produced by the mnlfa*/lmnlfa* scripts, so figures look
## consistent across the cross-sectional and longitudinal analyses.
## Source this near the top of a script, after loading ggplot2.

sunset_anchors <- c(
  slate    = "#335C67",
  cream    = "#FFF3B0",
  honey    = "#E09F3E",
  brick    = "#9E2A2B",
  bordeaux = "#540B0E"
)

# n evenly-spaced colors interpolated along the full gradient -- for
# sequential series (e.g. one color per wave/timepoint)
sunset_ramp <- function(n) {
  unname(colorRampPalette(sunset_anchors)(n))
}

# single accent color: scatter points, individual trajectory lines
pal_primary <- unname(sunset_anchors["slate"])
# translucent fill paired with pal_primary (e.g. under a mean-trajectory line)
pal_primary_fill <- unname(sunset_anchors["honey"])

# two-group contrasts (Parent vs Youth, Stage A vs Stage C, etc.)
pal_two <- unname(sunset_anchors[c("slate", "brick")])

# three-group contrasts with a neutral reference category (e.g. grand mean +
# two race categories)
pal_three_ref <- c("grey45", unname(sunset_anchors[c("slate", "brick")]))

# 4 well-separated colors for MCMC chain trace/density plots -- the anchor
# hues themselves (skipping the pale "cream", too low-contrast on white)
pal_chains <- unname(sunset_anchors[c("slate", "honey", "brick", "bordeaux")])

# ---------------------------------------------------------------------------
# Shape/linetype companions -- every categorical colour mapping above should
# be paired with one of these on the same grouping variable, so plots don't
# rely on hue alone (colorblind accessibility). Points get pal_shapes_*,
# lines get pal_linetypes_*, sized to match the corresponding pal_* above.
# ---------------------------------------------------------------------------
pal_linetypes_two   <- c("solid", "dashed")
pal_linetypes_three <- c("solid", "dashed", "dotted")
pal_linetypes_chains <- c("solid", "dashed", "dotted", "dotdash")

pal_shapes_two   <- c(16, 17)
pal_shapes_three <- c(16, 17, 15)
pal_shapes_chains <- c(16, 17, 15, 18)

# ---------------------------------------------------------------------------
# Study wave (bl, fu1..fu6) -- shared across every script that plots by
# wave, so "Baseline"/"Year 2"/etc. always get the same color+shape.
# Shapes: solid glyphs for the 4 MRI waves (bl, fu2, fu4, fu6 -- ABCD
# collects imaging every other year), hollow counterparts for the 3
# annual-only waves (fu1, fu3, fu5), so filled-vs-open is a bonus visual
# grouping on top of shape's primary job of wave identity.
# ---------------------------------------------------------------------------
wave_codes <- c("bl", "fu1", "fu2", "fu3", "fu4", "fu5", "fu6")
wave_display_labels <- setNames(
  c("Baseline", paste("Year", 1:6)),
  wave_codes
)
pal_waves <- setNames(sunset_ramp(7), wave_codes)
pal_shapes_waves <- setNames(
  c(bl = 16, fu1 = 0, fu2 = 17, fu3 = 5, fu4 = 15, fu5 = 2, fu6 = 18),
  wave_codes
)

# diverging fill scale (negative -> zero -> positive), built from the same
# palette family, for heatmaps/tiles of signed values (e.g. beta weights).
# Requires ggplot2 already loaded by the calling script.
sunset_diverging <- function(midpoint = 0, ...) {
  ggplot2::scale_fill_gradient2(
    low = unname(sunset_anchors["slate"]),
    mid = unname(sunset_anchors["cream"]),
    high = unname(sunset_anchors["bordeaux"]),
    midpoint = midpoint,
    ...
  )
}
