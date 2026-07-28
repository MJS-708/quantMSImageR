#' quantMSImageR: processing and quantification of targeted mass spectrometry
#' imaging data
#'
#' Tools (extending Cardinal) for processing and quantifying targeted DESI-MRM
#' mass spectrometry imaging data: signal-to-noise filtering, tissue/background
#' separation, per-feature ion images, quantile heatmaps and standards-based
#' quantification.
#'
#' @keywords internal
#' @importFrom methods as is new setValidity validObject
#' @importFrom stats median na.omit quantile sd setNames
#' @importFrom utils read.csv read.table write.csv write.table
#' @importFrom grDevices dev.new rainbow
#' @importFrom grid unit
#' @importFrom dplyr across all_of any_of arrange bind_rows distinct filter
#' @importFrom dplyr group_by left_join mutate pull right_join row_number select
#' @importFrom dplyr summarise ungroup
#' @importFrom ggplot2 aes discrete_scale element_blank element_line element_rect
#' @importFrom ggplot2 element_text facet_grid geom_point geom_tile ggplot labs
#' @importFrom ggplot2 coord_equal geom_histogram geom_line geom_vline
#' @importFrom ggplot2 coord_cartesian facet_wrap scale_fill_manual
#' @importFrom ggplot2 rel theme theme_minimal
#' @importFrom circlize colorRamp2
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom grDevices colorRampPalette
#' @importFrom patchwork wrap_plots
#' @importFrom yaml read_yaml
#' @importFrom viridis scale_fill_viridis
#' @importFrom chemCal inverse.predict
#' @importFrom tibble tibble column_to_rownames rownames_to_column
#' @importFrom tidyr pivot_wider
#' @importFrom ProtGenerics spectraData spectraData<-
"_PACKAGE"

# Non-standard-evaluation variables (dplyr/subset/aes) so R CMD check does not
# flag them as undefined globals.
utils::globalVariables(c(
  "Polarity", "ROI", "Type", "all_noise", "amount_pg", "analyte",
  "collision_eV", "cone_V", "duplicate_mrms", "feature", "label",
  "pg_perpixel", "response_perpixel", "sd_response",
  "pixel_ind", "precursor_mz", "product_mz", "response", "sample_name",
  "transition_id_int",
  "transition_id_name", "x", "x_loci", "y", "y_loci"
))

# Background-pixel label aliases.
#
# Non-tissue pixels are now labelled "background_pixels" (or "Background"),
# which describes them more accurately than the historical "noise_pixels" /
# "Noise". Both spellings are accepted everywhere a background label is matched,
# so masks and acquisitions created with earlier versions keep working.
# Centred rolling median over `k` consecutive values, truncated at both ends
# rather than padded or wrapped. Used by int2response(mode = "window") to
# smooth the internal standard along an acquisition line without borrowing
# pixels from the neighbouring line.
.roll_median <- function(v, k) {
  n <- length(v)
  if (n == 0L) return(v)
  k <- max(1L, as.integer(k))
  if (k >= n) return(rep(stats::median(v, na.rm = TRUE), n))

  half <- k %/% 2L
  vapply(seq_len(n), function(i) {
    lo <- max(1L, i - half)
    hi <- min(n,  i + half)
    stats::median(v[lo:hi], na.rm = TRUE)
  }, numeric(1))
}

.bg_labels <- function(x) {
  alias <- c(background_pixels = "noise_pixels",
             noise_pixels      = "background_pixels",
             Background        = "Noise",
             Noise             = "Background")
  unique(c(x, unname(alias[x][!is.na(alias[x])])))
}

# Feature-metadata column holding the analytical class of each feature -- "IS",
# "Analyte" and so on, taken from the ion library's `Type` column.
#
# It used to be stored as `analyte`, which collided with the other meaning of
# that word: `cal_metadata$analyte` names the *compound* being quantified,
# whereas this column names its *role*. Objects built by earlier versions carry
# the old spelling, so both are read.
.FEAT_TYPE <- "feature_type"

# "Not set", spelled the way a YAML config spells it.
#
# The configs use the literal string "None" for an option that is switched off
# (is_name, features$exclude), so an argument fed straight from YAML can arrive
# as "None" rather than NULL. Normalising once here keeps every caller from
# having to know that.
.none <- function(x) {
  if (is.null(x)) return(NULL)
  x <- as.character(x)
  if (length(x) != 1L || is.na(x) || !nzchar(trimws(x)) ||
      tolower(trimws(x)) %in% c("none", "null", "na")) return(NULL)
  trimws(x)
}

.feature_col <- function(obj, column) {
  column <- .none(column)
  if (is.null(column)) return(NULL)
  fd <- fData(obj)
  if (!column %in% names(fd)) return(NULL)
  fd[[column]]
}

.feature_type <- function(obj) {
  fd <- fData(obj)
  for (nm in c(.FEAT_TYPE, "analyte")) if (nm %in% names(fd)) return(fd[[nm]])
  NULL
}

# Which reference transition, if any, each query transition is.
#
# One definition of "same transition", used by both read_mrm() (measured
# transitions -> ion library) and build_feature_meta() (features -> ion
# library). They used to build string keys independently -- one rounding to
# whole numbers, the other not -- so the two could disagree about the same pair
# of acquisitions.
#
# A precursor/product pair is a transition's identity, so a query matching
# several reference rows is genuinely ambiguous and cannot be resolved from m/z
# alone. `ambiguity` says what to do:
#
#   "combine" (default) keeps the closest row's annotation but names the
#             feature for EVERY entry it matched, joined with " || ", e.g.
#             "LTC4 || 14_15-LTC4". Two isomers really can share a transition:
#             LTC4 (626.3 -> 308.17) and 14_15-LTC4 (626.3 -> 308.2) are 0.03 Da
#             apart and no acquisition can separate them. The joined name says
#             the measurement is one of these, which is what the data supports;
#             naming it for one of them asserts something it does not. The
#             separator matches the convention already used for merged
#             transitions elsewhere in the package.
#   "error"   refuses. Right when a collision means the library is wrong rather
#             than the chemistry ambiguous.
#   "warn"    takes the closest and says so.
#   "nearest" takes the closest quietly.
#
# A joined name is not a substitute for a clean library: two entries for the
# SAME compound are redundancy, and "combine" says so in its message rather
# than papering over it.
#
# Returns an integer vector of reference row indices (the closest match), NA
# where nothing matched, carrying attr "matches": a list holding, for each
# ambiguous query, every reference row it matched. NULL for unambiguous ones.
.match_transitions <- function(prec, prod, ref_prec, ref_prod,
                               tolerance = 0.05,
                               ambiguity = c("combine", "error", "warn",
                                             "nearest"),
                               labels = NULL, ref_labels = NULL,
                               context = "match_transitions") {
  ambiguity <- match.arg(ambiguity)
  prec <- suppressWarnings(as.numeric(prec))
  prod <- suppressWarnings(as.numeric(prod))
  ref_prec <- suppressWarnings(as.numeric(ref_prec))
  ref_prod <- suppressWarnings(as.numeric(ref_prod))

  if (is.null(labels))     labels     <- as.character(seq_along(prec))
  if (is.null(ref_labels)) ref_labels <- as.character(seq_along(ref_prec))

  out <- rep(NA_integer_, length(prec))
  hits <- vector("list", length(prec))
  amb <- character(0)

  for (i in seq_along(prec)) {
    if (!is.finite(prec[i]) || !is.finite(prod[i])) next
    hit <- which(abs(ref_prec - prec[i]) <= tolerance &
                 abs(ref_prod - prod[i]) <= tolerance)
    if (length(hit) == 0L) next
    if (length(hit) > 1L) {
      hits[[i]] <- hit
      amb <- c(amb, sprintf("  %s (%g -> %g) matches: %s", labels[i],
                            prec[i], prod[i],
                            paste(ref_labels[hit], collapse = ", ")))
      # Closest by summed absolute deviation, so every mode gets the best
      # available annotation rather than whichever row came first.
      hit <- hit[which.min(abs(ref_prec[hit] - prec[i]) +
                           abs(ref_prod[hit] - prod[i]))]
    }
    out[i] <- hit
  }
  attr(out, "matches") <- hits

  if (length(amb)) {
    if (ambiguity == "combine") {
      message(context, ": ", length(amb), " transition(s) match more than one ",
              "ion-library entry within ", tolerance, " Da; each is named for ",
              "all of them:\n", paste(amb, collapse = "\n"),
              "\nIf any pair above is the same compound entered twice, that is ",
              "library redundancy - de-duplicate it rather than shipping a ",
              "joined name.")
    } else if (ambiguity != "nearest") {
      msg <- paste0(context, ": ", length(amb),
                    " transition(s) match more than one ion-library entry within ",
                    tolerance, " Da:\n", paste(amb, collapse = "\n"),
                    "\nm/z alone cannot say which. Tighten mz_tolerance, give ",
                    "the colliding entries distinct product ions, or set ",
                    "ambiguity = \"combine\" to name it for both.")
      if (ambiguity == "error") stop(msg, call. = FALSE)
      warning(msg, call. = FALSE)
    }
  }
  out
}

# Regression weights for a calibration fit.
#
# 1/x and 1/x^2 are undefined at a zero (blank) level. Substituting a tiny
# amount and taking its reciprocal gives that one point a weight of order 1e9,
# which is why this is handled explicitly: a blank is given the same weight as
# the lowest positive standard, so it contributes once rather than deciding the
# fit. Negative amounts (which standard addition can produce below the
# x-intercept) are treated the same way.
.cal_weights <- function(x, weighting = c("1/x", "none", "1/x2")) {
  weighting <- match.arg(weighting)
  if (weighting == "none") return(rep(1, length(x)))

  pow <- if (weighting == "1/x") 1 else 2
  pos <- is.finite(x) & x > 0
  if (!any(pos)) return(rep(1, length(x)))

  w <- rep(NA_real_, length(x))
  w[pos] <- 1 / x[pos]^pow
  w[!pos] <- min(w[pos], na.rm = TRUE)
  w
}

# Calibrated-amount slot names.
#
# int2conc() writes `pg_pixel` / `pg_mm2`. Objects calibrated with earlier
# versions carry the old `conc - pg/pixel` / `conc - pg/mm2` names, and the
# values are amount estimates rather than concentrations, so the current names
# say so. Readers accept both spellings.
.SLOT_PG  <- "pg_pixel"
.SLOT_MM2 <- "pg_mm2"

.conc_slots <- function(x) {
  alias <- c(pg_pixel               = "conc - pg/pixel",
             `conc - pg/pixel`      = "pg_pixel",
             pg_mm2                 = "conc - pg/mm2",
             `conc - pg/mm2`        = "pg_mm2")
  unique(c(x, unname(alias[x][!is.na(alias[x])])))
}

# First calibrated-amount slot actually present on `obj`, preferring areal
# amount (comparable between acquisitions of different pixel size) and falling
# back to per-pixel. Returns NA_character_ when the object is not calibrated.
.pick_conc_slot <- function(obj, prefer = c("mm2", "pixel")) {
  prefer  <- match.arg(prefer)
  present <- names(spectraData(obj))
  wanted  <- if (prefer == "mm2") c(.SLOT_MM2, .SLOT_PG) else c(.SLOT_PG, .SLOT_MM2)
  for (w in wanted) {
    hit <- .conc_slots(w)[.conc_slots(w) %in% present]
    if (length(hit)) return(hit[1])
  }
  NA_character_
}

# Display unit for a calibrated-amount slot, either spelling.
.conc_unit <- function(slot) {
  if (isTRUE(slot %in% .conc_slots(.SLOT_MM2))) "pg/mm2" else "pg/pixel"
}
