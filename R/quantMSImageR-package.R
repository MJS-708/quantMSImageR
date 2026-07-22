#' quantMSImageR: quantification of DESI-MRM mass spectrometry imaging data
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
