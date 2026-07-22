#' Calibration metadata for an imaging experiment
#'
#' Holds everything needed to turn a measured response into an estimated amount:
#' the calibration design, the summarised response at each level, the fitted
#' per-analyte models, and their fit quality. It is carried in the
#' `calibrationInfo` slot of a [quant_MSImagingExperiment].
#'
#' @details
#' The slots are filled in pipeline order and are empty (zero-row data frames,
#' empty list) until the corresponding step has run:
#'
#' \describe{
#'   \item{`cal_metadata`}{`data.frame`. The calibration design supplied by the
#'     user: one row per (spot, analyte), with columns `identifier`, `analyte`,
#'     `amount_pg` and `level`. Set by [summarise_cal_levels()].}
#'   \item{`cal_response_data`}{`data.frame`. Summarised response per
#'     calibration level, with columns `analyte`, `pg_perpixel`,
#'     `response_perpixel` and `level`. Set by [summarise_cal_levels()].}
#'   \item{`cal_list`}{`list`. One `stats::lm` per analyte, of
#'     `response_perpixel ~ pg_perpixel`, **named by analyte** so that
#'     [int2conc()] can match models to features by `fData(x)$name`. Set by
#'     [create_cal_curve()].}
#'   \item{`r2_df`}{`data.frame`. One row per analyte, with `feature` and `r2`;
#'     [int2conc()] adds an `out_of_range` column giving the proportion of
#'     pixels extrapolated beyond the calibrated range. Set by
#'     [create_cal_curve()].}
#' }
#'
#' Invariants checked by the validity method: `cal_list` must be named if it is
#' non-empty, and every name must appear in `r2_df$feature` when `r2_df` is
#' populated. An unnamed `cal_list` is accepted only when empty -- matching
#' models to features positionally would silently apply one analyte's curve to
#' another.
#'
#' Use the accessors ([calibrationModels()], [calibrationDiagnostics()],
#' [calibrationLevels()], [calibrationMetadata()]) rather than `@`.
#'
#' @slot cal_metadata `data.frame` of the calibration design.
#' @slot cal_list Named `list` of fitted `lm` models, one per analyte.
#' @slot cal_response_data `data.frame` of summarised response per level.
#' @slot r2_df `data.frame` of fit quality per analyte.
#'
#' @return An object of class `calibrationInfo`.
#'
#' @examples
#' # Empty, as attached to a freshly coerced object
#' ci <- calibrationInfo()
#' ci
#'
#' # Populated by the calibration pipeline
#' cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
#' cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
#'           "quant_MSImagingExperiment")
#' cal_metadata <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
#' cal <- summarise_cal_levels(cal, cal_metadata, val_slot = "intensity",
#'                             cal_label = "Cal", id = "identifier")
#' cal <- create_cal_curve(cal, cal_type = "cal")
#' calibrationDiagnostics(cal)
#'
#' @seealso [quant_MSImagingExperiment], [create_cal_curve()],
#'   [quantMSImageR-accessors]
#' @family calibration
#' @name calibrationInfo-class
#' @aliases calibrationInfo-class
#' @export calibrationInfo
calibrationInfo = setClass("calibrationInfo",
         slots = c(
           cal_metadata = "data.frame",
           cal_list = "list",
           cal_response_data = "data.frame",
           r2_df = "data.frame"
         )
)

setValidity("calibrationInfo", function(object) {
  msg <- character()

  # A named cal_list is what lets int2conc() match models to features by name.
  # Allowing an unnamed non-empty list would reintroduce positional matching.
  if (length(object@cal_list) > 0 && is.null(names(object@cal_list)))
    msg <- c(msg, "cal_list must be named by analyte when non-empty")

  if (length(object@cal_list) > 0 && nrow(object@r2_df) > 0 &&
      "feature" %in% names(object@r2_df)) {
    missing <- setdiff(names(object@cal_list), object@r2_df$feature)
    if (length(missing) > 0)
      msg <- c(msg, paste0("cal_list names absent from r2_df$feature: ",
                           paste(missing, collapse = ", ")))
  }

  if (length(msg)) msg else TRUE
})


#' Tissue-level summaries for an imaging experiment
#'
#' Holds the tabular views of an imaging experiment produced by
#' [createMSIDatamatrix()]: one row per pixel, and optionally one row per region
#' of interest, together with the accompanying metadata. It is carried in the
#' `tissueInfo` slot of a [quant_MSImagingExperiment].
#'
#' @details
#' \describe{
#'   \item{`all_pixel_matrix`}{`data.frame`, one row per pixel and one column
#'     per feature.}
#'   \item{`roi_average_matrix`}{`data.frame`, one row per region of interest,
#'     populated only when `createMSIDatamatrix(roi_header=)` is given.}
#'   \item{`sample_metadata`}{`data.frame` describing the rows of the matrices
#'     above (sample, run, ROI identifier).}
#'   \item{`feature_metadata`}{`data.frame` describing their columns.}
#' }
#'
#' All four are zero-row data frames until [createMSIDatamatrix()] has run.
#' Access them with [tissueMatrix()] and [tissueData()] rather than `@`.
#'
#' @slot roi_average_matrix `data.frame` of per-ROI averages.
#' @slot all_pixel_matrix `data.frame` of per-pixel values.
#' @slot sample_metadata `data.frame` describing the matrix rows.
#' @slot feature_metadata `data.frame` describing the matrix columns.
#'
#' @return An object of class `tissueInfo`.
#'
#' @examples
#' ti <- tissueInfo()
#'
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' obj <- createMSIDatamatrix(obj, val_slot = "intensity", roi_header = NA)
#' dim(tissueMatrix(obj, which = "pixel"))
#'
#' @seealso [quant_MSImagingExperiment], [createMSIDatamatrix()],
#'   [quantMSImageR-accessors]
#' @name tissueInfo-class
#' @aliases tissueInfo-class
#' @export tissueInfo
tissueInfo = setClass("tissueInfo",
                           slots = c(
                             roi_average_matrix = "data.frame",
                             all_pixel_matrix = "data.frame",
                             sample_metadata = "data.frame",
                             feature_metadata = "data.frame"
                           )
)


#' Quantifiable MS imaging experiment
#'
#' Extends `Cardinal::MSImagingExperiment` with the metadata `quantMSImageR`
#' needs to filter against background pixels and to convert a measured response
#' into an estimated amount. Nearly every function in the package takes one of
#' these, so a Cardinal object is normally coerced on the way in.
#'
#' @details
#' The object is **features x pixels**: one row per MRM transition, one column
#' per pixel. Derived layers are added as additional spectra slots rather than
#' overwriting `intensity`, so the raw values survive the whole pipeline:
#'
#' \describe{
#'   \item{`intensity`}{Raw measured response, as read from the acquisition.}
#'   \item{`response`}{Internal-standard-normalised, added by [int2response()].}
#'   \item{`snr`}{Background-referenced signal-to-noise, added by [int2snr()].}
#'   \item{`pg_pixel`, `pg_mm2`}{Calibrated amount estimates, added by
#'     [int2conc()].}
#' }
#'
#' Inherited structure comes from Cardinal: `fData()` for feature metadata
#' (`mz`, `name`, transition m/z, ion-library columns), `pData()` for pixel
#' metadata (`x`, `y`, `run`, and the tissue mask in `sample_name`) and
#' `experimentData()` for acquisition metadata -- `pixelSize` in particular,
#' without which [int2conc()] cannot produce `pg_mm2`.
#'
#' Several acquisitions are held in a single object, distinguished by the `run`
#' column of `pData()`; this is what [combine_MSIs()] produces and what the
#' per-sample summaries group on.
#'
#' @slot calibrationInfo A [calibrationInfo-class] object.
#' @slot tissueInfo A [tissueInfo-class] object.
#'
#' @return An object of class `quant_MSImagingExperiment`, extending
#'   `Cardinal::MSImagingExperiment` with `calibrationInfo` and `tissueInfo`
#'   slots.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#'
#' # Cardinal object in, quant_MSImagingExperiment out
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#'
#' dim(obj)                      # features x pixels
#' names(spectraData(obj))       # available layers
#' head(as.data.frame(pData(obj)), 3)
#'
#' # Derived layers accumulate alongside intensity
#' obj <- int2snr(obj, snr_thresh = 3)
#' names(spectraData(obj))
#'
#' @seealso [calibrationInfo-class], [tissueInfo-class],
#'   [quantMSImageR-accessors], [combine_MSIs()]
#'
#' @import Cardinal
#' @name quant_MSImagingExperiment
#' @aliases quant_MSImagingExperiment-class
#' @export
quant_MSImagingExperiment = setClass("quant_MSImagingExperiment",
         contains = 'MSImagingExperiment',
         slots = c(
           calibrationInfo = "calibrationInfo",
           tissueInfo = "tissueInfo"
         )
)
