#' Accessors for calibration and tissue metadata
#'
#' Extract and replace the contents of a `quant_MSImagingExperiment`'s
#' `calibrationInfo` and `tissueInfo` slots. Prefer these to reaching into the
#' slots with `@`: the internal representation may change, but these functions
#' will not.
#'
#' Note the naming: `calibrationInfo()` and `tissueInfo()` are the *class
#' constructors*, so the slot accessors are `calibrationData()` and
#' `tissueData()`.
#'
#' @details
#' The calibration accessors follow the workflow. [summariseCalLevels()]
#' populates `calibrationLevels()`, [createCalCurve()] fills
#' `calibrationModels()` and `calibrationR2()`, and [int2conc()] reads the
#' models back out. `calibrationData()` moves the whole block at once, which is
#' what carrying fitted models from a standards acquisition to a study object
#' requires.
#'
#' @param x A `quant_MSImagingExperiment` object.
#' @param value Replacement value.
#' @param which For `tissueMatrix()`, whether to return the per-ROI average
#'   matrix (`"roi"`, the default) or the per-pixel matrix (`"pixel"`).
#'
#' @return
#' `calibrationModels()` a named list of `lm` objects, one per calibrated
#' analyte. `calibrationDiagnostics()` a data frame of fit quality, one row per
#' analyte: `feature`, `r2`, and an `out_of_range` column once [int2conc()] has
#' run. `calibrationR2()` the `feature`/`r2` columns of that table only.
#' `calibrationLevels()` the summarised response per calibration level.
#' `calibrationMetadata()` the calibration design supplied by the user.
#' `calibrationData()` the whole `calibrationInfo` object.
#' `tissueMatrix()` a data frame of features by pixels or by ROI.
#' The replacement forms return the modified object.
#'
#' @examples
#' cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
#' cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
#'           "quant_MSImagingExperiment")
#' cal_metadata <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
#'
#' cal <- summariseCalLevels(cal, cal_metadata, val_slot = "intensity",
#'                             cal_label = "Cal", id = "identifier")
#' cal <- createCalCurve(cal, cal_type = "cal")
#'
#' names(calibrationModels(cal))
#' calibrationDiagnostics(cal)
#' head(calibrationLevels(cal))
#'
#' # Carry the fitted models onto a study acquisition
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' study <- as(readRDS(p), "quant_MSImagingExperiment")
#' calibrationData(study) <- calibrationData(cal)
#'
#' @import Cardinal
#' @include setClasses.R
#' @name quantMSImageR-accessors
#' @aliases calibrationData calibrationModels calibrationR2 calibrationDiagnostics calibrationLevels calibrationMetadata tissueData tissueMatrix
#' @family accessors
NULL

# ---- calibrationData --------------------------------------------------------

#' @rdname quantMSImageR-accessors
#' @export
setGeneric("calibrationData", function(x) standardGeneric("calibrationData"))

#' @rdname quantMSImageR-accessors
#' @export
setMethod("calibrationData", "quant_MSImagingExperiment",
          function(x) x@calibrationInfo)

#' @rdname quantMSImageR-accessors
#' @export
setGeneric("calibrationData<-", function(x, value) standardGeneric("calibrationData<-"))

#' @rdname quantMSImageR-accessors
#' @export
setReplaceMethod("calibrationData", "quant_MSImagingExperiment",
                 function(x, value) {
                   # Validate the calibrationData itself: validObject() on the
                   # parent does not recurse into slots unless complete = TRUE.
                   validObject(value)
                   x@calibrationInfo <- value
                   x
                 })

# ---- calibrationModels ------------------------------------------------------

#' @rdname quantMSImageR-accessors
#' @export
setGeneric("calibrationModels", function(x) standardGeneric("calibrationModels"))

#' @rdname quantMSImageR-accessors
#' @export
setMethod("calibrationModels", "quant_MSImagingExperiment",
          function(x) x@calibrationInfo@cal_list)

#' @rdname quantMSImageR-accessors
#' @export
setGeneric("calibrationModels<-", function(x, value) standardGeneric("calibrationModels<-"))

#' @rdname quantMSImageR-accessors
#' @export
setReplaceMethod("calibrationModels", "quant_MSImagingExperiment",
                 function(x, value) {
                   ci <- x@calibrationInfo
                   ci@cal_list <- value
                   validObject(ci)
                   x@calibrationInfo <- ci
                   x
                 })

# ---- calibration diagnostics ------------------------------------------------

#' @rdname quantMSImageR-accessors
#' @export
setGeneric("calibrationDiagnostics", function(x) standardGeneric("calibrationDiagnostics"))

#' @rdname quantMSImageR-accessors
#' @export
setMethod("calibrationDiagnostics", "quant_MSImagingExperiment",
          function(x) x@calibrationInfo@r2_df)

#' @rdname quantMSImageR-accessors
#' @export
setGeneric("calibrationR2", function(x) standardGeneric("calibrationR2"))

#' @rdname quantMSImageR-accessors
#' @export
setMethod("calibrationR2", "quant_MSImagingExperiment",
          function(x) {
            # Kept for compatibility, and now narrowed to what its name claims:
            # the full table (R2 plus the out_of_range column int2conc() adds)
            # is calibrationDiagnostics().
            d <- x@calibrationInfo@r2_df
            d[, intersect(c("feature", "r2"), names(d)), drop = FALSE]
          })

# ---- calibrationLevels ------------------------------------------------------

#' @rdname quantMSImageR-accessors
#' @export
setGeneric("calibrationLevels", function(x) standardGeneric("calibrationLevels"))

#' @rdname quantMSImageR-accessors
#' @export
setMethod("calibrationLevels", "quant_MSImagingExperiment",
          function(x) x@calibrationInfo@cal_response_data)

# ---- calibrationMetadata ----------------------------------------------------

#' @rdname quantMSImageR-accessors
#' @export
setGeneric("calibrationMetadata", function(x) standardGeneric("calibrationMetadata"))

#' @rdname quantMSImageR-accessors
#' @export
setMethod("calibrationMetadata", "quant_MSImagingExperiment",
          function(x) x@calibrationInfo@cal_metadata)

# ---- tissueInfo / tissueMatrix ----------------------------------------------

# Note: `tissueInfo` is also the generator function for the tissueInfo class
# (see setClasses.R). Extracting the slot is spelled `tissueData()` to avoid
# masking the constructor.

#' @rdname quantMSImageR-accessors
#' @export
setGeneric("tissueData", function(x) standardGeneric("tissueData"))

#' @rdname quantMSImageR-accessors
#' @export
setMethod("tissueData", "quant_MSImagingExperiment",
          function(x) x@tissueInfo)

#' @rdname quantMSImageR-accessors
#' @export
setGeneric("tissueMatrix", function(x, which = c("roi", "pixel"))
  standardGeneric("tissueMatrix"))

#' @rdname quantMSImageR-accessors
#' @export
setMethod("tissueMatrix", "quant_MSImagingExperiment",
          function(x, which = c("roi", "pixel")) {
            which <- match.arg(which)
            if (which == "roi") x@tissueInfo@roi_average_matrix
            else                x@tissueInfo@all_pixel_matrix
          })
