#' calibrationInfo
#'
#' Class to store information about the calibration spots
#'
#' @return An object of class `calibrationInfo`.
#'
#' @examples
#' ci <- calibrationInfo()
#'
#' @name calibrationInfo
#' @aliases calibrationInfo-class
#' @export
calibrationInfo = setClass("calibrationInfo",
         slots = c(
           cal_metadata = "data.frame",
           cal_list = "list",
           cal_response_data = "data.frame",
           r2_df = "data.frame"
         )
)

#' tissueInfo
#'
#' Class to store information about the tissue pixels and ROIs
#'
#' @return An object of class `tissueInfo`.
#'
#' @examples
#' ti <- tissueInfo()
#'
#' @name tissueInfo
#' @aliases tissueInfo-class
#' @export
tissueInfo = setClass("tissueInfo",
                           slots = c(
                             roi_average_matrix = "data.frame",
                             all_pixel_matrix = "data.frame",
                             sample_metadata = "data.frame",
                             feature_metadata = "data.frame"
                           )
)


#' quant_MSImagingExperiment
#'
#' Class containting calibration metadata and MSImaging experiment
#'
#' @return An object of class `quant_MSImagingExperiment`, extending
#'   `Cardinal::MSImagingExperiment` with `calibrationInfo` and `tissueInfo` slots.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
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
