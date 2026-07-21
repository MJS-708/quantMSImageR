setGeneric("applySNR", function(MSIobject, ...) standardGeneric("applySNR"))

#' Apply SNR mask to intensity values
#'
#' Sets pixels to `NA` in `val_slot` wherever the `snr` spectra slot is `NA`.
#' Typically called after `int2snr()`: pixels that did not pass the SNR threshold
#' (stored as `NA` in the `snr` slot) are suppressed in the intensity slot so
#' they are excluded from downstream analysis and visualisation.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` object containing both
#'   `val_slot` and an `snr` spectra slot (populated by `int2snr()`).
#' @param val_slot Character. Name of the intensity slot to mask (default
#'   `"intensity"`).
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object with sub-threshold pixels set to `NA` in `val_slot`.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' obj <- int2snr(obj, val_slot = "intensity", sample_type = "sample_name",
#'                background = "background_pixels", tissue = "tissue_pixels",
#'                snr_thresh = 3)
#' obj <- applySNR(obj, val_slot = "intensity")
#'
#' @seealso [int2snr()]
#' @family filtering
#' @aliases applySNR
#' @export
setMethod("applySNR", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity", ...){

            # Iterate over features in study
            for(mz_ind in seq_len(nrow(fData(MSIobject)))){

              na_pixels = which(is.na(spectraData(MSIobject)[["snr"]][mz_ind, ]))

              # Save background response vector
              spectraData(MSIobject)[[val_slot]][mz_ind, na_pixels] = NA
            }

            return(MSIobject)
          })
