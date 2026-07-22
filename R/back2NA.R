setGeneric("back2NA", function(MSIobject, ...) standardGeneric("back2NA"))

#' Set background pixel intensities to NA
#'
#' Replaces the spectral values of background pixels with `NA` in a chosen
#' spectra slot. Tissue pixels are left unchanged.
#'
#' @section Destructive:
#' This modifies `val_slot` in place rather than adding a new layer. Keep a copy
#' of the object if the original values are needed afterwards.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` object whose `pData()` contains
#'   a column identifying pixel type.
#' @param val_slot Character. Name of the spectra slot to modify (default `"intensity"`).
#' @param pixel_header Character. Column of `pData()` holding the pixel-type
#'   labels (default `"sample_name"`, the imaging convention).
#' @param background Character. Value in `pixel_header` that labels background
#'   pixels (default `"background_pixels"`). The historical `"noise_pixels"` /
#'   `"Noise"` labels are matched too, so masks made with earlier versions keep
#'   working.
#' @param sample_type Deprecated alias for `pixel_header`, kept for backward
#'   compatibility. When supplied it overrides `pixel_header`.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input `quant_MSImagingExperiment` with background pixels set to
#'   `NA` in `val_slot`. Returns the object unchanged if no background pixels are
#'   found.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' obj <- back2NA(obj)
#'
#' @family filtering
#' @aliases back2NA
#' @export
setMethod("back2NA", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity",
                   pixel_header = "sample_name",
                   background = "background_pixels",
                   sample_type = NULL, ...){

            # Backward compatibility: `sample_type` was the previous name.
            if (!is.null(sample_type)) pixel_header <- sample_type
            .bg <- .bg_labels(background)

            if(!any(pData(MSIobject)[[pixel_header]] %in% .bg)){
              message("No background pixels. Return same values")
              return(MSIobject)
            }

            #Set background pixels
            background_pixels = which(pData(MSIobject)[[pixel_header]] %in% .bg)

            # Iterate over features in study
            for(mz_ind in seq_len(nrow(fData(MSIobject)))){

              # Save background response vector
              spectraData(MSIobject)[[val_slot]][mz_ind, background_pixels] = NA

            }

            return(MSIobject)
          })
