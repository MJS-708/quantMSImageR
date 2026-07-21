setGeneric("zero2na", function(MSIobject, ...) standardGeneric("zero2na"))

#' Replace zero intensities with NA
#'
#' Zeros in an MRM acquisition mean "nothing recorded" rather than "measured as
#' zero", and left in place they bias means and compress colour scales. This
#' converts them to `NA` so they are excluded from summaries.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` object.
#' @param val_slot Character. Spectra slot to modify (default `"intensity"`).
#' @return The input object with zero values in `val_slot` replaced by `NA`.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' obj <- zero2na(obj, val_slot = "intensity")
#'
#' @family filtering
#' @aliases zero2na
#' @export
setMethod("zero2na", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity"){

            MSIobject = as(MSIobject, "quant_MSImagingExperiment")

            for(i in seq_len(nrow(fData(MSIobject)))){

              ints = spectraData(MSIobject)[[val_slot]][i,]

              if(all(is.na(ints))){
                message(sprintf("all intensities are NA for m/z %s. Doing nothing.", i))
              } else if(sum(ints, na.rm = TRUE) == 0){
                message(sprintf("all intensities are 0 for m/z %s. Making NA.", i))
                spectra(MSIobject)[i, ] = NA
              } else{
                ints[which(ints == 0)] = NA
                spectraData(MSIobject)[[val_slot]][i, ] <- ints

              }
            }

            return(MSIobject)

          })
