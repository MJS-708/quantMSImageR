setGeneric("remove_blank_mzs", function(MSIobject) standardGeneric("remove_blank_mzs"))

#' Remove features without observed signal
#'
#' Drops any feature whose intensities are entirely zero or `NA` -- typically a
#' transition present in the acquisition method but never detected, or a panel
#' row padded in by [setCommonAxis()]. Reporting these alongside real
#' measurements is how empty transitions reach a figure.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` object.
#' @return The input object with features carrying no data removed.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' obj <- remove_blank_mzs(obj)
#'
#' @family acquisition
#' @aliases remove_blank_mzs
#' @export
setMethod("remove_blank_mzs", "quant_MSImagingExperiment",
          function(MSIobject){

            remove_inds = c()
            for(mz_ind in seq_len(nrow(fData(MSIobject)))){
              if(all(is.na(spectra(MSIobject)[mz_ind, ]))){
                remove_inds = c(remove_inds, mz_ind)
              }
            }

            if(length(remove_inds) > 0){
              MSIobject = MSIobject[-remove_inds, ]
            }

            return(MSIobject)
          })
