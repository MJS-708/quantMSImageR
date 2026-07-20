setGeneric("combine_MSIs", function(MSIobject, ...) standardGeneric("combine_MSIs"))

#' Function to combine MSImagingExperiment objects.
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject MSImagingExperiment object from Cardinal
#' @return quant_MSImagingExperiment object with intensity values replaced with response
#' @param ... additional MSImagingExperiment object to combine - must have matching fData() and same columns form pData().
#'
#' @examples
#' p1 <- system.file("extdata", "example.raw", "section01.RDS",
#'                   package = "quantMSImageR")
#' p2 <- system.file("extdata", "example.raw", "section02.RDS",
#'                   package = "quantMSImageR")
#' combined <- combine_MSIs(readRDS(p1), readRDS(p2))
#'
#' @aliases combine_MSIs
#' @export
setMethod("combine_MSIs", "MSImagingExperiment",
          function(MSIobject, ...){

            objects <- c(as.list(environment()), list(...))

            f_data = fData(MSIobject)

            for(ind in seq_along(objects)[-1]){

              MSIobject = as( cbind(MSIobject, objects[[ind]]), 'MSImagingExperiment')

              # Cardinal's cbind() concatenates the two objects' featureData
              # columns, so the accumulated fData grows and would not match the
              # next section's fData on the following iteration. The features are
              # identical across sections, so restore the shared fData each step.
              fData(MSIobject) = f_data

            }

            MSIobject = as(MSIobject, "quant_MSImagingExperiment")


            fData(MSIobject) = f_data

            return(MSIobject)

          })
