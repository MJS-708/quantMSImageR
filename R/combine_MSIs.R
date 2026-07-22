setGeneric("combine_MSIs", function(MSIobject, ...) standardGeneric("combine_MSIs"))

#' Combine MSI experiments across acquisitions
#'
#' Concatenates two or more acquisitions pixel-wise into a single object, so
#' that a whole study can be filtered, summarised and plotted together. Each
#' input keeps its own `run`, which is what downstream per-sample summaries
#' group on.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject An `MSImagingExperiment`; the first acquisition, whose
#'   `fData()` defines the shared feature axis.
#' @param ... Further `MSImagingExperiment` objects to combine. All must share
#'   the feature axis and the `pData()` columns of `MSIobject` -- use
#'   [align_features()] first if they do not.
#' @return A single `quant_MSImagingExperiment` holding every input's pixels,
#'   with the shared `fData()` restored and one `run` level per acquisition.
#'
#' @examples
#' p1 <- system.file("extdata", "example.raw", "section01.RDS",
#'                   package = "quantMSImageR")
#' p2 <- system.file("extdata", "example.raw", "section02.RDS",
#'                   package = "quantMSImageR")
#' combined <- combine_MSIs(readRDS(p1), readRDS(p2))
#'
#' @family combining acquisitions
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
