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

            objects <- c(list(MSIobject), list(...))

            f_data = fData(MSIobject)

            # Validate before cbind(). Cardinal's cbind fails on mismatched
            # feature keys with a message about XDataFrames that says nothing
            # about which acquisition is the odd one out, and it will happily
            # combine objects whose run identifiers collide -- which silently
            # merges two sections into one sample everywhere downstream.
            .nm  <- function(o) as.character(fData(o)$name)
            .lay <- function(o) sort(names(spectraData(o)))
            for(ind in seq_along(objects)[-1]){
              if(!identical(.nm(objects[[ind]]), .nm(objects[[1]])))
                stop("combine_MSIs: object ", ind, " has different features ",
                     "from object 1. Align them with align_features() first, ",
                     "after checking they really are the same transitions.",
                     call. = FALSE)
              if(!identical(.lay(objects[[ind]]), .lay(objects[[1]])))
                stop("combine_MSIs: object ", ind, " has spectra layers (",
                     paste(.lay(objects[[ind]]), collapse = ", "),
                     ") that differ from object 1 (",
                     paste(.lay(objects[[1]]), collapse = ", "),
                     "). Run the same processing steps on every acquisition ",
                     "before combining.", call. = FALSE)
            }

            .runs <- unlist(lapply(objects, function(o)
                       unique(as.character(pData(o)$run))))
            if(anyDuplicated(.runs))
              stop("combine_MSIs: run identifiers must be unique across ",
                   "acquisitions, but these repeat: ",
                   paste(unique(.runs[duplicated(.runs)]), collapse = ", "),
                   ".", call. = FALSE)

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
