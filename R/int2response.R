setGeneric("int2response", function(MSIobject, ...) standardGeneric("int2response"))

#' Normalise pixel intensities to internal-standard response
#'
#' Divides each feature's intensity by the internal standard measured in the
#' same pixel, line or sample, which suppresses drift and much of the local
#' variation in ionisation efficiency. Only a single internal standard is
#' supported.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` object.
#' @param val_slot Character. Spectra slot to normalise (default
#'   `"intensity"`).
#' @param IS_name Character. Value identifying the internal standard in the
#'   **`analyte` column** of `fData(MSIobject)` -- the ion library's `Type`
#'   column, typically `"IS"`. This is a type label, not a transition name:
#'   passing a feature name will not match. `"None"` (the default) performs
#'   **within-feature normalisation** instead: each feature is divided by its own
#'   summary at the chosen `mode`. That is not internal-standard normalisation,
#'   and it is not the same as skipping normalisation -- it removes
#'   between-line or between-sample scale differences within each feature. To
#'   leave values untouched, simply do not call this function.
#' @param mode Character. Level at which the standard is summarised:
#'   \describe{
#'     \item{`"line"`}{Median across a whole acquisition line (the default). A
#'       *line* is one horizontal raster row -- constant `y`, varying `x` --
#'       which is the order DESI acquires in, so it is also the axis along which
#'       source drift accumulates.}
#'     \item{`"sample"`}{Median across the whole acquisition.}
#'     \item{`"pixel"`}{The standard in that pixel alone -- responsive to local
#'       suppression, but carries the standard's own shot noise.}
#'     \item{`"window"`}{Rolling median over `window` consecutive pixels along
#'       the same horizontal row, ordered by `x`. A compromise between `"pixel"`
#'       and `"line"`: it smooths the standard's own shot noise while still
#'       tracking drift across the row. The window is truncated at the ends of a
#'       row rather than wrapping, so pixels from different rows are never
#'       mixed.}
#'   }
#' @param window Integer. Number of consecutive pixels averaged when
#'   `mode = "window"` (default `15`). Ignored for the other modes.
#' @param remove_IS Logical. Drop the internal-standard feature from the
#'   returned object (default `TRUE`).
#' @param ... Additional arguments (currently unused).
#' @return The input object with a `response` spectra slot holding the
#'   internal-standard-normalised values.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' # normalise each feature to itself per line (no internal standard)
#' obj <- int2response(obj, val_slot = "intensity", IS_name = "None")
#'
#' @family filtering
#' @aliases int2response
#' @export
setMethod("int2response", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity", IS_name = "None",
                   mode = c("line", "sample", "pixel", "window"),
                   window = 15, remove_IS = TRUE, ...){

            mode <- match.arg(mode)

            if(IS_name == "None"){
              IS_ind = NULL
            } else if(!any(fData(MSIobject)$analyte == IS_name)){
              message("No IS in this study so normalise to individual lipids")

              IS_ind = NULL
            } else{
              IS_ind = which(fData(MSIobject)$analyte == IS_name)
            }

            spectra(MSIobject, "response") = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))

            # Iterate over samples in study
            for(sample in unique(pData(MSIobject)$run)){

              # Select IS mz and pixels for specific sample
              sample_pixels = which(pData(MSIobject)$run == sample)

              tempMSIobject = MSIobject[, sample_pixels]

              # Save IS intensity vector
              IS_vec = spectraData(tempMSIobject)[[val_slot]][IS_ind, ]

              for(mz_ind in seq_len(nrow(fData(tempMSIobject)))){

                ints = spectraData(tempMSIobject)[[val_slot]][mz_ind, ]

                # Denominator vector: the internal standard where one was
                # found, otherwise the feature itself (within-feature
                # normalisation).
                denom_src = if(!is.null(IS_ind)) IS_vec else ints

                # Build `response` by index rather than by concatenation. The
                # per-line and per-window modes group pixels, and groups are not
                # guaranteed to be contiguous or in pixel order, so appending
                # results group-by-group would misalign them on assignment.
                response = rep(NA_real_, length(ints))

                if(mode == "pixel"){
                  response = ints / denom_src

                } else if(mode == "sample"){
                  response = ints / median(denom_src, na.rm = TRUE)

                } else if(mode == "line"){
                  for(line in unique(pData(tempMSIobject)$y)){
                    lp = which(pData(tempMSIobject)$y == line)
                    response[lp] = ints[lp] / median(denom_src[lp], na.rm = TRUE)
                  }

                } else if(mode == "window"){
                  # Rolling median of the standard over `window` consecutive
                  # pixels along each acquisition line, ordered by x. The window
                  # is truncated at the ends of a line rather than wrapping, so
                  # pixels from different lines are never mixed.
                  for(line in unique(pData(tempMSIobject)$y)){
                    lp = which(pData(tempMSIobject)$y == line)
                    lp = lp[order(pData(tempMSIobject)$x[lp])]
                    d  = .roll_median(denom_src[lp], window)
                    response[lp] = ints[lp] / d
                  }
                }

                spectra(MSIobject, "response")[mz_ind, sample_pixels] = response

              }
            }

            # Remove IS m/z (only when an internal standard feature was found;
            # guards against MSIobject[-NULL, ] emptying the object)
            if(remove_IS == TRUE && !is.null(IS_ind)){
              MSIobject = MSIobject[-IS_ind, ]
            }

            return(MSIobject)
          })
