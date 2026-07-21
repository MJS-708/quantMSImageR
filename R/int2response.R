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
#' @param IS_name Character. Name of the internal standard in
#'   `fData(MSIobject)` under the analyte header. `"None"` (the default)
#'   normalises each feature to itself.
#' @param mode Character. Level at which the standard is summarised:
#'   `"sample"` (median internal-standard intensity per sample), `"line"`
#'   (per acquisition line, the default) or `"pixel"` (per pixel).
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
          function(MSIobject, val_slot = "intensity", IS_name = "None", mode = "line", remove_IS = TRUE, ...){

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

                if(mode == "pixel"){
                  if(!is.null(IS_ind)){
                    response = ints / IS_vec
                  } else{
                    response = ints / ints
                  }
                }
                if(mode == "sample"){
                  if(!is.null(IS_ind)){
                    response = ints / median(IS_vec, na.rm = TRUE)
                  } else{
                    response = ints / median(ints, na.rm = TRUE)
                  }
                }
                if(mode == "line"){
                  response = c()
                  for(line in unique(pData(tempMSIobject)$y)){

                    line_pixels = which(pData(tempMSIobject)$y == line)

                    if(!is.null(IS_ind)){
                      response = c(response, (ints[line_pixels] / median(IS_vec[line_pixels], na.rm = TRUE)))
                    } else{
                      response = c(response, (ints[line_pixels] / median(ints[line_pixels], na.rm = TRUE)))
                    }
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
