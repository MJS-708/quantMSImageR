setGeneric("int2snr", function(MSIobject, ...) standardGeneric("int2snr"))

#' Function to convert the intensity values to SNR per pixel based on same transitions in noise/background pixels.Run after IS normalization.
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject quant_MSImagingExperiment - including pData(MSIobject)$sample_type == Noise
#' @param noise character in pData(MSIobject)$sample_type which indicates background / noise pixels
#' @param tissue character in pData(MSIobject)$sample_type which indicates tissue pixels to calculate SNR for
#' @param val_slot character defining slot name to normalise - takes "intensity" as default
#' @param snr_thresh Global minimum SNR to accept (below this value SNR = NA).
#'   Applied to every feature unless overridden by `snr_overrides`.
#' @param snr_overrides Optional named numeric vector mapping feature name
#'   (as in `fData(MSIobject)$name`) to a feature-specific SNR threshold. A
#'   feature listed here uses its own threshold instead of `snr_thresh`;
#'   features not listed fall back to `snr_thresh`. Default `NULL` (all
#'   features use the global threshold).
#' @param sample_type character column in pData(MSIobject) holding the pixel-type
#'   labels (default "sample_type").
#' @param average character, "mean" or "median": statistic used to summarise the
#'   noise-pixel vector per feature.
#' @param ... Additional arguments (currently unused).
#' @return MSIobject with intensity values replaced with SNR values
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' obj <- int2snr(obj, val_slot = "intensity", sample_type = "sample_name",
#'                noise = "noise_pixels", tissue = "tissue_pixels", snr_thresh = 3)
#'
#' @aliases int2snr
#' @export
setMethod("int2snr", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "response", noise = "Noise", tissue = "Tissue", snr_thresh = 3,
                   sample_type = "sample_type", average = c("mean", "median"),
                   snr_overrides = NULL, ...){
            average <- match.arg(average)

            if(!any(pData(MSIobject)[[sample_type]] == noise)){
              # No noise pixels to compute SNR against -- fall back to a
              # non-filtering snr slot (copy of intensity). Adding the slot
              # here is required so that downstream applySNR() and combine_MSIs()
              # (cbind) see a consistent set of spectra arrays across sections.
              warning("int2snr: no '", noise, "' pixels in sample_type='",
                      sample_type, "'. Returning intensity as snr (no SNR filtering).",
                      call. = FALSE)
              spectra(MSIobject, "snr") <- spectra(MSIobject, val_slot)
              return(MSIobject)
            }

            #Set noise and tissue pixels
            noise_pixels = which(pData(MSIobject)[[sample_type]] == noise)
            tissue_pixels = which(pData(MSIobject)[[sample_type]] == tissue)

            spectra(MSIobject, "snr") = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))

            # Iterate over features in study
            for(mz_ind in seq_len(nrow(fData(MSIobject)))){

              # Save noise response vector
              noise_vec = spectraData(MSIobject)[[val_slot]][mz_ind, noise_pixels]

              # Skip features with no usable noise signal -- typical of panel-
              # padded rows (a transition that exists in the other acquisition's
              # MRM panel but not this one). snr stays NA for these features.
              if (all(is.na(noise_vec) | noise_vec == 0)) next

              #Deal with 0 values in noise vector - 10% of lowest
              noise_vec[which(noise_vec ==0)] = NA
              mv_impute = min(noise_vec, na.rm = TRUE) / 10
              noise_vec[which(is.na(noise_vec))] = mv_impute

              noise_level <- if (average == "mean") mean(noise_vec) else median(noise_vec)

              # Feature-specific SNR threshold if provided, else the global one.
              feat_nm <- as.character(fData(MSIobject)$name[mz_ind])
              thr_i <- if (!is.null(snr_overrides) &&
                           feat_nm %in% names(snr_overrides))
                         snr_overrides[[feat_nm]] else snr_thresh

              # Calculate S/N of tissue pixels
              vec = spectraData(MSIobject)[[val_slot]][mz_ind, ]
              snr =  vec / noise_level
              snr =  ifelse(snr < thr_i, NA, snr)
              snr[which(!seq_along(vec) %in% tissue_pixels)] = NA
              spectra(MSIobject, "snr")[mz_ind, ] = snr

            }

            return(MSIobject)
          })
