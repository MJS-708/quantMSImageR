setGeneric("int2SNR", function(MSIobject, ...) standardGeneric("int2SNR"))

#' Calculate background-referenced signal-to-noise ratios
#'
#' Expresses each tissue pixel's response as a ratio to the same feature's
#' response in the background pixels of the same acquisition.
#'
#' @details
#' For feature `f` and tissue pixel `p`, with background-pixel set `B`:
#'
#' \deqn{SNR_{f,p} = I_{f,p} / \mathrm{summary}_{b \in B}(I_{f,b})}
#'
#' where `summary` is the mean or median selected by `average`. This is a
#' background-referenced signal ratio, not the classical analytical definition
#' based on the standard deviation of the noise.
#'
#' `median` is the default because a background region that clips part of the
#' tissue, or contains a spot of carryover, shifts a mean far more than a
#' median. Zero background values are treated as missing and imputed at one
#' tenth of the smallest non-zero background value for that feature, so a
#' feature whose background is entirely zero or `NA` is skipped and its `snr`
#' stays `NA`.
#'
#' When **more than half** of a feature's background pixels are zero -- common
#' in DESI-MRM, where many transitions read exactly zero off tissue -- the
#' median *is* that imputed value. The noise level is then set by the
#' imputation rule rather than by the background, and `snr_thresh = t` passes
#' any pixel above `t / 10` times the dimmest non-zero background value: at
#' `t = 10`, anything at least as bright as the dimmest background pixel.
#' `int2SNR()` warns when this happens, naming the features (the warning has
#' class `quantMSImageR_imputed_noise`). `average = "mean"` is computed from
#' the background's own values, to which the imputed zeros add almost nothing,
#' and gives a stricter noise level for such features. Results are unchanged
#' by the warning; it says how much the filter is doing.
#'
#' Run this **after** internal-standard normalisation ([int2response()]) when
#' one is used. A threshold of 3 is a detection-level criterion; a stricter
#' threshold is advisable before converting a feature to calibrated amounts.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` whose `pData()` carries a
#'   background label in the `pixel_header` column.
#' @param pixel_header Character. Column of `pData()` holding the pixel-type
#'   labels (default `"sample_name"`, the imaging convention).
#' @param background Character. Value in `pixel_header` marking the background
#'   (non-tissue) pixels used to estimate the noise level (default
#'   `"background_pixels"`). The historical `"noise_pixels"` / `"Noise"` labels
#'   are matched too, so masks made with earlier versions keep working.
#' @param tissue Character. Value in `pixel_header` marking the tissue pixels
#'   SNR is calculated for (default `"tissue_pixels"`).
#' @param noise Deprecated alias for `background`, kept for backward
#'   compatibility. When supplied it overrides `background`.
#' @param sample_type Deprecated alias for `pixel_header`, kept for backward
#'   compatibility. When supplied it overrides `pixel_header`.
#' @param val_slot Character. Spectra slot holding the measured response
#'   (default `"intensity"`; use `"response"` after [int2response()]).
#' @param snr_thresh Global minimum SNR to accept (below this value SNR = NA).
#'   Applied to every feature unless overridden by `snr_overrides`.
#' @param no_background Character. What to do when `pixel_header` contains no
#'   background pixels, so there is nothing to reference the signal against:
#'   `"error"` (default), `"all_na"` (create the slot, all `NA`) or `"copy"`
#'   (copy `val_slot` into the `snr` slot -- note the slot then holds measured
#'   response, not a ratio). `"error"` and `"all_na"` apply only when
#'   `snr_thresh > 0`; with `snr_thresh = 0` no filtering was requested, so the
#'   values are copied through regardless.
#' @param snr_overrides Optional named numeric vector mapping feature name
#'   (as in `fData(MSIobject)$name`) to a feature-specific SNR threshold. A
#'   feature listed here uses its own threshold instead of `snr_thresh`;
#'   features not listed fall back to `snr_thresh`. Default `NULL` (all
#'   features use the global threshold).
#' @param average Character, `"median"` (default) or `"mean"`: statistic used to
#'   summarise the background-pixel vector per feature.
#' @param ... Additional arguments (currently unused).
#' @return The input object with an additional `snr` spectra slot. Values below
#'   the selected threshold, and values outside the pixels labelled `tissue`,
#'   are stored as `NA` in that slot. No existing slot is modified.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#'
#' # The defaults match the imaging convention used throughout the package
#' obj <- int2SNR(obj, snr_thresh = 3)
#' names(spectraData(obj))
#'
#' @family filtering
#' @aliases int2SNR
#' @export
setMethod("int2SNR", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity",
                   pixel_header = "sample_name",
                   background = "background_pixels", tissue = "tissue_pixels",
                   snr_thresh = 3,
                   average = c("median", "mean"), snr_overrides = NULL,
                   no_background = c("error", "all_na", "copy"),
                   noise = NULL, sample_type = NULL, ...){
            average <- match.arg(average)

            # Backward compatibility: `noise` and `sample_type` were the
            # previous argument names.
            if (!is.null(noise)) background <- noise
            if (!is.null(sample_type)) pixel_header <- sample_type
            .bg <- .bg_labels(background)

            no_background <- match.arg(no_background)

            if(!any(pData(MSIobject)[[pixel_header]] %in% .bg)){
              # Nothing to reference the signal against. The slot still has to
              # be created either way, because applySNR() and combineMSIs()
              # (cbind) need the same set of spectra arrays on every section.
              msg <- paste0("int2SNR: no '", background, "' pixels in ",
                            "pixel_header='", pixel_header, "'. ")

              # Only fatal when filtering was actually asked for. snr_thresh = 0
              # means "run unfiltered", so a missing background is not a failure
              # to do something the caller wanted.
              if(no_background == "error" && snr_thresh > 0)
                stop(msg, "Without background pixels there is nothing to ",
                     "reference the signal against. Label background pixels ",
                     "with selectTissuePixels(), or pass ",
                     "no_background = \"all_na\" / \"copy\" to say what should ",
                     "happen instead.", call. = FALSE)

              if(no_background == "all_na" && snr_thresh > 0){
                warning(msg, "snr set to NA for every pixel.", call. = FALSE)
                spectra(MSIobject, "snr") <- matrix(NA_real_,
                                                     nrow = nrow(MSIobject),
                                                     ncol = ncol(MSIobject))
              } else {
                # Kept for workflows that deliberately run unfiltered, but the
                # snr layer then holds response, not a ratio.
                warning(msg, "Returning ", val_slot, " as snr: the snr slot ",
                        "will hold measured response, not a ratio.",
                        call. = FALSE)
                spectra(MSIobject, "snr") <- spectra(MSIobject, val_slot)
              }
              return(MSIobject)
            }

            #Set background and tissue pixels
            bg_pixels = which(pData(MSIobject)[[pixel_header]] %in% .bg)
            tissue_pixels = which(pData(MSIobject)[[pixel_header]] %in% tissue)

            spectra(MSIobject, "snr") = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))

            # Features whose median noise level is just the imputed value.
            imputed_feats <- character()

            # Iterate over features in study
            for(mz_ind in seq_len(nrow(fData(MSIobject)))){

              # Save background response vector
              noise_vec = spectraData(MSIobject)[[val_slot]][mz_ind, bg_pixels]

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

              # More than half the background zero: the median IS the imputed
              # value, so the imputation rule sets the threshold, not the
              # background. Only the median can land here -- the mean of a
              # vector holding a non-zero value is above min/10.
              if (noise_level == mv_impute)
                imputed_feats <- c(imputed_feats, feat_nm)
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

            # Nothing is filtered at 0, so the noise level would not matter.
            if (length(imputed_feats) && snr_thresh > 0)
              .warn_imputed_noise(MSIobject, imputed_feats, snr_thresh)

            return(MSIobject)
          })

#' Warn that a median noise level is the imputed value
#'
#' A classed warning, so [generateTxtImages()] can collect one per acquisition
#' and report them once for the study rather than once per SNR threshold.
#'
#' @param MSIobject The object `int2SNR()` was given.
#' @param feats Names of the features affected.
#' @param thr The global SNR threshold, for the message.
#' @keywords internal
#' @noRd
.warn_imputed_noise <- function(MSIobject, feats, thr) {
  run <- if ("run" %in% names(pData(MSIobject)))
           unique(as.character(pData(MSIobject)$run)) else character()
  n <- nrow(MSIobject)
  shown <- paste(utils::head(feats, 5), collapse = ", ")
  if (length(feats) > 5) shown <- paste0(shown, ", ...")
  msg <- sprintf(paste0(
    "int2SNR%s: in %d of %d features (%s) more than half the background ",
    "pixels are zero or missing, so the median noise level is the imputed value -- one ",
    "tenth of the smallest non-zero background value. An SNR of %g there ",
    "passes anything above %g x the dimmest non-zero background pixel. ",
    "average = \"mean\" (parameters: average_method: mean) is set by the ",
    "background's own values instead. See ?int2SNR."),
    if (length(run)) paste0(" (", paste(run, collapse = ", "), ")") else "",
    length(feats), n, shown, thr, thr / 10)
  cond <- structure(
    class = c("quantMSImageR_imputed_noise", "warning", "condition"),
    list(message = msg, call = NULL, run = run, features = feats,
         n_features = n))
  warning(cond)
}
