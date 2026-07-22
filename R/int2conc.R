setGeneric("int2conc", function(MSIobject, ...) standardGeneric("int2conc"))

#' Convert response to calibrated amount estimates
#'
#' Inverts each analyte's calibration model to turn a measured response into an
#' **estimated amount per pixel**. Features without a calibration model are
#' dropped.
#'
#' @details
#' Mass spectrometry imaging measures analyte response at each pixel rather than
#' analyte amount directly. The values produced here are calibrated *amount*
#' estimates (pg per pixel, and pg per mm-squared where the pixel size is
#' known), not volumetric concentrations. They are conditional on the
#' calibration design, matrix matching, acquisition conditions and the validity
#' of the fitted model -- see [plot_cal_coverage()] for the coverage check that
#' should accompany them.
#'
#' @section Areal conversion:
#' `pg_pixel` is inverted directly from the calibration model. `pg_mm2` is then
#' derived from the pixel size recorded in the acquisition metadata:
#'
#' \deqn{pg/mm^2 = (pg/pixel) \times (1000 / s)^2}
#'
#' where `s` is `experimentData(MSIobject)$pixelSize`, **expected in
#' micrometres**, and the factor 1000 converts micrometres to millimetres. A
#' single scalar is assumed, i.e. square pixels of side `s`; rectangular pixels
#' are not currently supported, and supplying `s` in millimetres would inflate
#' the result by a factor of 10^6. When `pixelSize` is missing or not a positive
#' number the `pg_mm2` slot is omitted and a warning is issued.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` carrying calibration models
#'   fitted by [create_cal_curve()].
#' @param pixels Character. Label(s) in `pixel_header` marking the pixels to
#'   quantify (default `"tissue_pixels"`, the imaging convention). Pass
#'   `pixels = "Tissue"` for a calibration acquisition.
#' @param pixel_header Character. Column of `pData()` holding the pixel-type
#'   labels (default `"sample_name"`). Pass `pixel_header = "sample_type"` for
#'   a calibration acquisition.
#' @param val_slot Character. Spectra slot holding the measured response
#'   (default `"intensity"`).
#' @param max_out_of_range Numeric in `[0, 1]`. Warn when more than this
#'   proportion of a feature's pixels fall outside the calibrated range, i.e.
#'   are extrapolated rather than interpolated (default `0.1`, so at least 90%
#'   of pixels must be interpolated between real standards). Any extrapolation
#'   at all is reported by `message()`; this controls only the escalation to a
#'   warning. Set to `1` to silence it.
#' @return The object subset to `pixels`, with two spectra slots added:
#'   `pg_pixel` (estimated amount per pixel) and, when
#'   `experimentData(MSIobject)$pixelSize` is available, `pg_mm2` (estimated
#'   areal amount). Features with no calibration model are removed, and
#'   `r2_df` gains an `out_of_range` column giving the proportion of each
#'   feature's pixels that were extrapolated beyond the calibrated range.
#'
#' @examples
#' cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
#' cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
#'           "quant_MSImagingExperiment")
#' cal_metadata <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
#' cal <- summarise_cal_levels(cal, cal_metadata, val_slot = "intensity",
#'                             cal_label = "Cal", id = "identifier")
#' cal <- create_cal_curve(cal, cal_type = "cal")
#'
#' # A calibration acquisition labels its pixels in `sample_type`, so the
#' # imaging defaults are overridden here
#' cal <- int2conc(cal, pixel_header = "sample_type", pixels = "Tissue")
#'
#' @family calibration
#' @aliases int2conc
#' @export
setMethod("int2conc", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity", pixel_header = "sample_name",
                   pixels = "tissue_pixels", max_out_of_range = 0.1){

            cal_list = MSIobject@calibrationInfo@cal_list
            pixel_inds = which(pData(MSIobject)[[pixel_header]] %in% pixels)

            if(length(pixel_inds) == 0)
              stop("int2conc: no pixels matched pixels = '",
                   paste(pixels, collapse = "', '"), "' in pixel_header = '",
                   pixel_header, "'. Labels present: ",
                   paste(unique(as.character(pData(MSIobject)[[pixel_header]])),
                         collapse = ", "), call. = FALSE)

            MSIobject = MSIobject[, pixel_inds]

            # Pixel size drives the pg/pixel -> pg/mm2 conversion. It comes from
            # the acquisition metadata, which some objects (e.g. ones assembled
            # by hand) never set -- convert once here so a missing value is
            # reported clearly instead of silently producing a zero-length
            # replacement inside the loop.
            px_size = suppressWarnings(as.numeric(experimentData(MSIobject)$pixelSize))
            has_px  = length(px_size) == 1L && !is.na(px_size) && px_size > 0
            if(!has_px)
              warning("int2conc: experimentData(MSIobject)$pixelSize is missing or ",
                      "not a positive number, so the '", .SLOT_MM2, "' slot cannot ",
                      "be computed. Only '", .SLOT_PG, "' will be added.",
                      call. = FALSE)

            spectra(MSIobject, .SLOT_PG) = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))
            if(has_px)
              spectra(MSIobject, .SLOT_MM2) = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))

            no_cal_indices = c()
            oor_frac       = list()

            for(i in seq_len(nrow(fData(MSIobject)))){

              # Look up the calibration model by feature name (cal_list is keyed
              # by name in create_cal_curve). Only fall back to the positional
              # entry when cal_list is unnamed (older objects) -- never index a
              # named list positionally, which would silently apply a different
              # analyte's curve. Features with no curve are dropped below, so a
              # study may carry more features than were calibrated.
              feat_nm = fData(MSIobject)$name[i]
              eqn = if(!is.null(names(cal_list))){
                      if(feat_nm %in% names(cal_list)) cal_list[[feat_nm]] else NULL
                    } else if(i <= length(cal_list)) cal_list[[i]] else NULL

              if(is.null(eqn) || typeof(eqn) != "list"){

                message(sprintf("No calibration curve for feature %s", fData(MSIobject)$name[i]))
                no_cal_indices = c(no_cal_indices, i)

              } else{

                pixel_ints = spectraData(MSIobject)[[val_slot]][i,]
                pixel_concs = vapply(pixel_ints,
                                     function(y) chemCal::inverse.predict(eqn, y)$Prediction,
                                     numeric(1))

                # Inverting the model outside the calibrated range is
                # extrapolation, not quantification. Values are kept (dropping
                # them silently would be worse) but the proportion affected is
                # recorded per feature and reported below.
                cal_rng = range(stats::model.frame(eqn)[["pg_perpixel"]], na.rm = TRUE)
                finite  = is.finite(pixel_concs)
                oor     = finite & (pixel_concs < cal_rng[1] | pixel_concs > cal_rng[2])
                oor_frac[[feat_nm]] = if(any(finite)) sum(oor) / sum(finite) else NA_real_

                spectra(MSIobject, .SLOT_PG)[i,] = pixel_concs

                if(has_px){
                  mm2_scalar = (1000 / px_size) ^2
                  spectra(MSIobject, .SLOT_MM2)[i,] = pixel_concs * mm2_scalar
                }

              }
            }

            if(length(no_cal_indices) > 0){
              MSIobject = MSIobject[-no_cal_indices, ]
            }

            # Report extrapolation. Anything out of range is worth a message;
            # a feature that is mostly out of range is not being quantified at
            # all, so that escalates to a warning unless explicitly allowed.
            oor_vec = unlist(oor_frac)
            if(length(oor_vec) > 0){

              any_oor = oor_vec[!is.na(oor_vec) & oor_vec > 0]
              if(length(any_oor) > 0)
                message("int2conc: pixels outside the calibrated range (",
                        "extrapolated, not interpolated):\n",
                        paste(sprintf("  %-14s %.1f%%", names(any_oor),
                                      100 * any_oor), collapse = "\n"))

              bad = oor_vec[!is.na(oor_vec) & oor_vec > max_out_of_range]
              if(length(bad) > 0 && max_out_of_range < 1)
                warning("int2conc: more than ", round(100 * max_out_of_range),
                        "% of pixels fall outside the calibrated range for: ",
                        paste(sprintf("%s (%.1f%%)", names(bad), 100 * bad),
                              collapse = ", "),
                        ". These amounts are extrapolated and should not be ",
                        "treated as quantitative -- check plot_cal_coverage(). ",
                        "Set max_out_of_range = 1 to silence this.",
                        call. = FALSE)
            }

            MSIobject@calibrationInfo@r2_df$out_of_range =
              unname(oor_vec[match(MSIobject@calibrationInfo@r2_df$feature,
                                   names(oor_vec))])

            return(MSIobject)
          })
