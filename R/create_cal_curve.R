setGeneric("create_cal_curve", function(MSIobject, ...) standardGeneric("create_cal_curve"))

#' Fit per-analyte calibration models
#'
#' Fits a linear response-versus-amount model for each standard, which
#' [int2conc()] later inverts to turn a measured response into an estimated
#' amount per pixel.
#'
#' @details
#' `cal_type = "cal"` fits the summarised response directly against the amount
#' deposited, appropriate for standards spotted onto the slide beside the
#' tissue (on-slide / external calibration).
#'
#' `cal_type = "std_addition"` first estimates the background contribution and
#' subtracts it from the deposited amounts before fitting, appropriate for
#' standards deposited onto tissue. Note that this is
#' **background-corrected on-tissue calibration** rather than classical standard
#' addition, which would estimate the endogenous amount from the x-intercept
#' across added amounts. The argument value is retained for compatibility.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` whose `cal_response_data` has
#'   been populated by [summarise_cal_levels()].
#' @param cal_type Character. `"cal"` for standards deposited onto the slide,
#'   `"std_addition"` for background-corrected on-tissue calibration.
#' @param level Character. Column of `cal_response_data` holding the level
#'   labels (default `"level"`).
#' @param background Character. Value of `level` marking the background level
#'   subtracted when `cal_type = "std_addition"`.
#' @return The input object with `cal_list` (one `lm` per standard, response
#'   versus amount in pg per pixel), `r2_df` (fit R-squared per standard) and
#'   the calibration metadata populated.
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
#' @family calibration
#' @aliases create_cal_curve
#' @export
setMethod("create_cal_curve", "quant_MSImagingExperiment",
          function(MSIobject, cal_type = "std_addition", level = "level", background = "background"){

            cal_data = MSIobject@calibrationInfo@cal_response_data

            features = unique(cal_data$lipid)

            # Set outputs
            cal_list = list()
            r2_df = tibble::tibble(feature= features, r2 = NA)

            for(i in seq_along(features)){

              feature = features[i]
              cal_subset = subset(cal_data, lipid == feature)

              if(cal_type == "std_addition"){

                eqn = stats::lm(response_perpixel~pg_perpixel, data = cal_subset, na.action = stats::na.exclude)

                # Calculate background conc
                background_conc = chemCal::inverse.predict(eqn, 0)$Prediction

                # Update conc values (original conc - background)
                cal_subset$pg_perpixel = cal_subset$pg_perpixel - background_conc


                cal_subset = cal_subset[cal_subset[[level]] != background, ]

              }

              cal_subset = dplyr::mutate(cal_subset, pg_perpixel = ifelse(pg_perpixel == 0, yes= 1e-9, no = pg_perpixel))

              # Update equation
              eqn = stats::lm(response_perpixel~pg_perpixel, data = cal_subset, na.action = stats::na.exclude, weights = (1/pg_perpixel))

              r2_df$r2[i] = summary(eqn)[["r.squared"]]
              cal_list[[feature]] = eqn
            }

            MSIobject@calibrationInfo@cal_list = cal_list
            MSIobject@calibrationInfo@r2_df = r2_df

            return(MSIobject)
          })
