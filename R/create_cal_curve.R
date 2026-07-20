setGeneric("create_cal_curve", function(MSIobject, ...) standardGeneric("create_cal_curve"))

#' Function to create calibration curves (response v concentration, where concentration is pg/pixel)
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` whose `@calibrationInfo@cal_response_data`
#'   slot has been populated by [summarise_cal_levels()].
#' @param cal_type string of approach to generate claibration curve - 'std_addition' if standards are on tissue and 'cal' if direct onto glass slide.
#' @param level Column header to find background label from 'MSIobject@calibrationInfo@cal_response_data'
#' @param background string referring to background level from "level" label in 'MSIobject@calibrationInfo@cal_response_data'.
#' @return MSIobject with slots updated for i) cal_list - List of linear models for each m/z (response v concentration, where concentration is pg/pixel) and ii) r2 values for each calibration iii) calibration metadata
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
