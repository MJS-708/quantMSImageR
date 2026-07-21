setGeneric("summarise_cal_levels", function(MSIobject, ...) standardGeneric("summarise_cal_levels"))

#' Summarise the response at each calibration level
#'
#' Averages the signal across the pixels of each calibration spot, giving one
#' response value per standard per level. This is the input
#' [create_cal_curve()] fits its models to.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` whose `pData()` carries the
#'   calibration labels and spot identifiers.
#' @param cal_metadata Data frame of calibration metadata, with columns
#'   `lipid` (feature name as in `fData()$name`), `identifier` (maps to
#'   `pData()`), `amount_pg` (amount of standard deposited at each spot) and
#'   `level` (dilution level).
#' @param val_slot Character. Spectra slot to summarise (default `"intensity"`).
#' @param cal_header Character. Column of `pData()` used to select calibration
#'   pixels (default `"sample_type"`).
#' @param cal_label Character. Value of `cal_header` marking calibration pixels
#'   (default `"Cal"`).
#' @param id Character. Column present in both `cal_metadata` and `pData()`
#'   that identifies a unique calibration spot (default `"identifier"`).
#' @return The input object with `cal_response_data` populated: mean response
#'   per pixel for every (standard, calibration level) pair, together with the
#'   per-level pixel counts. Amounts are expressed in **pg per pixel**.
#'
#' @examples
#' cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
#' cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
#'           "quant_MSImagingExperiment")
#' cal_metadata <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
#' cal <- summarise_cal_levels(cal, cal_metadata, val_slot = "intensity",
#'                             cal_label = "Cal", id = "identifier")
#'
#' @family calibration
#' @aliases summarise_cal_levels
#' @export
setMethod("summarise_cal_levels", "quant_MSImagingExperiment",
          function(MSIobject, cal_metadata, val_slot = "intensity", cal_header = "sample_type", cal_label = "Cal", id = "identifier"){

            MSIobject@calibrationInfo@cal_metadata = cal_metadata

            # create pixel data to associate pixel indices to cal levels.
            # Selection honours `cal_header` (the pixel-type column) so the
            # calibration label need not live in a column named "sample_type".
            pixel_data = data.frame(pData(MSIobject))
            pixel_data$pixel_ind = seq_len(nrow(pixel_data))
            pixel_data = pixel_data[pixel_data[[cal_header]] == cal_label &
                                      !is.na(pixel_data[[id]]), ]

            # Create output response df
            response_df = tibble::tibble(cal_spot = unique(pixel_data[[id]]),
                                 response_perpixel = NA,
                                 pixels = NA) |>
              dplyr::left_join(MSIobject@calibrationInfo@cal_metadata, by=c("cal_spot" = id)) |>
              dplyr::select(dplyr::any_of(c("cal_spot", "response_perpixel", "pixels", "level", "lipid", "amount_pg")))

            for(i in seq_len(nrow(response_df))){

              # Select feature
              lipid_n = response_df$lipid[i]
              lipid_ind = which(fData(MSIobject)$name == lipid_n)

              # select pixels
              cal_n = response_df$cal_spot[i]
              inds = which(pixel_data[[id]] == cal_n)
              pixels = as.numeric(pixel_data$pixel_ind[inds])


              response_df$pixels[i] = length(pixels)

                ints = spectraData(MSIobject)[[val_slot]][lipid_ind, pixels]
              ints = replace(ints, ints ==0, NA)

              response_df$response_perpixel[i] = mean(ints, na.rm = TRUE)

            }

            response_df = dplyr::mutate(response_df, pg_perpixel = amount_pg / pixels)


            MSIobject@calibrationInfo@cal_response_data = response_df


            return(MSIobject)
          })
