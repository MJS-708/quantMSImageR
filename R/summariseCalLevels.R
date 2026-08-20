setGeneric("summariseCalLevels", function(MSIobject, ...) standardGeneric("summariseCalLevels"))

#' Summarise the response at each calibration level
#'
#' Averages the signal across the pixels of each calibration spot, giving one
#' response value per standard per level. This is the input
#' [createCalCurve()] fits its models to.
#'
#' @details
#' Every spot is also summarised for dispersion, because a mean and a pixel
#' count cannot show whether a calibration level was measured reproducibly. The
#' returned table carries `sd_response` and `cv_response` alongside the mean, so
#' an imprecise level is visible rather than averaged away.
#'
#' Problems in the calibration design are errors rather than quiet `NA`s: an
#' analyte that is not in the object, an identifier that matches no pixel, a
#' spot with no usable response, and duplicated identifiers that disagree about
#' amount or level all stop the function. Each of these would otherwise produce
#' a curve that looks fitted but is not anchored to the data it claims.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` whose `pData()` carries the
#'   calibration labels and spot identifiers.
#' @param cal_metadata Data frame of calibration metadata, with columns
#'   `analyte` (feature name as in `fData()$name`), `identifier` (maps to
#'   `pData()`), `amount_pg` (amount of standard deposited at each spot) and
#'   `level` (dilution level). A legacy `lipid` column is read as `analyte`.
#' @param val_slot Character. Spectra slot to summarise (default `"intensity"`).
#' @param cal_header Character. Column of `pData()` used to select calibration
#'   pixels (default `"sample_type"`).
#' @param cal_label Character. Value of `cal_header` marking calibration pixels
#'   (default `"Cal"`).
#' @param id Character. Column present in both `cal_metadata` and `pData()`
#'   that identifies a unique calibration spot (default `"identifier"`).
#' @return The input object with `cal_response_data` populated: one row per
#'   (standard, calibration spot) pair holding `response_perpixel` (the mean),
#'   `median_response`, `sd_response`, `cv_response`, `pixels`,
#'   `n_nonmissing` and `pg_perpixel`. Amounts are expressed in **pg per
#'   pixel**.
#'
#' @examples
#' cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
#' cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
#'           "quant_MSImagingExperiment")
#' cal_metadata <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
#' cal <- summariseCalLevels(cal, cal_metadata, val_slot = "intensity",
#'                             cal_label = "Cal", id = "identifier")
#'
#' @family calibration
#' @aliases summariseCalLevels
#' @export
setMethod("summariseCalLevels", "quant_MSImagingExperiment",
          function(MSIobject, cal_metadata, val_slot = "intensity", cal_header = "sample_type", cal_label = "Cal", id = "identifier"){

            # `lipid` was the original column name. It described the assay
            # rather than the role, so metadata for a non-lipid panel read
            # oddly; both spellings are accepted and `analyte` is emitted.
            cal_metadata <- as.data.frame(cal_metadata)
            if(!"analyte" %in% names(cal_metadata) &&
               "lipid" %in% names(cal_metadata))
              names(cal_metadata)[names(cal_metadata) == "lipid"] <- "analyte"

            .need <- c(id, "analyte", "amount_pg", "level")
            .miss <- setdiff(.need, names(cal_metadata))
            if(length(.miss))
              stop("summariseCalLevels: cal_metadata is missing column(s): ",
                   paste(.miss, collapse = ", "), ". Columns present: ",
                   paste(names(cal_metadata), collapse = ", "), ".",
                   call. = FALSE)

            # One identifier must describe one spot. Repeated rows for the same
            # (identifier, analyte) that disagree would make the summary depend
            # on row order.
            .key <- paste(cal_metadata[[id]], cal_metadata$analyte, sep = "\r")
            .dup <- .key[duplicated(.key)]
            if(length(.dup)){
              .cmp <- cal_metadata[.key %in% .dup, c("amount_pg", "level")]
              if(nrow(unique(cbind(key = .key[.key %in% .dup], .cmp))) !=
                 length(unique(.dup)))
                stop("summariseCalLevels: cal_metadata has duplicated ",
                     id, "/analyte rows that disagree about amount_pg or ",
                     "level: ", paste(unique(sub("\r", " / ", .dup)),
                                      collapse = ", "), ".", call. = FALSE)
            }

            MSIobject@calibrationInfo@cal_metadata = cal_metadata

            # create pixel data to associate pixel indices to cal levels.
            # Selection honours `cal_header` (the pixel-type column) so the
            # calibration label need not live in a column named "sample_type".
            pixel_data = data.frame(pData(MSIobject))
            if(!cal_header %in% names(pixel_data))
              stop("summariseCalLevels: pData() has no column '", cal_header,
                   "'. Columns present: ",
                   paste(names(pixel_data), collapse = ", "), ".", call. = FALSE)
            pixel_data$pixel_ind = seq_len(nrow(pixel_data))
            pixel_data = pixel_data[pixel_data[[cal_header]] == cal_label &
                                      !is.na(pixel_data[[id]]), ]

            if(nrow(pixel_data) == 0)
              stop("summariseCalLevels: no pixels are labelled '", cal_label,
                   "' in pData()$", cal_header, ".", call. = FALSE)

            # Every analyte named in the metadata must exist in the object,
            # otherwise its curve would be fitted from nothing.
            .fnames <- as.character(fData(MSIobject)$name)
            .orphan <- setdiff(unique(as.character(cal_metadata$analyte)), .fnames)
            if(length(.orphan))
              stop("summariseCalLevels: cal_metadata names analyte(s) that ",
                   "are not features of this object: ",
                   paste(.orphan, collapse = ", "), ".", call. = FALSE)

            # Create output response df
            response_df = tibble::tibble(cal_spot = unique(pixel_data[[id]])) |>
              dplyr::left_join(MSIobject@calibrationInfo@cal_metadata, by=c("cal_spot" = id)) |>
              dplyr::select(dplyr::any_of(c("cal_spot", "level", "analyte", "amount_pg")))

            .n <- nrow(response_df)
            response_df$response_perpixel <- NA_real_
            response_df$median_response   <- NA_real_
            response_df$sd_response       <- NA_real_
            response_df$cv_response       <- NA_real_
            response_df$pixels            <- NA_integer_
            response_df$n_nonmissing      <- NA_integer_

            for(i in seq_len(.n)){

              # Select feature
              analyte_n = response_df$analyte[i]
              analyte_ind = which(.fnames == analyte_n)

              # select pixels
              cal_n = response_df$cal_spot[i]
              inds = which(pixel_data[[id]] == cal_n)
              pixels = as.numeric(pixel_data$pixel_ind[inds])

              if(length(pixels) == 0)
                stop("summariseCalLevels: calibration spot '", cal_n,
                     "' matches no pixel in pData()$", id, ".", call. = FALSE)

              response_df$pixels[i] = length(pixels)

              ints = spectraData(MSIobject)[[val_slot]][analyte_ind, pixels]
              # A recorded zero in MRM means "nothing recorded", so it is not a
              # measured response and must not pull the spot mean down.
              ints = replace(ints, ints == 0, NA)

              .ok <- is.finite(ints)
              response_df$n_nonmissing[i] = sum(.ok)

              if(!any(.ok))
                stop("summariseCalLevels: calibration spot '", cal_n,
                     "' has no usable response for '", analyte_n,
                     "' (all pixels are zero or NA in slot '", val_slot,
                     "'). Drop the spot from cal_metadata or check the ROI.",
                     call. = FALSE)

              response_df$response_perpixel[i] = mean(ints, na.rm = TRUE)
              response_df$median_response[i]   = stats::median(ints, na.rm = TRUE)
              response_df$sd_response[i]       = stats::sd(ints, na.rm = TRUE)
            }

            response_df = dplyr::mutate(response_df,
                                        pg_perpixel = amount_pg / pixels,
                                        cv_response = sd_response / response_perpixel)

            MSIobject@calibrationInfo@cal_response_data = response_df

            return(MSIobject)
          })
