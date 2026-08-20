setGeneric("createCalCurve", function(MSIobject, ...) standardGeneric("createCalCurve"))

#' Fit per-analyte calibration models
#'
#' Fits a linear response-versus-amount model for each standard, which
#' [int2conc()] later inverts to turn a measured response into an estimated
#' amount per pixel.
#'
#' @details
#' `cal_type` has no default and must be chosen explicitly, because the two
#' modes apply materially different processing.
#'
#' `cal_type = "cal"` fits the summarised response directly against the amount
#' deposited, appropriate for standards spotted onto the slide beside the
#' tissue (on-slide / external calibration).
#'
#' `cal_type = "std_addition"` is classical **standard addition**, for standards
#' deposited onto tissue where endogenous analyte is already present. It:
#' \enumerate{
#'   \item fits an unweighted `response ~ deposited amount` across all levels;
#'   \item takes the x-intercept of that line -- the amount at which predicted
#'     response is zero -- as an estimate of the endogenous amount already in
#'     the tissue (this value is negative when endogenous signal is present);
#'   \item subtracts it from every deposited amount, which shifts the x-axis so
#'     it expresses **total** amount present (endogenous + added);
#'   \item drops the rows at the `background` level and refits on that shifted
#'     axis.
#' }
#'
#' @section Weighting:
#' Calibration response in MSI is usually heteroscedastic -- absolute scatter
#' grows with amount -- so an unweighted fit lets the top standards dominate and
#' the low end, where most tissue pixels sit, is fitted worst. `weighting`
#' therefore defaults to `"1/x"`; `"1/x2"` weights the low end harder still, and
#' `"none"` is ordinary least squares.
#'
#' A zero (blank) level has no `1/x` weight. It is given the weight of the
#' lowest positive level rather than a weight derived from a substituted tiny
#' amount, which would let a single blank determine the whole regression.
#'
#' Weighting is an empirical model choice, not a property of the data. Check it
#' against residual structure, back-calculated accuracy at each level and
#' replicate precision rather than against R-squared -- see
#' [calibrationDiagnostics()] and [plotCalCoverage()].
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` whose `cal_response_data` has
#'   been populated by [summariseCalLevels()].
#' @param cal_type Character, required. `"cal"` for standards deposited onto the
#'   slide, `"std_addition"` for standard addition on tissue. See Details.
#' @param level Character. Column of `cal_response_data` holding the level
#'   labels (default `"level"`).
#' @param background Character. Value of `level` marking the background level
#'   subtracted when `cal_type = "std_addition"`.
#' @param weighting Character. Regression weights: `"1/x"` (default), `"1/x2"`
#'   or `"none"`. See the Weighting section.
#' @return The input object with `cal_list` (one `lm` per standard, response
#'   versus amount in pg per pixel), `r2_df` (fit R-squared per standard) and
#'   the calibration metadata populated.
#'
#' @examples
#' cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
#' cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
#'           "quant_MSImagingExperiment")
#' cal_metadata <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
#' cal <- summariseCalLevels(cal, cal_metadata, val_slot = "intensity",
#'                             cal_label = "Cal", id = "identifier")
#' cal <- createCalCurve(cal, cal_type = "cal")
#'
#' @family calibration
#' @aliases createCalCurve
#' @export
setMethod("createCalCurve", "quant_MSImagingExperiment",
          function(MSIobject, cal_type, level = "level", background = "background",
                   weighting = c("1/x", "none", "1/x2")){

            weighting <- match.arg(weighting)

            # No default: "cal" and "std_addition" apply materially different
            # processing, and silently assuming one of them is a quantitative
            # risk rather than a convenience.
            if (missing(cal_type))
              stop("createCalCurve: `cal_type` must be given explicitly -- ",
                   "\"cal\" for standards deposited on the slide, or ",
                   "\"std_addition\" for standard addition on tissue.",
                   call. = FALSE)
            cal_type <- match.arg(cal_type, c("cal", "std_addition"))

            cal_data = MSIobject@calibrationInfo@cal_response_data

            features = unique(cal_data$analyte)

            # Set outputs
            cal_list = list()
            r2_df = tibble::tibble(feature= features, r2 = NA)

            for(i in seq_along(features)){

              feature = features[i]
              cal_subset = subset(cal_data, analyte == feature)

              if(cal_type == "std_addition"){

                eqn = stats::lm(response_perpixel~pg_perpixel, data = cal_subset, na.action = stats::na.exclude)

                # Calculate background conc
                background_conc = chemCal::inverse.predict(eqn, 0)$Prediction

                # Update conc values (original conc - background)
                cal_subset$pg_perpixel = cal_subset$pg_perpixel - background_conc


                cal_subset = cal_subset[cal_subset[[level]] != background, ]

              }

              # Fit weights. A zero (blank) level has no 1/x weight, and the
              # obvious workaround -- substituting a tiny amount -- is not
              # harmless: 1/1e-9 is a weight of a billion, which lets the blank
              # dictate the entire regression. Give it the weight of the lowest
              # positive level instead, so it counts once rather than
              # overwhelmingly.
              w <- .cal_weights(cal_subset$pg_perpixel, weighting)

              # Update equation
              eqn = stats::lm(response_perpixel~pg_perpixel, data = cal_subset,
                              na.action = stats::na.exclude, weights = w)

              r2_df$r2[i] = summary(eqn)[["r.squared"]]
              cal_list[[feature]] = eqn
            }

            MSIobject@calibrationInfo@cal_list = cal_list
            MSIobject@calibrationInfo@r2_df = r2_df

            return(MSIobject)
          })
