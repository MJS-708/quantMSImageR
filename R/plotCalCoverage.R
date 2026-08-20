setGeneric("plotCalCoverage", function(MSIobject, ...) standardGeneric("plotCalCoverage"))

#' Check that measured pixels fall within the calibrated range
#'
#' Converting a response into a calibrated amount is only trustworthy where the
#' standards actually constrain the fit. A calibration whose levels sit above or
#' below the signal measured in tissue is extrapolation, not quantification.
#' [int2conc()] reports the proportion of extrapolated pixels per feature, but
#' without an explicit coverage check that is still easy to overlook -- and the
#' plot shows *where* the mismatch is, not just how much of it there is.
#'
#' This function plots the two together on log axes so the overlap can be
#' judged directly. The lower panel is the calibration itself: the standard
#' spots and the fitted model that `int2conc()` inverts. The upper panel is the
#' distribution of the object's own pixels on the same x-axis, with the
#' calibration levels marked. When the histogram sits inside the span of those
#' marks, every converted pixel is interpolated between real standards.
#'
#' @details
#' The standards are read back from each model's own model frame
#' (`stats::model.frame`) rather than from `cal_response_data`. This matters for
#' `cal_type = "std_addition"`, where [createCalCurve()] shifts the amount
#' axis by the estimated endogenous amount before fitting: the stored response
#' data is on the unshifted axis, so plotting it against the fitted line would
#' put the points and the model on different scales. The model frame always
#' holds the points the line was actually fitted to.
#'
#' The fitted model is linear in amount, so on log axes it is drawn as a curve
#' evaluated across the plotted range rather than as a straight line. Any part
#' of it predicting a non-positive response has no log-scale counterpart and is
#' dropped.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` that has been through
#'   [createCalCurve()] and [int2conc()], so it carries both the fitted models
#'   (see [calibrationModels()]) and a calibrated-amount slot.
#' @param features Features to plot: a character vector of names as in
#'   `fData(MSIobject)$name`, or numeric row indices. Defaults to every feature
#'   that has a calibration model.
#' @param val_slot Character. Spectra slot holding the measured response
#'   (default `"intensity"`).
#' @param conc_slot Character. Spectra slot holding the calibrated amount
#'   estimate written by [int2conc()]. Defaults to `"pg_pixel"`; objects
#'   calibrated with earlier versions carry `"conc - pg/pixel"`, which is
#'   matched automatically.
#' @param log_base Numeric. Base for both axes (default `2`).
#' @param bins Integer. Number of histogram bins for the pixel distribution
#'   (default `40`).
#' @param cal_colour Colour used for the standards and the level markers
#'   (default `"red"`).
#' @param cal_alpha Numeric. Opacity of the calibration-level marker lines
#'   (default `0.4`), kept faint so they read as a reference grid rather than
#'   competing with the data.
#' @param cal_linetype Line type for the calibration-level markers (default
#'   `"dotted"`).
#' @param ... Additional arguments (currently unused).
#'
#' @return A `patchwork` object: the pixel histogram stacked above the
#'   calibration curve, sharing an x-axis. Printing it draws the figure; it can
#'   also be further modified with `ggplot2` layers.
#'
#' @examples
#' cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
#' cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
#'           "quant_MSImagingExperiment")
#' cal_metadata <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
#' cal <- summariseCalLevels(cal, cal_metadata, val_slot = "intensity",
#'                             cal_label = "Cal", id = "identifier")
#' cal <- createCalCurve(cal, cal_type = "cal")
#' cal <- int2conc(cal, pixel_header = "sample_type", pixels = "Tissue")
#'
#' plotCalCoverage(cal, val_slot = "intensity")
#'
#' @family calibration
#' @aliases plotCalCoverage
#' @export
setMethod("plotCalCoverage", "quant_MSImagingExperiment",
          function(MSIobject, features = NULL, val_slot = "intensity",
                   conc_slot = "pg_pixel", log_base = 2, bins = 40,
                   cal_colour = "red", cal_alpha = 0.4,
                   cal_linetype = "dotted", ...){

            cal_list = calibrationModels(MSIobject)
            if(length(cal_list) == 0)
              stop("plotCalCoverage: no calibration models found. Run ",
                   "createCalCurve() first.", call. = FALSE)

            slots_present = names(spectraData(MSIobject))

            # Accept either spelling of the calibrated-amount slot, so objects
            # written by earlier versions still plot.
            .cs = .conc_slots(conc_slot)
            .cs = .cs[.cs %in% slots_present]
            if(length(.cs) == 0)
              stop("plotCalCoverage: no '", conc_slot, "' slot. Run int2conc() ",
                   "first. Slots present: ", paste(slots_present, collapse = ", "),
                   call. = FALSE)
            conc_slot = .cs[1]
            if(!val_slot %in% slots_present)
              stop("plotCalCoverage: no '", val_slot, "' slot. Slots present: ",
                   paste(slots_present, collapse = ", "), call. = FALSE)

            feat_nms  = as.character(fData(MSIobject)$name)
            available = intersect(feat_nms, names(cal_list))

            if(is.null(features)){
              features = available
            } else {
              if(is.numeric(features)) features = feat_nms[features]
              features = as.character(features)
              unknown  = setdiff(features, available)
              if(length(unknown) > 0)
                stop("plotCalCoverage: no calibration model for feature(s) ",
                     paste(unknown, collapse = ", "), ". Calibrated features: ",
                     paste(available, collapse = ", "), call. = FALSE)
            }

            if(length(features) == 0)
              stop("plotCalCoverage: none of the object's features have a ",
                   "calibration model.", call. = FALSE)

            lg = function(v) log(v, base = log_base)

            px_l = list(); st_l = list(); ln_l = list()

            for(nm in features){

              i    = match(nm, feat_nms)
              eqn  = cal_list[[nm]]

              conc = as.numeric(spectra(MSIobject, conc_slot)[i, ])
              resp = as.numeric(spectra(MSIobject, val_slot)[i, ])
              keep = is.finite(conc) & conc > 0 & is.finite(resp) & resp > 0
              if(!any(keep)){
                warning("plotCalCoverage: no positive finite pixels for '", nm,
                        "'; skipping.", call. = FALSE)
                next
              }
              px_l[[nm]] = data.frame(feature = nm, x = lg(conc[keep]),
                                      y = lg(resp[keep]))

              # Points the model was actually fitted to (see @details).
              md    = stats::model.frame(eqn)
              s_x   = md[["pg_perpixel"]]
              s_y   = md[["response_perpixel"]]
              s_ok  = is.finite(s_x) & s_x > 0 & is.finite(s_y) & s_y > 0
              st_l[[nm]] = data.frame(feature = nm, x = lg(s_x[s_ok]),
                                      y = lg(s_y[s_ok]))

              rng  = range(c(px_l[[nm]]$x, st_l[[nm]]$x))
              xx   = seq(rng[1], rng[2], length.out = 200)
              pred = stats::coef(eqn)[1] + stats::coef(eqn)[2] * log_base ^ xx
              ln_l[[nm]] = data.frame(feature = nm, x = xx[pred > 0],
                                      y = lg(pred[pred > 0]))
            }

            if(length(px_l) == 0)
              stop("plotCalCoverage: no feature had usable pixel values.",
                   call. = FALSE)

            pixels = do.call(rbind, px_l)
            stds   = do.call(rbind, st_l)
            line   = do.call(rbind, ln_l)

            # A single shared x range keeps the two panels (and every facet)
            # aligned, which is the whole point of the figure.
            xlim = range(c(pixels$x, stds$x))

            ax_lab = function(txt) bquote(log[.(log_base)] ~ .(txt))

            # Publication styling: no grid in the panel, solid black axes, and
            # calibration levels as faint dotted references behind the data.
            cov_theme = ggplot2::theme_minimal() +
              ggplot2::theme(
                panel.grid   = ggplot2::element_blank(),
                panel.border = ggplot2::element_blank(),
                axis.line    = ggplot2::element_line(colour = "black",
                                                     linewidth = 0.7),
                axis.ticks   = ggplot2::element_line(colour = "black",
                                                     linewidth = 0.5),
                strip.text   = ggplot2::element_text(face = "bold")
              )

            # geom_vline must be data-driven, not a bare xintercept vector, or
            # every facet would show every analyte's calibration levels.
            cal_lines = ggplot2::geom_vline(
              data = stds, ggplot2::aes(xintercept = x), colour = cal_colour,
              linewidth = 0.4, linetype = cal_linetype, alpha = cal_alpha,
              inherit.aes = FALSE)

            top = ggplot2::ggplot(pixels, ggplot2::aes(x = x)) +
              cal_lines +
              ggplot2::geom_histogram(bins = bins, fill = "grey55",
                                      colour = NA) +
              ggplot2::coord_cartesian(xlim = xlim) +
              ggplot2::labs(x = NULL, y = "pixel count") +
              cov_theme +
              ggplot2::theme(axis.text.x = ggplot2::element_blank())

            bottom = ggplot2::ggplot(stds, ggplot2::aes(x = x, y = y)) +
              cal_lines +
              ggplot2::geom_line(data = line, linetype = "dashed") +
              ggplot2::geom_point(colour = cal_colour, size = 2) +
              ggplot2::coord_cartesian(xlim = xlim) +
              ggplot2::labs(x = ax_lab("(pg / pixel)"),
                            y = ax_lab("(response)")) +
              cov_theme

            # Free the y scale per analyte but keep x shared, so panels stay
            # comparable across features and aligned between the two plots.
            if(length(unique(pixels$feature)) > 1){
              top    = top    + ggplot2::facet_wrap(~ feature, scales = "free_y",
                                                    nrow = 1)
              bottom = bottom + ggplot2::facet_wrap(~ feature, scales = "free_y",
                                                    nrow = 1)
            }

            patchwork::wrap_plots(top, bottom, ncol = 1, heights = c(1, 3))
          })
