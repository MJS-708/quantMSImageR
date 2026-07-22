#' Run the quantMSImageR example study
#'
#' Executes the full DESI-MRM pipeline on a small synthetic dataset bundled
#' with the package.  The dataset contains two samples -- `SampleA` (circular
#' tissue, sections 01-03) and `SampleB` (square tissue, sections 04-06) --
#' acquired on a panel of eight oxylipin features (seven analytes + one
#' internal standard, all negative-ion mode).
#'
#' Distinct tissue shapes per sample make the rendered ion images visually
#' easy to distinguish, and the bundled ion library carries a `Met-1` column
#' (`COX`, `LOX`, `IS`) that drives the row-split colour bar in the report's
#' main heatmap and the diagonal colour blocks in the correlation triangles.
#'
#' All file paths are resolved automatically via [system.file()], so the
#' example works on any machine without editing a YAML.
#'
#' @param render_report Logical. Render and open the HTML report in the RStudio
#'   viewer or default browser?  Requires pandoc. Default `TRUE`.
#' @param output_txt Logical. Write per-feature `.txt` image matrices to a
#'   temporary directory?  Default `FALSE`.
#' @param snr_thresh Numeric. SNR threshold passed to [generate_txt_images()].
#'   Default `3`.
#' @param shapes Character vector of tissue shapes to include -- any of
#'   `"circle"`, `"square"`. `NULL` (the default) uses both, i.e. the full
#'   two-group example (group A circle + group B square).
#' @param calibrate Logical. Also run the quantification demo on the
#'   bundled synthetic calibration standards (`summarise_cal_levels()` ->
#'   `create_cal_curve()` -> `int2conc()`)? In a real study this step is driven
#'   by the YAML `calibration:` block. Default `TRUE`.
#'
#' @return Invisibly returns the list produced by [generate_txt_images()], with
#'   an added `calibrated` element when `calibrate = TRUE`: the study sections
#'   with their tissue pixels converted to `pg_pixel` and `pg_mm2` layers,
#'   restricted to the analytes that have a standard.
#'   When `render_report = TRUE` a `report` element gives the path to the
#'   rendered HTML.
#'
#' @examples
#' # The pipeline itself, without rendering or opening the HTML report
#' \donttest{
#' res <- run_example(render_report = FALSE, shapes = "circle")
#' names(res)
#' }
#'
#' # The full demo: renders the report and opens it in the viewer/browser
#' if (interactive()) {
#'   res <- run_example()
#'   res$report        # path to the rendered HTML
#' }
#'
#' @family workflow
#' @export
run_example <- function(render_report = TRUE,
                        output_txt    = FALSE,
                        snr_thresh    = 3,
                        shapes        = NULL,
                        calibrate     = TRUE) {

  available_shapes <- c("circle", "square")

  # ---- Shape selection -----------------------------------------------------
  # Default to the full two-group example (both shapes: group A circle + group B
  # square). Pass `shapes = "circle"` or `"square"` to restrict to one group.
  if (is.null(shapes)) shapes <- available_shapes
  shapes <- intersect(shapes, available_shapes)
  if (!length(shapes))
    stop("`shapes` must contain at least one of: ",
         paste(available_shapes, collapse = ", "))
  message("Including shape(s): ", paste(shapes, collapse = ", "))

  tmp_dir  <- tempfile("quantMSImageR_example_")
  dir.create(tmp_dir, recursive = TRUE)

  extdata <- system.file("extdata", package = "quantMSImageR",
                         mustWork = TRUE)
  rmd     <- system.file("quantMSImageR___general_heatmap.Rmd",
                         package = "quantMSImageR", mustWork = TRUE)

  lib_ion_path <- file.path(extdata, "example_ion_library.csv")
  data_path    <- extdata          # example.raw/ lives inside extdata/
  image_dir    <- file.path(tmp_dir, "images", "Example")

  # ---- Section catalogue, filtered by selected shape(s) --------------------
  all_sections <- list(
    list(section = "section01", label = "SampleA_1", shape = "circle"),
    list(section = "section02", label = "SampleA_2", shape = "circle"),
    list(section = "section03", label = "SampleA_3", shape = "circle"),
    list(section = "section04", label = "SampleB_1", shape = "square"),
    list(section = "section05", label = "SampleB_2", shape = "square"),
    list(section = "section06", label = "SampleB_3", shape = "square")
  )
  selected <- Filter(function(s) s$shape %in% shapes, all_sections)

  fns <- lapply(selected, function(s)
    list(neg = "example", section = s$section, label = s$label))

  heatmap_order  <- vapply(selected, `[[`, character(1), "label")
  heatmap_labs   <- sub("_[0-9]+$", "", heatmap_order)
  baseline_label <- heatmap_labs[1]   # whichever shape comes first

  sample_map <- data.frame(
    run_id    = heatmap_order,
    group     = heatmap_labs,
    neg_files = "example",
    section   = vapply(selected, `[[`, character(1), "section"),
    shape     = vapply(selected, `[[`, character(1), "shape"),
    stringsAsFactors = FALSE
  )

  # ---- Demo: tissue selection ---------------------------------------------
  message("\n--- Tissue pixel selection demo ---")
  message("In a real study you would run select_tissue_pixels() before the")
  message("YAML pipeline to interactively draw a tissue ROI and save")
  message("tissue_pixels.csv.  Here the mask is already embedded in the RDS.")
  message("Showing the tissue vs background mask for the first section...")

  demo_section <- selected[[1]]$section
  demo_rds <- readRDS(file.path(extdata, "example.raw",
                                  paste0(demo_section, ".RDS")))
  .pd <- as.data.frame(pData(demo_rds))
  print(
    ggplot2::ggplot(.pd, ggplot2::aes(x, -y, fill = sample_name)) +
      ggplot2::geom_tile() +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = sprintf("Tissue vs background mask: %s", demo_section),
                    fill = NULL) +
      ggplot2::theme_minimal()
  )

  message("------------------------------------------------------------------\n")

  # Demonstrate a per-analyte SNR override: give the internal standard(s) half
  # the general threshold (internal standards are abundant, so a lower SNR cut
  # is appropriate). In a real study these come from parameters$snr_overrides in
  # the YAML; here we derive them from the ion library's Type == "IS" rows.
  .ion_lib      <- read.csv(lib_ion_path, check.names = FALSE)
  .is_names     <- .ion_lib$transition_id[.ion_lib$Type == "IS"]
  snr_overrides <- if (length(.is_names))
    stats::setNames(rep(snr_thresh / 2, length(.is_names)), .is_names) else NULL

  result <- generate_txt_images(
    fns            = fns,
    data_path      = data_path,
    image_dir      = image_dir,
    lib_ion_path   = lib_ion_path,
    snr_thresh     = snr_thresh,
    snr_overrides  = snr_overrides,
    thresh         = 20,
    perc           = 97,
    rot_clockwise  = 0,
    average_method = "median",
    output_txt     = output_txt
  )

  # ---- Quantification demo (calibration) -----------------------------------
  # Runs the calibration chain on the bundled synthetic standards data. In a
  # real study this is driven by the YAML `calibration:` block (which loads the
  # standards acquisition via read_mrm); the synthetic standards ship as an RDS,
  # so here we load them directly. Runs BEFORE the report render so the report
  # can pick `calibrated` up from this frame and add its calibration section.
  calibrated <- NULL
  if (calibrate) {
    cal_dir <- file.path(extdata, "cal_example.raw")
    cal_rds <- file.path(cal_dir, "cal_MSI.RDS")
    if (file.exists(cal_rds)) {
      message("\n--- Quantification demo (calibration) ---")
      cal_obj  <- as(readRDS(cal_rds), "quant_MSImagingExperiment")
      cal_meta <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))

      cal_obj <- summarise_cal_levels(cal_obj, cal_meta, val_slot = "intensity",
                                      cal_label = "Cal", id = "identifier")
      cal_obj <- create_cal_curve(cal_obj, cal_type = "cal")

      # Apply the curves to the STUDY sections, exactly as run_study.R does:
      # carry the whole calibrationInfo across, then convert the tissue pixels.
      # Features without a standard are dropped by int2conc(), so the calibrated
      # object holds only the three analytes that were calibrated.
      cal_study <- result$combined
      calibrationData(cal_study) <- calibrationData(cal_obj)
      cal_study <- int2conc(cal_study, val_slot = "intensity",
                            pixel_header = "sample_name",
                            pixels = "tissue_pixels")

      cal_out <- file.path(tmp_dir, "Example_calibrated.RDS")
      saveRDS(cal_study, cal_out)

      .r2 <- calibrationR2(cal_obj)
      message("Per-lipid calibration curve R^2:")
      message(paste(sprintf("  %-14s R2 = %.3f", .r2$feature, .r2$r2),
                    collapse = "\n"))
      message("Calibrated object (pg_pixel, pg_mm2) saved to: ", cal_out)

      calibrated        <- cal_study
      result$calibrated <- cal_study
    } else {
      message("Calibration data not found in extdata; skipping quantification demo.")
    }
  }

  if (render_report) {
    if (!nzchar(rmd) || !file.exists(rmd))
      stop("HTML report template not found -- reinstall quantMSImageR.")

    # Variables the Rmd looks up in its parent env
    combined          <- result$combined_snr
    # Per-feature SNR threshold used, for the SNR_used column in report 1.2: the
    # general threshold, with the IS override (half) applied. Mirrors run_study.R.
    snr_used          <- stats::setNames(rep(snr_thresh, nrow(fData(combined))),
                                          fData(combined)$name)
    if (!is.null(snr_overrides)) {
      .ov <- intersect(names(snr_overrides), names(snr_used))
      snr_used[.ov] <- snr_overrides[.ov]
    }
    ion_lib_meta      <- .ion_lib
    feature_meta      <- build_feature_meta(combined, ion_lib_meta)
    heatmap_row_split <- "Met-1"
    # Mirrors the YAML `colours:` block, so the example report is coloured the
    # same way a real study would be.
    pal_ion           <- "heatmap0"
    pal_heatmap       <- "heatmap2"
    pal_group         <- "hat"
    pal_feature       <- "reading"
    hm_cell_border    <- "white"
    out_path          <- tmp_dir
    report_fn         <- "Example"
    ratios_cfg        <- NULL

    html_file <- file.path(tmp_dir,
                            paste0("Example_SNR", snr_thresh,
                                    "_response_SNRfiltered.html"))

    rmarkdown::render(
      rmd,
      output_file       = html_file,
      intermediates_dir = tmp_dir,
      knit_root_dir     = tmp_dir,
      quiet             = TRUE
    )

    # Returned so the self-contained HTML can be copied somewhere permanent --
    # it is written to a session tempdir and would otherwise be hard to find.
    result$report <- html_file
    message("Report written to: ", html_file)

    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      rstudioapi::viewer(html_file)
    } else {
      utils::browseURL(html_file)
    }
  }

  invisible(result)
}
