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
#'   Default `1.5` (relaxed for synthetic data).
#' @param shapes Character vector of tissue shapes to include -- any of
#'   `"circle"`, `"square"`. `NULL` (the default) prompts interactively when
#'   `interactive()` is `TRUE`, otherwise both shapes are used.
#'
#' @return Invisibly returns the list produced by [generate_txt_images()].
#'
#' @examples
#' \dontrun{
#'   run_example()                       # interactive shape picker
#'   run_example(shapes = "circle")      # SampleA only
#'   run_example(shapes = c("circle","square"))   # both, no prompt
#' }
#'
#' @export
run_example <- function(render_report = TRUE,
                        output_txt    = FALSE,
                        snr_thresh    = 1.5,
                        shapes        = NULL) {

  available_shapes <- c("circle", "square")

  # ---- Shape selection (interactive default) -------------------------------
  if (is.null(shapes)) {
    if (interactive()) {
      message("\nSelect tissue shape(s) to include:")
      idx <- utils::menu(c("Circle  (SampleA, sections 01-03)",
                            "Square  (SampleB, sections 04-06)",
                            "Both    (full example)"),
                          title = "quantMSImageR example -- shape filter")
      shapes <- switch(as.character(idx),
                        "1" = "circle",
                        "2" = "square",
                        "3" = available_shapes,
                        available_shapes)  # 0 = quit -> default to both
    } else {
      shapes <- available_shapes
    }
  }
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
  message("Displaying the first selected section / feature 1 as an example...")

  demo_section <- selected[[1]]$section
  demo_rds <- readRDS(file.path(extdata, "example.raw",
                                  paste0(demo_section, ".RDS")))
  dev.new()
  print(image(demo_rds, enhance = "histogram", i = 1L,
              main = sprintf("Example: %s / feature 1 (%s)\n(select_tissue_pixels() opens this window interactively)",
                              demo_section, selected[[1]]$shape)))

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
    tiss_fc        = 0.6,
    thresh         = 20,
    perc           = 97,
    rot_clockwise  = 0,
    average_method = "median",
    output_txt     = output_txt
  )

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

    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      rstudioapi::viewer(html_file)
    } else {
      utils::browseURL(html_file)
    }
  }

  invisible(result)
}
