#' Run the quantMSImageR example study
#'
#' Executes the full DESI-MRM pipeline on a small synthetic dataset bundled
#' with the package.  The dataset contains two samples (`SampleA`, `SampleB`)
#' each represented by three 20×20-pixel sections, and eight oxylipin features
#' (seven analytes + one internal standard, all negative-ion mode).
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
#'
#' @return Invisibly returns the list produced by [generate_txt_images()].
#'
#' @examples
#' \dontrun{
#'   run_example()
#' }
#'
#' @export
run_example <- function(render_report = TRUE,
                        output_txt    = FALSE,
                        snr_thresh    = 1.5) {

  tmp_dir  <- tempfile("quantMSImageR_example_")
  dir.create(tmp_dir, recursive = TRUE)

  extdata <- system.file("extdata", package = "quantMSImageR",
                         mustWork = TRUE)
  rmd     <- system.file("quantMSImageR___general_heatmap.Rmd",
                         package = "quantMSImageR", mustWork = TRUE)

  lib_ion_path <- file.path(extdata, "example_ion_library.csv")
  data_path    <- extdata          # example.raw/ lives inside extdata/
  image_dir    <- file.path(tmp_dir, "images", "Example")

  heatmap_order  <- c("SampleA_1", "SampleA_2", "SampleA_3",
                      "SampleB_1", "SampleB_2", "SampleB_3")
  heatmap_labs   <- c("SampleA", "SampleA", "SampleA",
                      "SampleB", "SampleB", "SampleB")
  baseline_label <- "SampleA"

  fns <- list(
    list(neg = "example", section = "section01", label = "SampleA_1"),
    list(neg = "example", section = "section02", label = "SampleA_2"),
    list(neg = "example", section = "section03", label = "SampleA_3"),
    list(neg = "example", section = "section04", label = "SampleB_1"),
    list(neg = "example", section = "section05", label = "SampleB_2"),
    list(neg = "example", section = "section06", label = "SampleB_3")
  )

  sample_map <- data.frame(
    run_id    = heatmap_order,
    group     = heatmap_labs,
    neg_files = "example",
    section   = c("section01", "section02", "section03",
                  "section04", "section05", "section06"),
    stringsAsFactors = FALSE
  )

  # ---- Demo: show what tissue selection looks like -------------------------
  # Load one section RDS and display feature 1 (12-HHTrE) so users can see
  # the tissue-vs-background contrast that select_tissue_pixels() uses.
  message("\n--- Tissue pixel selection demo ---")
  message("In a real study you would run select_tissue_pixels() before the")
  message("YAML pipeline to interactively draw a tissue ROI and save")
  message("tissue_pixels.csv.  Here the mask is already embedded in the RDS.")
  message("Displaying section01 / feature '12-HHTrE' as an example...")

  demo_rds <- readRDS(file.path(extdata, "example.raw", "section01.RDS"))
  dev.new()
  print(image(demo_rds, enhance = "histogram", i = 1L,
              main = "Example: 12-HHTrE — use this to draw tissue ROI\n(select_tissue_pixels() opens this window interactively)"))

  message("------------------------------------------------------------------\n")

  result <- generate_txt_images(
    fns            = fns,
    data_path      = data_path,
    image_dir      = image_dir,
    lib_ion_path   = lib_ion_path,
    snr_thresh     = snr_thresh,
    tiss_fc        = 0.6,
    thresh         = 20,
    perc           = 97,
    rot_clockwise  = 0,
    average_method = "median",
    output_txt     = output_txt
  )

  if (render_report) {
    if (!nzchar(rmd) || !file.exists(rmd))
      stop("HTML report template not found — reinstall quantMSImageR.")

    combined  <- result$combined_snr
    html_file <- file.path(tmp_dir, "Example_SNR1.5_response_SNRfiltered.html")

    rmarkdown::render(
      rmd,
      output_file       = html_file,
      intermediates_dir = tmp_dir,
      knit_root_dir     = tmp_dir,
      quiet             = TRUE
    )

    # Open in RStudio viewer pane if available, otherwise default browser
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable()) {
      rstudioapi::viewer(html_file)
    } else {
      utils::browseURL(html_file)
    }
  }

  invisible(result)
}
