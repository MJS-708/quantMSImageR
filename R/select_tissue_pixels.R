#' Interactively select tissue pixels for an acquisition
#'
#' Loads a DESI-MRM acquisition, displays ion images of all features so the
#' user can identify which best discriminates tissue from background, then
#' opens an interactive ROI-selection window on a chosen feature.  The
#' resulting logical mask is saved as `tissue_pixels.csv` inside the
#' acquisition's `.raw` folder, where [generate_txt_images()] will find it
#' automatically.
#'
#' This function must be run **before** the YAML pipeline for each acquisition
#' that does not yet have a `tissue_pixels.csv`.  It requires an interactive R
#' session with a graphics device (i.e. not inside `Rscript --no-save`).
#'
#' @import Cardinal
#'
#' @param name Character. Acquisition name, without the `.raw` suffix
#'   (e.g. `"14Dec_AntibodyStudy_1"`).
#' @param data_path Character. Path to the folder that contains the `.raw`
#'   acquisition directory.
#' @param lib_ion_path Character. Full path to the MRM ion-library CSV.
#' @param feature Integer index or character name of the feature to use for
#'   tissue selection.  When `NULL` (default), all features are displayed
#'   first and the user is prompted to choose one.
#' @param enhance Character. Contrast enhancement passed to `Cardinal::image()`.
#'   Default `"histogram"`.
#' @param overwrite Logical. If `TRUE`, an existing `tissue_pixels.csv` will
#'   be overwritten.  Default `FALSE`.
#' @param preview Logical. Draw the saved mask (tissue versus background over
#'   the pixel coordinates) once selection finishes, so it can be checked before
#'   the pipeline uses it. Default `TRUE`.
#'
#' @return Invisibly returns a data frame with columns `x`, `y`,
#'   `tissue_pixels` and `background_pixels` (logical), one row per pixel, in
#'   acquisition order.
#'
#' @examples
#' # Requires a graphics device and user input, so it only runs interactively.
#' if (interactive()) {
#'   # Step 1 (once per acquisition, before running the study):
#'   select_tissue_pixels(
#'     name         = "my_acquisition",
#'     data_path    = "path/to/raw",
#'     lib_ion_path = "path/to/ion_library.csv"
#'   )
#'
#'   # Step 2: run the full pipeline from the study config
#'   run_study("path/to/study.yaml")
#' }
#'
#' @family acquisition
#' @export
select_tissue_pixels <- function(name,
                                 data_path,
                                 lib_ion_path,
                                 feature   = NULL,
                                 enhance   = "histogram",
                                 overwrite = FALSE,
                                 preview   = TRUE) {

  out_csv <- file.path(data_path, paste0(name, ".raw"), "tissue_pixels.csv")

  if (file.exists(out_csv) && !overwrite)
    stop("tissue_pixels.csv already exists for '", name, "'.\n",
         "  Path: ", out_csv, "\n",
         "  Use overwrite = TRUE to replace it.")

  message("Loading acquisition: ", name)
  obj <- read_mrm(name = name, folder = data_path, lib_ion_path = lib_ion_path)

  n_feat     <- nrow(fData(obj))
  feat_names <- fData(obj)$name

  # ---- Choose feature -------------------------------------------------------
  if (is.null(feature)) {
    message("\nDisplaying all ", n_feat, " features -- identify which best",
            " separates tissue from background, then close the window.")
    dev.new()
    print(image(obj, enhance = enhance, i = seq_len(n_feat), free = "xy"))

    message("\nAvailable features:")
    for (i in seq_len(n_feat))
      message(sprintf("  [%2d]  %s", i, feat_names[i]))

    raw_input <- readline(
      prompt = "Enter feature index (or name) for tissue selection: "
    )

    # Accept numeric index or partial/full name
    feat_idx <- suppressWarnings(as.integer(raw_input))
    if (is.na(feat_idx)) {
      feat_idx <- grep(raw_input, feat_names, ignore.case = TRUE, fixed = FALSE)
      if (length(feat_idx) == 0)
        stop("No feature matching '", raw_input, "' found.")
      if (length(feat_idx) > 1) {
        message("Multiple matches: ", paste(feat_names[feat_idx], collapse = ", "))
        stop("Be more specific.")
      }
    }

  } else if (is.character(feature)) {
    feat_idx <- which(feat_names == feature)
    if (length(feat_idx) == 0)
      stop("Feature '", feature, "' not found.\n",
           "  Available: ", paste(feat_names, collapse = ", "))
  } else {
    feat_idx <- as.integer(feature)
    if (feat_idx < 1L || feat_idx > n_feat)
      stop("feature index ", feat_idx, " out of range (1-", n_feat, ")")
  }

  message("\nUsing feature [", feat_idx, "]: ", feat_names[feat_idx])
  message("Draw a region of interest around the TISSUE (not background).",
          " Close the window when done.")

  # ---- Interactive ROI selection --------------------------------------------
  dev.new()
  tissue_pixels <- selectROI(obj, enhance = enhance, i = feat_idx)

  # ---- Build and save CSV ---------------------------------------------------
  # Coordinates make the mask portable across MRM panels of the same physical
  # sample -- pixel counts can differ between panels because of cycle timing,
  # but (x, y) positions identify the same tissue region.
  tpdf <- data.frame(
    x                 = Cardinal::pData(obj)$x,
    y                 = Cardinal::pData(obj)$y,
    tissue_pixels     = as.logical(tissue_pixels),
    background_pixels = !as.logical(tissue_pixels)
  )

  write.csv(tpdf, out_csv, row.names = FALSE)

  n_tiss <- sum(tpdf$tissue_pixels)
  n_bg   <- sum(tpdf$background_pixels)

  message("\nSaved: ", out_csv)
  message(sprintf("  Tissue pixels     : %d  (%.1f%%)", n_tiss,
                  100 * n_tiss / nrow(tpdf)))
  message(sprintf("  Background pixels : %d  (%.1f%%)", n_bg,
                  100 * n_bg / nrow(tpdf)))

  # ---- Mask preview ---------------------------------------------------------
  # The ROI window closes as soon as selection finishes, so draw the saved mask
  # to confirm what was actually written before it is used downstream.
  if (isTRUE(preview)) {
    tpdf$label <- ifelse(tpdf$tissue_pixels, "tissue_pixels", "background_pixels")
    print(
      ggplot2::ggplot(tpdf, ggplot2::aes(x = x, y = -y, fill = label)) +
        ggplot2::geom_tile() +
        ggplot2::coord_equal() +
        ggplot2::scale_fill_manual(
          values = c(tissue_pixels = "#1b7837", background_pixels = "grey85")) +
        ggplot2::labs(x = NULL, y = NULL, fill = NULL,
                      title = paste0(name, ": saved tissue mask")) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text = ggplot2::element_blank(),
                       panel.grid = ggplot2::element_blank())
    )
    tpdf$label <- NULL
  }

  invisible(tpdf)
}
