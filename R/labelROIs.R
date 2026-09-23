#' Label regions of interest within the tissue of an acquisition
#'
#' Run after [selectTissuePixels()]. Shows the ion image of a chosen transition
#' across the acquisition's **tissue pixels only**, lets you draw a region --
#' an airway, a vessel, parenchyma -- and asks what it is. Repeat for as many
#' regions as you need; the study report then compares them (see the `roi:`
#' block of the config template).
#'
#' ```
#' #> Draw a region on the image, then close the window.
#' #> Label for this region (used in this study: airway, vessel): airway
#' #> Saved as airway_03.
#' #> Add another region? [y/n]: y
#' ```
#'
#' @section Labels and numbers:
#'
#' Several regions can share a label: a section has many airways. Each is
#' numbered as it is saved (`airway_01`, `airway_02`, ...), so the report can
#' treat them as independent regions or pool them per sample. The prompt lists
#' the labels already used anywhere in the study (every `roi_labels.csv` beside
#' this one in `data_path`), so the same thing is called the same name in every
#' acquisition. Type the label only; the number is added for you.
#'
#' A pixel belongs to at most one region. A region drawn over part of an
#' earlier one takes those pixels, and an earlier region left with no pixels is
#' dropped. Tissue outside every region is `unassigned`, which the report can
#' keep or leave out.
#'
#' @section The file:
#'
#' Regions are saved as `roi_labels.csv` in the acquisition's `.raw` folder,
#' one per acquisition, beside `tissue_pixels.csv` -- which is not touched. It
#' lists the labelled pixels by `x`, `y`, `roi_label` and `roi_id`, and is
#' written after every region, so nothing is lost if the session ends early.
#' Running `labelROIs()` again adds regions to the file; `overwrite = TRUE`
#' starts again, keeping the old file as a dated backup.
#'
#' Panels acquired over the same area share their pixel grid, so the regions
#' drawn on one apply to the others: `copy_to` writes the same file into their
#' folders. Regions are matched by pixel position, as the tissue mask is.
#'
#' @import Cardinal
#'
#' @param name Character. Acquisition name, without the `.raw` suffix.
#' @param data_path Character. Path to the folder that contains the `.raw`
#'   acquisition directory.
#' @param lib_ion_path Character. Full path to the MRM ion-library CSV.
#' @param feature Integer index or character name of the transition whose ion
#'   image to draw on. When `NULL` (default), all features are displayed first
#'   and you are prompted to choose one.
#' @param enhance Character. Contrast enhancement passed to `Cardinal::image()`.
#'   Default `"histogram"`.
#' @param overwrite Logical. `TRUE` starts a new `roi_labels.csv`, keeping the
#'   existing one as `roi_labels_<date-time>.csv.bak`. `FALSE` (default) adds to
#'   it.
#' @param copy_to Character. Other acquisitions over the same area (the other
#'   panels of this sample) to copy the finished regions to. An acquisition that
#'   already has regions is skipped unless `overwrite = TRUE`. Default `NULL`.
#' @param preview Logical. Draw the saved regions over the tissue when done.
#'   Default `TRUE`.
#'
#' @return Invisibly, the saved regions: a data frame with columns `x`, `y`,
#'   `roi_label` and `roi_id`, one row per labelled pixel.
#'
#' @examples
#' # Requires a graphics device and user input, so it only runs interactively.
#' if (interactive()) {
#'   # Once per acquisition, after its tissue mask:
#'   selectTissuePixels("my_acquisition", data_path = "path/to/raw",
#'                      lib_ion_path = "path/to/ion_library.csv")
#'   labelROIs("my_acquisition", data_path = "path/to/raw",
#'             lib_ion_path = "path/to/ion_library.csv",
#'             copy_to = c("my_acquisition_pos01", "my_acquisition_pos02"))
#' }
#'
#' @seealso [selectTissuePixels()], [runStudy()]
#' @family acquisition
#' @export
labelROIs <- function(name,
                      data_path,
                      lib_ion_path,
                      feature   = NULL,
                      enhance   = "histogram",
                      overwrite = FALSE,
                      copy_to   = NULL,
                      preview   = TRUE) {

  if (!interactive())
    stop("labelROIs() draws regions by hand, so it needs an interactive R ",
         "session with a graphics device.", call. = FALSE)

  raw_dir  <- file.path(data_path, paste0(name, ".raw"))
  mask_csv <- file.path(raw_dir, "tissue_pixels.csv")
  roi_csv  <- file.path(raw_dir, "roi_labels.csv")

  if (!file.exists(mask_csv))
    stop("No tissue_pixels.csv for '", name, "'. Regions are drawn within the ",
         "tissue, so draw the tissue mask first:\n",
         "  selectTissuePixels(name = \"", name, "\", data_path = \"",
         data_path, "\", lib_ion_path = ...)", call. = FALSE)

  message("Loading acquisition: ", name)
  obj  <- readMRM(name = name, folder = data_path, lib_ion_path = lib_ion_path)
  tiss <- obj[, .mask_tissue(obj, mask_csv, name)]
  if (!ncol(tiss))
    stop("tissue_pixels.csv for '", name, "' marks no pixel as tissue.",
         call. = FALSE)
  key_t <- paste(pData(tiss)$x, pData(tiss)$y, sep = "_")

  # ---- Regions so far ---------------------------------------------------------
  rois <- data.frame(x = integer(), y = integer(), roi_label = character(),
                     roi_id = character(), stringsAsFactors = FALSE)
  if (file.exists(roi_csv)) {
    if (isTRUE(overwrite)) {
      bak <- file.path(raw_dir, sprintf("roi_labels_%s.csv.bak",
                                        format(Sys.time(), "%Y%m%d-%H%M%S")))
      file.copy(roi_csv, bak)
      message("Starting again; the previous regions are kept in ", basename(bak))
    } else {
      rois <- .read_roi_csv(roi_csv)
      message(sprintf("Adding to %d existing region(s): %s", length(unique(rois$roi_id)),
                      paste(unique(rois$roi_id), collapse = ", ")))
    }
  }

  study_labels <- .study_roi_labels(data_path, exclude = roi_csv)

  feat_idx <- .choose_feature(obj, feature, enhance,
                              "shows the regions you want to label")
  message("\nUsing feature [", feat_idx, "]: ", fData(obj)$name[feat_idx],
          " (tissue pixels only)")

  # ---- Draw, label, repeat ------------------------------------------------------
  repeat {
    message("\nDraw a region on the image, then close the window.")
    dev.new()
    sel <- as.logical(selectROI(tiss, enhance = enhance, i = feat_idx))

    if (!any(sel, na.rm = TRUE)) {
      message("No pixels selected.")
    } else {
      label <- .ask_roi_label(unique(c(study_labels, rois$roi_label)))
      used  <- suppressWarnings(as.integer(
        sub(paste0("^", label, "_(\\d+)$"), "\\1",
            unique(rois$roi_id[rois$roi_label == label]))))
      id    <- sprintf("%s_%02d", label, max(c(0L, used), na.rm = TRUE) + 1L)

      new_k  <- key_t[which(sel)]
      old_k  <- paste(rois$x, rois$y, sep = "_")
      taken  <- old_k %in% new_k
      if (any(taken)) {
        from <- table(rois$roi_id[taken])
        message(sprintf("  Takes %d pixel(s) from %s.", sum(taken),
                        paste(sprintf("%s (%d)", names(from), from), collapse = ", ")))
        gone <- setdiff(unique(rois$roi_id[taken]), rois$roi_id[!taken])
        if (length(gone))
          message("  ", paste(gone, collapse = ", "),
                  " had no pixels left and is removed.")
        rois <- rois[!taken, , drop = FALSE]
      }
      rois <- rbind(rois, data.frame(
        x = pData(tiss)$x[which(sel)], y = pData(tiss)$y[which(sel)],
        roi_label = label, roi_id = id, stringsAsFactors = FALSE))

      write.csv(rois, roi_csv, row.names = FALSE)
      message(sprintf("Saved as %s (%d pixels).", id, sum(sel)))
    }

    again <- tolower(trimws(readline("Add another region? [y/n]: ")))
    if (!again %in% c("y", "yes")) break
  }

  message(sprintf("\n%s: %d region(s), %d of %d tissue pixels labelled -> %s",
                  name, length(unique(rois$roi_id)), nrow(rois), ncol(tiss),
                  roi_csv))

  # ---- Same regions for the other panels over this area --------------------------
  for (target in as.character(copy_to)) {
    t_dir <- file.path(data_path, paste0(target, ".raw"))
    t_csv <- file.path(t_dir, "roi_labels.csv")
    if (!dir.exists(t_dir)) {
      message("copy_to: no acquisition folder ", t_dir, " -- skipped.")
    } else if (file.exists(t_csv) && !isTRUE(overwrite)) {
      message("copy_to: '", target, "' already has roi_labels.csv -- skipped ",
              "(overwrite = TRUE to replace it).")
    } else {
      if (file.exists(t_csv))
        file.copy(t_csv, file.path(t_dir, sprintf("roi_labels_%s.csv.bak",
                                   format(Sys.time(), "%Y%m%d-%H%M%S"))))
      write.csv(rois, t_csv, row.names = FALSE)
      message("copy_to: regions written to '", target, "'.")
    }
  }

  if (isTRUE(preview) && nrow(rois)) print(.roi_preview(tiss, rois, name))

  invisible(rois)
}

# Logical tissue flag per pixel of `obj`, from its tissue_pixels.csv: matched on
# (x, y) when the mask has coordinates, by row order for a legacy mask of the
# same length. The same rules generateTxtImages() applies.
.mask_tissue <- function(obj, mask_csv, name) {
  tpdf <- read.csv(mask_csv)
  if (all(c("x", "y") %in% names(tpdf))) {
    m <- match(paste(pData(obj)$x, pData(obj)$y, sep = "_"),
               paste(tpdf$x, tpdf$y, sep = "_"))
    out <- rep(FALSE, ncol(obj))
    out[!is.na(m)] <- as.logical(tpdf$tissue_pixels[m[!is.na(m)]])
    return(out & !is.na(out))
  }
  if (nrow(tpdf) == ncol(obj)) return(as.logical(tpdf$tissue_pixels))
  stop(sprintf(paste0(
    "tissue_pixels.csv for '%s' has %d rows but the acquisition has %d pixels, ",
    "and the CSV has no x/y columns to match on. Redraw it with ",
    "selectTissuePixels(overwrite = TRUE)."), name, nrow(tpdf), ncol(obj)),
    call. = FALSE)
}

.read_roi_csv <- function(path) {
  r <- read.csv(path, stringsAsFactors = FALSE)
  need <- c("x", "y", "roi_label", "roi_id")
  if (!all(need %in% names(r)))
    stop(path, " needs columns ", paste(need, collapse = ", "), ".",
         call. = FALSE)
  r[, need]
}

# Labels in use anywhere in the study: every roi_labels.csv beside this one.
.study_roi_labels <- function(data_path, exclude = NULL) {
  files <- Sys.glob(file.path(data_path, "*.raw", "roi_labels.csv"))
  files <- setdiff(normalizePath(files, winslash = "/", mustWork = FALSE),
                   normalizePath(exclude, winslash = "/", mustWork = FALSE))
  labs <- unlist(lapply(files, function(f)
    tryCatch(unique(read.csv(f, stringsAsFactors = FALSE)$roi_label),
             error = function(e) NULL)))
  sort(unique(as.character(labs)))
}

# Ask until the answer is a usable label. Spaces become underscores; the
# region's number is added afterwards, so a typed number would be doubled.
.ask_roi_label <- function(known) {
  repeat {
    shown <- if (length(known)) paste(known, collapse = ", ") else "none yet"
    lab <- trimws(readline(sprintf(
      "Label for this region (used in this study: %s): ", shown)))
    lab <- gsub("\\s+", "_", lab)
    if (!nzchar(lab)) {
      message("  A label is needed.")
    } else if (!grepl("^[A-Za-z][A-Za-z0-9_-]*$", lab)) {
      message("  Use letters, digits, '_' or '-', starting with a letter.")
    } else if (tolower(lab) == "unassigned") {
      message("  'unassigned' is reserved for tissue outside every region.")
    } else if (grepl("_\\d+$", lab)) {
      message("  Leave the number off -- it is added for you (airway, not airway_03).")
    } else {
      # Match an existing label regardless of case, so Airway and airway do
      # not become two regions types.
      hit <- known[tolower(known) == tolower(lab)]
      return(if (length(hit)) hit[1] else lab)
    }
  }
}

.roi_preview <- function(tiss, rois, name) {
  df <- data.frame(x = pData(tiss)$x, y = pData(tiss)$y)
  m  <- match(paste(df$x, df$y, sep = "_"), paste(rois$x, rois$y, sep = "_"))
  df$roi_label <- ifelse(is.na(m), "unassigned", rois$roi_label[m])
  df$roi_id    <- ifelse(is.na(m), NA_character_, rois$roi_id[m])
  labs <- setdiff(sort(unique(df$roi_label)), "unassigned")
  # The palette the report's region maps use.
  cols <- c(.anno_cols(labs, "Set 2"), unassigned = "grey75")
  cen  <- stats::aggregate(cbind(x, y) ~ roi_id, data = df[!is.na(df$roi_id), ],
                           FUN = stats::median)
  ggplot2::ggplot(df, ggplot2::aes(x = x, y = -y, fill = roi_label)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(data = cen, ggplot2::aes(x = x, y = -y, label = roi_id),
                       inherit.aes = FALSE, size = 3) +
    ggplot2::coord_equal() +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = NULL, y = NULL, fill = NULL,
                  title = paste0(name, ": saved regions")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
}
