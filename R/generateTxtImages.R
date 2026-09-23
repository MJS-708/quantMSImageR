#' Load, process and optionally export per-feature text-image matrices
#'
#' Reads one or more DESI-MRM acquisitions, applies SNR filtering and
#' tissue/background separation, and writes per-feature text matrices that
#' can be imported by imaging software (e.g. MassLynx QuanOptimise overlay
#' tools). Setting `output_txt = FALSE` runs all processing but skips file
#' writing, which is useful when the caller only needs the processed MSI
#' objects (e.g. for HTML report rendering).
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param fns Character vector of acquisition names (without `.raw` suffix)
#'   for single-polarity studies. For dual-polarity studies pass a list where
#'   each element is either a plain string (single acquisition) or a named list
#'   with fields `pos`, `neg`, `label`, and optionally `prefix_pos`/`prefix_neg`
#'   (see [bindPanels()]). `runStudy.R` builds this list automatically
#'   from the YAML `samples` section.
#'
#'   A tissue acquired in pieces (top and bottom, say) is one element with a
#'   `pieces` field instead of `pos`/`neg`: a list of areas, each with its own
#'   `pos` and/or `neg`. Each piece's panels are bound, and the pieces are then
#'   placed by stage position with [stitchAcquisitions()] and become one run.
#' @param data_path Path to the folder that contains the `.raw` acquisition
#'   directories.
#' @param image_dir Root output directory. One sub-directory per acquisition
#'   is created inside it.
#' @param lib_ion_path Full path to the MRM ion-library CSV. Must contain
#'   columns: `transition_id`, `precursor_mz`, `product_mz`, `collision_eV`,
#'   `cone_V`, `Polarity`, `Type`.
#' @param snr_thresh Numeric scalar **or vector**. One or more minimum
#'   signal-to-noise thresholds; pixels below each threshold are set to `NA`.
#'   When a vector is supplied, one complete output directory tree is created
#'   per threshold value. A value of `0` skips SNR filtering entirely *and*
#'   removes the requirement for a `tissue_pixels.csv` mask -- handy for a
#'   smoke-test pass before ROIs are drawn (default `0`).
#' @param thresh Numeric. Cold-spot percentile passed to `imageR()` as
#'   `threshold` (default `20`).
#' @param perc Numeric. Hot-spot percentile passed to `imageR()` as
#'   `percentile` (default `97`).
#' @param rot_clockwise Integer (0-3). Number of 90 deg clockwise rotations
#'   applied via `pracma::rot90()` (default `0`).
#' @param average_method Character, `"mean"` or `"median"`. Statistic used to
#'   summarise the noise pixel vector in `int2SNR()`. Applied consistently to
#'   every acquisition (default `"median"`).
#' @param output_txt Logical. When `FALSE`, processing runs but no files are
#'   written to disk (default `TRUE`).
#' @param exclude Character vector of feature names to drop before processing.
#'   Names must match `fData()$name` exactly (default `NULL` -- keep all).
#' @param rename Named character vector or list mapping old feature names to new
#'   display names, e.g. `c("old name" = "new name")`. Applied to all combined
#'   objects after loading; the new names appear in file names and the heatmap
#'   (default `NULL` -- no renaming).
#' @param ratios Optional list of metabolite ratio pairs. Each element must be
#'   a named list with fields:
#'   \describe{
#'     \item{`num`}{Name of the numerator feature (matches `fData()$name`
#'       **after** any `rename` overrides are applied).}
#'     \item{`den`}{Name of the denominator feature.}
#'     \item{`label`}{*(optional)* Output filename prefix. Defaults to
#'       `"<num>_over_<den>"`.}
#'   }
#'   Ratio images (numerator / denominator, pixel-wise) are written to the
#'   same SNR-filtered subdirectories as the individual ion images, immediately
#'   after them. Pixels where the denominator is zero or either ion is `NA` are
#'   set to `NA` (default `NULL` -- no ratios).
#' @param output_ratios Logical. When `FALSE`, ratio txt files are not written
#'   even if `ratios` pairs are defined. Has no effect when `output_txt` is
#'   `FALSE` (default `TRUE`).
#' @param snr_overrides Optional named map of feature name to a feature-specific
#'   SNR threshold. Keys may be original ion-library names **or** post-`rename`
#'   display names. Listed features use their own threshold in every report;
#'   all other features use the global `snr_thresh`. Ignored for any report
#'   whose global threshold is `0` (the mask-free smoke-test pass). Default
#'   `NULL` (all features use `snr_thresh`).
#' @param is_name Character. Value identifying the internal standard in
#'   `fData()$feature_type` -- that is, the ion library's `Type` column, typically
#'   `"IS"`. It is **not** a transition name.
#'   When supplied, each acquisition is normalised with [int2response()] before
#'   SNR filtering, and every downstream step -- SNR, background masking, images
#'   and ratios -- works on the resulting `response` layer instead of
#'   `intensity`. `NULL` or `"None"` (the default) skips normalisation.
#' @param is_norm_header Character. Ion-library column mapping each analyte to
#'   the standard that normalises it, needed only when `is_name` matches more
#'   than one feature (default `"IS_norm"`). See [int2response()].
#' @param is_mode Character. Level at which the internal standard is
#'   summarised: `"line"` (default), `"sample"`, `"pixel"` or `"window"` (a
#'   rolling median over `is_window` consecutive pixels along the acquisition
#'   line). See [int2response()].
#' @param is_window Integer. Number of consecutive pixels averaged when
#'   `is_mode = "window"` (default `15`).
#' @param remove_IS Logical. Drop the internal-standard feature after
#'   normalising (default `TRUE`).
#' @param type_header Character. Ion-library column holding the feature type,
#'   passed to [readMRM()] and matched by `is_name` (default `"Type"`).
#' @param rois Regions of interest drawn with [labelROIs()]. `TRUE` reads
#'   `roi_labels.csv` from each acquisition's `.raw` folder and adds `roi_label`
#'   and `roi_id` to `pData()`: the region's name and its numbered identifier
#'   (`airway_01`), `"unassigned"` for tissue outside every region and `NA` for
#'   background. Acquisitions without the file are entirely `"unassigned"`.
#'   `"auto"` is `TRUE` when any acquisition has the file; `FALSE` (default)
#'   ignores them. Regions already on an object -- a `section:` RDS saved with
#'   `roi_label`/`roi_id`, as the bundled example sections are -- are kept as
#'   they are, but `"auto"` cannot see them before loading, so switch regions on
#'   with `TRUE` for a study built from such sections.
#'
#' @return Invisibly returns a named list:
#'   \describe{
#'     \item{`combined`}{Raw combined object (background pixels retained).}
#'     \item{`combined_snr`}{SNR-filtered intensity for the **first** threshold
#'       in `snr_thresh`; sub-threshold pixels `NA`. Provided for
#'       backward compatibility.}
#'     \item{`combined_snr_list`}{Named list of SNR-filtered objects, one per
#'       threshold in `snr_thresh`, named `"snr<value>"` (e.g. `"snr3"`).}
#'     \item{`combined_NAbackground`}{Raw intensity with background set to `NA`.}
#'     \item{`rois`}{Logical: whether region-of-interest labels were attached.}
#'   }
#'
#' @seealso [int2SNR()], [applySNR()], [back2NA()], [imageR()]
#'
#' @examples
#' \donttest{
#' # Section-mode run on the bundled synthetic data (no tissue mask needed at
#' # snr_thresh = 0), returning the processed objects without writing files.
#' fns <- list(list(neg = "example", section = "section01", label = "A"))
#' res <- generateTxtImages(
#'   fns          = fns,
#'   data_path    = system.file("extdata", package = "quantMSImageR"),
#'   image_dir    = tempfile(),
#'   lib_ion_path = system.file("extdata", "example_ion_library.csv",
#'                              package = "quantMSImageR"),
#'   snr_thresh   = 0,
#'   output_txt   = FALSE
#' )
#' }
#'
#' @family workflow
#' @export
generateTxtImages <- function(
  fns,
  data_path,
  image_dir,
  lib_ion_path,
  snr_thresh     = 0,
  thresh         = 20,
  perc           = 97,
  rot_clockwise  = 0,
  average_method = "median",
  output_txt     = TRUE,
  exclude        = NULL,
  rename         = NULL,
  ratios         = NULL,
  output_ratios  = TRUE,
  snr_overrides  = NULL,
  is_name        = NULL,
  is_norm_header = "IS_norm",
  is_mode        = "line",
  is_window      = 15,
  remove_IS      = TRUE,
  type_header    = "Type",
  rois           = FALSE
) {

  average_method <- match.arg(average_method, c("mean", "median"))

  # Internal-standard normalisation. int2response() adds a `response` layer
  # rather than overwriting `intensity`, so when it runs everything downstream
  # -- SNR, background masking, images -- must work on `response` instead.
  .use_is <- !is.null(is_name) && nzchar(is_name) && !identical(is_name, "None")
  .val    <- if (.use_is) "response" else "intensity"
  if (.use_is)
    message("Internal-standard normalisation: dividing by '", is_name,
            "' per ", is_mode, "; downstream steps use the 'response' layer.")
  snr_thresh_vec <- sort(unique(as.numeric(unlist(snr_thresh))))

  # Normalise snr_overrides (a name -> threshold map, possibly a YAML list)
  # into a named numeric vector. Keys may be display (post-rename) names, so
  # we translate them to original ion-library names below, since int2SNR runs
  # before apply_renames().
  if (!is.null(snr_overrides) && length(snr_overrides) > 0) {
    snr_overrides <- vapply(snr_overrides, as.numeric, numeric(1))
    message("Per-analyte SNR overrides in effect: ",
            paste(sprintf("%s=%s", names(snr_overrides), snr_overrides),
                  collapse = ", "),
            " (all other features use the global snr_thresh).")
  } else {
    snr_overrides <- NULL
  }

  # Given the features present on an object (original ion-library names),
  # return a named numeric vector of per-feature thresholds for any feature
  # that has an override. An override key matches either the feature's
  # original name or its post-rename display name (so users can write either).
  resolve_overrides <- function(orig_names) {
    if (is.null(snr_overrides)) return(NULL)
    out <- vapply(orig_names, function(nm) {
      disp <- if (!is.null(rename) && nm %in% names(rename))
                as.character(rename[[nm]]) else nm
      if (disp %in% names(snr_overrides)) return(snr_overrides[[disp]])
      if (nm   %in% names(snr_overrides)) return(snr_overrides[[nm]])
      NA_real_
    }, numeric(1))
    out <- out[!is.na(out)]
    if (length(out) == 0) NULL else out
  }

  # If every requested threshold is 0, the workflow can skip the tissue mask
  # entirely (no SNR filtering, no tissue/background separation). A single
  # tissue_pixels.csv-free smoke-test pass is then possible, useful for
  # checking that acquisitions load cleanly before any ROI is drawn.
  .needs_mask <- any(snr_thresh_vec > 0)
  if (!.needs_mask)
    message("snr_thresh = 0: running without tissue masks (raw intensities).")

  # ----- Internal helpers ------------------------------------------------

  make_txt_mat <- function(MSIobject, feat_ind, val_slot, value_label,
                            threshold, percentile) {
    # When every pixel is NA (analyte absent in this sample), return a
    # correctly-dimensioned zero matrix so downstream tools (spatialData /
    # tissuUmaps) receive a file with the right pixel grid.
    ivals <- as.numeric(spectraData(MSIobject[feat_ind, ])[[val_slot]])
    if (all(is.na(ivals))) {
      nx <- length(unique(coord(MSIobject)$x))
      ny <- length(unique(coord(MSIobject)$y))
      return(pracma::rot90(matrix(0, nrow = ny, ncol = nx), k = rot_clockwise))
    }

    tryCatch({
      result <- imageR(
        MSIobject  = MSIobject, val_slot   = val_slot, value      = value_label,
        scale      = "suppress", threshold  = threshold, sample_lab = "run",
        pixels     = NA, percentile = percentile,
        feat_ind   = feat_ind, blank_back  = FALSE, text_image  = TRUE
      )
      pracma::rot90(as.matrix(result), k = rot_clockwise)
    }, error = function(e) {
      # A warning, not a message: this handler exists for per-feature data
      # problems, but it also swallows programming errors, and an empty matrix
      # is silently skipped by write_if_nonempty() -- so a broken call here
      # looks exactly like a successful run that wrote nothing.
      warning("imageR failed (", value_label, "): ", e$message, call. = FALSE)
      matrix(NA_real_, 0, 0)
    })
  }

  write_if_nonempty <- function(mat, path) {
    if (nrow(mat) > 0)
      write.table(mat, path, sep = "\t", row.names = FALSE, col.names = FALSE)
  }

  normalize_0_100 <- function(mat) {
    if (nrow(mat) == 0) return(mat)
    rng <- range(mat, na.rm = TRUE, finite = TRUE)
    if (!all(is.finite(rng)) || diff(rng) == 0)
      return(matrix(0, nrow(mat), ncol(mat)))
    (mat - rng[1]) / diff(rng) * 100
  }

  # Apply hot-spot suppression and cold-spot removal to a plain matrix
  # (mirrors the imageR scale="suppress" + threshold logic for derived matrices
  # such as ratios that are not stored in an MSI slot).
  suppress_mat <- function(mat, cold_pct, hot_pct) {
    if (nrow(mat) == 0) return(mat)
    vals <- mat
    cap <- quantile(vals, hot_pct / 100, na.rm = TRUE)
    if (is.finite(cap))
      vals[!is.na(vals) & vals > cap] <- cap
    floor_val <- quantile(vals, cold_pct / 100, na.rm = TRUE)
    if (is.finite(floor_val))
      vals[!is.na(vals) & vals < floor_val] <- 0
    vals
  }

  safe_feat_name <- function(x) {
    x <- gsub("\\|\\|", "_or_", x)
    x <- gsub(":",      "_",    x)
    x <- gsub("/",      ".",    x)
    x
  }

  apply_renames <- function(obj, rename_map) {
    if (is.null(rename_map) || length(rename_map) == 0) return(obj)
    old_nms <- names(rename_map)
    for (i in seq_along(old_nms)) {
      ind <- fData(obj)$name == old_nms[i]
      if (any(ind)) fData(obj)$name[ind] <- rename_map[[i]]
    }
    featureNames(obj) <- fData(obj)$name
    obj
  }

  `%||%` <- function(a, b) if (is.null(a)) b else a

  # Load a single acquisition, attach tissue/noise labels, trim, exclude.
  # When `section` is non-NULL, load a pre-saved per-section RDS file from
  # inside the .raw folder (e.g. slide1_brain01.RDS) instead of reading the
  # Waters raw data. The RDS is assumed to already carry pData$sample_name
  # and ion-library-matched fData (saved by an earlier preprocessing step).
  # Regions of interest are attached last, once the tissue labels are known.
  load_and_prep_acq <- function(fn_name, section = NULL) {
    attach_rois(.load_acq(fn_name, section), fn_name, section)
  }

  .load_acq <- function(fn_name, section = NULL) {
    if (!is.null(section) && nzchar(section)) {
      rds_path <- file.path(data_path, paste0(fn_name, ".raw"),
                            paste0(section, ".RDS"))
      if (!file.exists(rds_path))
        stop("Section RDS not found: ", rds_path)
      obj <- readRDS(rds_path)
      if (!"sample_name" %in% names(pData(obj)))
        stop("RDS for section '", section, "' lacks pData$sample_name; ",
             "regenerate RDS with tissue_pixels/background_pixels labels attached.")
      if (!is(obj, "quant_MSImagingExperiment"))
        obj <- as(obj, "quant_MSImagingExperiment")
    } else {
      obj  <- readMRM(name = fn_name, folder = data_path, lib_ion_path = lib_ion_path,
                       type_header = type_header, is_norm_header = is_norm_header)

      # Skip tissue mask entirely when no SNR > 0 is requested. Every pixel
      # is labelled `tissue_pixels` so downstream code that reads sample_name
      # (e.g. report quantiles) treats the whole acquisition as tissue.
      if (!.needs_mask) {
        pData(obj)$sample_name <- factor(
          rep("tissue_pixels", ncol(obj)),
          levels = c("tissue_pixels", "background_pixels")
        )
        obj <- as(obj, "quant_MSImagingExperiment")
        obj <- trimMSI(MSI_data = obj)
        if (!is.null(exclude) && length(exclude) > 0) {
          keep <- !fData(obj)$name %in% exclude
          if (!all(keep)) {
            message(sprintf("  Excluding %d transition(s): %s",
                            sum(!keep),
                            paste(fData(obj)$name[!keep], collapse = ", ")))
            obj <- obj[keep, ]
          }
        }
        return(obj)
      }

      tpdf_path <- sprintf("%s/%s.raw/tissue_pixels.csv", data_path, fn_name)
      tpdf <- read.csv(tpdf_path)

      has_xy <- all(c("x", "y") %in% names(tpdf))

      if (has_xy) {
        # Coordinate-based matching -- portable across MRM panels of the same
        # physical sample. Pixels in the MSI without a CSV match are labelled
        # background (conservative; they are excluded from tissue quantiles).
        key_obj <- paste(pData(obj)$x, pData(obj)$y, sep = "_")
        key_csv <- paste(tpdf$x,        tpdf$y,        sep = "_")
        midx    <- match(key_obj, key_csv)
        is_tiss <- rep(FALSE, ncol(obj))
        ok      <- !is.na(midx)
        is_tiss[ok] <- as.logical(tpdf$tissue_pixels[midx[ok]])

        n_unmatched <- sum(!ok)
        if (n_unmatched > 0)
          message(sprintf(
            "  '%s': %d / %d MSI pixels not present in tissue_pixels.csv -- labelled as background.",
            fn_name, n_unmatched, ncol(obj)))

        pData(obj)$sample_name <- factor(
          ifelse(is_tiss, "tissue_pixels", "background_pixels"),
          levels = c("tissue_pixels", "background_pixels")
        )
      } else if (nrow(tpdf) == ncol(obj)) {
        # Legacy mask -- no coordinates, row-order alignment with current MSI.
        # Masks written before the rename carry a `noise_pixels` column, so
        # accept either spelling.
        bg_col <- if ("background_pixels" %in% names(tpdf))
                    "background_pixels" else "noise_pixels"
        pData(obj)$sample_name <- makeFactor(
          tissue_pixels     = tpdf[["tissue_pixels"]],
          background_pixels = tpdf[[bg_col]]
        )
      } else {
        stop(sprintf(
          paste0("tissue_pixels.csv for '%s' has %d rows but the loaded MSI ",
                 "has %d pixels, and the CSV has no x/y columns to match on.\n",
                 "  CSV: %s\n",
                 "  Regenerate the mask with the current package version (which ",
                 "writes x/y columns) so it remains portable across MRM panels:\n",
                 "    selectTissuePixels(name = \"%s\", data_path = \"%s\", overwrite = TRUE)"),
          fn_name, nrow(tpdf), ncol(obj), tpdf_path, fn_name, data_path
        ))
      }

      obj <- as(obj, "quant_MSImagingExperiment")
    }
    # trimMSI() drops fully-background border rows/columns. When the tissue fills
    # a perfect rectangle every background pixel lies in such a border, so trimming
    # would remove them all and leave int2SNR() with nothing to reference against.
    # Keep the untrimmed object in that case.
    .has_bg  <- function(o) any(as.character(pData(o)$sample_name) %in%
                                  .bg_labels("background_pixels"))
    .pre_trim <- obj
    obj <- trimMSI(MSI_data = obj)
    if (.has_bg(.pre_trim) && !.has_bg(obj)) obj <- .pre_trim
    if (!is.null(exclude) && length(exclude) > 0) {
      keep <- !fData(obj)$name %in% exclude
      if (!all(keep)) {
        message(sprintf("  Excluding %d transition(s): %s",
                        sum(!keep),
                        paste(fData(obj)$name[!keep], collapse = ", ")))
        obj <- obj[keep, ]
      }
    }
    obj
  }

  # Load one or more acquisitions of the same polarity and merge them into a
  # single MSI by coordinate-matched rbind. bindPanels() handles the typical
  # case where two MRM panels of the same tissue sample at slightly different
  # rates and therefore produce different pixel grids -- pixels are matched on
  # (x, y), features are stacked, so each matched pixel ends up with both
  # panels' transitions.
  # `mode` says what several .raw files under one sample mean:
  #   "panels" (default) different transitions over the same pixels
  #   "stitch"           pieces of one tissue, acquired in separate passes
  # They are different operations -- bindPanels() intersects pixels and would
  # keep nothing from two halves that do not overlap -- so the YAML has to say
  # which, rather than the code guessing from the data.
  load_and_prep_multiple <- function(fns_vec, label, mode = "panels") {
    fns_vec <- as.character(unlist(fns_vec))
    objs    <- lapply(fns_vec, load_and_prep_acq)
    if (length(objs) == 1L) {
      obj <- objs[[1]]
    } else if (identical(mode, "stitch")) {
      obj <- stitchAcquisitions(objs, label = label)
    } else {
      obj <- objs[[1]]
      for (i in seq_along(objs)[-1])
        obj <- bind_panels(obj, objs[[i]], label = label)
    }
    pData(obj)$run <- factor(rep(label, ncol(obj)))
    obj
  }

  # One area of tissue: every panel acquired over it, both polarities, bound
  # into one object. A sample is one area, or -- with `pieces:` -- several,
  # stitched together afterwards. `combine` describes the several files WITHIN
  # a polarity; the two polarities are always panels of one another.
  load_area <- function(pos_fns, neg_fns, label, mode = "panels") {
    if (!is.null(pos_fns) && !is.null(neg_fns)) {
      pos_obj <- load_and_prep_multiple(pos_fns, label = label, mode = mode)
      neg_obj <- load_and_prep_multiple(neg_fns, label = label, mode = mode)
      bind_panels(pos_obj, neg_obj, label = label)
    } else {
      load_and_prep_multiple(pos_fns %||% neg_fns, label = label, mode = mode)
    }
  }

  # bindPanels() keeps obj1's pixel metadata, so regions drawn on obj2's
  # acquisition would be lost whenever obj1 was not labelled. Carry them over
  # where obj1 has none; panels of one area share their pixel grid.
  bind_panels <- function(obj1, obj2, label) {
    out <- bindPanels(obj1, obj2, label = label)
    if (!.use_roi) return(out)
    pd2  <- pData(obj2)
    k2   <- paste(pd2$x, pd2$y, sep = "_")
    kout <- paste(pData(out)$x, pData(out)$y, sep = "_")
    m    <- match(kout, k2)
    lab  <- as.character(pData(out)$roi_label)
    id   <- as.character(pData(out)$roi_id)
    lab2 <- as.character(pd2$roi_label)[m]
    id2  <- as.character(pd2$roi_id)[m]
    has2 <- !is.na(lab2) & lab2 != "unassigned"
    fill <- has2 & !is.na(lab) & lab == "unassigned"
    clash <- sum(has2 & !is.na(lab) & lab != "unassigned" & lab != lab2)
    if (clash > 0)
      message(sprintf(paste0(
        "  '%s': %d pixel(s) carry a different region label in each panel; ",
        "the first panel's label is kept."), label, clash))
    if (any(fill)) {
      lab[fill] <- lab2[fill]
      id[fill]  <- id2[fill]
      pData(out)$roi_label <- lab
      pData(out)$roi_id    <- id
    }
    out
  }

  # Regions of interest drawn with labelROIs(), read from roi_labels.csv in the
  # acquisition's .raw folder and matched on (x, y) like the tissue mask.
  # Background pixels belong to no region (NA); tissue outside every region is
  # "unassigned". With regions switched on, every acquisition gets both columns
  # -- unlabelled ones entirely "unassigned" -- so samples still combine.
  attach_rois <- function(obj, fn_name, section = NULL) {
    if (!.use_roi) return(obj)
    pd <- pData(obj)
    # Regions the object already carries -- a section RDS saved with them, as
    # the bundled example sections are -- are what they say they are.
    if (all(c("roi_label", "roi_id") %in% names(pd)) &&
        any(!is.na(pd$roi_label) & pd$roi_label != "unassigned"))
      return(obj)
    is_bg <- as.character(pData(obj)$sample_name) %in%
               .bg_labels("background_pixels")
    lab <- rep("unassigned", ncol(obj))
    id  <- lab
    roi_path <- file.path(data_path, paste0(fn_name, ".raw"), "roi_labels.csv")
    if ((is.null(section) || !nzchar(section)) && file.exists(roi_path)) {
      rl <- read.csv(roi_path, stringsAsFactors = FALSE)
      if (!all(c("x", "y", "roi_label", "roi_id") %in% names(rl)))
        stop("roi_labels.csv for '", fn_name, "' needs columns x, y, ",
             "roi_label and roi_id; redraw it with labelROIs().", call. = FALSE)
      m  <- match(paste(pData(obj)$x, pData(obj)$y, sep = "_"),
                  paste(rl$x, rl$y, sep = "_"))
      ok <- !is.na(m)
      lab[ok] <- as.character(rl$roi_label[m[ok]])
      id[ok]  <- as.character(rl$roi_id[m[ok]])
      n_bg <- sum(ok & is_bg & lab != "unassigned")
      if (n_bg > 0)
        message(sprintf(paste0(
          "  '%s': %d region pixel(s) are background in tissue_pixels.csv and ",
          "belong to no region."), fn_name, n_bg))
      message(sprintf("  '%s': %d region(s) from roi_labels.csv.", fn_name,
                      length(unique(id[ok & !is_bg & lab != "unassigned"]))))
    }
    lab[is_bg] <- NA_character_
    id[is_bg]  <- NA_character_
    pData(obj)$roi_label <- lab
    pData(obj)$roi_id    <- id
    obj
  }


  # ----- Normalise fns -> fn_list / fn_labels ----------------------------
  # fns: character vector (backward-compat) OR list where each element is
  # a string (single acq) or named list with pos/neg/label fields.
  fn_list   <- as.list(fns)
  fn_labels <- vapply(fn_list, function(e) {
    if (!is.list(e)) return(as.character(e))
    first <- if (!is.null(e$pieces)) e$pieces[[1]] else e
    as.character(e$label %||% (unlist(first$pos)[1] %||% unlist(first$neg)[1]))
  }, character(1))

  # Regions of interest: "auto" switches them on when any acquisition in the
  # study has a roi_labels.csv.
  .acq_names <- unique(unlist(lapply(fn_list, function(e) {
    if (!is.list(e)) return(as.character(e))
    areas <- if (!is.null(e$pieces)) e$pieces else list(e)
    unlist(lapply(areas, function(a) c(unlist(a$pos), unlist(a$neg))))
  })))
  .use_roi <- if (identical(tolower(as.character(rois)), "auto")) {
    any(file.exists(file.path(data_path, paste0(.acq_names, ".raw"),
                              "roi_labels.csv")))
  } else isTRUE(as.logical(rois))
  if (.use_roi)
    message("Regions of interest: on (roi_labels.csv read where present).")

  # ----- Load and process acquisitions -----------------------------------

  combined          <- NULL
  combined_snr_list <- vector("list", length(snr_thresh_vec))

  # Internal-standard normalisation, per acquisition and before SNR, since
  # int2SNR() references the background of whichever layer it is given.
  normalise_is <- function(obj, what) {
    if (!.use_is) return(obj)
    # int2response() matches the standard on the feature-type column -- the
    # ion library's Type column, typically "IS" -- not on the transition
    # name. Checked here as well so the message names the acquisition.
    .ft <- as.character(.feature_type(obj))
    if (!is_name %in% .ft)
      stop("generateTxtImages: no feature typed '", is_name,
           "' in '", what, "'. `is_name` is the value of the ion ",
           "library's Type column (e.g. \"IS\"), not a transition name. ",
           "Values present: ",
           paste(unique(.ft), collapse = ", "),
           call. = FALSE)
    int2response(obj, val_slot = "intensity",
                 normalisation = "internal_standard",
                 IS_name = is_name, is_norm_header = is_norm_header,
                 mode = is_mode, window = is_window,
                 remove_IS = remove_IS)
  }

  # SNR filtering at one threshold. A threshold of 0 short-circuits to the raw
  # object (no SNR computation) -- the smoke-test path -- so overrides are
  # ignored for a 0-valued global report.
  snr_filter <- function(obj, thr, ov) {
    if (thr == 0) return(obj)
    tmp <- int2SNR(
      MSIobject = obj, val_slot = .val, pixel_header = "sample_name",
      background = "background_pixels", tissue = "tissue_pixels",
      snr_thresh = thr, average = average_method,
      snr_overrides = ov
    )
    applySNR(MSIobject = tmp, val_slot = .val)
  }

  for (ind in seq_along(fn_list)) {
    fn_entry <- fn_list[[ind]]
    fn_label <- fn_labels[ind]
    message(sprintf("Loading %s (%d/%d)", fn_label, ind, length(fn_list)))

    # `parts` are normalised and SNR-filtered on their own: the sample, or each
    # piece of a sample acquired in pieces. Pieces are separate acquisitions --
    # often on different days -- so each is referenced to its own internal
    # standard and its own background before they are stitched into one.
    if (is.list(fn_entry)) {
      pos_fns <- if (!is.null(fn_entry$pos)) unlist(fn_entry$pos) else NULL
      neg_fns <- if (!is.null(fn_entry$neg)) unlist(fn_entry$neg) else NULL
      section <- if (!is.null(fn_entry$section)) as.character(fn_entry$section) else NULL
      combine_mode <- if (!is.null(fn_entry$combine))
                        as.character(fn_entry$combine) else "panels"

      if (!is.null(section) && nzchar(section)) {
        # RDS-section mode: one section of one acquisition
        if (!is.null(pos_fns) && !is.null(neg_fns))
          stop("`section:` with dual polarity (pos + neg) is not supported")
        raw_name <- (neg_fns %||% pos_fns)[1]
        parts    <- list(load_and_prep_acq(raw_name, section = section))
      } else if (!is.null(fn_entry$pieces)) {
        # Pieces of one tissue acquired separately (top and bottom, say). Each
        # piece is one area, its panels bound as usual; the pieces are placed
        # by stage position below and become one sample.
        parts <- lapply(fn_entry$pieces, function(pc) {
          p_fns <- if (!is.null(pc$pos)) unlist(pc$pos) else NULL
          n_fns <- if (!is.null(pc$neg)) unlist(pc$neg) else NULL
          load_area(p_fns, n_fns, label = as.character((n_fns %||% p_fns)[1]))
        })
      } else {
        parts <- list(load_area(pos_fns, neg_fns, label = fn_label,
                                mode = combine_mode))
      }
    } else {
      # Plain string: backward-compatible single acquisition
      parts <- list(load_and_prep_acq(fn_entry))
    }

    parts <- lapply(parts, normalise_is, what = fn_label)

    # One sample from its parts. Stitching is repeated for every SNR threshold;
    # it places the pieces identically each time, so only the first says so.
    assemble <- function(objs, quiet = FALSE) {
      o <- if (length(objs) == 1L) objs[[1]]
           else if (quiet) suppressMessages(stitchAcquisitions(objs, label = fn_label))
           else stitchAcquisitions(objs, label = fn_label)
      pData(o)$run <- factor(rep(fn_label, ncol(o)))
      o
    }
    tissue <- assemble(parts)

    # Per-feature overrides resolved against THIS acquisition's features.
    .ov <- resolve_overrides(as.character(fData(tissue)$name))

    # One SNR-filtered object per requested threshold.
    for (si in seq_along(snr_thresh_vec)) {
      thr <- snr_thresh_vec[si]
      tissue_snr <- if (thr == 0) tissue
                    else assemble(lapply(parts, snr_filter, thr = thr, ov = .ov),
                                  quiet = TRUE)

      if (is.null(combined_snr_list[[si]])) {
        combined_snr_list[[si]] <- tissue_snr
      } else {
        al <- alignFeatures(combined_snr_list[[si]], tissue_snr)
        combined_snr_list[[si]] <- combineMSIs(al$obj1, al$obj2)
      }
    }

    if (is.null(combined)) {
      combined <- tissue
    } else {
      al       <- alignFeatures(combined, tissue)
      combined <- combineMSIs(al$obj1, al$obj2)
    }
  }

  combined_NAbackground <- if (.needs_mask) back2NA(
    combined, val_slot = .val,
    background = "background_pixels", pixel_header = "sample_name"
  ) else combined

  # Apply display-name overrides to all objects
  if (!is.null(rename) && length(rename) > 0) {
    combined              <- apply_renames(combined,              rename)
    combined_snr_list     <- lapply(combined_snr_list, function(o) apply_renames(o, rename))
    combined_NAbackground <- apply_renames(combined_NAbackground, rename)
  }

  out <- list(
    combined              = combined,
    combined_snr          = combined_snr_list[[1]],
    combined_snr_list     = setNames(combined_snr_list, paste0("snr", snr_thresh_vec)),
    combined_NAbackground = combined_NAbackground,
    rois                  = .use_roi
  )

  if (!output_txt) return(invisible(out))

  # ----- Generate text-image files ---------------------------------------

  for (snr_i in seq_along(snr_thresh_vec)) {
    snr_t          <- snr_thresh_vec[snr_i]
    combined_snr_i <- combined_snr_list[[snr_i]]

    for (fn_label in fn_labels) {

      combined_snr_tmp  <- combined_snr_i[,        pData(combined_snr_i)$run        == fn_label]
      combined_back_tmp <- combined_NAbackground[,  pData(combined_NAbackground)$run == fn_label]

      image_path <- file.path(image_dir, fn_label)
      dirs <- list(
        snr_filt  = file.path(image_path, sprintf("intensity_SNRfiltered%s",                   snr_t)),
        snr_norm  = file.path(image_path, sprintf("response_SNRfiltered%s_NORM",               snr_t)),
        raw       = file.path(image_path, "intensity_raw"),
        hs_filt   = file.path(image_path, sprintf("intensity_SNRfiltered%s_hs%s_cs%s_removal", snr_t, perc, thresh)),
        hs_norm   = file.path(image_path, sprintf("response_SNRfiltered%s_hs%s_cs%s_NORM",     snr_t, perc, thresh)),
        combined  = file.path(image_path, sprintf("response_SNRfiltered%s_hs%s_cs%s_COMBINED", snr_t, perc, thresh))
      )
      for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

      feat_names_all <- fData(combined_snr_i)$name
      n_features     <- length(feat_names_all)

      # --- Ion images ---
      for (feat_ind in seq_len(n_features)) {

        feat_name <- safe_feat_name(feat_names_all[feat_ind])

        # SNR-filtered raw intensities
        mat_snr <- make_txt_mat(combined_snr_tmp, feat_ind, .val,
                                 "DESI-MRM response - S/N filtered", 0, 100)
        write_if_nonempty(mat_snr, file.path(dirs$snr_filt, paste0(feat_name, ".txt")))

        # 0-100 normalised SNR
        scaled_snr <- normalize_0_100(mat_snr)
        write_if_nonempty(scaled_snr, file.path(dirs$snr_norm, paste0(feat_name, ".txt")))

        # COMBINED: embed global max in [1,1] for cross-sample colour scaling
        if (nrow(scaled_snr) > 0) {
          global_max <- max(
            as.numeric(spectraData(combined_snr_i[feat_ind, ])[[.val]]),
            na.rm = TRUE
          )
          scaled_snr_combined       <- scaled_snr
          scaled_snr_combined[1, 1] <- global_max
          write_if_nonempty(scaled_snr_combined,
                            file.path(dirs$combined, paste0(feat_name, ".txt")))
        }

        # raw doesn't vary with SNR threshold -- write only on first pass
        if (snr_i == 1L) {
          mat_raw <- make_txt_mat(combined_back_tmp, feat_ind, .val,
                                   "DESI-MRM response", 0, 100)
          write_if_nonempty(mat_raw, file.path(dirs$raw, paste0(feat_name, ".txt")))
        }

        # SNR-filtered + hot/cold-spot removal
        mat_hs <- make_txt_mat(combined_snr_tmp, feat_ind, .val,
                                "DESI-MRM response - S/N filtered", thresh, perc)
        write_if_nonempty(mat_hs, file.path(dirs$hs_filt, paste0(feat_name, ".txt")))

        # 0-100 normalised hs/cs version
        scaled_hs <- normalize_0_100(mat_hs)
        write_if_nonempty(scaled_hs, file.path(dirs$hs_norm, paste0(feat_name, ".txt")))
      }

      # --- Ratio images (written right after ion images) ---
      if (output_ratios && !is.null(ratios) && length(ratios) > 0) {
        for (rp in ratios) {
          num_name    <- as.character(rp$num)
          den_name    <- as.character(rp$den)
          ratio_label <- if (!is.null(rp$label) && nzchar(as.character(rp$label)))
                           safe_feat_name(as.character(rp$label))
                         else
                           paste0(safe_feat_name(num_name), "_over_",
                                  safe_feat_name(den_name))

          num_idx <- which(feat_names_all == num_name)
          den_idx <- which(feat_names_all == den_name)

          if (length(num_idx) == 0) {
            message(sprintf("  Ratio '%s': numerator '%s' not found -- skipping",
                            ratio_label, num_name))
            next
          }
          if (length(den_idx) == 0) {
            message(sprintf("  Ratio '%s': denominator '%s' not found -- skipping",
                            ratio_label, den_name))
            next
          }

          # Extract raw pixel matrices with no hs/cs scaling (percentile=100,
          # threshold=0) so the ratio reflects true signal proportions.
          num_raw <- make_txt_mat(combined_snr_tmp, num_idx[1], .val,
                                   "ratio_num", 0, 100)
          den_raw <- make_txt_mat(combined_snr_tmp, den_idx[1], .val,
                                   "ratio_den", 0, 100)

          if (nrow(num_raw) == 0 || nrow(den_raw) == 0) next

          ratio_mat <- num_raw / den_raw
          ratio_mat[!is.finite(ratio_mat)] <- NA

          # SNR-filtered ratio (raw)
          write_if_nonempty(ratio_mat,
                            file.path(dirs$snr_filt, paste0(ratio_label, ".txt")))

          # 0-100 normalised
          write_if_nonempty(normalize_0_100(ratio_mat),
                            file.path(dirs$snr_norm, paste0(ratio_label, ".txt")))

          # hs/cs-suppressed ratio and its normalised version
          ratio_hs <- suppress_mat(ratio_mat, thresh, perc)
          write_if_nonempty(ratio_hs,
                            file.path(dirs$hs_filt, paste0(ratio_label, ".txt")))
          write_if_nonempty(normalize_0_100(ratio_hs),
                            file.path(dirs$hs_norm, paste0(ratio_label, ".txt")))
        }
      }
    }
  }

  invisible(out)
}
