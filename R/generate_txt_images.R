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
#'   (see [bind_panels()]). `run_study.R` builds this list automatically
#'   from the YAML `samples` section.
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
#' @param tiss_fc Numeric. SNR threshold used for the tissue fold-change layer
#'   (default `0.6`).
#' @param thresh Numeric. Cold-spot percentile passed to `imageR()` as
#'   `threshold` (default `20`).
#' @param perc Numeric. Hot-spot percentile passed to `imageR()` as
#'   `percentile` (default `97`).
#' @param rot_clockwise Integer (0-3). Number of 90 deg clockwise rotations
#'   applied via `pracma::rot90()` (default `0`).
#' @param average_method Character, `"mean"` or `"median"`. Statistic used to
#'   summarise the noise pixel vector in `int2snr()`. Applied consistently to
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
#'
#' @return Invisibly returns a named list:
#'   \describe{
#'     \item{`combined`}{Raw combined object (background pixels retained).}
#'     \item{`combined_snr`}{SNR-filtered intensity for the **first** threshold
#'       in `snr_thresh`; sub-threshold pixels `NA`. Provided for
#'       backward compatibility.}
#'     \item{`combined_snr_list`}{Named list of SNR-filtered objects, one per
#'       threshold in `snr_thresh`, named `"snr<value>"` (e.g. `"snr3"`).}
#'     \item{`combined_FC`}{Tissue fold-change layer (`snr` slot).}
#'     \item{`combined_NAbackground`}{Raw intensity with background set to `NA`.}
#'   }
#'
#' @seealso [int2snr()], [applySNR()], [back2NA()], [imageR()]
#'
#' @examples
#' \donttest{
#' # Section-mode run on the bundled synthetic data (no tissue mask needed at
#' # snr_thresh = 0), returning the processed objects without writing files.
#' fns <- list(list(neg = "example", section = "section01", label = "A"))
#' res <- generate_txt_images(
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
#' @export
generate_txt_images <- function(
  fns,
  data_path,
  image_dir,
  lib_ion_path,
  snr_thresh     = 0,
  tiss_fc        = 0.6,
  thresh         = 20,
  perc           = 97,
  rot_clockwise  = 0,
  average_method = "median",
  output_txt     = TRUE,
  exclude        = NULL,
  rename         = NULL,
  ratios         = NULL,
  output_ratios  = TRUE,
  snr_overrides  = NULL
) {

  average_method <- match.arg(average_method, c("mean", "median"))
  snr_thresh_vec <- sort(unique(as.numeric(unlist(snr_thresh))))

  # Normalise snr_overrides (a name -> threshold map, possibly a YAML list)
  # into a named numeric vector. Keys may be display (post-rename) names, so
  # we translate them to original ion-library names below, since int2snr runs
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

  # If every requested threshold is 0, the pipeline can skip the tissue mask
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
        pixels     = NA, percentile = percentile, overlay    = FALSE,
        feat_ind   = feat_ind, blank_back  = FALSE, text_image  = TRUE
      )
      pracma::rot90(as.matrix(result), k = rot_clockwise)
    }, error = function(e) {
      message("imageR failed (", value_label, "): ", e$message)
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
  load_and_prep_acq <- function(fn_name, section = NULL) {
    if (!is.null(section) && nzchar(section)) {
      rds_path <- file.path(data_path, paste0(fn_name, ".raw"),
                            paste0(section, ".RDS"))
      if (!file.exists(rds_path))
        stop("Section RDS not found: ", rds_path)
      obj <- readRDS(rds_path)
      if (!"sample_name" %in% names(pData(obj)))
        stop("RDS for section '", section, "' lacks pData$sample_name; ",
             "regenerate RDS with tissue_pixels/noise_pixels labels attached.")
      if (!is(obj, "quant_MSImagingExperiment"))
        obj <- as(obj, "quant_MSImagingExperiment")
    } else {
      obj  <- read_mrm(name = fn_name, folder = data_path, lib_ion_path = lib_ion_path)

      # Skip tissue mask entirely when no SNR > 0 is requested. Every pixel
      # is labelled `tissue_pixels` so downstream code that reads sample_name
      # (e.g. report quantiles) treats the whole acquisition as tissue.
      if (!.needs_mask) {
        pData(obj)$sample_name <- factor(
          rep("tissue_pixels", ncol(obj)),
          levels = c("tissue_pixels", "noise_pixels")
        )
        obj <- as(obj, "quant_MSImagingExperiment")
        obj <- trim_MSI(MSI_data = obj)
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
        # noise (conservative; they are excluded from tissue quantiles).
        key_obj <- paste(pData(obj)$x, pData(obj)$y, sep = "_")
        key_csv <- paste(tpdf$x,        tpdf$y,        sep = "_")
        midx    <- match(key_obj, key_csv)
        is_tiss <- rep(FALSE, ncol(obj))
        ok      <- !is.na(midx)
        is_tiss[ok] <- as.logical(tpdf$tissue_pixels[midx[ok]])

        n_unmatched <- sum(!ok)
        if (n_unmatched > 0)
          message(sprintf(
            "  '%s': %d / %d MSI pixels not present in tissue_pixels.csv -- labelled as noise.",
            fn_name, n_unmatched, ncol(obj)))

        pData(obj)$sample_name <- factor(
          ifelse(is_tiss, "tissue_pixels", "noise_pixels"),
          levels = c("tissue_pixels", "noise_pixels")
        )
      } else if (nrow(tpdf) == ncol(obj)) {
        # Legacy mask -- no coordinates, row-order alignment with current MSI
        pData(obj)$sample_name <- makeFactor(
          tissue_pixels = tpdf[["tissue_pixels"]],
          noise_pixels  = tpdf[["noise_pixels"]]
        )
      } else {
        stop(sprintf(
          paste0("tissue_pixels.csv for '%s' has %d rows but the loaded MSI ",
                 "has %d pixels, and the CSV has no x/y columns to match on.\n",
                 "  CSV: %s\n",
                 "  Regenerate the mask with the current package version (which ",
                 "writes x/y columns) so it remains portable across MRM panels:\n",
                 "    select_tissue_pixels(name = \"%s\", data_path = \"%s\", overwrite = TRUE)"),
          fn_name, nrow(tpdf), ncol(obj), tpdf_path, fn_name, data_path
        ))
      }

      obj <- as(obj, "quant_MSImagingExperiment")
    }
    obj <- trim_MSI(MSI_data = obj)
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
  # single MSI by coordinate-matched rbind. bind_panels() handles the typical
  # case where two MRM panels of the same tissue sample at slightly different
  # rates and therefore produce different pixel grids -- pixels are matched on
  # (x, y), features are stacked, so each matched pixel ends up with both
  # panels' transitions.
  load_and_prep_multiple <- function(fns_vec, label) {
    fns_vec <- as.character(unlist(fns_vec))
    objs    <- lapply(fns_vec, load_and_prep_acq)
    obj     <- objs[[1]]
    for (i in seq_along(objs)[-1])
      obj <- bind_panels(obj, objs[[i]], label = label)
    pData(obj)$run <- factor(rep(label, ncol(obj)))
    obj
  }


  # ----- Normalise fns -> fn_list / fn_labels ----------------------------
  # fns: character vector (backward-compat) OR list where each element is
  # a string (single acq) or named list with pos/neg/label fields.
  fn_list   <- as.list(fns)
  fn_labels <- vapply(fn_list, function(e)
    if (is.list(e)) e$label %||% (e$pos[1] %||% e$neg[1]) else as.character(e),
    character(1))

  # ----- Load and process acquisitions -----------------------------------

  combined          <- NULL
  combined_FC       <- NULL
  combined_snr_list <- vector("list", length(snr_thresh_vec))

  for (ind in seq_along(fn_list)) {
    fn_entry <- fn_list[[ind]]
    fn_label <- fn_labels[ind]
    message(sprintf("Loading %s (%d/%d)", fn_label, ind, length(fn_list)))

    if (is.list(fn_entry)) {
      pos_fns <- if (!is.null(fn_entry$pos)) unlist(fn_entry$pos) else NULL
      neg_fns <- if (!is.null(fn_entry$neg)) unlist(fn_entry$neg) else NULL
      section <- if (!is.null(fn_entry$section)) as.character(fn_entry$section) else NULL

      if (!is.null(section) && nzchar(section)) {
        # RDS-section mode: one section of one acquisition
        if (!is.null(pos_fns) && !is.null(neg_fns))
          stop("`section:` with dual polarity (pos + neg) is not supported")
        raw_name <- (neg_fns %||% pos_fns)[1]
        tissue   <- load_and_prep_acq(raw_name, section = section)
        pData(tissue)$run <- factor(rep(fn_label, ncol(tissue)))
      } else if (!is.null(pos_fns) && !is.null(neg_fns)) {
        # Both polarities: combine within each polarity then bind across
        pos_obj <- load_and_prep_multiple(pos_fns, label = fn_label)
        neg_obj <- load_and_prep_multiple(neg_fns, label = fn_label)
        tissue  <- bind_panels(pos_obj, neg_obj, label = fn_label)
      } else {
        # Single polarity (pos: OR neg: only)
        tissue <- load_and_prep_multiple(pos_fns %||% neg_fns, label = fn_label)
      }
    } else {
      # Plain string: backward-compatible single acquisition
      tissue <- load_and_prep_acq(fn_entry)
    }

    tissue_fc <- if (.needs_mask) int2snr(
      MSIobject = tissue, val_slot = "intensity", sample_type = "sample_name",
      noise = "tissue_pixels", tissue = "tissue_pixels",
      snr_thresh = tiss_fc, average = average_method
    ) else tissue

    # Per-feature overrides resolved against THIS acquisition's features.
    .ov <- resolve_overrides(as.character(fData(tissue)$name))

    # Compute one SNR-filtered object per requested threshold. A threshold of
    # 0 short-circuits to the raw `tissue` object (no SNR computation) -- the
    # smoke-test path -- so overrides are ignored for a 0-valued global report.
    for (si in seq_along(snr_thresh_vec)) {
      thr <- snr_thresh_vec[si]
      tissue_snr <- if (thr == 0) {
        tissue
      } else {
        tmp <- int2snr(
          MSIobject = tissue, val_slot = "intensity", sample_type = "sample_name",
          noise = "noise_pixels", tissue = "tissue_pixels",
          snr_thresh = thr, average = average_method,
          snr_overrides = .ov
        )
        applySNR(MSIobject = tmp, val_slot = "intensity")
      }

      if (is.null(combined_snr_list[[si]])) {
        combined_snr_list[[si]] <- tissue_snr
      } else {
        al <- align_features(combined_snr_list[[si]], tissue_snr)
        combined_snr_list[[si]] <- combine_MSIs(al$obj1, al$obj2)
      }
    }

    if (is.null(combined)) {
      combined    <- tissue
      combined_FC <- tissue_fc
    } else {
      al          <- align_features(combined,    tissue)
      combined    <- combine_MSIs(al$obj1, al$obj2)
      al_fc       <- align_features(combined_FC, tissue_fc)
      combined_FC <- combine_MSIs(al_fc$obj1, al_fc$obj2)
    }
  }

  combined_NAbackground <- if (.needs_mask) back2NA(
    combined, val_slot = "intensity",
    background = "noise_pixels", tissue = "tissue_pixels",
    sample_type = "sample_name"
  ) else combined

  # Apply display-name overrides to all objects
  if (!is.null(rename) && length(rename) > 0) {
    combined              <- apply_renames(combined,              rename)
    combined_snr_list     <- lapply(combined_snr_list, function(o) apply_renames(o, rename))
    combined_FC           <- apply_renames(combined_FC,           rename)
    combined_NAbackground <- apply_renames(combined_NAbackground, rename)
  }

  out <- list(
    combined              = combined,
    combined_snr          = combined_snr_list[[1]],
    combined_snr_list     = setNames(combined_snr_list, paste0("snr", snr_thresh_vec)),
    combined_FC           = combined_FC,
    combined_NAbackground = combined_NAbackground
  )

  if (!output_txt) return(invisible(out))

  # ----- Generate text-image files ---------------------------------------

  for (snr_i in seq_along(snr_thresh_vec)) {
    snr_t          <- snr_thresh_vec[snr_i]
    combined_snr_i <- combined_snr_list[[snr_i]]

    for (fn_label in fn_labels) {

      combined_snr_tmp  <- combined_snr_i[,        pData(combined_snr_i)$run        == fn_label]
      combined_FC_tmp   <- combined_FC[,            pData(combined_FC)$run           == fn_label]
      combined_back_tmp <- combined_NAbackground[,  pData(combined_NAbackground)$run == fn_label]

      image_path <- file.path(image_dir, fn_label)
      dirs <- list(
        snr_filt  = file.path(image_path, sprintf("intensity_SNRfiltered%s",                   snr_t)),
        snr_norm  = file.path(image_path, sprintf("response_SNRfiltered%s_NORM",               snr_t)),
        tissue_fc = file.path(image_path, sprintf("tissue-FC%s",                               tiss_fc)),
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
        mat_snr <- make_txt_mat(combined_snr_tmp, feat_ind, "intensity",
                                 "DESI-MRM response - S/N filtered", 0, 100)
        write_if_nonempty(mat_snr, file.path(dirs$snr_filt, paste0(feat_name, ".txt")))

        # 0-100 normalised SNR
        scaled_snr <- normalize_0_100(mat_snr)
        write_if_nonempty(scaled_snr, file.path(dirs$snr_norm, paste0(feat_name, ".txt")))

        # COMBINED: embed global max in [1,1] for cross-sample colour scaling
        if (nrow(scaled_snr) > 0) {
          global_max <- max(
            as.numeric(spectraData(combined_snr_i[feat_ind, ])[["intensity"]]),
            na.rm = TRUE
          )
          scaled_snr_combined       <- scaled_snr
          scaled_snr_combined[1, 1] <- global_max
          write_if_nonempty(scaled_snr_combined,
                            file.path(dirs$combined, paste0(feat_name, ".txt")))
        }

        # tissue_fc and raw don't vary with SNR threshold -- write only on first pass
        if (snr_i == 1L) {
          mat_fc <- make_txt_mat(combined_FC_tmp, feat_ind, "snr",
                                  "Ratio to tissue", 0, 100)
          write_if_nonempty(mat_fc, file.path(dirs$tissue_fc, paste0(feat_name, ".txt")))

          mat_raw <- make_txt_mat(combined_back_tmp, feat_ind, "intensity",
                                   "DESI-MRM response", 0, 100)
          write_if_nonempty(mat_raw, file.path(dirs$raw, paste0(feat_name, ".txt")))
        }

        # SNR-filtered + hot/cold-spot removal
        mat_hs <- make_txt_mat(combined_snr_tmp, feat_ind, "intensity",
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
          num_raw <- make_txt_mat(combined_snr_tmp, num_idx[1], "intensity",
                                   "ratio_num", 0, 100)
          den_raw <- make_txt_mat(combined_snr_tmp, den_idx[1], "intensity",
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
