#' Run a full DESI-MRM study from a YAML configuration
#'
#' Executes the whole workflow for a study described by a single YAML file:
#' loads every acquisition, applies the tissue masks, runs SNR filtering at each
#' requested threshold, optionally builds calibration models and converts tissue
#' pixels to estimated amounts, writes per-feature `.txt` images, and renders one
#' HTML report per SNR threshold.
#'
#' @details
#' The configuration format is documented by the template shipped with the
#' package:
#'
#' ```r
#' file.show(system.file("config_template.yaml", package = "quantMSImageR"))
#' ```
#'
#' A study can also be run from the command line, which uses this function
#' internally:
#'
#' ```
#' Rscript -e 'quantMSImageR::run_study("config.yaml")'
#' ```
#'
#' @param config_file Path to the study YAML configuration.
#'
#' @return Invisibly, the list produced by [generate_txt_images()] for the
#'   study, with a `calibrated` element added when the configuration's
#'   `calibration:` block is enabled. Called for its side effects: `.txt`
#'   images, the calibrated `.RDS` and the HTML reports written under the
#'   configured output paths.
#'
#' @examples
#' # A minimal study driven from a config, using the bundled synthetic data.
#' # snr_thresh = 0 skips SNR filtering and the tissue-mask requirement, and
#' # both outputs are switched off, so this runs without writing a report.
#' \donttest{
#' out <- file.path(tempdir(), "qmsi_study")
#' dir.create(out, showWarnings = FALSE)
#'
#' cfg <- list(
#'   study = "example_study",
#'   paths = list(
#'     data_path    = system.file("extdata", package = "quantMSImageR"),
#'     out_path     = out,
#'     image_dir    = file.path(out, "images"),
#'     lib_ion_path = system.file("extdata", "example_ion_library.csv",
#'                                package = "quantMSImageR")
#'   ),
#'   samples = list(
#'     list(neg = "example", sections = list("section01"),
#'          run_id = "S1", label = "A")
#'   ),
#'   parameters = list(snr_thresh = 0),
#'   output = list(render_report = FALSE, output_txt = FALSE)
#' )
#'
#' cfg_file <- file.path(out, "config.yaml")
#' yaml::write_yaml(cfg, cfg_file)
#'
#' res <- run_study(cfg_file)
#' names(res)
#' }
#'
#' @family workflow
#' @export
run_study <- function(config_file) {

  `%||%` <- function(a, b) if (is.null(a)) b else a

  if (!file.exists(config_file))
    stop("Config file not found: ", config_file, call. = FALSE)

  cfg <- yaml::read_yaml(config_file)

  # Validate before anything is read or written. Every check runs, so a config
  # with several problems reports all of them once rather than failing on each
  # in turn across successive runs.
  .v <- validate_config(cfg)
  if (length(.v$warnings))
    for (w in .v$warnings) warning("Config: ", w, call. = FALSE)
  if (length(.v$errors)) {
    print(.v)
    stop("Config validation failed with ", length(.v$errors),
         " error(s); see above. Nothing has been processed.", call. = FALSE)
  }


  # ---------------------------------------------------------------------------
  # Expand samples: when a YAML entry has `sections:` (list of pre-saved RDS
  # basenames inside the `.raw` folder), emit one internal sample per section.
  # All expanded entries share the same `label`, so the replicate-ID logic
  # below assigns unique run IDs (e.g. 1a_1, 1a_2, ...).
  # ---------------------------------------------------------------------------
  .expanded <- list()
  for (s in cfg$samples) {
    if (!is.null(s$sections) && length(s$sections) > 0) {
      for (sec_entry in s$sections) {
        # Per-section labels: `sec_entry` is a list with `file` (+ optional `label`).
        # Shared label: `sec_entry` is a plain string -> inherit parent `label`.
        if (is.list(sec_entry)) {
          sec_file  <- as.character(sec_entry$file  %||% sec_entry[[1]])
          sec_label <- if (!is.null(sec_entry$label))
                         as.character(sec_entry$label) else s$label
        } else {
          sec_file  <- as.character(sec_entry)
          sec_label <- s$label
        }
        e <- s
        e$sections <- NULL
        e$section  <- sec_file
        e$label    <- sec_label
        .expanded[[length(.expanded) + 1L]] <- e
      }
    } else {
      .expanded[[length(.expanded) + 1L]] <- s
    }
  }
  cfg$samples <- .expanded
  rm(.expanded)

  # ---------------------------------------------------------------------------
  # Unpack config
  # ---------------------------------------------------------------------------
  data_path     <- cfg$paths$data_path
  out_path      <- cfg$paths$out_path
  image_dir     <- file.path(cfg$paths$image_dir, cfg$study)
  lib_ion_path  <- cfg$paths$lib_ion_path
  markdown_rmd  <- cfg$paths$markdown_rmd

  # Ion library metadata for optional row annotations in heatmap.
  # Force UTF-8 on all character columns -- Excel-on-Windows often saves
  # Windows-1252, which trips trimws()/sub() in the report.
  ion_lib_meta <- if (!is.null(lib_ion_path) && nzchar(lib_ion_path %||% "") &&
                      file.exists(lib_ion_path))
    read.csv(lib_ion_path, check.names = FALSE) else NULL
  if (!is.null(ion_lib_meta)) {
    .chr <- vapply(ion_lib_meta, is.character, logical(1))
    ion_lib_meta[.chr] <- lapply(ion_lib_meta[.chr], function(x)
      iconv(x, from = "", to = "UTF-8", sub = ""))
  }

  heatmap_labs  <- trimws(vapply(cfg$samples, `[[`, character(1), "label"))

  # Warn if any labels look like near-duplicates (differ only in case / _ vs space).
  # Biological replicates legitimately share a label, so we only flag a collision
  # when MORE THAN ONE distinct original string normalises to the same value --
  # i.e. true near-duplicates like "Ctrl_M" vs "Ctrl M".
  .norm   <- tolower(gsub("[ _]+", "_", heatmap_labs))
  .by_n   <- split(heatmap_labs, .norm)
  .typos  <- .by_n[vapply(.by_n, function(x) length(unique(x)) > 1L, logical(1))]
  if (length(.typos)) {
    .near <- sort(unique(unlist(.typos)))
    warning(sprintf(
      "Possible label typos -- these labels differ only in case/underscores/spaces: %s\n  Verify your YAML.",
      paste(.near, collapse = ", ")
    ))
  }
  rm(.norm, .by_n, .typos)

  # Unique run ID per sample. Honour an explicit `run_id` field per sample if
  # provided; otherwise fall back to label, appending _1, _2, ... when a label
  # appears more than once (biological replicates).
  .explicit_id <- vapply(cfg$samples, function(s) {
    if (!is.null(s$run_id) && nzchar(trimws(as.character(s$run_id))))
      trimws(as.character(s$run_id)) else NA_character_
  }, character(1))

  .seen <- list()
  heatmap_order <- vapply(seq_along(heatmap_labs), function(i) {
    if (!is.na(.explicit_id[i])) return(.explicit_id[i])
    lab <- heatmap_labs[i]
    if (sum(heatmap_labs == lab) == 1L) return(lab)
    .seen[[lab]] <<- (.seen[[lab]] %||% 0L) + 1L
    paste0(lab, "_", .seen[[lab]])
  }, character(1))
  rm(.seen, .explicit_id)

  # Run IDs must be unique -- otherwise pData$run can't disambiguate samples.
  if (anyDuplicated(heatmap_order)) {
    .dup <- unique(heatmap_order[duplicated(heatmap_order)])
    stop("Duplicate run_id values: ", paste(.dup, collapse = ", "),
         ". Provide unique `run_id:` entries in the YAML.")
  }

  # Build fns list: each element has pos (string or list), neg (string or list),
  # and label.  pos: / neg: may be a single string or a YAML sequence.
  # Use the unique run ID as the label so pData$run is unique per sample.
  fns <- lapply(seq_along(cfg$samples), function(i) {
    s <- cfg$samples[[i]]
    list(
      pos     = if (!is.null(s$pos)) as.character(unlist(s$pos)) else NULL,
      neg     = if (!is.null(s$neg)) as.character(unlist(s$neg)) else NULL,
      section = if (!is.null(s$section)) as.character(s$section) else NULL,
      label   = heatmap_order[i]
    )
  })

  # Sample mapping table (original filenames <-> run IDs <-> group labels) for the report
  sample_map <- data.frame(
    run_id     = heatmap_order,
    group      = heatmap_labs,
    pos_files  = vapply(fns, function(f)
      if (!is.null(f$pos)) paste(f$pos, collapse = "; ") else "", character(1)),
    neg_files  = vapply(fns, function(f)
      if (!is.null(f$neg)) paste(f$neg, collapse = "; ") else "", character(1)),
    section    = vapply(fns, function(f)
      if (!is.null(f$section)) f$section else "", character(1)),
    stringsAsFactors = FALSE
  )
  # Drop section column if unused across the whole study
  if (all(sample_map$section == "")) sample_map$section <- NULL

  snr_thresh     <- as.numeric(unlist(cfg$parameters$snr_thresh %||% 0))
  # Optional per-analyte SNR overrides: named map (analyte -> SNR). These
  # analytes use their own threshold in every report; all others use the
  # global snr_thresh above. Keys may be display (post-rename) or original
  # ion-library names.
  snr_overrides  <- cfg$parameters$snr_overrides %||% NULL
  thresh         <- cfg$parameters$thresh         %||% 20
  perc           <- cfg$parameters$perc           %||% 97
  rot_clockwise  <- cfg$parameters$rot_clockwise  %||% 0
  average_method <- cfg$parameters$average_method %||% "median"
  is_name        <- cfg$parameters$is_name        %||% NULL
  is_norm_header  <- cfg$parameters$is_norm_header  %||% "IS_norm"
  is_mode        <- cfg$parameters$is_mode        %||% "line"
  is_window      <- cfg$parameters$is_window      %||% 15
  type_header    <- cfg$parameters$type_header    %||% "Type"
  remove_IS      <- cfg$parameters$remove_IS      %||% TRUE
  baseline_label    <- trimws(cfg$parameters$baseline_label %||% heatmap_labs[1])
  heatmap_row_split <- cfg$parameters$heatmap_row_split %||% NULL

  # Report palettes (see the `colours:` block of config_template.yaml). Read
  # here so the report picks them up from this frame at render time.
  pal_ion     <- cfg$colours$ion_image %||% "heatmap0"
  pal_heatmap <- cfg$colours$heatmap   %||% "heatmap2"
  pal_group   <- cfg$colours$group     %||% "hat"
  pal_feature <- cfg$colours$feature   %||% "reading"
  hm_cell_border <- cfg$colours$cell_border %||% "white"

  fig_dpi        <- cfg$output$fig_dpi        %||% 300
  render_report  <- cfg$output$render_report  %||% TRUE
  output_txt     <- cfg$output$output_txt     %||% TRUE
  output_ratios  <- cfg$output$output_ratios  %||% TRUE
  # report_fn is set per-SNR inside the render loop (see below)

  # ---- Feature control (optional) -------------------------------------------
  # exclude: "None"/null/empty -> keep all; otherwise a list of ion-library names.
  feat_exclude <- cfg$features$exclude
  if (is.null(feat_exclude) || length(feat_exclude) == 0 ||
      (length(feat_exclude) == 1 &&
       tolower(trimws(as.character(feat_exclude))) %in% c("none", "")))
    feat_exclude <- NULL else
    feat_exclude <- as.character(unlist(feat_exclude))

  # rename: gated by `enabled:`. New style: rename: {enabled: T/F, "old":"new"}.
  # Back-compat: a rename block with no toggle key is treated as on.
  feat_rename <- NULL
  .rn <- cfg$features$rename
  if (!is.null(.rn)) {
    .rn_flag <- .rn$enabled %||% .rn$execute      # accept either key
    .do_rn <- if (!is.null(.rn_flag)) isTRUE(.rn_flag) else TRUE
    .rn$enabled <- NULL; .rn$execute <- NULL
    if (.do_rn && length(.rn) > 0) feat_rename <- as.list(unlist(.rn))
  }

  # ratios: gated by `enabled:`. New style: ratios: {enabled: T/F, entries: [..]}.
  # Back-compat: a bare top-level sequence (no toggle/entries) is used as-is.
  feat_ratios <- NULL
  .rt <- cfg$ratios
  if (!is.null(.rt)) {
    .rt_has_flag <- any(c("enabled", "execute", "entries") %in% names(.rt))
    if (.rt_has_flag) {
      if (isTRUE(.rt$enabled %||% .rt$execute)) feat_ratios <- .rt$entries %||% NULL
    } else {
      feat_ratios <- .rt
    }
  }

  # ---------------------------------------------------------------------------
  # Process acquisitions (always runs; controls output via flags)
  # ---------------------------------------------------------------------------
  result <- generate_txt_images(
    fns            = fns,
    data_path      = data_path,
    image_dir      = image_dir,
    lib_ion_path   = lib_ion_path,
    snr_thresh     = snr_thresh,
    thresh         = thresh,
    perc           = perc,
    rot_clockwise  = rot_clockwise,
    average_method = average_method,
    output_txt     = output_txt,
    exclude        = feat_exclude,
    rename         = feat_rename,
    ratios         = feat_ratios,
    output_ratios  = output_ratios,
    snr_overrides  = snr_overrides,
    is_name        = is_name,
    is_norm_header  = is_norm_header,
    is_mode        = is_mode,
    is_window      = is_window,
    remove_IS      = remove_IS,
    type_header    = type_header
  )

  # ---------------------------------------------------------------------------
  # Calibration (optional): gated by calibration.enabled. Builds a response-vs-
  # amount curve per analyte from a standards acquisition, then converts the
  # study's tissue-pixel intensities to pg/pixel + pg/mm2. Writes a calibrated
  # RDS alongside the reports.
  # ---------------------------------------------------------------------------
  .cal <- cfg$calibration
  if (!is.null(.cal) && isTRUE(.cal$enabled %||% .cal$execute)) {
    message("Calibration enabled: building curves from '", .cal$cal_acquisition, "'.")

    .use_is  <- !is.null(is_name) && nzchar(is_name) && !identical(is_name, "None")
  # The study is quantified from whichever layer the workflow produced. Using
  # un-normalised standards against an IS-normalised tissue would mix units and
  # silently produce wrong amounts, so the default follows the workflow.
  cal_val  <- .cal$val_slot         %||% (if (.use_is) "response" else "intensity")

    # No fallback: "cal" and "std_addition" do materially different things, so
    # the config must say which.
    cal_type      <- .cal$cal_type
    cal_weighting <- .cal$weighting %||% "1/x"
    cal_max_oor   <- .cal$max_out_of_range %||% 0.1
    if (is.null(cal_type) || !nzchar(cal_type))
      stop("Calibration is enabled but `calibration: cal_type:` is not set. ",
           "Use \"cal\" for standards on the slide, or \"std_addition\" for ",
           "standard addition on tissue.", call. = FALSE)
    bg_level <- .cal$background_level %||% "background"
    q_pixels <- as.character(unlist(.cal$quantify_pixels %||% "tissue_pixels"))

    # 1. Load the calibration acquisition and label its Cal ROIs.
    cal_obj <- read_mrm(name = .cal$cal_acquisition, folder = data_path,
                        lib_ion_path = lib_ion_path, overwrite = FALSE,
                        type_header = type_header,
                        is_norm_header = is_norm_header)
    cal_obj <- as(cal_obj, "quant_MSImagingExperiment")

  # Normalise the standards the same way as the study.
  if (.use_is && identical(cal_val, "response"))
    cal_obj <- int2response(cal_obj, val_slot = "intensity",
                            normalisation = "internal_standard",
                            IS_name = is_name, is_norm_header = is_norm_header,
                            mode = is_mode, window = is_window,
                            remove_IS = remove_IS)

    roi <- read.csv(.cal$cal_roi_csv)
    if (all(c("x", "y") %in% names(roi))) {
      .m <- match(paste(pData(cal_obj)$x, pData(cal_obj)$y, sep = "_"),
                  paste(roi$x, roi$y, sep = "_"))
      pData(cal_obj)$sample_type <- roi$sample_type[.m]
      pData(cal_obj)$identifier  <- roi$identifier[.m]
    } else {
      pData(cal_obj)$sample_type <- roi$sample_type
      pData(cal_obj)$identifier  <- roi$identifier
    }

    cal_metadata <- read.csv(.cal$cal_metadata)

    # 2. Summarise calibration levels -> fit curves.
    cal_obj <- summarise_cal_levels(cal_obj, cal_metadata, val_slot = cal_val,
                                    cal_label = "Cal", id = "identifier")
    cal_obj <- create_cal_curve(cal_obj, cal_type = cal_type,
                                background = bg_level, weighting = cal_weighting)

    # 3. Apply curves to the study's tissue pixels. The imaging workflow labels
    #    pixels in `sample_name` (tissue_pixels/background_pixels), so bridge via
    #    pixel_header = "sample_name".
    combined_cal <- result$combined
    # Carry the whole calibrationInfo across, not just cal_list: the report's
    # calibration section also reports r2_df, and keeping cal_response_data with
    # the calibrated object preserves its provenance.
    calibrationData(combined_cal) <- calibrationData(cal_obj)
    combined_cal <- int2conc(combined_cal, val_slot = cal_val,
                             pixel_header = "sample_name", pixels = q_pixels,
                             max_out_of_range = cal_max_oor)

    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    .cal_out <- file.path(out_path, paste0(cfg$study, "_calibrated.RDS"))
    saveRDS(combined_cal, .cal_out)
    message("Calibration written to: ", .cal_out)

    # Picked up by the report's section 6 (see quantMSImageR___general_heatmap.Rmd).
    calibrated <- combined_cal
  }

  # ---------------------------------------------------------------------------
  # Render HTML report
  # ---------------------------------------------------------------------------
  if (render_report) {

    # Fall back to the Rmd shipped with the installed package when the YAML
    # doesn't set paths$markdown_rmd. This keeps the single source of truth in
    # inst/ and avoids stale user-side copies drifting from the package.
    if (is.null(markdown_rmd) || !nzchar(markdown_rmd)) {
      markdown_rmd <- system.file("quantMSImageR___general_heatmap.Rmd",
                                   package = "quantMSImageR")
    }
    if (!nzchar(markdown_rmd) || !file.exists(markdown_rmd))
      stop("Rmd not found: ", markdown_rmd,
           "\n  Set paths$markdown_rmd in the YAML, or reinstall the package.")

    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

    # Expose objects and metadata expected by the Rmd
    ratios_cfg <- feat_ratios   # ratio pairs for section 6
    # heatmap_order, heatmap_labs, baseline_label already in scope

    # Render one report per SNR threshold. Each report's `combined` object is
    # the SNR-filtered MSI for that threshold; `report_fn` is updated per
    # iteration so the embedded xlsx (Fold_Change, Cor_*) lands beside the
    # matching HTML.
    snr_objs    <- result$combined_snr_list
    user_pinned <- !is.null(cfg$output$report_fn) && nzchar(cfg$output$report_fn)

    # The list is named "snr<value>" in sorted-threshold order -- the authoritative
    # threshold per report (snr_thresh from the YAML may be in a different order).
    snr_report_vals <- as.numeric(sub("^snr", "", names(snr_objs)))

    # Normalise overrides to a named numeric vector once (keys as written in YAML)
    snr_ov_vec <- if (!is.null(snr_overrides) && length(snr_overrides) > 0)
      vapply(snr_overrides, as.numeric, numeric(1)) else NULL

    for (.i in seq_along(snr_objs)) {

      combined     <- snr_objs[[.i]]
      feature_meta <- build_feature_meta(combined, ion_lib_meta)
      snr_report   <- snr_report_vals[.i]   # global SNR for this report

      # Effective per-feature SNR for the report's feature table. combined has
      # already been through rename, so fData names are DISPLAY names. An override
      # keyed by an original library name is mapped to its display name here.
      .disp <- as.character(fData(combined)$name)
      snr_used <- stats::setNames(rep(snr_report, length(.disp)), .disp)
      if (snr_report > 0 && !is.null(snr_ov_vec)) {
        for (k in names(snr_ov_vec)) {
          disp_k <- if (!is.null(feat_rename) && k %in% names(feat_rename))
                      as.character(feat_rename[[k]]) else k
          if (disp_k %in% names(snr_used)) snr_used[disp_k] <- snr_ov_vec[[k]]
        }
      }

      # File-stem for this threshold. If the user explicitly pinned report_fn
      # AND only one threshold is requested, honour that; otherwise auto-name
      # per SNR so multiple runs don't clobber each other.
      report_fn <- if (user_pinned && length(snr_objs) == 1L)
        as.character(cfg$output$report_fn)[1] else
        paste0(cfg$study, "_SNR", snr_report)

      out_file <- file.path(out_path,
                            paste0(report_fn, "_response_SNRfiltered.html"))

      message(sprintf("Rendering report %d/%d (SNR = %s) -> %s",
                       .i, length(snr_objs), snr_report, out_file))

      rmarkdown::render(
        markdown_rmd,
        output_file       = out_file,
        intermediates_dir = out_path,
        knit_root_dir     = out_path
      )
    }
  }

  message("Done.")

  invisible(result)
}
