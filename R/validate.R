#' Collect validation findings
#'
#' A tiny accumulator so each check reports itself the same way and the caller
#' gets one consolidated answer instead of the first failure.
#'
#' @return An environment with `errors`, `warnings` and `rows`.
#' @keywords internal
#' @noRd
.chk_new = function(){
  e = new.env(parent = emptyenv())
  e$errors = character(); e$warnings = character(); e$rows = list()
  e
}

#' @param env Accumulator from `.chk_new()`.
#' @param check Short check name.
#' @param status `"ok"`, `"warning"` or `"error"`.
#' @param message Explanation; `NA` for a passing check.
#' @keywords internal
#' @noRd
.chk_add = function(env, check, status, message = NA_character_){
  env$rows[[length(env$rows) + 1L]] = data.frame(
    check = check, status = status, message = message,
    stringsAsFactors = FALSE)
  if(identical(status, "error"))   env$errors   = c(env$errors, message)
  if(identical(status, "warning")) env$warnings = c(env$warnings, message)
  invisible(env)
}

#' @keywords internal
#' @noRd
.chk_out = function(env){
  structure(list(errors   = env$errors,
                 warnings = env$warnings,
                 summary  = do.call(rbind, env$rows)),
            class = "quantValidation")
}

#' Print a validation result
#'
#' @param x A `quantValidation` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @family workflow
#' @export
print.quantValidation = function(x, ...){
  n_e = length(x$errors); n_w = length(x$warnings)
  cat(sprintf("quantMSImageR validation: %d error(s), %d warning(s), %d check(s)\n",
              n_e, n_w, nrow(x$summary)))
  if(n_e) cat("\nErrors:\n",   paste0("  - ", x$errors,   collapse = "\n"), "\n")
  if(n_w) cat("\nWarnings:\n", paste0("  - ", x$warnings, collapse = "\n"), "\n")
  if(!n_e && !n_w) cat("All checks passed.\n")
  invisible(x)
}

#' Validate a study YAML configuration
#'
#' Structural and cross-file checks on a study config before anything is
#' processed: required blocks present, referenced files and ion-library columns
#' actually exist, enumerated values spelled correctly, and the calibration and
#' internal-standard settings coherent with each other.
#'
#' Every check runs, so one call reports every problem rather than stopping at
#' the first. [runStudy()] calls this first and refuses to start if there are
#' errors -- these are mistakes that would otherwise surface as an obscure
#' failure halfway through a long run, or worse, as a plausible-looking result
#' computed from the wrong column.
#'
#' @param config A config list, or a path to a study YAML file.
#' @param check_paths Logical. Also check that the referenced data directory,
#'   ion library and calibration files exist on disk (default `TRUE`). Set
#'   `FALSE` to validate a config's structure away from its data.
#'
#' @return A `quantValidation` object with `errors`, `warnings` and a `summary`
#'   data frame of every check.
#'
#' @examples
#' # The shipped template is structurally valid; its placeholder paths are not
#' # meant to exist, so the path checks are skipped here.
#' cfg <- system.file("config_template.yaml", package = "quantMSImageR")
#' validateConfig(cfg, check_paths = FALSE)
#'
#' @family workflow
#' @export
validateConfig = function(config, check_paths = TRUE){

  `%||%` = function(a, b) if (is.null(a)) b else a

  if(is.character(config) && length(config) == 1L){
    if(!file.exists(config))
      stop("validateConfig: config file not found: ", config, call. = FALSE)
    config = yaml::read_yaml(config)
  }

  chk = .chk_new()

  # ---- required blocks -----------------------------------------------------
  for(blk in c("paths", "samples")){
    if(is.null(config[[blk]]))
      .chk_add(chk, sprintf("%s block", blk), "error",
               sprintf("config has no '%s:' block.", blk))
    else .chk_add(chk, sprintf("%s block", blk), "ok")
  }

  p   = config$paths
  par = config$parameters

  # ---- paths ---------------------------------------------------------------
  if(!is.null(p)){
    for(k in c("data_path", "out_path", "image_dir", "lib_ion_path")){
      if(is.null(p[[k]]) || !nzchar(p[[k]]))
        .chk_add(chk, sprintf("paths$%s", k), "error",
                 sprintf("paths$%s is not set.", k))
      else .chk_add(chk, sprintf("paths$%s", k), "ok")
    }

    if(isTRUE(check_paths)){
      for(k in c("data_path", "lib_ion_path")){
        v = p[[k]]
        if(!is.null(v) && nzchar(v) && !file.exists(v))
          .chk_add(chk, sprintf("%s exists", k), "error",
                   sprintf("paths$%s does not exist: %s", k, v))
        else if(!is.null(v) && nzchar(v))
          .chk_add(chk, sprintf("%s exists", k), "ok")
      }
    }
  }

  # ---- ion library ---------------------------------------------------------
  type_header = par$type_header %||% "Type"
  ion_lib = NULL
  if(isTRUE(check_paths) && !is.null(p$lib_ion_path) &&
     nzchar(p$lib_ion_path) && file.exists(p$lib_ion_path)){

    ion_lib = tryCatch(read.csv(p$lib_ion_path, check.names = FALSE),
                       error = function(e) NULL)
    if(is.null(ion_lib)){
      .chk_add(chk, "ion library readable", "error",
               sprintf("could not read the ion library: %s", p$lib_ion_path))
    } else {
      .chk_add(chk, "ion library readable", "ok")

      need = c("transition_id", "precursor_mz", "product_mz")
      miss = setdiff(need, colnames(ion_lib))
      if(length(miss))
        .chk_add(chk, "ion library columns", "error", sprintf(
          "ion library is missing required column(s): %s. Present: %s.",
          paste(miss, collapse = ", "),
          paste(colnames(ion_lib), collapse = ", ")))
      else .chk_add(chk, "ion library columns", "ok")

      if(!type_header %in% colnames(ion_lib))
        .chk_add(chk, "type_header column", "error", sprintf(
          "parameters$type_header = '%s' is not an ion-library column. Present: %s.",
          type_header, paste(colnames(ion_lib), collapse = ", ")))
      else .chk_add(chk, "type_header column", "ok")

      # A typo here would otherwise fall through to the default and quietly
      # draw the other heatmap, which is hard to notice in a long report.
      hs = par$heatmap_style
      if(!is.null(hs) && !hs %in% c("auto", "per_sample", "contribution"))
        .chk_add(chk, "heatmap_style", "error", sprintf(
          "parameters$heatmap_style = '%s'; use 'auto', 'per_sample' or 'contribution'.",
          hs))
      else .chk_add(chk, "heatmap_style", "ok")

      hsc = par$heatmap_scale
      if(!is.null(hsc) && !as.character(hsc) %in% c("shared", "feature"))
        .chk_add(chk, "heatmap_scale", "error", sprintf(
          "parameters$heatmap_scale = '%s'; use 'shared' or 'feature'.",
          as.character(hsc)))
      else .chk_add(chk, "heatmap_scale", "ok")

      hrs = par$heatmap_row_split
      if(!is.null(hrs) && nzchar(hrs) && !hrs %in% colnames(ion_lib))
        .chk_add(chk, "heatmap_row_split column", "warning", sprintf(
          "parameters$heatmap_row_split = '%s' is not an ion-library column; heatmap rows will not be split.",
          hrs))
      else .chk_add(chk, "heatmap_row_split column", "ok")
    }
  }

  # ---- samples -------------------------------------------------------------
  s = config$samples

  # A typo in `combine:` silently falls back to "panels", which for two spatial
  # halves of one tissue either errors deep inside bindPanels() or keeps only
  # the overlapping pixels. Catch it here, where the message can name it.
  #
  # Checked outside the ion-library block on purpose: this is a per-sample key,
  # so gating it on a readable library would skip it whenever check_paths is
  # FALSE - which is exactly when a config is being checked before its data
  # exists.
  .cmb = unlist(lapply(config$samples, function(x) x$combine))
  .bad = setdiff(unique(as.character(.cmb)), c("panels", "stitch"))
  if(length(.bad))
    .chk_add(chk, "sample combine mode", "error", sprintf(
      "samples$combine = %s; use 'panels' (different transitions over the same pixels) or 'stitch' (pieces of one tissue).",
      paste(sQuote(.bad), collapse = ", ")))
  else .chk_add(chk, "sample combine mode", "ok")
  if(!is.null(s)){
    if(!length(s)){
      .chk_add(chk, "samples listed", "error", "samples: is empty.")
    } else {
      .chk_add(chk, "samples listed", "ok")

      no_src = vapply(s, function(e)
        is.null(e$pos) && is.null(e$neg) && is.null(e$sections), logical(1))
      if(any(no_src))
        .chk_add(chk, "sample acquisitions", "error", sprintf(
          "%d sample(s) have neither pos:, neg: nor sections:.", sum(no_src)))
      else .chk_add(chk, "sample acquisitions", "ok")

      labs = vapply(s, function(e) as.character(e$label %||% NA_character_), character(1))
      if(anyNA(labs) || any(!nzchar(trimws(labs))))
        .chk_add(chk, "sample labels", "error",
                 "every sample needs a label: (it is the group name).")
      else .chk_add(chk, "sample labels", "ok")

      ids = vapply(s, function(e) as.character(e$run_id %||% NA_character_), character(1))
      dup = unique(ids[!is.na(ids) & duplicated(ids)])
      if(length(dup))
        .chk_add(chk, "unique run_id", "error", sprintf(
          "duplicate run_id value(s): %s. Each acquisition needs its own identifier.",
          paste(dup, collapse = ", ")))
      else .chk_add(chk, "unique run_id", "ok")

      bl = par$baseline_label
      if(!is.null(bl) && nzchar(bl) && !bl %in% labs)
        .chk_add(chk, "baseline_label", "warning", sprintf(
          "parameters$baseline_label = '%s' is not one of the sample labels (%s).",
          bl, paste(unique(labs), collapse = ", ")))
      else .chk_add(chk, "baseline_label", "ok")
    }
  }

  # ---- processing parameters ----------------------------------------------
  if(!is.null(par)){
    snr = suppressWarnings(as.numeric(unlist(par$snr_thresh)))
    if(!length(snr) || anyNA(snr) || any(snr < 0))
      .chk_add(chk, "snr_thresh", "error",
               "parameters$snr_thresh must be one or more non-negative numbers.")
    else .chk_add(chk, "snr_thresh", "ok")

    am = as.character(par$average_method %||% "median")
    if(!am %in% c("mean", "median"))
      .chk_add(chk, "average_method", "error", sprintf(
        "parameters$average_method = '%s'; use 'mean' or 'median'.", am))
    else .chk_add(chk, "average_method", "ok")

    # Internal standard
    is_name = par$is_name
    use_is  = !is.null(is_name) && nzchar(is_name) && !identical(is_name, "None")

    im = as.character(par$is_mode %||% "line")
    if(!im %in% c("line", "sample", "pixel", "window"))
      .chk_add(chk, "is_mode", "error", sprintf(
        "parameters$is_mode = '%s'; use line, window, sample or pixel.", im))
    else .chk_add(chk, "is_mode", "ok")

    if(identical(im, "window")){
      w = suppressWarnings(as.numeric(par$is_window %||% 15))
      if(is.na(w) || w < 2)
        .chk_add(chk, "is_window", "error",
                 "parameters$is_window must be a number of at least 2 when is_mode: window.")
      else .chk_add(chk, "is_window", "ok")
    }

    # The standard is matched on the ion library's type column, which is the
    # easiest thing to get wrong -- a transition name silently does nothing.
    if(use_is && !is.null(ion_lib) && type_header %in% colnames(ion_lib)){
      vals = unique(as.character(ion_lib[[type_header]]))
      if(!is_name %in% vals)
        .chk_add(chk, "is_name in library", "error", sprintf(
          "parameters$is_name = '%s' is not a value of the ion library's '%s' column. Values present: %s. Note this is a type label, not a transition name.",
          is_name, type_header, paste(vals, collapse = ", ")))
      else .chk_add(chk, "is_name in library", "ok")

      # The analyte-to-standard map is only needed when the panel actually has
      # several standards, and the library says how many. Checking here means a
      # multi-standard panel is caught before any data is read, and a
      # single-standard one is never asked for a column it does not need.
      n_is = sum(as.character(ion_lib[[type_header]]) == is_name, na.rm = TRUE)
      hdr  = .none(par$is_norm_header %||% "IS_norm")

      if(n_is > 1 && is.null(hdr))
        .chk_add(chk, "is_norm_header", "error", sprintf(
          "%d features are typed '%s', so parameters$is_norm_header cannot be None: each analyte must name the standard it uses.",
          n_is, is_name))
      else if(n_is > 1 && !hdr %in% colnames(ion_lib))
        .chk_add(chk, "is_norm_header", "error", sprintf(
          "%d features are typed '%s', but the ion library has no '%s' column mapping each analyte to its standard. Columns present: %s.",
          n_is, is_name, hdr, paste(colnames(ion_lib), collapse = ", ")))
      else if(n_is > 1){
        # Every analyte must map to something that is itself a standard.
        rows    = as.character(ion_lib[[type_header]]) != is_name
        tgt     = trimws(as.character(ion_lib[[hdr]])[rows])
        stds    = as.character(ion_lib$transition_id)[!rows]
        blank   = as.character(ion_lib$transition_id)[rows][!nzchar(tgt)]
        unknown = setdiff(tgt[nzchar(tgt)], stds)
        if(length(blank) || length(unknown))
          .chk_add(chk, "is_norm_header", "error", sprintf(
            "ion library column '%s' is incomplete: %s.", hdr,
            paste(c(if(length(blank))
                      sprintf("no standard named for %s",
                              paste(utils::head(blank, 5), collapse = ", ")),
                    if(length(unknown))
                      sprintf("%s is not a feature typed '%s'",
                              paste(utils::head(unknown, 5), collapse = ", "),
                              is_name)),
                  collapse = "; ")))
        else .chk_add(chk, "is_norm_header", "ok")
      }
      else .chk_add(chk, "is_norm_header", "ok")
    }
  }

  # ---- colours -------------------------------------------------------------
  if(!is.null(config$colours)){
    known = names(quantPalettes())
    for(k in names(config$colours)){
      v = as.character(config$colours[[k]])

      # cell_border is a single colour, not a palette: anything grDevices can
      # turn into a colour, or "none" for no border at all.
      if(identical(k, "cell_border")){
        ok = length(v) == 1L &&
             (tolower(v) %in% c("none", "") ||
              !inherits(try(grDevices::col2rgb(v), silent = TRUE), "try-error"))
        if(!ok)
          .chk_add(chk, "colours$cell_border", "error", sprintf(
            "colours$cell_border = '%s' is not a colour. Use a name or hex (e.g. 'white', '#FFFFFF'), or 'none'.",
            paste(v, collapse = ", ")))
        else .chk_add(chk, "colours$cell_border", "ok")
        next
      }

      ok = if(identical(k, "ion_image")) v %in% c(known, "viridis") else v %in% known
      if(!ok)
        .chk_add(chk, sprintf("colours$%s", k), "error", sprintf(
          "colours$%s = '%s' is not a known palette. Available: %s.",
          k, v, paste(c(known, "viridis"), collapse = ", ")))
      else .chk_add(chk, sprintf("colours$%s", k), "ok")
    }
  }

  # ---- output --------------------------------------------------------------
  if(!is.null(config$output$fig_dpi)){
    d = suppressWarnings(as.numeric(config$output$fig_dpi))
    if(length(d) != 1L || is.na(d) || d < 72)
      .chk_add(chk, "output$fig_dpi", "error", sprintf(
        "output$fig_dpi = '%s' must be a number of at least 72 (300 is the usual floor for a publication figure).",
        paste(config$output$fig_dpi, collapse = ", ")))
    else if(d > 600)
      .chk_add(chk, "output$fig_dpi", "warning", sprintf(
        "output$fig_dpi = %g will make a very large HTML report; 300 is usually enough.", d))
    else .chk_add(chk, "output$fig_dpi", "ok")
  }

  # ---- calibration ---------------------------------------------------------
  cal = config$calibration
  if(isTRUE(cal$enabled %||% cal$execute)){

    ct = as.character(cal$cal_type %||% "")
    if(!nzchar(ct))
      .chk_add(chk, "calibration cal_type", "error",
               "calibration is enabled but cal_type: is not set. Use 'cal' (on-slide) or 'std_addition' (on-tissue).")
    else if(!ct %in% c("cal", "std_addition"))
      .chk_add(chk, "calibration cal_type", "error", sprintf(
        "calibration$cal_type = '%s'; use 'cal' or 'std_addition'.", ct))
    else .chk_add(chk, "calibration cal_type", "ok")

    if(is.null(cal$cal_acquisition) || !nzchar(cal$cal_acquisition))
      .chk_add(chk, "calibration acquisition", "error",
               "calibration is enabled but cal_acquisition: is not set.")
    else .chk_add(chk, "calibration acquisition", "ok")

    if(isTRUE(check_paths)){
      for(k in c("cal_roi_csv", "cal_metadata")){
        v = cal[[k]]
        if(is.null(v) || !nzchar(v))
          .chk_add(chk, sprintf("calibration %s", k), "error",
                   sprintf("calibration is enabled but %s: is not set.", k))
        else if(!file.exists(v))
          .chk_add(chk, sprintf("calibration %s", k), "error",
                   sprintf("calibration$%s does not exist: %s", k, v))
        else .chk_add(chk, sprintf("calibration %s", k), "ok")
      }

      if(!is.null(cal$cal_metadata) && nzchar(cal$cal_metadata) &&
         file.exists(cal$cal_metadata)){
        cm   = tryCatch(read.csv(cal$cal_metadata, check.names = FALSE),
                        error = function(e) NULL)
        # `analyte` replaced `lipid`, which named the assay rather than the
        # role; both are accepted.
        a_col = if(is.null(cm)) NA_character_
                else if("analyte" %in% colnames(cm)) "analyte"
                else if("lipid" %in% colnames(cm)) "lipid"
                else NA_character_
        need = c("identifier", "amount_pg", "level")
        miss = if(is.null(cm)) c(need, "analyte")
               else c(setdiff(need, colnames(cm)),
                      if(is.na(a_col)) "analyte" else character(0))
        if(length(miss))
          .chk_add(chk, "calibration metadata columns", "error", sprintf(
            "calibration metadata is missing column(s): %s.",
            paste(miss, collapse = ", ")))
        else .chk_add(chk, "calibration metadata columns", "ok")

        # Calibration selects analytes by name, so an analyte that matches no
        # transition simply produces no curve -- silently, unless flagged.
        if(!is.null(cm) && !is.na(a_col) && !is.null(ion_lib) &&
           "transition_id" %in% colnames(ion_lib)){
          orphan = setdiff(unique(as.character(cm[[a_col]])),
                           as.character(ion_lib$transition_id))
          if(length(orphan))
            .chk_add(chk, "calibration analytes in library", "warning", sprintf(
              "%d calibration analyte(s) match no transition_id and will not be quantified: %s.",
              length(orphan), paste(utils::head(orphan, 5), collapse = ", ")))
          else .chk_add(chk, "calibration analytes in library", "ok")
        }
      }
    }
  }

  # ---- ratios --------------------------------------------------------------
  rt = config$ratios
  if(isTRUE(rt$enabled %||% rt$execute)){
    ent = rt$entries
    if(!length(ent)){
      .chk_add(chk, "ratio entries", "error",
               "ratios are enabled but no entries: are defined.")
    } else {
      bad = vapply(ent, function(e)
        is.null(e$num) || is.null(e$den), logical(1))
      if(any(bad))
        .chk_add(chk, "ratio entries", "error", sprintf(
          "%d ratio entry/entries lack num: or den:.", sum(bad)))
      else .chk_add(chk, "ratio entries", "ok")
    }
  }

  .chk_out(chk)
}
