#!/usr/bin/env Rscript
# =============================================================================
# run_study.R  —  quantMSImageR DESI-MRM study runner
# =============================================================================
#
# Usage (command line):
#   Rscript run_study.R path/to/config.yaml
#
# Usage (interactive):
#   CONFIG_FILE <- "path/to/config.yaml"
#   source(system.file("run_study.R", package = "quantMSImageR"))
#
# The YAML must follow the structure in inst/config_template.yaml.
# =============================================================================

library(quantMSImageR)
library(Cardinal)
library(yaml)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------------------------------------------------------------------------
# Resolve config path
# ---------------------------------------------------------------------------
if (!exists("CONFIG_FILE")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0)
    stop("Provide a YAML config path:\n  Rscript run_study.R config.yaml")
  CONFIG_FILE <- args[1]
}

if (!file.exists(CONFIG_FILE))
  stop("Config file not found: ", CONFIG_FILE)

cfg <- yaml::read_yaml(CONFIG_FILE)

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
# Force UTF-8 on all character columns — Excel-on-Windows often saves
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
# These would form separate groups and silently break baseline matching.
.norm <- tolower(gsub("[ _]+", "_", heatmap_labs))
.near <- heatmap_labs[duplicated(.norm) | duplicated(.norm, fromLast = TRUE)]
if (length(.near)) {
  warning(sprintf(
    "Possible label typos — these labels differ only in case/underscores/spaces: %s\n  Verify your YAML.",
    paste(sort(unique(.near)), collapse = ", ")
  ))
}
rm(.norm, .near)

# Unique run ID per sample. Honour an explicit `run_id` field per sample if
# provided; otherwise fall back to label, appending _1, _2, … when a label
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

# Run IDs must be unique — otherwise pData$run can't disambiguate samples.
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

# Sample mapping table (original filenames ↔ run IDs ↔ group labels) for the report
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
tiss_fc        <- cfg$parameters$tiss_fc        %||% 0.6
thresh         <- cfg$parameters$thresh         %||% 20
perc           <- cfg$parameters$perc           %||% 97
rot_clockwise  <- cfg$parameters$rot_clockwise  %||% 0
average_method <- cfg$parameters$average_method %||% "median"
baseline_label    <- trimws(cfg$parameters$baseline_label %||% heatmap_labs[1])
heatmap_row_split <- cfg$parameters$heatmap_row_split %||% NULL

render_report  <- cfg$output$render_report  %||% TRUE
output_txt     <- cfg$output$output_txt     %||% TRUE
output_ratios  <- cfg$output$output_ratios  %||% TRUE
# report_fn is set per-SNR inside the render loop (see below)

# Feature overrides (optional)
feat_exclude <- cfg$features$exclude %||% NULL
feat_rename  <- if (!is.null(cfg$features$rename))
                  as.list(unlist(cfg$features$rename)) else NULL

# Ratio pairs (optional): list of {num, den, label} entries from YAML
feat_ratios  <- cfg$ratios %||% NULL

# ---------------------------------------------------------------------------
# Process acquisitions (always runs; controls output via flags)
# ---------------------------------------------------------------------------
result <- generate_txt_images(
  fns            = fns,
  data_path      = data_path,
  image_dir      = image_dir,
  lib_ion_path   = lib_ion_path,
  snr_thresh     = snr_thresh,
  tiss_fc        = tiss_fc,
  thresh         = thresh,
  perc           = perc,
  rot_clockwise  = rot_clockwise,
  average_method = average_method,
  output_txt     = output_txt,
  exclude        = feat_exclude,
  rename         = feat_rename,
  ratios         = feat_ratios,
  output_ratios  = output_ratios
)

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

  for (.i in seq_along(snr_objs)) {

    combined     <- snr_objs[[.i]]
    feature_meta <- build_feature_meta(combined, ion_lib_meta)

    # File-stem for this threshold. If the user explicitly pinned report_fn
    # AND only one threshold is requested, honour that; otherwise auto-name
    # per SNR so multiple runs don't clobber each other.
    report_fn <- if (user_pinned && length(snr_objs) == 1L)
      as.character(cfg$output$report_fn)[1] else
      paste0(cfg$study, "_SNR", snr_thresh[.i])

    out_file <- file.path(out_path,
                          paste0(report_fn, "_response_SNRfiltered.html"))

    message(sprintf("Rendering report %d/%d (SNR = %s) -> %s",
                     .i, length(snr_objs), snr_thresh[.i], out_file))

    rmarkdown::render(
      markdown_rmd,
      output_file       = out_file,
      intermediates_dir = out_path,
      knit_root_dir     = out_path
    )
  }
}

message("Done.")
