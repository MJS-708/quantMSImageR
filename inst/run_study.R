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

heatmap_labs  <- vapply(cfg$samples, `[[`, character(1), "label")

# Unique run ID per sample: same as label when all labels are unique;
# append _1, _2, ... when a label appears more than once (biological replicates).
.seen <- list()
heatmap_order <- vapply(seq_along(heatmap_labs), function(i) {
  lab <- heatmap_labs[i]
  if (sum(heatmap_labs == lab) == 1L) return(lab)
  .seen[[lab]] <<- (.seen[[lab]] %||% 0L) + 1L
  paste0(lab, "_", .seen[[lab]])
}, character(1))
rm(.seen)

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

snr_thresh     <- cfg$parameters$snr_thresh     %||% 3
tiss_fc        <- cfg$parameters$tiss_fc        %||% 0.6
thresh         <- cfg$parameters$thresh         %||% 20
perc           <- cfg$parameters$perc           %||% 97
rot_clockwise  <- cfg$parameters$rot_clockwise  %||% 0
average_method <- cfg$parameters$average_method %||% "median"
baseline_label <- cfg$parameters$baseline_label %||% heatmap_labs[1]

render_report  <- cfg$output$render_report %||% TRUE
output_txt     <- cfg$output$output_txt    %||% TRUE
report_fn      <- cfg$output$report_fn     %||% paste0(cfg$study, "_SNR", snr_thresh)

# Feature overrides (optional)
feat_exclude <- cfg$features$exclude %||% NULL
feat_rename  <- if (!is.null(cfg$features$rename))
                  as.list(unlist(cfg$features$rename)) else NULL

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
  rename         = feat_rename
)

# ---------------------------------------------------------------------------
# Render HTML report
# ---------------------------------------------------------------------------
if (render_report) {

  if (is.null(markdown_rmd) || !nzchar(markdown_rmd))
    stop("render_report is TRUE but paths$markdown_rmd is not set in the YAML.")
  if (!file.exists(markdown_rmd))
    stop("Rmd not found: ", markdown_rmd)

  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

  # Expose objects and metadata expected by the Rmd
  combined       <- result$combined_snr
  # heatmap_order, heatmap_labs, baseline_label already in scope

  rmarkdown::render(
    markdown_rmd,
    output_file      = file.path(
      out_path,
      paste0(report_fn, "_response_SNRfiltered.html")
    ),
    intermediates_dir = out_path,
    knit_root_dir     = out_path
  )

  message("Report written to: ",
          file.path(out_path, paste0(report_fn, "_response_SNRfiltered.html")))
}

message("Done.")
