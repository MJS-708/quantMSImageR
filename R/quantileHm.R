#' Quantile heatmap of MSI features across samples
#'
#' For each feature, computes the nth quantile of pixel intensities within each
#' acquisition run, then z-scores the resulting per-feature profile across
#' samples (clipped to \[-1, 1\]) and renders a `ComplexHeatmap::Heatmap`.
#'
#' **Rows are samples and columns are features**, in the order given by
#' `heatmap_order`. Studies normally have more samples than features, so this
#' puts the long dimension vertically where there is room for it and keeps
#' sample names horizontally readable.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` object.
#' @param quant_val Numeric in (0, 1). Quantile to summarise per feature per
#'   sample (e.g. `0.95` for the 95th percentile).
#' @param heatmap_order Character vector of run names in the desired column order.
#'   Must match values in `pData(MSIobject)$run`. Defaults to `NA` (use
#'   discovery order).
#' @param heatmap_labs Character vector of display labels, one per entry in
#'   `heatmap_order`. Used to create column-split groups. Defaults to `NA`
#'   (no splitting).
#' @param feature_split Optional character/factor vector of length equal to the
#'   number of features, used to split the heatmap columns into labelled groups.
#'   Typically derived from an ion library metadata column (e.g. `Met-1`).
#'   Defaults to `NULL` (no splitting).
#' @param feature_split_name Display name for the feature-group annotation
#'   legend (e.g. `"Met-1"`, `"Pathway"`). Defaults to `"Pathway"`.
#' @param group_split_name Display name for the sample-group annotation legend
#'   (e.g. `"Treatment"`, `"Group"`). Defaults to `"Group"`.
#' @param palette Character. Diverging colour ramp for the z-scores:
#'   `"heatmap2"` (default, blue-white-red) or `"heatmap0"`. See
#'   [quantPalettes()].
#' @param group_palette,feature_palette Character. Qualitative palettes for the
#'   sample-group and feature-group colour bars (defaults `"hat"` and
#'   `"reading"`). Any [grDevices::hcl.colors()] palette name also works.
#' @param row_split,row_split_name,column_split_name Deprecated aliases for
#'   `feature_split`, `feature_split_name` and `group_split_name`. They were
#'   named for the axis each grouping landed on before samples and features
#'   swapped axes; supplying them still works.
#' @param cell_size Numeric. Side of a heatmap cell, in millimetres (default
#'   `8`). Giving the body an absolute size is what makes cells square; leaving
#'   it to `ComplexHeatmap` stretches them to fill the device, which produces
#'   very oblong cells when there are far more features than samples. Set to
#'   `NA` to restore the fill-the-device behaviour.
#' @param cell_border Colour for the line drawn around every cell (default
#'   `"white"`), which separates neighbouring cells instead of letting equal
#'   colours merge into one block. `NA` draws no border.
#' @param max_aspect Numeric >= 1. Largest cell width-to-height ratio allowed
#'   when one dimension has many more entries than the other (default `1.5`).
#'   Cells stay square until the counts differ by more than four-fold; beyond
#'   that the shorter dimension's cells are widened up to this ratio so the
#'   plotting area is not reduced to a sliver.
#' @param fontsize Numeric. Point size for the sample and feature labels and the
#'   annotation names (default `8`). `ComplexHeatmap` defaults to 12, which
#'   crowds the panel once feature names are long enough to need rotating.
#'
#' @return A `ComplexHeatmap::Heatmap` object (rows = samples, columns =
#'   features), with an absolutely-sized body unless `cell_size` is `NA`.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' hm <- quantileHm(obj, quant_val = 0.5,
#'                   heatmap_order = "section01", heatmap_labs = "A")
#'
#' @family visualisation
#' @export

quantileHm = function(MSIobject, quant_val, heatmap_order = NA, heatmap_labs = NA,
                        feature_split = NULL,
                        feature_split_name = "Pathway",
                        group_split_name = "Group",
                        cell_size = 8, max_aspect = 1.5,
                        cell_border = "white", fontsize = 8,
                        palette = c("heatmap2", "heatmap0"),
                        group_palette = "hat",
                        feature_palette = "reading",
                        row_split = NULL,
                        row_split_name = NULL,
                        column_split_name = NULL) {

  # Backward compatibility: these were named for where they appeared before
  # samples and features swapped axes. Feature grouping used to split rows and
  # sample grouping used to split columns; both now do the opposite, so the
  # arguments are named for what they group rather than where it lands.
  if (!is.null(row_split))         feature_split      <- row_split
  if (!is.null(row_split_name))    feature_split_name <- row_split_name
  if (!is.null(column_split_name)) group_split_name   <- column_split_name

  palette <- match.arg(palette)

  # Shared with contributionHm(), so the two views of one study cannot
  # disagree about what the underlying numbers are.
  #
  # This also fixes heatmap_order = NA, the documented default: the previous
  # code summarised every run and then indexed the result by `heatmap_order`,
  # so the default subscripted by NA and returned a column of NAs. Only the
  # report's always-supplied argument kept it working.
  out_matrix <- .quantile_matrix(MSIobject, quant_val, heatmap_order,
                                 context = "quantileHm")

  # Scale each row to the percentage of its maximum value, avoiding division by 0
  row_max <- matrixStats::rowMaxs(out_matrix, na.rm = TRUE)

  # Replace zeros in row_max with 1 to avoid NaN during division
  row_max[row_max == 0] <- 1

  # Scale the matrix
  scaled_out_matrix <- sweep(out_matrix, 1, row_max, FUN = "/")

  # Z-score scaling with clipping between -1 and 1.
  # Single-sample: sd = NA so z-scoring is undefined; use 0 (neutral) instead.
  if (ncol(out_matrix) == 1) {
    z_matrix <- matrix(0, nrow = nrow(out_matrix), ncol = 1,
                       dimnames = dimnames(out_matrix))
  } else {
    # t() ensures result is n_features x n_samples (apply over rows returns transposed)
    z_matrix <- t(apply(out_matrix, 1, function(x) {
      z <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      pmin(pmax(z, -1), 1)  # Clip to [-1, 1]
    }))
    # Ensure no NaN or Inf in z_matrix
    z_matrix[is.na(z_matrix)] <- -1
    z_matrix[is.infinite(z_matrix)] <- 1
  }

  # Orientation: samples are ROWS and features are COLUMNS. Studies normally
  # have more samples than features, so this keeps the long dimension vertical
  # where there is room for it, and keeps sample names horizontally readable.
  z_matrix <- t(z_matrix)

  # Feature grouping (Met-1 / Pathway) splits COLUMNS.
  fs <- if (!is.null(feature_split))
          factor(feature_split, levels = unique(feature_split)) else NULL

  # Sample grouping (Treatment / Group) splits ROWS.
  gs <- if (!all(is.na(heatmap_labs)))
          factor(heatmap_labs, levels = unique(heatmap_labs)) else NULL

  # Annotations, splits, cell sizing and fonts all come from one place, shared
  # with contributionHm(), so the two heatmaps of a study are the same panel
  # and only the cell fill differs between them.
  layout <- .hm_layout(n_row = nrow(z_matrix), n_col = ncol(z_matrix),
                       gs = gs, fs = fs,
                       group_split_name = group_split_name,
                       feature_split_name = feature_split_name,
                       group_palette = group_palette,
                       feature_palette = feature_palette,
                       cell_size = cell_size, max_aspect = max_aspect,
                       cell_border = cell_border, fontsize = fontsize)

  # Diverging ramp anchored at 0, since the matrix is a z-score clipped to
  # [-1, 1] and the midpoint is meaningful.
  .cols   <- quantPalettes(palette)
  .col_fn <- circlize::colorRamp2(
    seq(-1, 1, length.out = length(.cols)), .cols)

  hm <- do.call(ComplexHeatmap::Heatmap, c(list(
    z_matrix, name = "Z-score", col = .col_fn), layout))

  return(hm)
}


# The parts of a heatmap that are not the cells: group and feature colour
# bars, the splits they imply, cell sizing and label fonts.
#
# Extracted so quantileHm() and contributionHm() cannot drift apart. They
# are two readings of one matrix, so a reader has to be able to put them side
# by side and trust that a colour bar means the same thing in both.
#
# Returns the argument list to splice into ComplexHeatmap::Heatmap().
# Wrap a long legend label.
#
# ComplexHeatmap draws legend labels at full width, and an ion-library class
# such as "Membrane component || Lipid-mediated signalling" is 45 characters.
# A few of those widen the whole drawing past the device-width cap in the
# report template, where the excess is CROPPED rather than scaled -- which is
# how a legend ends up sliced off at the right edge. Breaking at the "||"
# separator first keeps the two halves of a compound name on their own lines;
# anything still too long falls back to ordinary word wrapping.
.wrap_label <- function(x, width = 28) {
  vapply(as.character(x), function(s) {
    if (is.na(s) || nchar(s) <= width) return(s)
    s2 <- gsub(" \\|\\| ", " ||\n", s)
    if (max(nchar(strsplit(s2, "\n", fixed = TRUE)[[1]])) > width)
      s2 <- paste(strwrap(s, width = width), collapse = "\n")
    s2
  }, character(1), USE.NAMES = FALSE)
}

.hm_layout <- function(n_row, n_col, gs, fs,
                       group_split_name = "Group",
                       feature_split_name = "Pathway",
                       group_palette = "hat", feature_palette = "reading",
                       cell_size = 8, max_aspect = 1.5,
                       cell_border = "white", fontsize = 8) {

  # One font size for every label on the panel, so the sample names, the
  # rotated feature names and the annotation titles stay in proportion.
  .lab_gp <- grid::gpar(fontsize = fontsize)

  # A hairline in the background colour between cells. Adjacent samples often
  # land in the same clipped z-score, and without a border they merge into one
  # block so the number of sections is no longer readable off the panel.
  .rect_gp <- if (length(cell_border) != 1L || is.na(cell_border))
                grid::gpar(col = NA)
              else grid::gpar(col = cell_border, lwd = 0.5)
  # Same treatment for the annotation bars, so they read as one tile per
  # sample / feature rather than a continuous stripe.
  .anno_gp <- .rect_gp

  # Coloured left annotation for the sample groups.
  left_anno <- NULL
  if (!is.null(gs)) {
    grp_cols <- .anno_cols(levels(gs), group_palette)
    .col_args <- list(); .col_args[[group_split_name]] <- grp_cols
    .anno_args <- list(); .anno_args[[group_split_name]] <- gs
    # Show the group colour key: the row slice titles are suppressed
    # (row_title = NULL), so without this legend the group colours are
    # unlabelled. Drawing several annotation legends together can trip a
    # ComplexHeatmap viewport bug ("depth applied to NULL"); the report draws
    # with merge_legends = TRUE, which packs the legends and avoids it.
    left_anno <- do.call(ComplexHeatmap::rowAnnotation,
      c(.anno_args, list(col = .col_args, show_legend = TRUE,
                          gp = .anno_gp,
                          show_annotation_name = TRUE,
                          annotation_name_gp = .lab_gp,
                          annotation_name_side = "bottom",
                          annotation_legend_param = list(
                            title_gp = grid::gpar(fontsize = fontsize,
                                                   fontface = "bold"),
                            at = levels(gs),
                            labels = .wrap_label(levels(gs)),
                            labels_gp = .lab_gp))))
  }

  # Coloured top annotation for the feature groups.
  top_anno <- NULL
  if (!is.null(fs)) {
    fs_cols <- .anno_cols(levels(fs), feature_palette)
    .col_args <- list(); .col_args[[feature_split_name]] <- fs_cols
    .anno_args <- list(); .anno_args[[feature_split_name]] <- fs
    top_anno <- do.call(ComplexHeatmap::HeatmapAnnotation,
      c(.anno_args, list(col = .col_args, gp = .anno_gp,
                          show_annotation_name = TRUE,
                          annotation_name_gp = .lab_gp,
                          annotation_legend_param = list(
                            title_gp = grid::gpar(fontsize = fontsize,
                                                   fontface = "bold"),
                            at = levels(fs),
                            labels = .wrap_label(levels(fs)),
                            labels_gp = .lab_gp),
                          annotation_name_side = "right")))
  }

  # Absolute cell sizing. Without width/height ComplexHeatmap stretches the
  # body to fill the device, which makes cells oblong whenever the two
  # dimensions differ. Sizing the body in millimetres makes cells square; when
  # one dimension has more than four times the entries of the other, cells on
  # the short dimension may stretch up to `max_aspect` so the body is not a
  # sliver.
  .size <- NULL
  if (!is.null(cell_size) && !all(is.na(cell_size))) {
    cell_w <- cell_size
    cell_h <- cell_size
    if (n_col > 0 && n_row > 0) {
      if (n_row / n_col > 4)      cell_w <- cell_size * max_aspect
      else if (n_col / n_row > 4) cell_h <- cell_size * max_aspect
    }
    .size <- list(width  = grid::unit(n_col * cell_w, "mm"),
                  height = grid::unit(n_row * cell_h, "mm"))
  }

  # row_title is suppressed: sample groups are labelled by the left colour
  # bar, so the group name is not repeated down the side.
  c(list(
    rect_gp = .rect_gp,
    row_names_gp = .lab_gp, column_names_gp = .lab_gp,
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
      labels_gp = .lab_gp),
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE,
    row_split = gs, row_title = NULL, cluster_row_slices = FALSE,
    column_split = fs, column_title = NULL, cluster_column_slices = FALSE,
    left_annotation = left_anno,
    top_annotation = top_anno), .size)
}
