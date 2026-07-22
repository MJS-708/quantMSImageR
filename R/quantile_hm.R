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
#'   [quant_palettes()].
#' @param group_palette,feature_palette Character. Qualitative palettes for the
#'   sample-group and feature-group colour bars (defaults `"hat"` and
#'   `"reading"`). Any [grDevices::hcl.colors()] palette name also works.
#' @param row_split,row_split_name,column_split_name Deprecated aliases for
#'   `feature_split`, `feature_split_name` and `group_split_name`. They were
#'   named for the axis each grouping landed on before samples and features
#'   swapped axes; supplying them still works.
#' @param cell_size Numeric. Side of a heatmap cell, in millimetres (default
#'   `6`). Giving the body an absolute size is what makes cells square; leaving
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
#'
#' @return A `ComplexHeatmap::Heatmap` object (rows = samples, columns =
#'   features), with an absolutely-sized body unless `cell_size` is `NA`.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' hm <- quantile_hm(obj, quant_val = 0.5,
#'                   heatmap_order = "section01", heatmap_labs = "A")
#'
#' @family visualisation
#' @export

quantile_hm = function(MSIobject, quant_val, heatmap_order = NA, heatmap_labs = NA,
                        feature_split = NULL,
                        feature_split_name = "Pathway",
                        group_split_name = "Group",
                        cell_size = 6, max_aspect = 1.5,
                        cell_border = "white",
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

  # A hairline in the background colour between cells. Adjacent samples often
  # land in the same clipped z-score, and without a border they merge into one
  # block so the number of sections is no longer readable off the panel.
  .rect_gp  <- if (is.na(cell_border)) grid::gpar(col = NA)
               else grid::gpar(col = cell_border, lwd = 0.5)
  # Same treatment for the annotation bars, so they read as one tile per
  # sample / feature rather than a continuous stripe.
  .anno_gp  <- .rect_gp

  # Prepare the output matrix
  sample_names = unique(pData(MSIobject)$run)

  if(!all(is.na(heatmap_order))){
    sample_names <- factor(sample_names, levels = heatmap_order)
  }

  featurenames = fData(MSIobject)$name

  out_matrix = matrix(NA, ncol = length(sample_names), nrow = length(featurenames))

  # Loop to fill the out_matrix with quantile values
  for(col in seq_along(sample_names)){
    subsetMSI = MSIobject[, which(pData(MSIobject)$run == sample_names[col])]

    intensity_data = spectraData(subsetMSI)[["intensity"]]
    out_matrix[,col] = apply(intensity_data, 1, quantile, probs = quant_val, na.rm = TRUE)
  }

  rownames(out_matrix) = featurenames
  colnames(out_matrix) = sample_names

  # Reorder the columns of out_matrix according to the custom sample order (heatmap_order)
  # drop=FALSE preserves matrix dimensions when heatmap_order has only one element
  out_matrix <- out_matrix[, heatmap_order, drop = FALSE]

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

  # Feature grouping (Met-1 / Pathway) now splits COLUMNS.
  fs <- if (!is.null(feature_split))
          factor(feature_split, levels = unique(feature_split)) else NULL

  # Sample grouping (Treatment / Group) now splits ROWS.
  gs <- if (!all(is.na(heatmap_labs)))
          factor(heatmap_labs, levels = unique(heatmap_labs)) else NULL

  # Coloured left annotation for the sample groups.
  left_anno <- NULL
  if (!is.null(gs)) {
    grp_levels <- levels(gs)
    grp_cols <- .anno_cols(grp_levels, group_palette)
    .col_args <- list(); .col_args[[group_split_name]] <- grp_cols
    .anno_args <- list(); .anno_args[[group_split_name]] <- gs
    # show_legend = FALSE: the group colour bar is already labelled by the row
    # titles, so its legend is redundant -- and drawing it alongside the
    # feature-group + Z-score legends triggers a ComplexHeatmap legend/viewport
    # bug ("depth applied to NULL") with multiple groups.
    left_anno <- do.call(ComplexHeatmap::rowAnnotation,
      c(.anno_args, list(col = .col_args, show_legend = FALSE,
                          gp = .anno_gp,
                          show_annotation_name = TRUE,
                          annotation_name_side = "bottom")))
  }

  # Coloured top annotation for the feature groups.
  top_anno <- NULL
  if (!is.null(fs)) {
    fs_levels <- levels(fs)
    fs_cols <- .anno_cols(fs_levels, feature_palette)
    .col_args <- list(); .col_args[[feature_split_name]] <- fs_cols
    .anno_args <- list(); .anno_args[[feature_split_name]] <- fs
    top_anno <- do.call(ComplexHeatmap::HeatmapAnnotation,
      c(.anno_args, list(col = .col_args, gp = .anno_gp,
                          show_annotation_name = TRUE,
                          annotation_name_side = "right")))
  }

  # Absolute cell sizing. Without width/height ComplexHeatmap stretches the body
  # to fill the device, which makes cells oblong whenever the two dimensions
  # differ. Sizing the body in millimetres makes cells square; when one
  # dimension has more than four times the entries of the other, cells on the
  # short dimension may stretch up to `max_aspect` so the body is not a sliver.
  .size <- NULL
  if (!is.null(cell_size) && !all(is.na(cell_size))) {
    n_row <- nrow(z_matrix)   # samples
    n_col <- ncol(z_matrix)   # features

    cell_w <- cell_size
    cell_h <- cell_size
    if (n_col > 0 && n_row > 0) {
      if (n_row / n_col > 4)      cell_w <- cell_size * max_aspect
      else if (n_col / n_row > 4) cell_h <- cell_size * max_aspect
    }

    .size <- list(width  = grid::unit(n_col * cell_w, "mm"),
                  height = grid::unit(n_row * cell_h, "mm"))
  }

  # row_title is suppressed: sample groups are labelled by the left colour bar,
  # so the group name is not repeated down the side.
  # Diverging ramp anchored at 0, since the matrix is a z-score clipped to
  # [-1, 1] and the midpoint is meaningful.
  .cols   <- quant_palettes(palette)
  .col_fn <- circlize::colorRamp2(
    seq(-1, 1, length.out = length(.cols)), .cols)

  hm <- do.call(ComplexHeatmap::Heatmap, c(list(
    z_matrix, name = "Z-score", col = .col_fn, rect_gp = .rect_gp,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_dend = FALSE,
    row_split = gs, row_title = NULL, cluster_row_slices = FALSE,
    column_split = fs, column_title = NULL, cluster_column_slices = FALSE,
    left_annotation = left_anno,
    top_annotation = top_anno), .size))

  return(hm)
}
