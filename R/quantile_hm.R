#' Quantile heatmap of MSI features across samples
#'
#' For each feature, computes the nth quantile of pixel intensities within each
#' acquisition run, then z-scores the resulting per-feature profile across samples
#' (clipped to \[-1, 1\]) and renders a `ComplexHeatmap::Heatmap`. Rows are
#' features; columns are samples in the order given by `heatmap_order`.
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
#' @param row_split Optional character/factor vector of length equal to the
#'   number of features, used to split heatmap rows into labelled groups.
#'   Typically derived from an ion library metadata column (e.g. `Met-1`).
#'   Defaults to `NULL` (no row splitting).
#' @param row_split_name Display name for the row-split annotation legend
#'   (e.g. `"Met-1"`, `"Pathway"`). Defaults to `"Pathway"`.
#' @param column_split_name Display name for the column-split annotation
#'   legend (e.g. `"Treatment"`, `"Group"`). Defaults to `"Group"`.
#'
#' @return A `ComplexHeatmap::Heatmap` object (rows = features, columns = samples).
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' hm <- quantile_hm(obj, quant_val = 0.5,
#'                   heatmap_order = "section01", heatmap_labs = "A")
#'
#' @export

quantile_hm = function(MSIobject, quant_val, heatmap_order = NA, heatmap_labs = NA,
                        row_split = NULL,
                        row_split_name = "Pathway",
                        column_split_name = "Group") {
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

  # Row split: preserve declared order, no re-clustering across slices
  rs <- if (!is.null(row_split)) factor(row_split, levels = unique(row_split)) else NULL

  # Coloured top annotation for column groups (Treatment / sample group)
  top_anno <- NULL
  cs <- NULL
  if (!all(is.na(heatmap_labs))) {
    cs <- factor(heatmap_labs, levels = unique(heatmap_labs))
    grp_levels <- levels(cs)
    grp_cols <- setNames(
      grDevices::hcl.colors(max(length(grp_levels), 2), palette = "Dark 3")[seq_along(grp_levels)],
      grp_levels
    )
    .col_args <- list(); .col_args[[column_split_name]] <- grp_cols
    .anno_args <- list(); .anno_args[[column_split_name]] <- cs
    top_anno <- do.call(ComplexHeatmap::HeatmapAnnotation,
      c(.anno_args, list(col = .col_args, show_annotation_name = TRUE,
                          annotation_name_side = "right")))
  }

  # Coloured right annotation for row groups (Met-1 / Pathway)
  right_anno <- NULL
  if (!is.null(rs)) {
    rs_levels <- levels(rs)
    rs_cols <- setNames(
      grDevices::hcl.colors(max(length(rs_levels), 2), palette = "Set 2")[seq_along(rs_levels)],
      rs_levels
    )
    .col_args <- list(); .col_args[[row_split_name]] <- rs_cols
    .anno_args <- list(); .anno_args[[row_split_name]] <- rs
    right_anno <- do.call(ComplexHeatmap::rowAnnotation,
      c(.anno_args, list(col = .col_args, show_annotation_name = TRUE,
                          annotation_name_side = "bottom")))
  }

  # Build the heatmap with optional split / colour-bar annotations.
  # - Column titles default to levels(cs) so each Group block is labelled in
  #   black text above its colour bar (e.g. "Ctrl_M", "HDM_F").
  # - row_title is suppressed: rows are grouped via the right colour bar only,
  #   so verbose Met-1 names don't crowd the LHS.
  hm <- ComplexHeatmap::Heatmap(z_matrix, name = "Z-score",
    cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = FALSE,
    column_split = cs,
    column_title_gp = grid::gpar(col = "black", fontsize = 11),
    row_split = rs, row_title = NULL, cluster_row_slices = FALSE,
    top_annotation = top_anno,
    right_annotation = right_anno)

  return(hm)
}
