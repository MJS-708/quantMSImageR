#' Quantile heatmap of MSI features across samples
#'
#' For each feature, computes the nth quantile of pixel intensities within each
#' acquisition run, then z-scores the resulting per-feature profile across samples
#' (clipped to \[-1, 1\]) and renders a `ComplexHeatmap::Heatmap`. Rows are
#' features; columns are samples in the order given by `heatmap_order`.
#'
#' @import Cardinal
#' @import ComplexHeatmap
#' @import matrixStats
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
#'
#' @return A `ComplexHeatmap::Heatmap` object (rows = features, columns = samples).
#'
#' @export

quantile_hm = function(MSIobject, quant_val, heatmap_order = NA, heatmap_labs=NA){
  # Prepare the output matrix
  sample_names = unique(pData(MSIobject)$run)

  if(!all(is.na(heatmap_order))){
    sample_names <- factor(sample_names, levels = heatmap_order)
  }

  featurenames = fData(MSIobject)$name

  out_matrix = matrix(NA, ncol = length(sample_names), nrow = length(featurenames))

  # Loop to fill the out_matrix with quantile values
  for(col in 1:length(sample_names)){
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

  # Creating the heatmap: rows = features, columns = samples
  if (!all(is.na(heatmap_labs))) {
    heatmap_labs <- factor(heatmap_labs, levels = unique(heatmap_labs))
    hm = ComplexHeatmap::Heatmap(z_matrix, name = "Z-score", cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = FALSE, column_split = heatmap_labs)
  } else {
    hm = ComplexHeatmap::Heatmap(z_matrix, name = "Z-score", cluster_rows = FALSE, cluster_columns = FALSE, show_column_dend = FALSE)
  }

  return(hm)
}
