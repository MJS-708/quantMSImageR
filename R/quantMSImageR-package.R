#' quantMSImageR: quantification of DESI-MRM mass spectrometry imaging data
#'
#' Tools (extending Cardinal) for processing and quantifying targeted DESI-MRM
#' mass spectrometry imaging data: signal-to-noise filtering, tissue/background
#' separation, per-feature ion images, quantile heatmaps and standards-based
#' absolute quantification.
#'
#' @keywords internal
#' @importFrom methods as is new
#' @importFrom stats median na.omit quantile sd setNames
#' @importFrom utils read.csv read.table write.csv write.table
#' @importFrom grDevices dev.new rainbow
#' @importFrom grid unit
#' @importFrom dplyr across all_of any_of arrange bind_rows distinct filter
#' @importFrom dplyr group_by left_join mutate pull right_join row_number select
#' @importFrom dplyr summarise ungroup
#' @importFrom ggplot2 aes discrete_scale element_blank element_line element_rect
#' @importFrom ggplot2 element_text facet_grid geom_point geom_tile ggplot labs
#' @importFrom ggplot2 rel theme theme_minimal
#' @importFrom viridis scale_fill_viridis
#' @importFrom chemCal inverse.predict
#' @importFrom tibble tibble column_to_rownames rownames_to_column
#' @importFrom tidyr pivot_wider
#' @importFrom ProtGenerics spectraData spectraData<-
"_PACKAGE"

# Non-standard-evaluation variables (dplyr/subset/aes) so R CMD check does not
# flag them as undefined globals.
utils::globalVariables(c(
  "Cell_type", "Polarity", "ROI", "Type", "all_noise", "amount_pg",
  "collision_eV", "cone_V", "duplicate_mrms", "lipid", "pg_perpixel",
  "pixel_ind", "precursor_mz", "product_mz", "response", "sample_name",
  "transformed_x", "transformed_y", "transition_id_int",
  "transition_id_name", "x", "x_loci", "y", "y_loci"
))
