#' Remove pure-background border rows and columns from an MSI object
#'
#' Drops any x-column or y-row in which **every** pixel is labelled as
#' `"background_pixels"` (or the historical `"noise_pixels"`) in
#' `pData(MSI_data)$sample_name`. This trims empty scan borders that arise from
#' the Waters acquisition geometry without affecting tissue or mixed-content
#' rows/columns.
#'
#' @import Cardinal
#'
#' @param MSI_data A `quant_MSImagingExperiment` object whose `pData()` contains
#'   columns `x`, `y`, and `sample_name` (populated by `makeFactor()`).
#'
#' @return A subset `quant_MSImagingExperiment` with pure-background border
#'   pixels removed.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' obj <- trimMSI(obj)
#'
#' @family acquisition
#' @export
trimMSI = function(MSI_data){
  if (!is(MSI_data, "MSImagingExperiment"))
    stop("trimMSI: MSI_data must be an MSImagingExperiment ",
         "(quant_MSImagingExperiment included); got ",
         class(MSI_data)[1L], ".", call. = FALSE)

  pd <- as.data.frame(pData(MSI_data))

  # Checked up front because the pixel columns are consumed inside dplyr
  # pipelines below, where a missing one surfaces as a bare
  # "object 'sample_name' not found" with no mention of trimMSI.
  .need <- c("x", "y", "sample_name")
  if (!all(.need %in% names(pd)))
    stop("trimMSI: pData() must carry ",
         paste(.need, collapse = ", "), "; missing ",
         paste(setdiff(.need, names(pd)), collapse = ", "), ".",
         call. = FALSE)

  # Identify x columns that are entirely noise
  bad_x <- pd |>
    dplyr::group_by(x) |>
    dplyr::summarise(all_noise = all(sample_name %in% .bg_labels("background_pixels"))) |>
    dplyr::filter(all_noise) |>
    dplyr::pull(x)

  # Identify y rows that are entirely noise
  bad_y <- pd |>
    dplyr::group_by(y) |>
    dplyr::summarise(all_noise = all(sample_name %in% .bg_labels("background_pixels"))) |>
    dplyr::filter(all_noise) |>
    dplyr::pull(y)

  # Pixels to KEEP
  keep_pixels <- which(
    !(pd$x %in% bad_x | pd$y %in% bad_y)
  )

  return(MSI_data[, keep_pixels])
}
