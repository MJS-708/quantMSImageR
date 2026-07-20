#' Merge positive- and negative-mode acquisitions from the same tissue area
#'
#' A thin convenience wrapper around [bind_panels()] for the common
#' positive/negative-mode pairing. `bind_panels()` matches pixels on `(x, y)`
#' rather than requiring identical pixel grids, which matters because pos/neg
#' pairs nearly always have slightly different pixel counts due to acquisition
#' timing.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param pos_obj A positive-mode `MSImagingExperiment`.
#' @param neg_obj A negative-mode `MSImagingExperiment` of the same tissue.
#' @param label Character. Optional. Written to `pData(result)$run`.
#'
#' @return Equivalent to `bind_panels(pos_obj, neg_obj, label)`.
#'
#' @seealso [bind_panels()]
#'
#' @examples
#' p1 <- system.file("extdata", "example.raw", "section01.RDS",
#'                   package = "quantMSImageR")
#' p2 <- system.file("extdata", "example.raw", "section02.RDS",
#'                   package = "quantMSImageR")
#' # convenience wrapper for bind_panels()
#' merged <- bind_polarities(readRDS(p1), readRDS(p2), label = "A")
#'
#' @export
bind_polarities <- function(pos_obj, neg_obj, label = NULL) {
  bind_panels(pos_obj, neg_obj, label = label)
}
