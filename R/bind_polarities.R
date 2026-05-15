#' Merge positive- and negative-mode acquisitions from the same tissue area
#'
#' Deprecated. Forwards to [bind_panels()], which generalises the operation
#' to any two acquisitions of the same tissue (different polarities, different
#' MRM panels) by matching pixels on `(x, y)` rather than requiring identical
#' pixel grids. The strict equality check the old implementation enforced
#' rarely holds in practice — pos/neg pairs nearly always have slightly
#' different pixel counts due to acquisition timing.
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
#' @export
bind_polarities <- function(pos_obj, neg_obj, label = NULL) {
  .Deprecated("bind_panels",
              msg = paste0("bind_polarities() is deprecated; use bind_panels(). ",
                            "bind_panels() handles pos+neg pairs with non-identical ",
                            "pixel grids by coordinate matching."))
  bind_panels(pos_obj, neg_obj, label = label)
}
