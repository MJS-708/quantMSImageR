#' Align two MSI objects to their common features
#'
#' Cardinal's \code{cbind} (used by \code{\link{combine_MSIs}}) requires the
#' \code{mz} key column to match exactly across objects.  When acquisitions
#' come from different instrument methods or ion-library versions, features
#' with the same display name can carry different \code{mz} values, causing
#' \code{cbind} to fail with a "non-matching key columns" error.
#'
#' \code{align_features} resolves this by:
#' \enumerate{
#'   \item Subsetting both objects to the \emph{intersection} of feature names.
#'   \item Reordering \code{obj2} to match \code{obj1}'s feature order.
#'   \item Forcing \code{obj2}'s \code{mz} values to equal \code{obj1}'s, which
#'     permits Cardinal's key check to pass after features have been matched by
#'     name.  \code{\link{combine_MSIs}} then restores \code{fData} from
#'     \code{obj1}, so all feature metadata originates from the first object.
#' }
#'
#' Matching is by display name alone.  Users must ensure that identically named
#' features represent the same analytical transition in both objects: two
#' methods can reuse a name while differing in precursor ion, product ion,
#' adduct, polarity, collision energy or transition definition, and rewriting
#' the \code{mz} key would then silently merge different measurements.  Check
#' \code{fData()} (precursor and product m/z in particular) before relying on
#' the result for anything quantitative.
#'
#' The function is called automatically inside
#' \code{\link{generate_txt_images}} before every cross-sample
#' \code{\link{combine_MSIs}}, so mixed-method studies (e.g. 1 Hz + 2 Hz
#' acquisitions) are handled transparently.  It is also exported so that it
#' can be used directly in custom scripts or Rmd reports whenever two
#' \code{quant_MSImagingExperiment} objects need to be combined.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param obj1 A \code{quant_MSImagingExperiment}.  Its feature set and
#'   \code{mz} values are used as the reference.
#' @param obj2 A \code{quant_MSImagingExperiment} to align against
#'   \code{obj1}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{obj1}}{Subset of \code{obj1} containing only common
#'       features.}
#'     \item{\code{obj2}}{Subset of \code{obj2} reordered and \code{mz}-keyed
#'       to match \code{obj1}, ready for \code{\link{combine_MSIs}}.}
#'   }
#'
#' @seealso \code{\link{combine_MSIs}}, \code{\link{generate_txt_images}}
#'
#' @examples
#' p1 <- system.file("extdata", "example.raw", "section01.RDS",
#'                   package = "quantMSImageR")
#' p2 <- system.file("extdata", "example.raw", "section02.RDS",
#'                   package = "quantMSImageR")
#' al <- align_features(readRDS(p1), readRDS(p2))
#'
#' @family combining acquisitions
#' @export
align_features <- function(obj1, obj2) {
  nms1   <- fData(obj1)$name
  nms2   <- fData(obj2)$name
  common <- intersect(nms1, nms2)

  if (length(common) == 0)
    stop("No common features between the two MSI objects -- ",
         "check that both acquisitions used a compatible ion library.")

  dropped1 <- setdiff(nms1, common)
  dropped2 <- setdiff(nms2, common)
  if (length(dropped1) > 0)
    message("  align_features: dropping from obj1 (", length(dropped1), "): ",
            paste(dropped1, collapse = ", "))
  if (length(dropped2) > 0)
    message("  align_features: dropping from obj2 (", length(dropped2), "): ",
            paste(dropped2, collapse = ", "))

  obj1 <- obj1[which(nms1 %in% common), ]
  obj2 <- obj2[which(nms2 %in% common), ]

  # Reorder obj2 so features are in the same order as obj1
  ord  <- match(fData(obj1)$name, fData(obj2)$name)
  obj2 <- obj2[ord, ]

  # Force obj2's mz key to exactly match obj1's so Cardinal cbind succeeds
  mz(obj2) <- mz(obj1)

  list(obj1 = obj1, obj2 = obj2)
}
