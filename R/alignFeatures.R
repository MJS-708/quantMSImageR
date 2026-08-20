# Reorder an object's features and give it a new mz key.
#
# The obvious `obj[ord, ]` only works while the permutation leaves mz ascending:
# Cardinal's MassDataFrame requires a sorted mz key, so two acquisitions that
# list the same features in a different order could not be aligned at all. The
# fast path is kept, with a rebuild behind it for the permutations it rejects.
.reorder_features <- function(obj, ord, new_mz) {

  out <- tryCatch({
    o <- obj[ord, ]
    mz(o) <- new_mz
    o
  }, error = function(e) NULL)
  if (!is.null(out)) return(out)

  # Pulled layer by layer: as.list() on the SpectraArrays of an already-subset
  # object fails, and obj2 always arrives here subset to the common features.
  .lyr <- names(spectraData(obj))
  sp <- lapply(stats::setNames(.lyr, .lyr),
               function(nm) as.matrix(spectra(obj, nm))[ord, , drop = FALSE])
  fd <- as.data.frame(fData(obj))[ord, , drop = FALSE]
  fd$mz <- new_mz

  out <- MSImagingExperiment(spectraData    = sp,
                             featureData    = do.call(MassDataFrame, as.list(fd)),
                             pixelData      = pixelData(obj),
                             experimentData = experimentData(obj))
  out <- as(out, "quant_MSImagingExperiment")

  # Coercion starts these empty, and losing a calibration silently would be
  # worse than the ordering problem this is working around.
  if (is(obj, "quant_MSImagingExperiment")) {
    out@calibrationInfo <- obj@calibrationInfo
    out@tissueInfo      <- obj@tissueInfo
  }
  featureNames(out) <- fData(out)$name
  out
}

#' Align two MSI objects to their common features
#'
#' Cardinal's \code{cbind} (used by \code{\link{combineMSIs}}) requires the
#' \code{mz} key column to match exactly across objects.  When acquisitions
#' come from different instrument methods or ion-library versions, features
#' with the same display name can carry different \code{mz} values, causing
#' \code{cbind} to fail with a "non-matching key columns" error.
#'
#' \code{alignFeatures} resolves this by:
#' \enumerate{
#'   \item Subsetting both objects to the \emph{intersection} of feature names.
#'   \item Reordering \code{obj2} to match \code{obj1}'s feature order.
#'   \item Forcing \code{obj2}'s \code{mz} values to equal \code{obj1}'s, which
#'     permits Cardinal's key check to pass after features have been matched by
#'     name.  \code{\link{combineMSIs}} then restores \code{fData} from
#'     \code{obj1}, so all feature metadata originates from the first object.
#' }
#'
#' Features are paired by display name, and by default that pairing is then
#' \emph{verified} against the transition each name refers to: precursor and
#' product m/z must agree to within \code{mz_tolerance}.  Two methods can reuse
#' a name while differing in precursor ion, product ion or transition
#' definition, and because this function rewrites the \code{mz} key, an
#' unverified pairing would silently merge different measurements into one
#' feature.  A disagreement is therefore an error naming the features involved.
#'
#' \code{feature_match = "name"} restores name-only matching for objects that
#' carry no precursor/product metadata.  It is not the default: it cannot
#' detect the failure above, and this function is called automatically inside
#' \code{\link{generateTxtImages}}.
#'
#' The function is called automatically inside
#' \code{\link{generateTxtImages}} before every cross-sample
#' \code{\link{combineMSIs}}, so mixed-method studies (e.g. 1 Hz + 2 Hz
#' acquisitions) are handled transparently.  It is also exported so that it
#' can be used directly in custom scripts or Rmd reports whenever two
#' \code{quant_MSImagingExperiment} objects need to be combined.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param obj1 An \code{MSImagingExperiment} -- a
#'   \code{quant_MSImagingExperiment} is one. Its feature set and \code{mz}
#'   values are used as the reference, and its \code{fData()} must carry a
#'   \code{name} column.
#' @param obj2 An \code{MSImagingExperiment} to align against \code{obj1},
#'   likewise carrying \code{fData()$name}.
#' @param feature_match Character. \code{"transition"} (default) verifies that
#'   same-named features share a precursor and product m/z; \code{"name"}
#'   matches on the display name alone, without checking what it refers to.
#' @param mz_tolerance Numeric. Half-width in Da within which two precursor or
#'   product m/z values are taken to be the same (default \code{0.4}, i.e.
#'   nominal mass). MRM selects Q1 and Q3 at unit resolution, so two panels can
#'   record the same channel as 308.17 and 308.2; at a tighter tolerance this
#'   check would reject them as different transitions. See \link{readMRM}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{obj1}}{Subset of \code{obj1} containing only common
#'       features.}
#'     \item{\code{obj2}}{Subset of \code{obj2} reordered and \code{mz}-keyed
#'       to match \code{obj1}, ready for \code{\link{combineMSIs}}.}
#'   }
#'
#' @seealso \code{\link{combineMSIs}}, \code{\link{generateTxtImages}}
#'
#' @examples
#' p1 <- system.file("extdata", "example.raw", "section01.RDS",
#'                   package = "quantMSImageR")
#' p2 <- system.file("extdata", "example.raw", "section02.RDS",
#'                   package = "quantMSImageR")
#' al <- alignFeatures(readRDS(p1), readRDS(p2))
#'
#' @family combining acquisitions
#' @export
alignFeatures <- function(obj1, obj2,
                           feature_match = c("transition", "name"),
                           mz_tolerance = 0.4) {
  feature_match <- match.arg(feature_match)

  # MSImagingExperiment, not quant_MSImagingExperiment: quant_ contains it, so
  # this accepts both, and the objects shipped in extdata -- including the ones
  # this function's own example reads -- are plain MSImagingExperiments.
  # Demanding the subclass would fail R CMD check on that example.
  if (!is(obj1, "MSImagingExperiment") || !is(obj2, "MSImagingExperiment"))
    stop("alignFeatures: obj1 and obj2 must both be MSImagingExperiment ",
         "objects (quant_MSImagingExperiment included); got ",
         class(obj1)[1L], " and ", class(obj2)[1L], ".", call. = FALSE)

  nms1   <- fData(obj1)$name
  nms2   <- fData(obj2)$name

  # Without this, a missing name column makes intersect(NULL, NULL) empty and
  # the error below blames the ion library for a structural problem.
  if (is.null(nms1) || is.null(nms2))
    stop("alignFeatures: both objects need a 'name' column in fData(); ",
         "features are paired by name. Objects read by readMRM() carry one.",
         call. = FALSE)

  common <- intersect(nms1, nms2)

  if (length(common) == 0)
    stop("No common features between the two MSI objects -- ",
         "check that both acquisitions used a compatible ion library.")

  dropped1 <- setdiff(nms1, common)
  dropped2 <- setdiff(nms2, common)
  if (length(dropped1) > 0)
    message("  alignFeatures: dropping from obj1 (", length(dropped1), "): ",
            paste(dropped1, collapse = ", "))
  if (length(dropped2) > 0)
    message("  alignFeatures: dropping from obj2 (", length(dropped2), "): ",
            paste(dropped2, collapse = ", "))

  obj1 <- obj1[which(nms1 %in% common), ]
  obj2 <- obj2[which(nms2 %in% common), ]

  # Reorder obj2 so features are in the same order as obj1, and take obj1's mz
  # key, which is what lets Cardinal's cbind accept the pair.
  ord  <- match(fData(obj1)$name, fData(obj2)$name)
  obj2 <- .reorder_features(obj2, ord, mz(obj1))

  # Confirm the paired names really are the same transition before the mz key
  # is rewritten, since rewriting it is what makes a mismatch unrecoverable.
  if (feature_match == "transition") {
    .num <- function(x) suppressWarnings(as.numeric(
      vapply(strsplit(as.character(x), " || ", fixed = TRUE), `[`,
             character(1), 1)))
    f1 <- fData(obj1); f2 <- fData(obj2)
    need <- c("precursor_mz", "product_mz")
    if (!all(need %in% names(f1)) || !all(need %in% names(f2)))
      stop("alignFeatures: feature_match = \"transition\" needs precursor_mz ",
           "and product_mz in fData() of both objects. Objects read by ",
           "readMRM() carry them. Pass feature_match = \"name\" to match on ",
           "the display name alone, accepting that identical names are then ",
           "assumed to be the same transition.", call. = FALSE)

    bad <- which(abs(.num(f1$precursor_mz) - .num(f2$precursor_mz)) > mz_tolerance |
                 abs(.num(f1$product_mz)   - .num(f2$product_mz))   > mz_tolerance)
    if (length(bad))
      stop("alignFeatures: ", length(bad), " feature(s) share a name but not ",
           "a transition:\n",
           paste(sprintf("  %s: %s -> %s vs %s -> %s",
                         as.character(f1$name)[bad],
                         .num(f1$precursor_mz)[bad], .num(f1$product_mz)[bad],
                         .num(f2$precursor_mz)[bad], .num(f2$product_mz)[bad]),
                 collapse = "\n"),
           "\nCombining these would merge different measurements into one ",
           "feature. Reconcile the ion libraries, or pass ",
           "feature_match = \"name\" if the names really are authoritative.",
           call. = FALSE)
  }

  list(obj1 = obj1, obj2 = obj2)
}
