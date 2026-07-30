#' Stitch acquisitions that are pieces of one tissue into a single sample
#'
#' A tissue larger than the stage range is sometimes acquired in two passes --
#' top then bottom -- producing two `.raw` folders that are *one* sample. This
#' concatenates their pixels into a single run.
#'
#' It is not the same operation as either neighbour, and picking the wrong one
#' is quiet rather than loud:
#'
#' \describe{
#'   \item{[bind_panels()]}{two MRM *panels* over one tissue: different
#'     transitions, shared pixel grid. It **intersects** pixels on `(x, y)`, so
#'     given two spatial halves it keeps only coordinates present in both --
#'     nothing, if the halves do not overlap.}
#'   \item{[combine_MSIs()]}{several *samples* into one object, each keeping
#'     its own run identifier. It refuses duplicate run identifiers, precisely
#'     to stop two sections silently becoming one sample.}
#'   \item{`stitch_acquisitions()`}{several *pieces of one sample* into one
#'     run. The union of pixels, one run identifier.}
#' }
#'
#' Why this matters beyond the picture: two halves of one tissue are not two
#' replicates. Left as separate samples they count twice in every group
#' summary -- the group mean, the box/violin plots, any downstream model -- and
#' because the halves agree with each other they look like a reproducible
#' effect rather than one measurement.
#'
#' @section Coordinates:
#'
#' No offset is applied. Waters records absolute stage coordinates, so the
#' pieces already sit in a common frame and shifting them would move the
#' tissue apart. Pixels appearing at the same `(x, y)` in more than one piece
#' are therefore an error, not something to silently average: it means the
#' acquisitions genuinely overlap, or that they are not pieces of one tissue.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param ... Two or more `MSImagingExperiment` objects, or a single list of
#'   them. All must carry the same features and the same spectra layers.
#' @param label Character. Run identifier for the stitched sample. Defaults to
#'   the first object's run.
#'
#' @return One object holding the union of the input pixels, with `run` set to
#'   `label` throughout.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- readRDS(p)
#' # Split by pixel, then stitch back: one sample, all pixels, one run.
#' half1 <- obj[, 1:50]
#' half2 <- obj[, 51:ncol(obj)]
#' whole <- stitch_acquisitions(half1, half2, label = "section01")
#' ncol(whole) == ncol(obj)
#'
#' @seealso [bind_panels()], [combine_MSIs()]
#' @family combining acquisitions
#' @export
stitch_acquisitions <- function(..., label = NULL) {

  objects <- list(...)
  # Tolerate a single list, which is how the YAML path arrives.
  if (length(objects) == 1L && is.list(objects[[1]]) &&
      !methods::is(objects[[1]], "MSImagingExperiment"))
    objects <- objects[[1]]

  if (length(objects) < 2L)
    stop("stitch_acquisitions: needs at least two acquisitions to stitch.",
         call. = FALSE)

  .nm  <- function(o) as.character(fData(o)$name)
  .lay <- function(o) sort(names(spectraData(o)))
  for (ind in seq_along(objects)[-1]) {
    if (!identical(.nm(objects[[ind]]), .nm(objects[[1]])))
      stop("stitch_acquisitions: piece ", ind, " has different features from ",
           "piece 1. Pieces of one acquisition should share a transition ",
           "list; if these are different panels of the same tissue you want ",
           "bind_panels() instead.", call. = FALSE)
    if (!identical(.lay(objects[[ind]]), .lay(objects[[1]])))
      stop("stitch_acquisitions: piece ", ind, " has spectra layers (",
           paste(.lay(objects[[ind]]), collapse = ", "),
           ") that differ from piece 1 (",
           paste(.lay(objects[[1]]), collapse = ", "),
           "). Run the same processing steps on every piece before ",
           "stitching.", call. = FALSE)
  }

  # Absolute stage coordinates, so a repeated (x, y) means the pieces really
  # do overlap. Averaging or silently keeping one would invent a pixel.
  keys <- lapply(objects, function(o)
    paste(pData(o)$x, pData(o)$y, sep = "_"))
  all_keys <- unlist(keys, use.names = FALSE)
  dup <- unique(all_keys[duplicated(all_keys)])
  if (length(dup))
    stop("stitch_acquisitions: ", length(dup), " pixel coordinate(s) appear ",
         "in more than one piece (e.g. ",
         paste(utils::head(sub("_", ", ", dup), 3), collapse = "; "),
         "). Pieces of one tissue should tile, not overlap -- check these are ",
         "not the same region acquired twice.", call. = FALSE)

  f_data <- fData(objects[[1]])
  out <- objects[[1]]
  for (ind in seq_along(objects)[-1]) {
    out <- as(cbind(out, objects[[ind]]), "MSImagingExperiment")
    # Cardinal's cbind() concatenates the two objects' featureData columns, so
    # the accumulated fData would not match the next piece's on the following
    # iteration. The features are identical by the check above, so restore the
    # shared fData each step.
    fData(out) <- f_data
  }

  run_lab <- if (!is.null(label) && nzchar(as.character(label)))
               as.character(label)
             else as.character(pData(objects[[1]])$run)[1]
  pData(out)$run <- factor(rep(run_lab, ncol(out)))
  out
}
