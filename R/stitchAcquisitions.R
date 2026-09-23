#' Stitch acquisitions that are pieces of one tissue into a single sample
#'
#' A tissue larger than the stage range is sometimes acquired in two passes --
#' top then bottom -- producing two `.raw` folders that are *one* sample. This
#' places their pixels on one grid and makes them a single run.
#'
#' It is not the same operation as either neighbour, and picking the wrong one
#' is quiet rather than loud:
#'
#' \describe{
#'   \item{[bindPanels()]}{two MRM *panels* over one tissue: different
#'     transitions, shared pixel grid. It **intersects** pixels on `(x, y)`, so
#'     given two spatial halves it keeps only coordinates present in both --
#'     nothing, if the halves do not overlap.}
#'   \item{[combineMSIs()]}{several *samples* into one object, each keeping
#'     its own run identifier. It refuses duplicate run identifiers, precisely
#'     to stop two sections silently becoming one sample.}
#'   \item{`stitchAcquisitions()`}{several *pieces of one sample* into one
#'     run. The union of pixels, one run identifier.}
#' }
#'
#' Why this matters beyond the picture: two halves of one tissue are not two
#' replicates. Left as separate samples they count twice in every group
#' summary -- the group mean, the box/violin plots, any downstream model -- and
#' because the halves agree with each other they look like a reproducible
#' effect rather than one measurement.
#'
#' @section Where the pieces go:
#'
#' Each acquisition numbers its pixels from its own first pixel, so every piece
#' starts at `(1, 1)` and the indices alone cannot say where one piece lies
#' relative to another. [readMRM()] also keeps the stage position of every
#' pixel (`x_stage`, `y_stage`, in mm), and that is what places the pieces: a
#' piece acquired below, above or beside another lands there, whatever order
#' the pieces are given in. The origin is the smallest stage position of any
#' piece.
#'
#' * **Same pixel size.** Pieces must share a pixel size on each axis (within
#'   2%); otherwise there is no common grid to put them on, and this stops.
#' * **Nearest grid position.** A piece whose origin falls between grid
#'   positions of the first -- the stage moved a fraction of a pixel -- is
#'   moved to the nearest one. The largest such shift is reported, and is at
#'   most half a pixel.
#' * **Gaps stay empty.** Stage space between two pieces is not filled in; no
#'   pixel is invented.
#' * **Overlap is an error.** Pixels from two pieces landing on the same grid
#'   position mean the acquisitions overlap, or are not pieces of one tissue.
#'   Averaging them or keeping one would invent a measurement, so this stops.
#'
#' Objects without stage positions -- cached by an earlier version of
#' [readMRM()], or built another way -- keep their own `(x, y)`, and any pixel
#' two of them share is the same error. Re-read the acquisitions with
#' `readMRM(overwrite = TRUE)` to place them by stage position instead.
#'
#' @section Regions of interest:
#'
#' Pieces labelled with [labelROIs()] each number their regions from `_01`, so
#' `airway_01` in the top piece and `airway_01` in the bottom piece are
#' different airways. They are renumbered across the stitched sample, in piece
#' order, so every region keeps an identifier of its own.
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
#'   `label` throughout and `x`, `y` on the common grid.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- readRDS(p)
#' # Split by pixel, then stitch back: one sample, all pixels, one run.
#' half1 <- obj[, 1:50]
#' half2 <- obj[, 51:ncol(obj)]
#' whole <- stitchAcquisitions(half1, half2, label = "section01")
#' ncol(whole) == ncol(obj)
#'
#' @seealso [bindPanels()], [combineMSIs()], [readMRM()]
#' @family combining acquisitions
#' @export
stitchAcquisitions <- function(..., label = NULL) {

  objects <- list(...)
  # Tolerate a single list, which is how the YAML path arrives.
  if (length(objects) == 1L && is.list(objects[[1]]) &&
      !methods::is(objects[[1]], "MSImagingExperiment"))
    objects <- objects[[1]]

  if (length(objects) < 2L)
    stop("stitchAcquisitions: needs at least two acquisitions to stitch.",
         call. = FALSE)

  .nm  <- function(o) as.character(fData(o)$name)
  .lay <- function(o) sort(names(spectraData(o)))
  for (ind in seq_along(objects)[-1]) {
    if (!identical(.nm(objects[[ind]]), .nm(objects[[1]])))
      stop("stitchAcquisitions: piece ", ind, " has different features from ",
           "piece 1. Pieces of one acquisition should share a transition ",
           "list; if these are different panels of the same tissue you want ",
           "bindPanels() instead.", call. = FALSE)
    if (!identical(.lay(objects[[ind]]), .lay(objects[[1]])))
      stop("stitchAcquisitions: piece ", ind, " has spectra layers (",
           paste(.lay(objects[[ind]]), collapse = ", "),
           ") that differ from piece 1 (",
           paste(.lay(objects[[1]]), collapse = ", "),
           "). Run the same processing steps on every piece before ",
           "stitching.", call. = FALSE)
  }

  .piece_nm <- function(ind) {
    r <- unique(as.character(pData(objects[[ind]])$run))
    if (length(r) == 1L && nzchar(r)) sprintf("'%s'", r) else paste("piece", ind)
  }

  .has_stage <- function(o) {
    pd <- pData(o)
    all(c("x_stage", "y_stage") %in% names(pd)) &&
      all(is.finite(pd$x_stage)) && all(is.finite(pd$y_stage))
  }
  by_stage <- all(vapply(objects, .has_stage, logical(1)))

  if (by_stage) {
    xy <- .stage_grid(objects, .piece_nm)
  } else {
    if (any(vapply(objects, .has_stage, logical(1))))
      message("stitchAcquisitions: only some pieces carry stage positions, ",
              "so none are placed by them. Re-read every piece with ",
              "readMRM(overwrite = TRUE).")
    xy <- lapply(objects, function(o)
      list(x = pData(o)$x, y = pData(o)$y))
  }

  # A repeated grid position means two pieces claim the same place on the
  # tissue. Averaging or silently keeping one would invent a pixel.
  keys <- lapply(xy, function(p) paste(p$x, p$y, sep = "_"))
  all_keys <- unlist(keys, use.names = FALSE)
  dup <- unique(all_keys[duplicated(all_keys)])
  if (length(dup)) {
    hint <- if (by_stage)
      "Pieces of one tissue should tile, not overlap -- check these are not the same region acquired twice."
    else
      "Without stage positions every piece starts at (1, 1); re-read the acquisitions with readMRM(overwrite = TRUE) so they are placed by where the stage recorded them."
    stop("stitchAcquisitions: ", length(dup), " pixel position(s) appear ",
         "in more than one piece (e.g. ",
         paste(utils::head(sub("_", ", ", dup), 3), collapse = "; "),
         "). ", hint, call. = FALSE)
  }

  objects <- .harmonise_pdata(objects)
  for (ind in seq_along(objects))
    objects[[ind]] <- .set_coords(objects[[ind]], xy[[ind]]$x, xy[[ind]]$y)
  objects <- .renumber_rois(objects)

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

# Grid spacing along one stage axis, in mm. Stage positions are rounded to
# 0.1 um first, so the tiny jitter in the recorded values does not count as a
# step, and steps much smaller than the typical one are ignored for the same
# reason. The estimate is then refined over the full extent, which is far more
# precise than any single step.
.stage_pitch <- function(v) {
  u <- sort(unique(round(v[is.finite(v)], 4)))
  if (length(u) < 2L) return(NA_real_)
  d <- diff(u)
  d <- d[d > 0.5 * stats::quantile(d, 0.9, names = FALSE)]
  p <- stats::median(d)
  span <- max(u) - min(u)
  span / max(1, round(span / p))
}

# Grid positions for every piece from its stage coordinates. Each piece keeps
# its own internal spacing -- positions within a piece are counted with its own
# pitch from its own origin -- and the piece as a whole is moved to the nearest
# position on the common grid. Returns one list(x, y) per piece.
.stage_grid <- function(objects, piece_nm) {
  axes <- c(x = "x_stage", y = "y_stage")

  pitch <- vapply(axes, function(col)
    vapply(objects, function(o) .stage_pitch(pData(o)[[col]]), numeric(1)),
    numeric(length(objects)))
  pitch <- matrix(pitch, nrow = length(objects),
                  dimnames = list(NULL, names(axes)))

  common <- apply(pitch, 2, function(p) {
    p <- p[is.finite(p)]
    if (!length(p)) NA_real_ else p[1]
  })
  if (anyNA(common))
    stop("stitchAcquisitions: cannot tell the pixel size along ",
         paste(names(common)[is.na(common)], collapse = " and "),
         " -- every piece is a single line of pixels on that axis.",
         call. = FALSE)

  for (ax in names(axes)) for (ind in seq_along(objects)) {
    p <- pitch[ind, ax]
    if (is.finite(p) && abs(p - common[[ax]]) / common[[ax]] > 0.02)
      stop(sprintf(paste0(
        "stitchAcquisitions: %s has a %.1f \u00b5m pixel along %s, but %s has ",
        "%.1f \u00b5m. Pieces of one tissue must be acquired at the same pixel ",
        "size to share a grid."),
        piece_nm(ind), p * 1000, ax, piece_nm(1), common[[ax]] * 1000),
        call. = FALSE)
    if (!is.finite(p)) pitch[ind, ax] <- common[[ax]]
  }

  origin <- vapply(axes, function(col)
    min(vapply(objects, function(o) min(pData(o)[[col]]), numeric(1))),
    numeric(1))

  placed <- lapply(seq_along(objects), function(ind) {
    pd  <- pData(objects[[ind]])
    res <- list()
    shift <- 0
    for (ax in names(axes)) {
      v     <- pd[[axes[[ax]]]]
      lo    <- min(v)
      local <- round((v - lo) / pitch[ind, ax])
      exact <- (lo - origin[[ax]]) / common[[ax]]
      off   <- round(exact)
      shift <- max(shift, abs(exact - off) * common[[ax]] * 1000)
      res[[ax]] <- local + off + 1
    }
    # Two pixels of one piece on one grid position means its own stage
    # positions are irregular, not that pieces overlap -- say which.
    k <- paste(res$x, res$y, sep = "_")
    if (anyDuplicated(k))
      stop("stitchAcquisitions: ", piece_nm(ind), " has ",
           sum(duplicated(k)), " pixel(s) whose stage positions fall on the ",
           "same grid position as another pixel of the same piece. Its stage ",
           "positions are not on a regular grid.", call. = FALSE)
    list(xy = res, shift = shift)
  })

  shifts   <- vapply(placed, `[[`, numeric(1), "shift")
  shift_at <- which.max(shifts)

  message(sprintf(paste0(
    "stitchAcquisitions: %d pieces placed by stage position on a %.1f x %.1f ",
    "\u00b5m grid; %s."),
    length(objects), common[["x"]] * 1000, common[["y"]] * 1000,
    if (max(shifts) < 0.05) "every piece already on it"
    else sprintf("largest shift onto it %.1f \u00b5m (%s)", max(shifts),
                 piece_nm(shift_at))))
  lapply(placed, `[[`, "xy")
}

# Region identifiers are numbered within each labelled acquisition, so after
# stitching, airway_01 from two pieces would be two different airways sharing
# a name. Renumber per label across pieces, in piece order.
.renumber_rois <- function(objects) {
  if (!all(vapply(objects, function(o)
        all(c("roi_label", "roi_id") %in% names(pData(o))), logical(1))))
    return(objects)
  counter <- list()
  for (ind in seq_along(objects)) {
    pd  <- pData(objects[[ind]])
    lab <- as.character(pd$roi_label)
    id  <- as.character(pd$roi_id)
    new <- id
    real <- !is.na(id) & !is.na(lab) & lab != "unassigned"
    for (old in unique(id[real])) {
      l <- lab[real & id == old][1]
      counter[[l]] <- if (is.null(counter[[l]])) 1L else counter[[l]] + 1L
      new[real & id == old] <- sprintf("%s_%02d", l, counter[[l]])
    }
    pData(objects[[ind]])$roi_id <- new
  }
  objects
}
