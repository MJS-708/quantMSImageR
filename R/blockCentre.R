#' Remove a blocking factor's offset from an MSI object
#'
#' Divides every pixel of a feature by a value computed within its own block --
#' the patient, donor, sex or slide the sample belongs to -- so that samples
#' from different blocks can be compared without the block's own level getting
#' in the way.
#'
#' The correction is multiplicative on intensity, which is additive in log
#' space: a block-level gain (one patient's sections ionising better, one
#' slide run on a different day) is removed, while genuine differences in how
#' far a treatment moves a feature survive. That survival is the point -- a
#' patient who responds twice as strongly as another should still look twice as
#' strong afterwards.
#'
#' Two references are available. `"center"` divides by the block's own
#' geometric mean, which needs no privileged sample and is steadier because it
#' pools every sample in the block. `"baseline"` divides by the geometric mean
#' of the block's `baseline` group, so each value reads directly as a fold
#' change against its own control; it is more interpretable and noisier, since
#' all of the reference's noise propagates into every other sample of that
#' block. A block with no `baseline` sample falls back to `"center"` for that
#' block alone, and says so.
#'
#' @param MSIobject A `MSImagingExperiment` (or `quant_MSImagingExperiment`).
#' @param block Character/factor vector, one entry per PIXEL, naming the block
#'   that pixel's sample belongs to.
#' @param mode `"center"` (default), `"baseline"` or `"off"`. `"off"` returns
#'   the object unchanged.
#' @param groups Character vector, one entry per pixel, giving each pixel's
#'   group. Required when `mode = "baseline"`.
#' @param baseline The value of `groups` that is the reference. Required when
#'   `mode = "baseline"`.
#' @param val_slot Name of the `spectraData` slot to correct (default
#'   `"intensity"`).
#' @param verbose Report what was done (default `TRUE`).
#'
#' @return `MSIobject` with the chosen slot divided through by its per-feature,
#'   per-block reference. Non-positive and missing values are left as they are:
#'   they carry no ratio.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- readRDS(p)
#'
#' # Two blocks, taken across the pixels so the example needs nothing else.
#' blk <- rep(c("A", "B"), length.out = ncol(obj))
#' centred <- blockCentre(obj, block = blk, mode = "center")
#'
#' # The reference is a geometric mean, so that is what lands on 1: each
#' # feature of each block is now expressed relative to its own level.
#' ints <- as.matrix(spectraData(centred)[["intensity"]])
#' geo  <- function(x) exp(mean(log(x[is.finite(x) & x > 0])))
#' round(apply(ints[, blk == "A"], 1L, geo), 2)
#'
#' @family workflow
#' @export
blockCentre <- function(MSIobject, block, mode = c("center", "baseline", "off"),
                        groups = NULL, baseline = NULL,
                        val_slot = "intensity", verbose = TRUE) {
  mode <- match.arg(mode)
  if (identical(mode, "off")) return(MSIobject)

  mat <- spectraData(MSIobject)[[val_slot]]
  if (is.null(mat))
    stop("blockCentre: no '", val_slot, "' slot in this object.", call. = FALSE)
  mat <- as.matrix(mat)

  block <- as.character(block)
  if (length(block) != ncol(mat))
    stop("blockCentre: block has ", length(block), " entries but there are ",
         ncol(mat), " pixels.", call. = FALSE)

  if (identical(mode, "baseline")) {
    if (is.null(groups) || is.null(baseline))
      stop("blockCentre: mode = \"baseline\" needs groups and baseline.",
           call. = FALSE)
    groups <- as.character(groups)
    if (length(groups) != ncol(mat))
      stop("blockCentre: groups has ", length(groups), " entries but there ",
           "are ", ncol(mat), " pixels.", call. = FALSE)
  }

  # Only positive values carry a ratio: a zero or a missing pixel has no
  # information about the block's level and must not drag the reference down.
  lmat <- mat
  lmat[!is.finite(lmat) | lmat <= 0] <- NA_real_
  lmat <- log2(lmat)

  lv <- unique(block[!is.na(block)])
  fell_back <- character(0)
  for (b in lv) {
    cols <- which(block == b)
    ref_cols <- cols
    if (identical(mode, "baseline")) {
      ref_cols <- cols[groups[cols] == baseline]
      if (!length(ref_cols)) {
        # Centring on a reference the block does not contain is undefined, and
        # dropping the block silently would be worse than using its own mean.
        ref_cols <- cols
        fell_back <- c(fell_back, b)
      }
    }
    # Geometric mean per feature, as an offset in log space.
    off <- rowMeans(lmat[, ref_cols, drop = FALSE], na.rm = TRUE)
    off[!is.finite(off)] <- 0          # nothing to correct by
    mat[, cols] <- mat[, cols] / 2^off
  }

  spectraData(MSIobject)[[val_slot]] <- mat

  if (verbose) {
    message(sprintf(
      "blockCentre: '%s' divided within %d block%s by its %s.",
      val_slot, length(lv), if (length(lv) == 1L) "" else "s",
      if (identical(mode, "baseline")) sprintf("'%s' group", baseline)
      else "own geometric mean"))
    if (length(fell_back))
      message(sprintf(
        "  %s had no '%s' sample, so %s centred on its own mean instead.",
        paste(fell_back, collapse = ", "), baseline,
        if (length(fell_back) == 1L) "it was" else "they were"))
  }

  MSIobject
}
