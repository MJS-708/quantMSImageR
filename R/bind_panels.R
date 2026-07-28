#' Merge two MSI objects of the same tissue by coordinate-matched rbind
#'
#' Combines two `MSImagingExperiment` objects acquired on the *same physical
#' sample* into one whose feature set spans both inputs. Pixels are matched
#' by `(x, y)` coordinates rather than by row order, so the two acquisitions
#' can have different pixel counts (which is the rule, not the exception:
#' different MRM panels or polarities almost always sample at slightly
#' different rates and therefore produce different grids).
#'
#' Pixels present in only one input are dropped -- a small percentage of edge
#' coverage in the more densely sampled acquisition typically. Pixels present
#' in both keep real intensities from both inputs, so cross-panel /
#' cross-polarity colocalisation is meaningful at the combined object.
#'
#' The common use is pairing a positive- and a negative-mode acquisition of the
#' same section, but any two panels of the same physical area work.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param obj1,obj2 `MSImagingExperiment` (or `quant_MSImagingExperiment`)
#'   objects of the same physical tissue. Both must carry `x` and `y` in
#'   `pData(.)`.
#' @param label Character. Optional. If supplied, written to
#'   `pData(result)$run` so the merged object combines cleanly downstream.
#' @param feature_match Character. How a feature name appearing in both inputs
#'   is treated: `"transition"` (default) keeps `obj1`'s copy only after
#'   confirming both refer to the same precursor and product m/z, and errors
#'   otherwise; `"name"` keeps `obj1`'s copy on the strength of the name alone.
#' @param mz_tolerance Numeric. Half-width in Da within which two precursor or
#'   product m/z values count as the same (default `0.4`, i.e. nominal mass --
#'   MRM selects Q1 and Q3 at unit resolution). See [read_mrm()].
#'
#' @return A `quant_MSImagingExperiment` with:
#'   \itemize{
#'     \item Features = obj1's features followed by obj2's features
#'           (feature names duplicated in both are kept once, with obj1's
#'           values).
#'     \item Pixels = intersection of (x, y) grids.
#'     \item Each pixel carries real intensities from both inputs.
#'     \item Every spectra layer both inputs carry -- `response`, `snr` and
#'           calibrated amounts as well as `intensity` -- together with
#'           `obj1`'s experiment metadata and calibration/tissue slots.
#'   }
#'
#' @section Processed inputs:
#' Both inputs must carry the same set of spectra layers, since a layer present
#' on only one side cannot be filled in for the other half of the features.
#' Bind panels at the same stage of processing: either both raw, or both after
#' the same steps.
#'
#' @seealso [combine_MSIs()], [generate_txt_images()]
#'
#' @examples
#' # Both inputs must be the SAME physical area measured twice, so the example
#' # splits one section into two disjoint "panels" rather than using two
#' # different sections -- binding unrelated sections by coordinate would
#' # silently pair pixels that are not the same piece of tissue.
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- readRDS(p)
#' panel_a <- obj[1:4, ]
#' panel_b <- obj[5:nrow(fData(obj)), ]
#'
#' merged <- bind_panels(panel_a, panel_b, label = "SampleA_1")
#' nrow(fData(merged))   # features from both panels
#'
#' @family combining acquisitions
#' @export
bind_panels <- function(obj1, obj2, label = NULL,
                        feature_match = c("transition", "name"),
                        mz_tolerance = 0.4) {

  feature_match <- match.arg(feature_match)

  pd1 <- pData(obj1); pd2 <- pData(obj2)
  if (!all(c("x", "y") %in% names(pd1)) ||
      !all(c("x", "y") %in% names(pd2)))
    stop("bind_panels: both inputs must carry x and y columns in pData().")

  # 1. Coordinate-match pixels
  k1 <- paste(pd1$x, pd1$y, sep = "_")
  k2 <- paste(pd2$x, pd2$y, sep = "_")
  common <- intersect(k1, k2)
  if (length(common) == 0)
    stop("bind_panels: no shared (x, y) coordinates between the two objects.")

  n1 <- length(k1); n2 <- length(k2); nC <- length(common)
  if (nC < n1 || nC < n2)
    message(sprintf(
      "  bind_panels: keeping %d pixels (intersection); dropping %d from obj1 and %d from obj2.",
      nC, n1 - nC, n2 - nC))

  obj1c <- obj1[, match(common, k1)]
  obj2c <- obj2[, match(common, k2)]

  # 2. Drop duplicate feature names from obj2 (obj1 wins)
  nms1 <- as.character(fData(obj1c)$name)
  nms2 <- as.character(fData(obj2c)$name)
  overlap <- intersect(nms1, nms2)
  if (length(overlap) > 0) {
    # obj2's copy is discarded, so it is worth knowing the two are the same
    # measurement before choosing one of them. A shared name across two panels
    # can easily be two different transitions.
    if (feature_match == "transition") {
      .num <- function(x) suppressWarnings(as.numeric(
        vapply(strsplit(as.character(x), " || ", fixed = TRUE), `[`,
               character(1), 1)))
      f1 <- fData(obj1c); f2 <- fData(obj2c)
      need <- c("precursor_mz", "product_mz")
      if (!all(need %in% names(f1)) || !all(need %in% names(f2)))
        stop("bind_panels: feature_match = \"transition\" needs precursor_mz ",
             "and product_mz in fData() of both objects. Pass ",
             "feature_match = \"name\" to drop duplicates on the name alone.",
             call. = FALSE)

      i1 <- match(overlap, nms1); i2 <- match(overlap, nms2)
      bad <- which(abs(.num(f1$precursor_mz)[i1] - .num(f2$precursor_mz)[i2]) > mz_tolerance |
                   abs(.num(f1$product_mz)[i1]   - .num(f2$product_mz)[i2])   > mz_tolerance)
      if (length(bad))
        stop("bind_panels: ", length(bad), " feature name(s) appear in both ",
             "panels over different transitions:
",
             paste(sprintf("  %s: %s -> %s vs %s -> %s", overlap[bad],
                           .num(f1$precursor_mz)[i1][bad], .num(f1$product_mz)[i1][bad],
                           .num(f2$precursor_mz)[i2][bad], .num(f2$product_mz)[i2][bad]),
                   collapse = "
"),
             "
Keeping one would discard a different measurement. Rename them ",
             "in the ion library, or pass feature_match = \"name\" if the names ",
             "really are authoritative.", call. = FALSE)
    }

    message(sprintf(
      "  bind_panels: dropping %d duplicate feature(s) from obj2: %s",
      length(overlap), paste(overlap, collapse = ", ")))
    obj2c <- obj2c[!(nms2 %in% overlap), ]
    nms2  <- as.character(fData(obj2c)$name)
  }

  # 3. rbind every spectra layer, not only intensity: binding a pair of
  #    processed panels used to silently return raw data, because response,
  #    snr and any calibrated layers were dropped here.
  lyr1 <- names(spectraData(obj1c)); lyr2 <- names(spectraData(obj2c))
  if (!identical(sort(lyr1), sort(lyr2)))
    stop("bind_panels: the two objects carry different spectra layers (",
         paste(sort(lyr1), collapse = ", "), " vs ",
         paste(sort(lyr2), collapse = ", "),
         "). A layer on only one side cannot be filled in for the other half ",
         "of the features -- process both panels the same way before binding.",
         call. = FALSE)

  sdata <- lapply(stats::setNames(lyr1, lyr1), function(nm)
    rbind(as.matrix(spectra(obj1c, nm)), as.matrix(spectra(obj2c, nm))))

  # 4. Combined fData with unique sequential mz keys (Cardinal needs unique mz)
  fd1 <- as.data.frame(fData(obj1c))
  fd2 <- as.data.frame(fData(obj2c))

  # Align columns (one side may have columns the other lacks; fill with NA)
  all_cols <- union(colnames(fd1), colnames(fd2))
  for (col in setdiff(all_cols, colnames(fd1))) fd1[[col]] <- NA
  for (col in setdiff(all_cols, colnames(fd2))) fd2[[col]] <- NA
  fd1 <- fd1[, all_cols, drop = FALSE]
  fd2 <- fd2[, all_cols, drop = FALSE]

  combined_fd     <- rbind(fd1, fd2)
  combined_fd$mz  <- seq_len(nrow(combined_fd))

  extra_cols <- setdiff(colnames(combined_fd), "mz")
  fd_mdf <- do.call(
    Cardinal::MassDataFrame,
    c(list(mz = combined_fd$mz),
      as.list(combined_fd[, extra_cols, drop = FALSE]))
  )

  # 5. Pixel metadata from obj1c (we trust obj1 since obj2 was matched to it)
  pdata <- pixelData(obj1c)
  if (!is.null(label))
    pdata$run <- factor(rep(label, ncol(obj1c)))

  # 6. Assemble and return
  out <- MSImagingExperiment(
    spectraData    = sdata,
    featureData    = fd_mdf,
    pixelData      = pdata,
    experimentData = experimentData(obj1c)
  )
  featureNames(out) <- combined_fd$name
  out <- as(out, "quant_MSImagingExperiment")

  # Coercion starts these empty. obj1 supplied the pixel grid and wins on
  # duplicate features, so its metadata is the one that still describes the
  # result; losing a calibration silently would be worse than not carrying it.
  if (is(obj1c, "quant_MSImagingExperiment")) {
    out@calibrationInfo <- obj1c@calibrationInfo
    out@tissueInfo      <- obj1c@tissueInfo
  }
  out
}
