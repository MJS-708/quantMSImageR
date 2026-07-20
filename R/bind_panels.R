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
#' This replaces the older `bind_polarities()` (kept as a deprecated alias)
#' and the per-acquisition `union_pad_features()` approach.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param obj1,obj2 `MSImagingExperiment` (or `quant_MSImagingExperiment`)
#'   objects of the same physical tissue. Both must carry `x` and `y` in
#'   `pData(.)`.
#' @param label Character. Optional. If supplied, written to
#'   `pData(result)$run` so the merged object combines cleanly downstream.
#'
#' @return A `quant_MSImagingExperiment` with:
#'   \itemize{
#'     \item Features = obj1's features followed by obj2's features
#'           (feature names duplicated in both are kept once, with obj1's
#'           values).
#'     \item Pixels = intersection of (x, y) grids.
#'     \item Each pixel carries real intensities from both inputs.
#'   }
#'
#' @seealso [combine_MSIs()], [generate_txt_images()]
#'
#' @examples
#' p1 <- system.file("extdata", "example.raw", "section01.RDS",
#'                   package = "quantMSImageR")
#' p2 <- system.file("extdata", "example.raw", "section02.RDS",
#'                   package = "quantMSImageR")
#' merged <- bind_panels(readRDS(p1), readRDS(p2), label = "A")
#'
#' @export
bind_panels <- function(obj1, obj2, label = NULL) {

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
    message(sprintf(
      "  bind_panels: dropping %d duplicate feature(s) from obj2: %s",
      length(overlap), paste(overlap, collapse = ", ")))
    obj2c <- obj2c[!(nms2 %in% overlap), ]
    nms2  <- as.character(fData(obj2c)$name)
  }

  # 3. rbind intensity matrices
  idata <- rbind(
    as.matrix(spectraData(obj1c)[["intensity"]]),
    as.matrix(spectraData(obj2c)[["intensity"]])
  )

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
    spectraData = idata,
    featureData = fd_mdf,
    pixelData   = pdata
  )
  featureNames(out) <- combined_fd$name
  as(out, "quant_MSImagingExperiment")
}
