require(testthat)
require(quantMSImageR)

# One pixel per sample, which is all the centring arithmetic needs: the
# correction is per feature and per block, not spatial.
.fake_msi <- function(m) {
  fdata <- MassDataFrame(mz = seq_len(nrow(m)),
                         name = paste0("feat_", seq_len(nrow(m))))
  pdata <- PositionDataFrame(run = paste0("s", seq_len(ncol(m))),
                             coord = data.frame(x = seq_len(ncol(m)), y = 1L))
  obj <- MSImagingExperiment(spectraData = m, featureData = fdata,
                             pixelData = pdata)
  spectraData(obj)[["intensity"]] <- m
  obj
}

test_that("centring removes a block's gain but keeps its response", {
  # Two blocks, two groups. Block B is measured at 4x the gain of block A, and
  # within each block the treated samples are genuinely 2x the control.
  base   <- c(100, 400)                       # two features
  design <- expand.grid(grp = c("Ctrl", "Trt"), blk = c("A", "B"),
                        stringsAsFactors = FALSE)
  gain   <- c(A = 1, B = 4)
  resp   <- c(Ctrl = 1, Trt = 2)
  m <- vapply(seq_len(nrow(design)), function(i)
    base * gain[[design$blk[i]]] * resp[[design$grp[i]]], numeric(2))

  obj <- .fake_msi(m)
  out <- blockCentre(obj, block = design$blk, mode = "center", verbose = FALSE)
  v <- as.matrix(spectraData(out)[["intensity"]])

  # The 4x block gain is gone: the two blocks now read alike.
  expect_equal(v[, design$blk == "A"], v[, design$blk == "B"],
               ignore_attr = TRUE)
  # And the 2x treatment response survives in both.
  for (b in c("A", "B")) {
    j <- which(design$blk == b)
    expect_equal(v[, j[design$grp[j] == "Trt"]] /
                   v[, j[design$grp[j] == "Ctrl"]],
                 c(2, 2), ignore_attr = TRUE)
  }
})

test_that("baseline centring reads as a fold change against the block's own control", {
  design <- expand.grid(grp = c("Ctrl", "Trt"), blk = c("A", "B"),
                        stringsAsFactors = FALSE)
  gain <- c(A = 1, B = 7)
  resp <- c(Ctrl = 1, Trt = 3)
  m <- vapply(seq_len(nrow(design)), function(i)
    c(100, 400) * gain[[design$blk[i]]] * resp[[design$grp[i]]], numeric(2))

  out <- blockCentre(.fake_msi(m), block = design$blk, mode = "baseline",
                     groups = design$grp, baseline = "Ctrl", verbose = FALSE)
  v <- as.matrix(spectraData(out)[["intensity"]])

  # Every control is exactly 1: it is its own reference.
  expect_equal(as.numeric(v[, design$grp == "Ctrl"]), rep(1, 4))
  # Every treated sample is the response, whatever its block's gain was.
  expect_equal(as.numeric(v[, design$grp == "Trt"]), rep(3, 4))
})

test_that("a block without the baseline falls back to its own mean", {
  design <- data.frame(grp = c("Ctrl", "Trt", "Trt"), blk = c("A", "A", "B"),
                       stringsAsFactors = FALSE)
  m <- matrix(c(100, 200, 300), nrow = 1)

  expect_message(
    out <- blockCentre(.fake_msi(m), block = design$blk, mode = "baseline",
                       groups = design$grp, baseline = "Ctrl"),
    "B had no 'Ctrl' sample")
  v <- as.numeric(spectraData(out)[["intensity"]])
  expect_equal(v[1:2], c(1, 2))     # block A against its control
  expect_equal(v[3], 1)             # block B against its own mean
})

test_that("zero and missing pixels do not pull the reference", {
  m <- matrix(c(100, 0, NA, 100), nrow = 1)
  out <- blockCentre(.fake_msi(m), block = rep("A", 4), mode = "center",
                     verbose = FALSE)
  v <- as.numeric(spectraData(out)[["intensity"]])
  # Reference is the geometric mean of the two real values, so both become 1.
  expect_equal(v[c(1, 4)], c(1, 1))
  # The others are carried through untouched rather than invented.
  expect_equal(v[2], 0)
  expect_true(is.na(v[3]))
})

test_that("off is a no-op and the arguments are checked", {
  m <- matrix(c(1, 2, 3, 4), nrow = 2)
  obj <- .fake_msi(m)
  expect_identical(blockCentre(obj, block = c("A", "B"), mode = "off"), obj)

  expect_error(blockCentre(obj, block = "A", mode = "center"),
               "block has 1 entries")
  expect_error(blockCentre(obj, block = c("A", "B"), mode = "baseline"),
               "needs groups and baseline")
})
