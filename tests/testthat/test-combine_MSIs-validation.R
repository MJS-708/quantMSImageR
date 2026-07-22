# combine_MSIs() cbinds acquisitions. Cardinal's own failure for mismatched
# features names neither object, and duplicate run identifiers are not caught at
# all -- two sections then merge into one sample everywhere downstream.

mk <- function(run = "s1", names = c("f1", "f2"), layers = character(0)) {
  n <- length(names)
  fdata <- MassDataFrame(mz = seq_len(n), name = names)
  pdata <- PositionDataFrame(run = rep(run, 4),
                             coord = expand.grid(x = 1:2, y = 1:2))
  obj <- as(MSImagingExperiment(spectraData = matrix(seq_len(n * 4), nrow = n),
                                 featureData = fdata, pixelData = pdata),
            "quant_MSImagingExperiment")
  for (l in layers) spectra(obj, l) <- matrix(1, nrow = n, ncol = 4)
  obj
}

test_that("acquisitions with the same features and layers combine", {
  out <- combine_MSIs(mk("s1"), mk("s2"))
  expect_equal(ncol(out), 8L)
  expect_setequal(unique(as.character(pData(out)$run)), c("s1", "s2"))
})

test_that("mismatched features name the offending object", {
  expect_error(combine_MSIs(mk("s1"), mk("s2", names = c("f1", "f9"))),
               regexp = "object 2 has different features")
})

test_that("mismatched spectra layers are caught before cbind", {
  expect_error(combine_MSIs(mk("s1"), mk("s2", layers = "snr")),
               regexp = "spectra layers")
})

test_that("duplicate run identifiers are rejected", {
  # Silently merging two sections into one sample is worse than refusing.
  expect_error(combine_MSIs(mk("s1"), mk("s1")),
               regexp = "run identifiers must be unique")
})

test_that("more than two objects still combine", {
  out <- combine_MSIs(mk("s1"), mk("s2"), mk("s3"))
  expect_equal(length(unique(as.character(pData(out)$run))), 3L)
})
