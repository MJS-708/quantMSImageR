test_that(".roll_median is centred and truncates at the ends", {
  rm_ <- quantMSImageR:::.roll_median
  v <- c(1, 2, 3, 100, 5, 6, 7)

  # k = 3, centred: ends use the shorter available window.
  expect_equal(rm_(v, 3),
               c(median(v[1:2]), median(v[1:3]), median(v[2:4]), median(v[3:5]),
                 median(v[4:6]), median(v[5:7]), median(v[6:7])))

  # A single spike is suppressed by the rolling median but not by no smoothing.
  expect_lt(rm_(v, 5)[4], v[4])

  # Window at least as long as the vector collapses to one median.
  expect_equal(rm_(v, 99), rep(median(v), length(v)))
  expect_equal(rm_(numeric(0), 3), numeric(0))
})

make_obj <- function() {
  # Two lines of five pixels, deliberately supplied in a scrambled pixel order
  # so that any group-order/pixel-order mismatch shows up.
  g <- expand.grid(x = 1:5, y = 1:2)
  ord <- c(7, 1, 9, 3, 5, 2, 8, 4, 10, 6)
  g <- g[ord, ]

  fdata <- MassDataFrame(mz = 1:2, name = c("feat", "std"),
                         analyte = c("Analyte", "IS"))
  pdata <- PositionDataFrame(run = rep("s1", nrow(g)), coord = g)
  ints  <- rbind(rep(10, nrow(g)), rep(2, nrow(g)))
  as(MSImagingExperiment(spectraData = ints, featureData = fdata,
                         pixelData = pdata), "quant_MSImagingExperiment")
}

test_that("window mode runs and respects pixel order", {
  obj <- make_obj()
  out <- int2response(obj, val_slot = "intensity", IS_name = "IS",
                      mode = "window", window = 3, remove_IS = TRUE)

  r <- as.numeric(spectra(out, "response"))
  # Constant data: every pixel normalises to 10 / 2 regardless of ordering.
  expect_true(all(abs(r - 5) < 1e-9))
})

test_that("line mode assigns by pixel index, not group order", {
  # Regression test: `response` used to be built by appending per-line results
  # and then assigned in pixel order, which misaligned whenever pixels were not
  # already sorted by y.
  obj <- make_obj()

  # Give each line a different standard level so misalignment would show.
  ints <- spectraData(obj)[["intensity"]]
  y <- pData(obj)$y
  ints[2, y == 2] <- 4          # standard is 2 on line 1, 4 on line 2
  spectraData(obj)[["intensity"]] <- ints

  out <- int2response(obj, val_slot = "intensity", IS_name = "IS",
                      mode = "line", remove_IS = TRUE)
  r <- as.numeric(spectra(out, "response"))

  expect_true(all(abs(r[y == 1] - 10 / 2) < 1e-9))
  expect_true(all(abs(r[y == 2] - 10 / 4) < 1e-9))
})

test_that("an invalid mode is rejected", {
  expect_error(int2response(make_obj(), IS_name = "IS", mode = "rolling"),
               "should be one of")
})
