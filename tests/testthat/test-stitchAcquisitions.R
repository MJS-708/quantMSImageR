# Two halves of one tissue are one sample. Getting this wrong is quiet: they
# agree with each other, so as two samples they look like a reproducible
# effect rather than one measurement.

obj <- function() {
  p <- system.file("extdata", "example.raw", "section01.RDS",
                   package = "quantMSImageR")
  readRDS(p)
}

test_that("stitching is the union of pixels under one run", {
  o <- obj()
  n <- ncol(o)
  a <- o[, 1:50]
  b <- o[, 51:n]
  w <- stitchAcquisitions(a, b, label = "whole")

  expect_equal(ncol(w), n)                       # union, not intersection
  expect_equal(nrow(fData(w)), nrow(fData(o)))   # features unchanged
  expect_equal(as.character(unique(pData(w)$run)), "whole")
  # Every original coordinate survives.
  key <- function(x) sort(paste(pData(x)$x, pData(x)$y, sep = "_"))
  expect_equal(key(w), key(o))
})

test_that("a list is accepted, as the YAML path supplies one", {
  o <- obj()
  w <- stitchAcquisitions(list(o[, 1:50], o[, 51:ncol(o)]), label = "whole")
  expect_equal(ncol(w), ncol(o))
})

test_that("overlapping pieces are refused, not silently merged", {
  # A repeated (x, y) means the pieces are not tiles of one tissue. Keeping
  # one of them, or averaging, would invent a pixel.
  o <- obj()
  expect_error(stitchAcquisitions(o[, 1:50], o[, 40:60], label = "whole"),
               regexp = "appear in more than one piece")
})

test_that("pieces must share features and layers", {
  o <- obj()
  a <- o[, 1:50]
  b <- o[2:nrow(fData(o)), 51:ncol(o)]
  expect_error(stitchAcquisitions(a, b, label = "whole"),
               regexp = "different features")
})

test_that("it refuses to stitch a single acquisition", {
  expect_error(stitchAcquisitions(obj(), label = "whole"),
               regexp = "at least two")
})

test_that("stitch differs from bindPanels on non-overlapping halves", {
  # bindPanels() intersects pixels, so on two halves that do not overlap it
  # has nothing to keep. That is the mistake stitchAcquisitions() exists to
  # prevent, and the YAML `combine:` key exists to choose between.
  o <- obj()
  a <- o[, 1:50]
  b <- o[, 51:ncol(o)]
  expect_error(bindPanels(a, b, label = "whole"),
               regexp = "no shared")
  expect_s4_class(stitchAcquisitions(a, b, label = "whole"),
                  "MSImagingExperiment")
})
