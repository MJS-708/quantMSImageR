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

# ---- Placement by stage position ------------------------------------------
# Every acquisition numbers its own pixels from (1, 1); readMRM() also keeps
# where each one was on the stage (mm). A piece is a small synthetic grid at a
# given stage origin, with one feature whose value encodes the pixel's stage
# position, so a test can tell whether data travelled with its pixel.

piece <- function(nx = 3, ny = 2, x0 = 0, y0 = 0, pitch = 0.038, run = "p",
                  roi_id = NULL) {
  g  <- expand.grid(x = seq_len(nx), y = seq_len(ny))
  xs <- x0 + (g$x - 1) * pitch
  ys <- y0 + (g$y - 1) * pitch
  pd <- PositionDataFrame(run = factor(rep(run, nrow(g))), coord = g,
                          x_stage = xs, y_stage = ys)
  if (!is.null(roi_id)) {
    pd$roi_label <- sub("_[0-9]+$", "", roi_id)
    pd$roi_id    <- roi_id
  }
  MSImagingExperiment(
    spectraData = matrix(round(xs * 1000) * 1e4 + round(ys * 1000), nrow = 1),
    featureData = MassDataFrame(mz = 1, name = "f1"),
    pixelData   = pd)
}
grid_of <- function(o) sort(paste(pData(o)$x, pData(o)$y, sep = "_"))

test_that("pieces are placed where the stage recorded them, gaps left empty", {
  p <- 0.038
  top <- piece(y0 = 0,     run = "top")
  btm <- piece(y0 = 5 * p, run = "btm")          # 3 empty rows between them
  w <- suppressMessages(stitchAcquisitions(top, btm, label = "whole"))

  expect_equal(ncol(w), ncol(top) + ncol(btm))
  expect_setequal(unique(pData(w)$y), c(1, 2, 6, 7))
  expect_setequal(unique(pData(w)$x), 1:3)
  expect_equal(as.character(unique(pData(w)$run)), "whole")
})

test_that("the order the pieces are given in does not matter", {
  p  <- 0.038
  a  <- piece(y0 = 0,     run = "a")
  b  <- piece(y0 = 5 * p, run = "b")
  ab <- suppressMessages(stitchAcquisitions(a, b, label = "w"))
  ba <- suppressMessages(stitchAcquisitions(b, a, label = "w"))

  expect_equal(grid_of(ab), grid_of(ba))
  # Each grid position holds the same measurement either way.
  val <- function(o) {
    v <- as.numeric(spectra(o, "intensity"))
    v[order(paste(pData(o)$x, pData(o)$y, sep = "_"))]
  }
  expect_equal(val(ab), val(ba))
})

test_that("pieces beside each other work too", {
  p <- 0.038
  w <- suppressMessages(stitchAcquisitions(
    piece(x0 = 0, run = "l"), piece(x0 = 3 * p, run = "r"), label = "w"))
  expect_setequal(unique(pData(w)$x), 1:6)
  expect_setequal(unique(pData(w)$y), 1:2)
})

test_that("a fractional offset moves to the nearest grid position and is reported", {
  p <- 0.038
  a <- piece(x0 = 0,         run = "a")
  b <- piece(x0 = 3.3 * p,   run = "b")          # 0.3 px off the grid of a
  expect_message(w <- stitchAcquisitions(a, b, label = "w"),
                 regexp = "largest shift onto it 11.4", fixed = TRUE)
  expect_setequal(unique(pData(w)$x), 1:6)
})

test_that("pieces that overlap on the stage are refused", {
  p <- 0.038
  expect_error(suppressMessages(stitchAcquisitions(
    piece(y0 = 0, run = "a"), piece(y0 = 1 * p, run = "b"), label = "w")),
    regexp = "appear in more than one piece")
})

test_that("pieces of different pixel size are refused", {
  expect_error(suppressMessages(stitchAcquisitions(
    piece(pitch = 0.038, run = "a"),
    piece(pitch = 0.050, y0 = 1, run = "b"), label = "w")),
    regexp = "pixel")
})

test_that("region numbers continue across pieces", {
  # Each piece numbers its own regions from _01; after stitching, the top
  # piece's airway_01 and the bottom piece's airway_01 are two airways.
  p   <- 0.038
  top <- piece(y0 = 0,     run = "top", roi_id = rep("airway_01", 6))
  btm <- piece(y0 = 5 * p, run = "btm", roi_id = rep("airway_01", 6))
  w <- suppressMessages(stitchAcquisitions(top, btm, label = "w"))
  expect_setequal(unique(pData(w)$roi_id), c("airway_01", "airway_02"))
  expect_equal(unique(pData(w)$roi_label), "airway")
})

test_that("without stage positions on every piece, none are placed by them", {
  p <- 0.038
  a <- piece(y0 = 0, run = "a")
  b <- piece(y0 = 5 * p, run = "b")
  pData(b)$x_stage <- NA_real_
  # Both then start at (1, 1), which is an overlap; the message says why.
  expect_error(suppressMessages(stitchAcquisitions(a, b, label = "w")),
               regexp = "readMRM(overwrite = TRUE)", fixed = TRUE)
})
