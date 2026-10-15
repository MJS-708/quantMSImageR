require(testthat)
require(quantMSImageR)

# Two pieces that share a column of positions, so the overlap path is exercised
# without depending on stage metadata.
.two_pieces <- function() {
  p <- system.file("extdata", "example.raw", "section01.RDS",
                   package = "quantMSImageR")
  obj <- readRDS(p)
  n <- ncol(obj)
  list(a = obj[, 1:60], b = obj[, 51:n])      # columns 51-60 are in both
}

test_that("overlapping pieces are refused by default, and the message says how", {
  pc <- .two_pieces()
  expect_error(stitchAcquisitions(pc$a, pc$b, label = "s"),
               "appear in more than one piece")
  expect_error(stitchAcquisitions(pc$a, pc$b, label = "s"),
               'overlap = "first"')
})

test_that("overlap = first keeps the earlier piece and drops the later one's", {
  pc <- .two_pieces()
  n_union <- 60L + (ncol(pc$b) - 10L)          # the shared 10 counted once

  expect_message(
    out <- stitchAcquisitions(pc$a, pc$b, label = "s", overlap = "first"),
    "overlapping pixel\\(s\\) resolved by keeping the first piece")
  expect_equal(ncol(out), n_union)

  # Keeping the later piece instead gives the same union, not a different size:
  # it is the same positions, resolved the other way.
  expect_message(
    rev_out <- stitchAcquisitions(pc$a, pc$b, label = "s", overlap = "last"),
    "keeping the last piece")
  expect_equal(ncol(rev_out), n_union)
})

test_that("pieces that tile are untouched by the overlap setting", {
  p <- system.file("extdata", "example.raw", "section01.RDS",
                   package = "quantMSImageR")
  obj <- readRDS(p)
  a <- obj[, 1:50]
  b <- obj[, 51:ncol(obj)]

  plain <- stitchAcquisitions(a, b, label = "s")
  first <- stitchAcquisitions(a, b, label = "s", overlap = "first")
  expect_equal(ncol(plain), ncol(obj))
  expect_equal(ncol(first), ncol(plain))
})

test_that("an unknown overlap value is refused", {
  pc <- .two_pieces()
  expect_error(stitchAcquisitions(pc$a, pc$b, label = "s", overlap = "average"))
})
