# Regression weights for calibration fits. The original implementation replaced
# a zero amount with 1e-9 and weighted by 1/amount, which handed the blank a
# weight of ~1e9 -- one point then determined the whole curve.

test_that("a zero level does not dominate the weighting", {
  x <- c(0, 1, 10, 100)

  w <- quantMSImageR:::.cal_weights(x, "1/x")
  expect_true(all(is.finite(w)))
  # The blank gets the lowest positive level's weight, not a billion.
  expect_equal(w[1], 1 / 100)
  expect_equal(w[2:4], 1 / x[2:4])
  expect_lt(max(w) / min(w), 1e3)

  w2 <- quantMSImageR:::.cal_weights(x, "1/x2")
  expect_equal(w2[1], 1 / 100^2)
  expect_equal(w2[2:4], 1 / x[2:4]^2)

  expect_equal(quantMSImageR:::.cal_weights(x, "none"), rep(1, 4))
})

test_that("negative amounts are handled like the blank", {
  # Standard addition shifts the amount axis, which can put a level below zero.
  w <- quantMSImageR:::.cal_weights(c(-5, 2, 20), "1/x")
  expect_true(all(is.finite(w)))
  expect_equal(w[1], 1 / 20)
})

test_that("the blank no longer decides the fit", {
  # True line is y = 10x + 5, measured exactly at the four standards. The blank
  # is contaminated and reads 40 instead of ~5 -- a real and common failure.
  # The standards-only fit is the reference: adding a blank should nudge the
  # line, not redraw it.
  d   <- data.frame(pg_perpixel       = c(0, 1, 2, 5, 10),
                    response_perpixel = c(40, 15, 25, 55, 105))
  std <- d[-1, ]

  fit <- function(dat, w) stats::lm(response_perpixel ~ pg_perpixel,
                                    data = dat, weights = w)

  ref <- fit(std, quantMSImageR:::.cal_weights(std$pg_perpixel, "1/x"))
  new <- fit(d,   quantMSImageR:::.cal_weights(d$pg_perpixel, "1/x"))
  old <- fit(d,   1 / ifelse(d$pg_perpixel == 0, 1e-9, d$pg_perpixel))

  expect_equal(unname(coef(ref)[2]), 10)          # reference recovers truth

  # Old behaviour: the blank's weight is ~1e9 against ~1, so the line is
  # dragged through it exactly and the slope collapses from 10 to about 2.
  expect_equal(unname(coef(old)[1]), 40, tolerance = 1e-6)
  expect_lt(unname(coef(old)[2]), 4)

  # New behaviour: the blank counts once. The slope stays within 10% of the
  # standards-only fit instead of being decided by one point.
  expect_equal(unname(coef(new)[2]), unname(coef(ref)[2]), tolerance = 0.1)
  expect_lt(abs(coef(new)[2] - coef(ref)[2]),
            abs(coef(old)[2] - coef(ref)[2]))

  # The mechanism, stated directly: share of total weight held by the blank.
  w_new <- quantMSImageR:::.cal_weights(d$pg_perpixel, "1/x")
  w_old <- 1 / ifelse(d$pg_perpixel == 0, 1e-9, d$pg_perpixel)
  expect_lt(w_new[1] / sum(w_new), 0.1)
  expect_gt(w_old[1] / sum(w_old), 0.99)
})

test_that("createCalCurve exposes weighting and defaults to 1/x", {
  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
            "quant_MSImagingExperiment")
  cm  <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
  cal <- summariseCalLevels(cal, cm, val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")

  w1 <- createCalCurve(cal, cal_type = "cal")
  w0 <- createCalCurve(cal, cal_type = "cal", weighting = "none")

  expect_s4_class(w1, "quant_MSImagingExperiment")
  expect_equal(length(calibrationModels(w1)), length(calibrationModels(w0)))
  # Weighting has to actually reach the fit.
  expect_false(isTRUE(all.equal(coef(calibrationModels(w1)[[1]]),
                                coef(calibrationModels(w0)[[1]]))))

  expect_error(createCalCurve(cal, cal_type = "cal", weighting = "1/sqrt(x)"))
})
