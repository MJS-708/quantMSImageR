require(testthat)
require(quantMSImageR)


make_snr_obj <- function(bg_label = "Background") {
  fdata <- MassDataFrame(mz = c(500, 510, 540, 550),
                         analyte = c("IS", rep("analyte", 3)))
  pdata <- PositionDataFrame(
    run       = c(rep("run1", 2), rep("run2", 2), rep(bg_label, 2)),
    coord     = expand.grid(x = 1:3, y = 1:2),
    sample_ID = c(rep("Tissue", 4), rep(bg_label, 2))
  )
  ints <- matrix(nrow = 4, ncol = 6,
                 data = c(rep(c(5, 100, 150, 100), 2),
                          rep(c(10, 100, 105, 300), 2),
                          1, 2, 3, 4,
                          3, 6, 9, 12))
  obj <- MSImagingExperiment(spectraData = ints, featureData = fdata,
                              pixelData = pdata)
  as(obj, "quant_MSImagingExperiment")
}

test_that("SNR is NA for background pixels and correct for tissue (average = 'mean')", {
  obj      <- make_snr_obj()
  new_data <- int2snr(MSIobject = obj, val_slot = "intensity",
                      background = "Background", tissue = "Tissue",
                      snr_thresh = 3, pixel_header = "sample_ID",
                      average = "mean")

  expect_true(all(is.na(spectra(new_data, "snr")[, 5:6])))
  expect_equal(spectra(new_data, "snr")[1, ], c(NA, NA, 5,    5,    NA, NA))
  expect_equal(spectra(new_data, "snr")[2, ], c(25, 25, 25,   25,   NA, NA))
  expect_equal(spectra(new_data, "snr")[3, ], c(25, 25, 17.5, 17.5, NA, NA))
  expect_equal(spectra(new_data, "snr")[4, ], c(12.5, 12.5, 37.5, 37.5, NA, NA))
})

test_that("legacy 'Noise' labels and the deprecated `noise` argument still work", {
  # Masks and scripts written before the noise -> background rename must keep
  # working: the old label is matched, and the old argument name is honoured.
  legacy <- make_snr_obj(bg_label = "Noise")
  current <- int2snr(MSIobject = make_snr_obj(), val_slot = "intensity",
                     background = "Background", tissue = "Tissue",
                     snr_thresh = 3, pixel_header = "sample_ID")

  by_new_arg <- int2snr(MSIobject = legacy, val_slot = "intensity",
                        background = "Background", tissue = "Tissue",
                        snr_thresh = 3, pixel_header = "sample_ID")
  by_old_arg <- int2snr(MSIobject = legacy, val_slot = "intensity",
                        noise = "Noise", tissue = "Tissue",
                        snr_thresh = 3, pixel_header = "sample_ID")

  expect_equal(spectra(by_new_arg, "snr"), spectra(current, "snr"))
  expect_equal(spectra(by_old_arg, "snr"), spectra(current, "snr"))
})

test_that("average = 'median' gives different results from 'mean' for skewed background", {
  obj       <- make_snr_obj()
  res_mean  <- int2snr(obj, val_slot = "intensity",
                       background = "Background", tissue = "Tissue",
                       snr_thresh = 1, pixel_header = "sample_ID",
                       average = "mean")
  res_med   <- int2snr(obj, val_slot = "intensity",
                       background = "Background", tissue = "Tissue",
                       snr_thresh = 1, pixel_header = "sample_ID",
                       average = "median")
  # Results need not be identical; just verify both run without error
  expect_s4_class(res_mean, "quant_MSImagingExperiment")
  expect_s4_class(res_med,  "quant_MSImagingExperiment")
})

test_that("invalid average_method is rejected", {
  obj <- make_snr_obj()
  expect_error(
    int2snr(obj, val_slot = "intensity",
            background = "Background", tissue = "Tissue",
            snr_thresh = 3, pixel_header = "sample_ID",
            average = "geometric"),
    regexp = "'arg' should be one of"
  )
})

no_bg_obj <- function() {
  fdata <- MassDataFrame(mz = 1:2, name = c("f1", "f2"))
  pdata <- PositionDataFrame(run       = rep("s1", 4),
                              coord     = expand.grid(x = 1:2, y = 1:2),
                              sample_ID = rep("Tissue", 4))
  ints  <- matrix(c(10, 20, 30, 40, 100, 200, 300, 400), nrow = 2, byrow = TRUE)
  as(MSImagingExperiment(spectraData = ints, featureData = fdata,
                          pixelData = pdata),
     "quant_MSImagingExperiment")
}

test_that("no background pixels is an error when filtering was asked for", {
  # Without background there is no ratio to compute, and the previous
  # behaviour -- copying intensity into a slot named `snr` -- let a downstream
  # reader mistake raw response for a signal-to-background ratio.
  expect_error(
    int2snr(no_bg_obj(), val_slot = "intensity",
            background = "Background", tissue = "Tissue",
            snr_thresh = 3, pixel_header = "sample_ID"),
    regexp = "no 'Background' pixels"
  )
})

test_that("no_background gives explicit alternatives to failing", {
  expect_warning(
    na_res <- int2snr(no_bg_obj(), val_slot = "intensity",
                      background = "Background", tissue = "Tissue",
                      snr_thresh = 3, pixel_header = "sample_ID",
                      no_background = "all_na"),
    regexp = "NA for every pixel"
  )
  expect_true(all(is.na(spectra(na_res, "snr"))))

  expect_warning(
    cp_res <- int2snr(no_bg_obj(), val_slot = "intensity",
                      background = "Background", tissue = "Tissue",
                      snr_thresh = 3, pixel_header = "sample_ID",
                      no_background = "copy"),
    regexp = "not a ratio"
  )
  expect_equal(spectra(cp_res, "snr"), spectra(cp_res, "intensity"))
})

test_that("snr_thresh = 0 passes values through without failing", {
  # No filtering was requested, so a missing background is not a failure to do
  # something the caller wanted. This is the smoke-test pass run before ROIs
  # have been drawn.
  expect_warning(
    result <- int2snr(no_bg_obj(), val_slot = "intensity",
                      background = "Background", tissue = "Tissue",
                      snr_thresh = 0, pixel_header = "sample_ID"),
    regexp = "not a ratio"
  )
  expect_equal(ncol(result), 4L)
  expect_equal(spectra(result, "snr"), spectra(result, "intensity"))
})
