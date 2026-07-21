test_that("summarise_cal_levels summarises calibration ROIs", {

  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  skip_if_not(file.exists(file.path(cal_dir, "cal_MSI.RDS")),
              "synthetic cal data not generated (run inst/scripts/generate_cal_data.R)")

  cal  <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
             "quant_MSImagingExperiment")
  meta <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))

  out <- summarise_cal_levels(cal, meta, val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")

  rd <- out@calibrationInfo@cal_response_data

  expect_equal(ncol(rd), 7L)
  expect_equal(nrow(rd), 45L)              # 15 spots x 3 lipids
  expect_true(all(rd$pixels == 4))         # 2x2 calibration spots
  expect_true(all(is.finite(rd$response_perpixel)))
  expect_true(all(rd$pg_perpixel > 0))
})
