test_that("createMSIDatamatrix summarises ROIs from the synthetic cal data", {

  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  skip_if_not(file.exists(file.path(cal_dir, "cal_MSI.RDS")),
              "synthetic cal data not generated (run inst/scripts/generate_cal_data.R)")

  cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
            "quant_MSImagingExperiment")

  out <- createMSIDatamatrix(cal, val_slot = "intensity",
                             roi_header = "identifier")

  ram <- out@tissueInfo@roi_average_matrix   # one row per calibration spot
  apm <- out@tissueInfo@all_pixel_matrix      # one row per calibration pixel

  expect_equal(nrow(ram), 15L)    # 5 levels x 3 reps
  expect_equal(nrow(apm), 60L)    # 15 spots x 4 pixels
  # Take the expected columns from the object itself rather than hard-coding
  # analyte names, which change whenever generate_cal_data.R is re-tuned.
  expect_true(all(fData(cal)$name %in% colnames(ram)))
  expect_length(fData(cal)$name, 3L)
})
