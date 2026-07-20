test_that("create_cal_curve fits a linear model per lipid", {

  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  skip_if_not(file.exists(file.path(cal_dir, "cal_MSI.RDS")),
              "synthetic cal data not generated (run inst/generate_cal_data.R)")

  cal  <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
             "quant_MSImagingExperiment")
  meta <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))

  cal <- summarise_cal_levels(cal, meta, val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")
  cal <- create_cal_curve(cal, cal_type = "cal")

  cl <- cal@calibrationInfo@cal_list
  r2 <- cal@calibrationInfo@r2_df

  expect_equal(length(cl), 3L)                                  # 3 lipids
  expect_true(all(vapply(cl, function(m) inherits(m, "lm"), logical(1))))
  expect_true(all(r2$r2 > 0.9))          # synthetic curves are near-linear
})
