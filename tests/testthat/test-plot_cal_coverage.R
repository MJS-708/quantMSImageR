calibrated_obj <- function() {
  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  skip_if_not(file.exists(file.path(cal_dir, "cal_MSI.RDS")),
              "synthetic cal data not generated (run inst/scripts/generate_cal_data.R)")

  cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
            "quant_MSImagingExperiment")
  cal_metadata <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))

  cal <- summarise_cal_levels(cal, cal_metadata, val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")
  cal <- create_cal_curve(cal, cal_type = "cal")
  int2conc(cal, val_slot = "intensity", pixels = "Tissue")
}

test_that("plot_cal_coverage returns a two-panel plot for all calibrated features", {
  obj <- calibrated_obj()
  p   <- plot_cal_coverage(obj, val_slot = "intensity")

  expect_s3_class(p, "patchwork")
  expect_length(p, 2L)   # histogram + calibration curve
})

test_that("features can be selected by name or by index", {
  obj  <- calibrated_obj()
  nm   <- fData(obj)$name[1]

  by_name  <- plot_cal_coverage(obj, features = nm, val_slot = "intensity")
  by_index <- plot_cal_coverage(obj, features = 1, val_slot = "intensity")

  expect_s3_class(by_name,  "patchwork")
  expect_equal(unique(by_index[[1]]$data$feature), nm)
})

test_that("standards are taken from the model frame, matching the fitted line", {
  # The plotted standards must be the points the model was actually fitted to,
  # otherwise std_addition curves would plot points and line on different axes.
  obj <- calibrated_obj()
  nm  <- fData(obj)$name[1]
  p   <- plot_cal_coverage(obj, features = nm, val_slot = "intensity")

  eqn      <- obj@calibrationInfo@cal_list[[nm]]
  expected <- sort(log2(stats::model.frame(eqn)[["pg_perpixel"]]))
  plotted  <- sort(p[[2]]$data$x)

  expect_equal(plotted, expected)
})

test_that("missing calibration or concentration slot fails informatively", {
  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  skip_if_not(file.exists(file.path(cal_dir, "cal_MSI.RDS")))

  raw <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
            "quant_MSImagingExperiment")

  # No models fitted yet
  expect_error(plot_cal_coverage(raw), "create_cal_curve")

  # Models fitted, but int2conc() not run
  cal_metadata <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
  fitted <- summarise_cal_levels(raw, cal_metadata, val_slot = "intensity",
                                 cal_label = "Cal", id = "identifier")
  fitted <- create_cal_curve(fitted, cal_type = "cal")
  expect_error(plot_cal_coverage(fitted, val_slot = "intensity"), "int2conc")
})

test_that("an uncalibrated feature name is rejected with the available list", {
  obj <- calibrated_obj()
  expect_error(
    plot_cal_coverage(obj, features = "not_an_analyte", val_slot = "intensity"),
    "no calibration model for feature"
  )
})
