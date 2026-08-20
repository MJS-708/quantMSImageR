fitted_cal <- function() {
  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  skip_if_not(file.exists(file.path(cal_dir, "cal_MSI.RDS")),
              "synthetic cal data not generated (run inst/scripts/generate_cal_data.R)")

  cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
            "quant_MSImagingExperiment")
  meta <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))

  cal <- summariseCalLevels(cal, meta, val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")
  createCalCurve(cal, cal_type = "cal")
}

test_that("calibration accessors return the underlying slots", {
  cal <- fitted_cal()

  expect_identical(calibrationModels(cal),   cal@calibrationInfo@cal_list)
  expect_identical(calibrationDiagnostics(cal), cal@calibrationInfo@r2_df)
  # calibrationR2() is narrowed to what its name promises
  expect_setequal(names(calibrationR2(cal)), c("feature", "r2"))
  expect_identical(calibrationLevels(cal),   cal@calibrationInfo@cal_response_data)
  expect_identical(calibrationMetadata(cal), cal@calibrationInfo@cal_metadata)
  expect_identical(calibrationData(cal),     cal@calibrationInfo)

  expect_type(calibrationModels(cal), "list")
  expect_true(all(vapply(calibrationModels(cal), inherits, logical(1), "lm")))
})

test_that("calibrationData<- carries fitted models onto another object", {
  cal <- fitted_cal()

  p <- system.file("extdata", "example.raw", "section01.RDS",
                   package = "quantMSImageR")
  study <- as(readRDS(p), "quant_MSImagingExperiment")
  expect_length(calibrationModels(study), 0L)

  calibrationData(study) <- calibrationData(cal)
  expect_identical(calibrationModels(study), calibrationModels(cal))
})

test_that("validity rejects an unnamed non-empty cal_list", {
  # Positional matching is what silently applied one analyte's curve to
  # another before; the validity method makes that unrepresentable.
  cal <- fitted_cal()
  bad <- calibrationModels(cal)
  names(bad) <- NULL

  expect_error(calibrationModels(cal) <- bad, "must be named")
})

test_that("validity rejects models absent from r2_df", {
  cal <- fitted_cal()
  bad <- calibrationModels(cal)
  names(bad)[1] <- "not_an_analyte"

  expect_error(calibrationModels(cal) <- bad, "absent from r2_df")
})

test_that("tissue accessors return the datamatrix slots", {
  p <- system.file("extdata", "example.raw", "section01.RDS",
                   package = "quantMSImageR")
  obj <- as(readRDS(p), "quant_MSImagingExperiment")
  obj <- createMSIDataMatrix(obj, val_slot = "intensity", roi_header = NA)

  expect_identical(tissueData(obj), obj@tissueInfo)
  expect_identical(tissueMatrix(obj, which = "pixel"),
                   obj@tissueInfo@all_pixel_matrix)
  expect_identical(tissueMatrix(obj, which = "roi"),
                   obj@tissueInfo@roi_average_matrix)
  expect_gt(nrow(tissueMatrix(obj, which = "pixel")), 0L)
})
