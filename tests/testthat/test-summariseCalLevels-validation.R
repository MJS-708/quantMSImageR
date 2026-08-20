# A calibration design problem must stop the run. Each of these cases used to
# produce a curve that looked fitted but was not anchored to the data it claimed.

cal_obj <- function() {
  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  as(readRDS(file.path(cal_dir, "cal_MSI.RDS")), "quant_MSImagingExperiment")
}
cal_md <- function() {
  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  read.csv(file.path(cal_dir, "calibration_metadata.csv"))
}

test_that("dispersion is reported alongside the spot mean", {
  # A mean and a pixel count cannot show whether a level was measured
  # reproducibly, which is the first thing to check before trusting a curve.
  out <- summariseCalLevels(cal_obj(), cal_md(), val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")
  rd  <- calibrationLevels(out)

  expect_true(all(c("response_perpixel", "median_response", "sd_response",
                    "cv_response", "pixels", "n_nonmissing", "pg_perpixel",
                    "analyte") %in% names(rd)))
  expect_true(all(rd$n_nonmissing > 0))
  expect_true(all(is.finite(rd$cv_response)))
})

test_that("a legacy lipid column is read as analyte", {
  md <- cal_md()
  names(md)[names(md) == "analyte"] <- "lipid"
  out <- summariseCalLevels(cal_obj(), md, val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")
  expect_true("analyte" %in% names(calibrationLevels(out)))
})

test_that("an analyte absent from the object is an error", {
  md <- cal_md()
  md$analyte[1] <- "not_a_feature"
  expect_error(
    summariseCalLevels(cal_obj(), md, val_slot = "intensity",
                         cal_label = "Cal", id = "identifier"),
    regexp = "not features of this object"
  )
})

test_that("missing metadata columns are named", {
  md <- cal_md()
  md$amount_pg <- NULL
  expect_error(
    summariseCalLevels(cal_obj(), md, val_slot = "intensity",
                         cal_label = "Cal", id = "identifier"),
    regexp = "missing column\\(s\\): amount_pg"
  )
})

test_that("contradictory duplicate identifiers are an error", {
  md  <- cal_md()
  bad <- md[1, ]
  bad$amount_pg <- bad$amount_pg * 2
  expect_error(
    summariseCalLevels(cal_obj(), rbind(md, bad), val_slot = "intensity",
                         cal_label = "Cal", id = "identifier"),
    regexp = "disagree about amount_pg or level"
  )
})

test_that("a calibration label that matches no pixel is an error", {
  expect_error(
    summariseCalLevels(cal_obj(), cal_md(), val_slot = "intensity",
                         cal_label = "NotALabel", id = "identifier"),
    regexp = "no pixels are labelled"
  )
})
