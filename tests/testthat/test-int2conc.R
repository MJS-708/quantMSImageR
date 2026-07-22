test_that("int2conc converts intensities to calibrated amount layers", {

  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  skip_if_not(file.exists(file.path(cal_dir, "cal_MSI.RDS")),
              "synthetic cal data not generated (run inst/scripts/generate_cal_data.R)")

  cal  <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
             "quant_MSImagingExperiment")
  meta <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))

  cal <- summarise_cal_levels(cal, meta, val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")
  cal <- create_cal_curve(cal, cal_type = "cal")

  out <- int2conc(cal, pixel_header = "sample_type", pixels = "Tissue")

  pp <- spectra(out, "pg_pixel")
  mm <- spectra(out, "pg_mm2")

  expect_true(any(is.finite(pp)))
  # pg_mm2 = pg_pixel * (1000 / pixelSize)^2, with pixelSize = 100 -> x 100
  fin <- which(is.finite(pp))
  expect_equal(as.numeric(mm[fin][1]),
               as.numeric(pp[fin][1]) * (1000 / 100)^2)
})

test_that("readers accept the pre-rename slot spellings", {
  # Objects calibrated with earlier versions carry `conc - pg/pixel`; the
  # helpers must still resolve them so old .RDS files keep working.
  expect_true("conc - pg/pixel" %in% quantMSImageR:::.conc_slots("pg_pixel"))
  expect_true("pg_mm2" %in% quantMSImageR:::.conc_slots("conc - pg/mm2"))
  expect_equal(quantMSImageR:::.conc_unit("conc - pg/mm2"), "pg/mm2")
  expect_equal(quantMSImageR:::.conc_unit("pg_pixel"), "pg/pixel")
})
