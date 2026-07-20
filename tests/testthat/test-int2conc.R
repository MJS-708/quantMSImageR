test_that("int2conc converts intensities to concentration layers", {

  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  skip_if_not(file.exists(file.path(cal_dir, "cal_MSI.RDS")),
              "synthetic cal data not generated (run inst/generate_cal_data.R)")

  cal  <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
             "quant_MSImagingExperiment")
  meta <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))

  cal <- summarise_cal_levels(cal, meta, val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")
  cal <- create_cal_curve(cal, cal_type = "cal")

  out <- int2conc(cal, val_slot = "intensity", pixels = "Tissue")

  pp <- spectra(out, "conc - pg/pixel")
  mm <- spectra(out, "conc - pg/mm2")

  expect_true(any(is.finite(pp)))
  # pg/mm2 = pg/pixel * (1000 / pixelSize)^2, with pixelSize = 100 -> x 100
  fin <- which(is.finite(pp))
  expect_equal(as.numeric(mm[fin][1]),
               as.numeric(pp[fin][1]) * (1000 / 100)^2)
})
