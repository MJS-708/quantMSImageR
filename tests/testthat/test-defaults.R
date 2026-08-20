# The documented defaults have drifted from the implementation before. These
# tests pin them so the help pages and the code cannot diverge silently again.

bundled_section <- function() {
  p <- system.file("extdata", "example.raw", "section01.RDS",
                   package = "quantMSImageR")
  as(readRDS(p), "quant_MSImagingExperiment")
}

test_that("int2SNR defaults to median, as documented", {
  obj <- bundled_section()

  default <- int2SNR(obj, snr_thresh = 3)
  median  <- int2SNR(obj, snr_thresh = 3, average = "median")
  mean    <- int2SNR(obj, snr_thresh = 3, average = "mean")

  expect_equal(spectra(default, "snr"), spectra(median, "snr"))
  # And the two summaries genuinely differ, so the test above has teeth.
  expect_false(isTRUE(all.equal(spectra(median, "snr"), spectra(mean, "snr"))))
})

test_that("filtering defaults match the imaging convention out of the box", {
  # sample_name / tissue_pixels / background_pixels, so the bundled example
  # needs no argument overrides.
  obj <- bundled_section()

  expect_silent(snr <- int2SNR(obj, snr_thresh = 3))
  expect_true("snr" %in% names(spectraData(snr)))
  expect_true(any(is.finite(spectra(snr, "snr"))))

  bg <- back2NA(obj)
  bg_px <- which(pData(obj)$sample_name == "background_pixels")
  expect_true(all(is.na(spectra(bg, "intensity")[, bg_px])))
})

test_that("int2SNR adds a slot rather than replacing intensity", {
  obj <- bundled_section()
  out <- int2SNR(obj, snr_thresh = 3)

  expect_equal(spectra(out, "intensity"), spectra(obj, "intensity"))
  expect_true("snr" %in% names(spectraData(out)))
})

test_that("createCalCurve requires cal_type explicitly", {
  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  skip_if_not(file.exists(file.path(cal_dir, "cal_MSI.RDS")))

  cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
            "quant_MSImagingExperiment")
  meta <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
  cal <- summariseCalLevels(cal, meta, val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")

  expect_error(createCalCurve(cal), "must be given explicitly")
  expect_error(createCalCurve(cal, cal_type = "nonsense"), "should be one of")
})

test_that("int2conc warns above 10 percent extrapolated pixels", {
  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  skip_if_not(file.exists(file.path(cal_dir, "cal_MSI.RDS")))

  cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
            "quant_MSImagingExperiment")
  meta <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
  cal <- summariseCalLevels(cal, meta, val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")
  cal <- createCalCurve(cal, cal_type = "cal")

  # S4 hides the real signature inside .local, so read the deparsed method.
  src <- paste(deparse(getMethod("int2conc", "quant_MSImagingExperiment")),
               collapse = " ")
  expect_match(src, "max_out_of_range\\s*=\\s*0\\.1")

  out <- suppressMessages(int2conc(cal, pixel_header = "sample_type",
                                   pixels = "Tissue"))
  diag <- calibrationDiagnostics(out)
  expect_true("out_of_range" %in% names(diag))

  # The bundled calibration brackets the bundled tissue signal, so nothing is
  # extrapolated. This doubles as a regression test on the example data: if a
  # future change to generate_cal_data.R breaks that coverage, the vignette's
  # plotCalCoverage() figure would stop demonstrating the good case.
  expect_true(all(diag$out_of_range == 0))
  expect_no_warning(
    suppressMessages(int2conc(cal, pixel_header = "sample_type",
                              pixels = "Tissue", max_out_of_range = 0)))
})

test_that("calibrationR2 returns only R-squared; diagnostics returns the rest", {
  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  skip_if_not(file.exists(file.path(cal_dir, "cal_MSI.RDS")))

  cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")),
            "quant_MSImagingExperiment")
  meta <- read.csv(file.path(cal_dir, "calibration_metadata.csv"))
  cal <- summariseCalLevels(cal, meta, val_slot = "intensity",
                              cal_label = "Cal", id = "identifier")
  cal <- createCalCurve(cal, cal_type = "cal")
  out <- suppressMessages(suppressWarnings(
    int2conc(cal, pixel_header = "sample_type", pixels = "Tissue")))

  expect_setequal(names(calibrationR2(out)), c("feature", "r2"))
  expect_true("out_of_range" %in% names(calibrationDiagnostics(out)))
})
