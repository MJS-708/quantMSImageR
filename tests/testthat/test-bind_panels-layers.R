# bind_panels() rebuilt the merged object from the intensity layer alone, so
# binding two processed panels silently returned raw data: response, snr and any
# calibrated layers were dropped, as were experimentData and the calibration
# slots. It also discarded obj2's copy of a shared feature name without
# checking the two were the same transition.

mk <- function(names, prec, prod, layers = character(0), px = 100) {
  n <- length(names)
  fd <- MassDataFrame(mz = seq_len(n), name = names,
                      precursor_mz = as.character(prec),
                      product_mz   = as.character(prod))
  pd <- PositionDataFrame(run = rep("s1", 4),
                          coord = expand.grid(x = 1:2, y = 1:2))
  o <- MSImagingExperiment(
    spectraData    = matrix(seq_len(n * 4), nrow = n),
    featureData    = fd, pixelData = pd,
    experimentData = CardinalIO::ImzMeta(pixelSize = px))
  o <- as(o, "quant_MSImagingExperiment")
  for (l in layers) spectra(o, l) <- matrix(seq_len(n * 4) * 2, nrow = n)
  o
}

test_that("every spectra layer survives the bind", {
  a <- mk(c("f1", "f2"), c(291.2, 355.2), c(193.1, 275.1), layers = c("response", "snr"))
  b <- mk(c("g1", "g2"), c(400.1, 410.1), c(100.1, 110.1), layers = c("response", "snr"))

  out <- bind_panels(a, b)

  expect_setequal(names(spectraData(out)), c("intensity", "response", "snr"))
  expect_equal(nrow(fData(out)), 4L)
  # Values must travel with their feature, obj1's first then obj2's.
  expect_equal(as.numeric(spectra(out, "response")[1, ]),
               as.numeric(spectra(a, "response")[1, ]))
  expect_equal(as.numeric(spectra(out, "response")[3, ]),
               as.numeric(spectra(b, "response")[1, ]))
})

test_that("mismatched layer sets are refused rather than silently reduced", {
  a <- mk(c("f1", "f2"), c(291.2, 355.2), c(193.1, 275.1), layers = "response")
  b <- mk(c("g1", "g2"), c(400.1, 410.1), c(100.1, 110.1))

  expect_error(bind_panels(a, b), regexp = "different spectra layers")
})

test_that("experiment metadata and calibration slots are carried over", {
  a <- mk(c("f1", "f2"), c(291.2, 355.2), c(193.1, 275.1))
  b <- mk(c("g1", "g2"), c(400.1, 410.1), c(100.1, 110.1))

  # pixelSize is what int2conc() needs for the pg/mm2 layer.
  out <- bind_panels(a, b)
  expect_equal(as.numeric(experimentData(out)$pixelSize), 100)

  cal_dir <- system.file("extdata", "cal_example.raw", package = "quantMSImageR")
  cal <- as(readRDS(file.path(cal_dir, "cal_MSI.RDS")), "quant_MSImagingExperiment")
  cal <- summarise_cal_levels(cal, read.csv(file.path(cal_dir, "calibration_metadata.csv")),
                              val_slot = "intensity", cal_label = "Cal",
                              id = "identifier")
  calibrationData(a) <- calibrationData(cal)
  expect_gt(nrow(calibrationLevels(bind_panels(a, b))), 0)
})

test_that("a shared name over a different transition is refused", {
  a <- mk(c("f1", "shared"), c(291.2, 355.2), c(193.1, 275.1))
  b <- mk(c("shared", "g2"), c(355.2, 410.1), c(300.0, 110.1))   # different product

  expect_error(bind_panels(a, b),
               regexp = "appear in both panels over different transitions")
  expect_error(bind_panels(a, b), regexp = "shared: 355.2 -> 275.1 vs 355.2 -> 300")

  # Opting out is possible, but has to be said.
  out <- suppressMessages(bind_panels(a, b, feature_match = "name"))
  expect_equal(as.character(fData(out)$name), c("f1", "shared", "g2"))
})

test_that("a genuinely shared transition is deduplicated quietly", {
  a <- mk(c("f1", "shared"), c(291.2, 355.2), c(193.1, 275.1))
  b <- mk(c("shared", "g2"), c(355.2, 410.1), c(275.1, 110.1))   # same transition

  out <- suppressMessages(bind_panels(a, b))
  expect_equal(as.character(fData(out)$name), c("f1", "shared", "g2"))
})
