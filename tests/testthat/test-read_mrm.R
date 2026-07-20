test_that("read_mrm returns a cached acquisition when imaging/ is absent", {

  extdata <- system.file("extdata", package = "quantMSImageR")
  lib     <- file.path(extdata, "ion_library_pos04.csv")

  # pos04_test.raw ships only the cached .rds (no imaging/ text files), so
  # read_mrm returns the cache.
  obj <- read_mrm("pos04_test", folder = extdata, lib_ion_path = lib,
                  overwrite = FALSE)

  expect_s4_class(obj, "MSImagingExperiment")
  expect_gt(ncol(obj), 0)
  expect_gt(nrow(obj), 0)
  expect_true("name" %in% names(fData(obj)))
})
