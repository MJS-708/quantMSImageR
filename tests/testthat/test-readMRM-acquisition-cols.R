require(testthat)
require(quantMSImageR)

# Two transitions can share a precursor and a product and differ only in
# collision energy, so the feature table has to carry the acquisition settings
# or there is no way to tell them apart short of opening the raw data.
test_that("the bundled ion library defines the acquisition settings", {
  lib <- read.csv(system.file("extdata", "example_ion_library.csv",
                              package = "quantMSImageR"), check.names = FALSE)
  expect_true(all(c("precursor_mz", "product_mz", "collision_eV", "cone_V")
                  %in% names(lib)))
  expect_true(is.numeric(lib$collision_eV))
  expect_true(is.numeric(lib$cone_V))
})

test_that("buildFeatureMeta carries them through to the feature table", {
  p <- system.file("extdata", "example.raw", "section01.RDS",
                   package = "quantMSImageR")
  obj <- readRDS(p)
  lib <- read.csv(system.file("extdata", "example_ion_library.csv",
                              package = "quantMSImageR"), check.names = FALSE)

  fm <- buildFeatureMeta(obj, lib)
  expect_true(all(c("collision_eV", "cone_V") %in% names(fm)))
  expect_equal(nrow(fm), nrow(fData(obj)))
  # Matched by m/z, so the values must be the library's own, not placeholders.
  expect_true(all(fm$collision_eV %in% lib$collision_eV))
})
