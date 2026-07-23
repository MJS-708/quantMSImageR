# A transition's identity is its precursor/product pair. read_mrm() and
# build_feature_meta() used to decide that separately -- one rounding m/z to
# whole numbers, the other not -- and neither noticed when one measured
# transition matched several library entries.

lib <- function() data.frame(
  transition_id = c("9-HOTE", "13-HOTE", "15-HOTE", "PGE2-d4"),
  precursor_mz  = c(291.2, 291.2, 291.2, 355.2),
  product_mz    = c(169.1, 193.1, 221.2, 275.1),
  collision_eV  = 12, cone_V = 25,
  Polarity      = "Negative",
  Type          = c("Analyte", "Analyte", "Analyte", "IS"),
  stringsAsFactors = FALSE)

trans <- function(prec, prod) data.frame(
  transition_id = seq_along(prec), precursor_mz = prec, product_mz = prod)

test_that("a measured transition is matched to its library entry", {
  j <- .join_ion_library(trans(c(291.2, 355.2), c(193.1, 275.1)), lib(),
                         polarity = "Negative")
  expect_equal(j$transition_id_name, c("13-HOTE", "PGE2-d4"))
  expect_equal(j$Type, c("Analyte", "IS"))
  expect_equal(j$new_transition_int, 1:2)
})

test_that("matching is to one decimal place by default", {
  # The instrument reports 291.23 where the library says 291.2: same
  # transition to the precision an MRM library is quoted at.
  j <- .join_ion_library(trans(291.23, 193.13), lib(), polarity = "Negative")
  expect_equal(j$transition_id_name, "13-HOTE")

  # A whole nominal unit away is a different transition, not a rounding
  # difference. Nominal-mass matching could not tell these apart.
  j0 <- .join_ion_library(trans(291.0, 193.0), lib(), polarity = "Negative")
  expect_equal(j0$transition_id_name, "1")
  expect_equal(j0$Type, "Unknown")

  # Widening to nominal mass brings it back.
  j1 <- .join_ion_library(trans(291.0, 193.0), lib(), polarity = "Negative",
                          mz_tolerance = 0.5)
  expect_equal(j1$transition_id_name, "13-HOTE")

  # A transition the library does not describe is kept, not dropped.
  j2 <- .join_ion_library(trans(500.0, 100.0), lib(), polarity = "Negative")
  expect_equal(nrow(j2), 1L)
  expect_equal(j2$transition_id_name, "1")
  expect_equal(j2$Type, "Unknown")
})

test_that("isomers separated only by product ion do not collide", {
  # All three share precursor 291.2; the product ion is what distinguishes
  # them, and at the default tolerance it does.
  j <- .join_ion_library(trans(rep(291.2, 3), c(169.1, 193.1, 221.2)), lib(),
                         polarity = "Negative")
  expect_equal(j$transition_id_name, c("9-HOTE", "13-HOTE", "15-HOTE"))
})

test_that("an ambiguous match is an error, and says why", {
  # Widening the tolerance to 30 Da makes the three 291.2 isomers
  # indistinguishable -- exactly the collapse nominal-mass matching risks.
  expect_error(
    .join_ion_library(trans(291.2, 193.1), lib(), polarity = "Negative",
                      mz_tolerance = 30),
    regexp = "match more than one ion-library entry")

  expect_warning(
    j <- .join_ion_library(trans(291.2, 193.1), lib(), polarity = "Negative",
                           mz_tolerance = 30, ambiguity = "warn"),
    regexp = "9-HOTE, 13-HOTE, 15-HOTE")
  # "warn" takes the closest, not the first.
  expect_equal(j$transition_id_name, "13-HOTE")

  expect_silent(
    j2 <- .join_ion_library(trans(291.2, 193.1), lib(), polarity = "Negative",
                            mz_tolerance = 30, ambiguity = "nearest"))
  expect_equal(j2$transition_id_name, "13-HOTE")
})

test_that("polarity is honoured and duplicate names are made unique", {
  l <- lib(); l$Polarity[1] <- "Positive"
  j <- .join_ion_library(trans(291.2, 169.1), l, polarity = "Negative")
  expect_equal(j$Type, "Unknown")          # the only match is positive-mode

  # The same transition acquired twice must not produce two identical names.
  j2 <- .join_ion_library(trans(c(291.2, 291.2), c(193.1, 193.1)), lib(),
                          polarity = "Negative")
  expect_equal(j2$transition_id_name, c("13-HOTE", "2:- 13-HOTE"))
})

test_that("build_feature_meta uses the same rule as the join", {
  p   <- system.file("extdata", "example.raw", "section01.RDS",
                     package = "quantMSImageR")
  obj <- readRDS(p)
  l   <- read.csv(system.file("extdata", "example_ion_library.csv",
                              package = "quantMSImageR"), check.names = FALSE)

  fm <- build_feature_meta(obj, l, verbose = FALSE)
  expect_equal(nrow(fm), nrow(fData(obj)))
  expect_equal(as.character(fm$transition_id), as.character(fData(obj)$name))

  # An ambiguous library must be refused here too, not resolved by position.
  expect_error(build_feature_meta(obj, l, mz_tolerance = 50, verbose = FALSE),
               regexp = "match more than one ion-library entry")
})
