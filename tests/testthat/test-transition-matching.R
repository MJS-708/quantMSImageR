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

test_that("matching is at unit resolution, because MRM is", {
  # A triple quadrupole running MRM sets Q1 and Q3 at unit resolution, so the
  # decimals in a library entry are a transcription detail, not a measurement.
  # 291.23 -> 193.13 and 291.0 -> 193.0 are the same channel as 291.2 -> 193.1.
  j <- .join_ion_library(trans(291.23, 193.13), lib(), polarity = "Negative")
  expect_equal(j$transition_id_name, "13-HOTE")

  j0 <- .join_ion_library(trans(291.0, 193.0), lib(), polarity = "Negative")
  expect_equal(j0$transition_id_name, "13-HOTE")

  # A whole nominal unit out on the product is a different channel.
  j1 <- .join_ion_library(trans(291.2, 194.3), lib(), polarity = "Negative")
  expect_equal(j1$transition_id_name, "1")
  expect_equal(j1$Type, "Unknown")

  # A transition the library does not describe is kept, not dropped.
  j2 <- .join_ion_library(trans(500.0, 100.0), lib(), polarity = "Negative")
  expect_equal(nrow(j2), 1L)
  expect_equal(j2$transition_id_name, "1")
  expect_equal(j2$Type, "Unknown")
})

test_that("precursor and product must BOTH match", {
  # Right precursor, wrong product: not a match. Matching on the precursor
  # alone would annotate a feature as a compound the product ion rules out.
  j <- .join_ion_library(trans(291.2, 300.0), lib(), polarity = "Negative")
  expect_equal(j$Type, "Unknown")

  # And the other way round.
  j2 <- .join_ion_library(trans(400.0, 193.1), lib(), polarity = "Negative")
  expect_equal(j2$Type, "Unknown")
})

test_that("isomers separated by more than a unit on the product do not collide", {
  # All three share precursor 291.2; the product ion is what distinguishes
  # them, and at unit resolution 169.1 / 193.1 / 221.2 still do.
  j <- .join_ion_library(trans(rep(291.2, 3), c(169.1, 193.1, 221.2)), lib(),
                         polarity = "Negative")
  expect_equal(j$transition_id_name, c("9-HOTE", "13-HOTE", "15-HOTE"))
})

test_that("an ambiguous match is named for every entry it matched", {
  # Widening to 30 Da makes the three 291.2 isomers indistinguishable. Two
  # isomers really can share a transition at unit resolution -- LTC4 and
  # 14_15-LTC4 are 0.03 Da apart on the product -- and naming the feature for
  # one of them would assert more than was measured.
  expect_message(
    j <- .join_ion_library(trans(291.2, 193.1), lib(), polarity = "Negative",
                           mz_tolerance = 30),
    regexp = "9-HOTE, 13-HOTE, 15-HOTE")
  expect_equal(j$transition_id_name, "9-HOTE || 13-HOTE || 15-HOTE")
  # transition_id must not disagree with the name it was built from.
  expect_equal(j$transition_id, "9-HOTE || 13-HOTE || 15-HOTE")
  # Annotation still comes from the closest entry.
  expect_equal(j$Type, "Analyte")
})

test_that("combining never merges an analyte with an internal standard", {
  # One feature cannot be both: it would silently break IS normalisation.
  # A library that puts them on one transition is wrong, not ambiguous.
  expect_error(
    .join_ion_library(trans(291.2, 193.1), lib(), polarity = "Negative",
                      mz_tolerance = 100),
    regexp = "differing Type")
})

test_that("the stricter ambiguity modes still behave", {
  expect_error(
    .join_ion_library(trans(291.2, 193.1), lib(), polarity = "Negative",
                      mz_tolerance = 30, ambiguity = "error"),
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

  # An ambiguous library must be refused here too when asked to be, not
  # resolved by position.
  expect_error(build_feature_meta(obj, l, mz_tolerance = 50, verbose = FALSE,
                                  ambiguity = "error"),
               regexp = "match more than one ion-library entry")
})
