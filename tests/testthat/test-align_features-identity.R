# align_features() rewrites obj2's mz key so Cardinal's cbind will accept it.
# That makes an incorrect pairing unrecoverable, so a shared name is no longer
# taken as proof that two features are the same transition.

mk <- function(names, prec, prod, run = "s1") {
  n <- length(names)
  as(MSImagingExperiment(
       spectraData = matrix(seq_len(n * 4), nrow = n),
       featureData = MassDataFrame(mz = seq_len(n), name = names,
                                   precursor_mz = as.character(prec),
                                   product_mz   = as.character(prod)),
       pixelData   = PositionDataFrame(run = rep(run, 4),
                                       coord = expand.grid(x = 1:2, y = 1:2))),
     "quant_MSImagingExperiment")
}

test_that("matching transitions align and keep obj1's mz key", {
  a <- mk(c("f1", "f2"), c(291.2, 355.2), c(193.1, 275.1))
  b <- mk(c("f2", "f1"), c(355.2, 291.2), c(275.1, 193.1), run = "s2")

  al <- align_features(a, b)
  expect_equal(as.character(fData(al$obj1)$name),
               as.character(fData(al$obj2)$name))
  expect_equal(mz(al$obj2), mz(al$obj1))
})

test_that("features listed in a different order are reordered, with their data", {
  # Cardinal's MassDataFrame needs a sorted mz key, so a permutation that
  # unsorts it cannot be expressed as a subset -- two acquisitions listing the
  # same features in a different order used to fail outright.
  a <- mk(c("f1", "f2"), c(291.2, 355.2), c(193.1, 275.1))
  b <- mk(c("f2", "f1"), c(355.2, 291.2), c(275.1, 193.1), run = "s2")
  spectra(b, "snr") <- matrix(seq_len(8) * 2, nrow = 2)

  before <- as.matrix(spectraData(b)[["intensity"]])
  al <- align_features(a, b)

  expect_equal(as.character(fData(al$obj2)$name), c("f1", "f2"))
  # The rows must travel with their feature, not stay put.
  expect_equal(as.numeric(spectraData(al$obj2)[["intensity"]][1, ]),
               as.numeric(before[2, ]))
  # Every layer is carried, and the acquisition is still identifiable.
  expect_setequal(names(spectraData(al$obj2)), c("intensity", "snr"))
  expect_equal(unique(as.character(pData(al$obj2)$run)), "s2")
})

test_that("a shared name over a different transition is refused", {
  # f2 is 355.2 -> 275.1 in one method and 355.2 -> 300.0 in the other: the
  # same label over a different product ion.
  a <- mk(c("f1", "f2"), c(291.2, 355.2), c(193.1, 275.1))
  b <- mk(c("f1", "f2"), c(291.2, 355.2), c(193.1, 300.0), run = "s2")

  expect_error(align_features(a, b),
               regexp = "share a name but not a transition")
  # The message must name the offender and show both transitions.
  expect_error(align_features(a, b), regexp = "f2: 355.2 -> 275.1 vs 355.2 -> 300")

  # Opting out is possible, but has to be said.
  expect_silent(al <- align_features(a, b, feature_match = "name"))
  expect_equal(mz(al$obj2), mz(al$obj1))
})

test_that("agreement is judged to one decimal place", {
  a <- mk(c("f1", "f2"), c(291.2, 355.2), c(193.1, 275.1))
  b <- mk(c("f1", "f2"), c(291.23, 355.24), c(193.13, 275.14), run = "s2")
  expect_silent(align_features(a, b))

  # A whole nominal unit apart is a different transition, and the default no
  # longer waves it through.
  c2 <- mk(c("f1", "f2"), c(291.0, 355.4), c(193.0, 275.3), run = "s2")
  expect_error(align_features(a, c2),
               regexp = "share a name but not a transition")
  expect_silent(align_features(a, c2, mz_tolerance = 0.5))
})

test_that("objects without transition metadata must opt in to name matching", {
  bare <- function(run) as(MSImagingExperiment(
    spectraData = matrix(1:8, nrow = 2),
    featureData = MassDataFrame(mz = 1:2, name = c("f1", "f2")),
    pixelData   = PositionDataFrame(run = rep(run, 4),
                                    coord = expand.grid(x = 1:2, y = 1:2))),
    "quant_MSImagingExperiment")

  expect_error(align_features(bare("s1"), bare("s2")),
               regexp = "needs precursor_mz and product_mz")
  expect_silent(align_features(bare("s1"), bare("s2"), feature_match = "name"))
})
