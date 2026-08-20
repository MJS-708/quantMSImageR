# zero2NA() modifies a chosen layer in place. Its all-zero branch used to write
# through spectra(), which targets the default layer -- so calling it with
# val_slot = "response" blanked the intensity data instead.

make_obj <- function() {
  fdata <- MassDataFrame(mz = 1:2, name = c("f1", "f2"))
  pdata <- PositionDataFrame(run = rep("s1", 4),
                             coord = expand.grid(x = 1:2, y = 1:2))
  ints  <- matrix(c(10, 20, 30, 40,      # f1: ordinary values
                    5,  0,  7,  0),      # f2: some zeros
                  nrow = 2, byrow = TRUE)
  obj <- as(MSImagingExperiment(spectraData = ints, featureData = fdata,
                                 pixelData = pdata),
            "quant_MSImagingExperiment")
  # A second layer whose second feature is entirely zero.
  resp <- matrix(c(1, 2, 3, 4,
                   0, 0, 0, 0), nrow = 2, byrow = TRUE)
  spectra(obj, "response") <- resp
  obj
}

test_that("an all-zero feature blanks the requested layer, not another one", {
  obj <- make_obj()
  before <- spectraData(obj)[["intensity"]]

  out <- suppressMessages(zero2NA(obj, val_slot = "response"))

  expect_true(all(is.na(spectraData(out)[["response"]][2, ])))
  # The layer that was not asked for must be untouched.
  expect_equal(spectraData(out)[["intensity"]], before)
})

test_that("scattered zeros in the requested layer become NA", {
  obj <- make_obj()
  out <- suppressMessages(zero2NA(obj, val_slot = "intensity"))

  expect_equal(as.numeric(spectraData(out)[["intensity"]][2, ]),
               c(5, NA, 7, NA))
  expect_equal(as.numeric(spectraData(out)[["intensity"]][1, ]),
               c(10, 20, 30, 40))
  # response was not the target, so its all-zero feature stays zero.
  expect_true(all(spectraData(out)[["response"]][2, ] == 0))
})
