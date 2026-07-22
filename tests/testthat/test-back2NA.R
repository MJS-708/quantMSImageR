require(testthat)
require(quantMSImageR)


make_obj <- function(sample_ids) {
  fdata <- MassDataFrame(mz = 1:2, name = c("feat_1", "feat_2"))
  n     <- length(sample_ids)
  coord <- data.frame(x = seq_len(n), y = rep(1L, n))
  pdata <- PositionDataFrame(run = rep("s1", n), coord = coord,
                             sample_ID = sample_ids)
  ints  <- matrix(seq_len(2 * n), nrow = 2)
  as(MSImagingExperiment(spectraData = ints, featureData = fdata,
                          pixelData = pdata),
     "quant_MSImagingExperiment")
}

test_that("background pixels become NA, tissue pixels are unchanged", {
  ids <- c("tissue_pixels", "tissue_pixels", "noise_pixels", "noise_pixels")
  obj <- make_obj(ids)

  result    <- back2NA(obj, val_slot = "intensity",
                       background = "noise_pixels",
                       pixel_header = "sample_ID")
  int_orig  <- spectraData(obj)[["intensity"]]
  int_new   <- spectraData(result)[["intensity"]]
  noise_px  <- which(ids == "noise_pixels")
  tissue_px <- which(ids == "tissue_pixels")

  expect_true(all(is.na(int_new[, noise_px])))
  expect_equal(int_new[, tissue_px], int_orig[, tissue_px])
})

test_that("returns object unchanged when no background pixels exist", {
  obj    <- make_obj(rep("tissue_pixels", 4))
  result <- back2NA(obj, val_slot = "intensity",
                    background = "noise_pixels",
                    pixel_header = "sample_ID")
  expect_equal(spectraData(result)[["intensity"]],
               spectraData(obj)[["intensity"]])
})

test_that("the deprecated sample_type alias still selects the column", {
  ids <- c("tissue_pixels", "tissue_pixels", "noise_pixels", "noise_pixels")
  obj <- make_obj(ids)

  new_arg <- back2NA(obj, val_slot = "intensity",
                     background = "noise_pixels", pixel_header = "sample_ID")
  old_arg <- back2NA(obj, val_slot = "intensity",
                     background = "noise_pixels", sample_type = "sample_ID")

  expect_equal(spectraData(old_arg)[["intensity"]],
               spectraData(new_arg)[["intensity"]])
})
