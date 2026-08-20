require(testthat)
require(quantMSImageR)
require(ComplexHeatmap)


make_obj <- function(n_features = 3, runs = c("s1", "s2"), n_px = 4,
                     seed = 42) {
  set.seed(seed)
  fdata <- MassDataFrame(mz   = seq_len(n_features),
                         name = paste0("lipid_", seq_len(n_features)))
  n_tot <- length(runs) * n_px
  pdata <- PositionDataFrame(
    run   = rep(runs, each = n_px),
    coord = do.call(rbind, lapply(seq_along(runs),
                                  function(i) expand.grid(x = 1:2, y = 1:2)))
  )
  ints <- matrix(runif(n_features * n_tot, 1, 100), nrow = n_features)
  as(MSImagingExperiment(spectraData = ints, featureData = fdata,
                          pixelData = pdata),
     "quant_MSImagingExperiment")
}

test_that("quantileHm returns a Heatmap object", {
  obj <- make_obj()
  hm  <- quantileHm(obj, quant_val = 0.5,
                     heatmap_order = c("s1", "s2"),
                     heatmap_labs  = c("A", "B"))
  expect_s4_class(hm, "Heatmap")
})

test_that("heatmap matrix is n_samples x n_features", {
  # Samples are rows and features are columns: studies usually have more
  # samples than features, so the long dimension runs vertically.
  n_feat <- 4; runs <- c("s1", "s2", "s3")
  obj <- make_obj(n_features = n_feat, runs = runs)
  hm  <- quantileHm(obj, quant_val = 0.75,
                     heatmap_order = runs,
                     heatmap_labs  = c("A", "B", "C"))
  mat <- hm@matrix
  expect_equal(nrow(mat), length(runs))
  expect_equal(ncol(mat), n_feat)
  expect_equal(rownames(mat), runs)
})

test_that("cells are square, stretching to at most max_aspect", {
  # Balanced panel -> square cells.
  obj <- make_obj(n_features = 3, runs = c("s1", "s2", "s3"))
  hm  <- quantileHm(obj, quant_val = 0.5, cell_size = 6)
  expect_equal(as.numeric(hm@matrix_param$width) / ncol(hm@matrix),
               as.numeric(hm@matrix_param$height) / nrow(hm@matrix))

  # Many features, few samples -> rows stretch, but no further than max_aspect.
  wide <- make_obj(n_features = 20, runs = c("s1", "s2"))
  hw   <- quantileHm(wide, quant_val = 0.5, cell_size = 6, max_aspect = 1.5)
  cw   <- as.numeric(hw@matrix_param$width)  / ncol(hw@matrix)
  ch   <- as.numeric(hw@matrix_param$height) / nrow(hw@matrix)
  expect_equal(ch / cw, 1.5)
})

test_that("cells carry a border so neighbours stay separable", {
  obj <- make_obj(n_features = 4, runs = c("s1", "s2"))

  hm <- quantileHm(obj, quant_val = 0.5)
  expect_equal(hm@matrix_param$gp$col, "white")

  none <- quantileHm(obj, quant_val = 0.5, cell_border = NA)
  expect_true(is.na(none@matrix_param$gp$col))
})

test_that("deprecated row_split arguments still work", {
  obj <- make_obj(n_features = 4, runs = c("s1", "s2"))
  fs  <- c("A", "A", "B", "B")

  new_arg <- quantileHm(obj, quant_val = 0.5, feature_split = fs,
                         feature_split_name = "Met-1")
  old_arg <- quantileHm(obj, quant_val = 0.5, row_split = fs,
                         row_split_name = "Met-1")

  expect_equal(old_arg@matrix, new_arg@matrix)
  expect_equal(levels(old_arg@column_order_list |> names()),
               levels(new_arg@column_order_list |> names()))
})

test_that("z-score values are clipped to [-1, 1]", {
  obj <- make_obj(n_features = 5, runs = c("s1", "s2", "s3", "s4"))
  hm  <- quantileHm(obj, quant_val = 0.5,
                     heatmap_order = c("s1", "s2", "s3", "s4"),
                     heatmap_labs  = c("A", "A", "B", "B"))
  mat <- hm@matrix
  expect_true(all(mat >= -1 & mat <= 1, na.rm = TRUE))
})

test_that("quantileHm works without heatmap_order / heatmap_labs", {
  obj <- make_obj()
  expect_s4_class(quantileHm(obj, quant_val = 0.5), "Heatmap")
})
