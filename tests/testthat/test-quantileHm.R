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

test_that("a feature with no variance is unscoreable, not drawn as low", {
  # A flat feature has sd 0, so every z-score is 0/0. An earlier version sent
  # those to -1, which painted an unvarying feature as uniformly depleted --
  # indistinguishable from one genuinely low in every sample.
  runs <- c("s1", "s2", "s3")
  obj  <- make_obj(n_features = 3, runs = runs)
  spectra(obj)[2, ] <- 50            # lipid_2 identical in every pixel

  hm  <- quantileHm(obj, quant_val = 0.5, heatmap_order = runs,
                    heatmap_labs = c("A", "B", "C"))
  flat <- hm@matrix[, "lipid_2"]

  expect_true(all(is.na(flat)))
  expect_false(any(flat == -1, na.rm = TRUE))
  # The features that do vary are still scored.
  expect_false(anyNA(hm@matrix[, "lipid_1"]))
  # apply() takes its row names from the first feature and drops them for the
  # whole matrix if another feature's differ, so the unscoreable row has to
  # carry the sample names too -- otherwise the panel loses its row labels,
  # and only in studies that happen to contain a flat feature.
  expect_equal(rownames(hm@matrix), runs)
})

test_that("unscoreable cells use na_col, matching contributionHm", {
  runs <- c("s1", "s2", "s3")
  obj  <- make_obj(n_features = 3, runs = runs)

  # ComplexHeatmap stores colours as #RRGGBBAA ("grey88" -> "#E0E0E0FF"), so
  # compare the colour, not the string.
  same_colour <- function(a, b) expect_equal(grDevices::col2rgb(a), grDevices::col2rgb(b))

  hm <- quantileHm(obj, quant_val = 0.5, heatmap_order = runs)
  same_colour(hm@matrix_color_mapping@na_col, "grey88")

  custom <- quantileHm(obj, quant_val = 0.5, heatmap_order = runs,
                       na_col = "black")
  same_colour(custom@matrix_color_mapping@na_col, "black")
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

# ---- Row blocks (regions of interest) --------------------------------------

test_that("without sample_block nothing about the panel changes", {
  # The region work added an argument to this function; a study with no regions
  # must get exactly the panel it got before.
  runs <- c("s1", "s2", "s3", "s4")
  obj  <- make_obj(n_features = 3, runs = runs)
  labs <- c("A", "A", "B", "B")

  hm <- quantileHm(obj, quant_val = 0.5, heatmap_order = runs,
                   heatmap_labs = labs)

  # Still the clipped z-score of each feature across ALL four samples, computed
  # here by hand so the expectation does not just restate the implementation.
  med <- vapply(runs, function(r)
    stats::quantile(as.numeric(spectra(obj)[1, pData(obj)$run == r]), 0.5),
    numeric(1))
  z <- pmin(pmax((med - mean(med)) / stats::sd(med), -1), 1)
  expect_equal(unname(hm@matrix[, "lipid_1"]), unname(z))
  # Rows are still split by group -- two slices, the A samples then the B ones
  # -- and the group bar is the only annotation.
  expect_length(hm@row_order_list, 2L)
  expect_equal(sort(hm@row_order_list[[1]]), 1:2)
  expect_equal(names(hm@left_annotation@anno_list), "Group")
})

test_that("sample_block scores each block against itself", {
  runs <- c("s1", "s2", "s3", "s4")
  obj  <- make_obj(n_features = 3, runs = runs)

  hm <- quantileHm(obj, quant_val = 0.5, heatmap_order = runs,
                   heatmap_labs  = c("A", "B", "A", "B"),
                   sample_block  = c("a", "a", "b", "b"))
  m <- hm@matrix

  # Within a block of two, the two rows are mirror images: each block was
  # centred on its own mean rather than on the panel's.
  expect_equal(unname(m["s1", ] + m["s2", ]), rep(0, ncol(m)), tolerance = 1e-8)
  expect_equal(unname(m["s3", ] + m["s4", ]), rep(0, ncol(m)), tolerance = 1e-8)
  # Rows split by REGION, not by group: the groups here alternate (A, B, A, B),
  # so a first slice of rows 1-2 can only have come from the block.
  expect_length(hm@row_order_list, 2L)
  expect_equal(sort(hm@row_order_list[[1]]), 1:2)
  # The group keeps its own colour bar beside the region one.
  expect_equal(names(hm@left_annotation@anno_list), c("Region", "Group"))
})

test_that("a block of one sample is unscoreable rather than mid-ramp", {
  # One airway on its own has nothing to be high or low against; drawing it at
  # the middle of the ramp would read as an ordinary result.
  runs <- c("s1", "s2", "s3")
  obj  <- make_obj(n_features = 3, runs = runs)

  hm <- quantileHm(obj, quant_val = 0.5, heatmap_order = runs,
                   sample_block = c("a", "a", "b"))
  expect_true(all(is.na(hm@matrix["s3", ])))
  expect_false(anyNA(hm@matrix["s1", ]))
})

test_that("sample_block must have one entry per sample", {
  obj <- make_obj(n_features = 3, runs = c("s1", "s2", "s3"))
  expect_error(quantileHm(obj, quant_val = 0.5,
                          heatmap_order = c("s1", "s2", "s3"),
                          sample_block  = c("a", "b")),
               regexp = "one entry per sample|sample_block has")
})
