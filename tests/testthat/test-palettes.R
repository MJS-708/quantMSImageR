test_that("quant_palettes returns the documented palettes", {
  p <- quant_palettes()
  expect_setequal(names(p), c("heatmap0", "heatmap2", "hat", "reading"))
  expect_true(all(vapply(p, function(x) all(grepl("^#[0-9A-Fa-f]{6}$", x)),
                         logical(1))))

  # Continuous ramps interpolate to n; qualitative sets recycle so hues stay
  # distinct rather than being blended into muddy intermediates.
  expect_length(quant_palettes("heatmap0", n = 20), 20)
  expect_length(quant_palettes("hat", n = 3), 3)
  expect_equal(quant_palettes("hat", n = 3), quant_palettes("hat")[1:3])
})

test_that("heatmap2 runs blue -> white -> red for z-scores", {
  # Low must be blue and high red: the ltc ordering is the other way round and
  # is reversed on import, so this guards that reversal.
  h <- quant_palettes("heatmap2")
  expect_equal(tolower(h[1]), "#0571b0")
  expect_equal(tolower(h[length(h)]), "#ca0020")
})

test_that("imageR palettes actually change the fill scale", {
  p <- system.file("extdata", "example.raw", "section01.RDS",
                   package = "quantMSImageR")
  obj <- as(readRDS(p), "quant_MSImagingExperiment")

  cols_of <- function(g) {
    toupper(substr(g$scales$scales[[1]]$palette(seq(0, 1, length.out = 5)), 1, 7))
  }
  vir <- cols_of(imageR(obj, feat_ind = 1, sample_lab = "run",
                        palette = "viridis"))
  h0  <- cols_of(imageR(obj, feat_ind = 1, sample_lab = "run",
                        palette = "heatmap0"))

  expect_false(identical(vir, h0))
  expect_equal(h0[1], toupper(quant_palettes("heatmap0")[1]))
})

test_that("quantile_hm maps the z-score range onto the chosen palette", {
  p <- system.file("extdata", "example.raw", "section01.RDS",
                   package = "quantMSImageR")
  obj <- as(readRDS(p), "quant_MSImagingExperiment")

  hm  <- quantile_hm(obj, quant_val = 0.5, palette = "heatmap2")
  got <- substr(hm@matrix_color_mapping@col_fun(c(-1, 0, 1)), 1, 7)
  expect_equal(toupper(got), toupper(quant_palettes("heatmap2")[c(1, 3, 5)]))
})

test_that("the shipped config template carries a colours block", {
  tpl <- yaml::read_yaml(system.file("config_template.yaml",
                                     package = "quantMSImageR"))
  expect_setequal(names(tpl$colours),
                  c("ion_image", "heatmap", "group", "feature"))
  expect_equal(tpl$colours$ion_image, "heatmap0")
  # Every default named in the template must be a palette we actually ship.
  expect_true(all(unlist(tpl$colours) %in% names(quant_palettes())))
})
