require(testthat)
require(quantMSImageR)

# A two-block, two-group config: the shape block_centre is for.
blocked_cfg <- function(out = tempfile()) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  samp <- function(id, blk, grp)
    list(neg = "example", sections = list("section01"),
         run_id = id, block = blk, group = grp)
  list(
    study = "v",
    paths = list(
      data_path    = system.file("extdata", package = "quantMSImageR"),
      out_path     = out,
      image_dir    = file.path(out, "images"),
      lib_ion_path = system.file("extdata", "example_ion_library.csv",
                                 package = "quantMSImageR")),
    samples = list(samp("S1", "P1", "Ctrl"), samp("S2", "P1", "Trt"),
                   samp("S3", "P2", "Ctrl"), samp("S4", "P2", "Trt")),
    parameters = list(snr_thresh = 0, baseline_group = "Ctrl"),
    output = list(render_report = FALSE, output_txt = FALSE))
}

test_that("block_centre takes off, center or baseline", {
  cfg <- blocked_cfg()
  for (v in list("off", "center", "baseline", "centre")) {
    cfg$parameters$block_centre <- v
    expect_length(validateConfig(cfg)$errors, 0)
  }
  cfg$parameters$block_centre <- "z"
  expect_true(any(grepl("block_centre", validateConfig(cfg)$errors)))
})

test_that("off/on are YAML booleans, so the value may arrive as a logical", {
  # `block_centre: off` parses as FALSE, not "off" -- the config still means
  # what it says.
  expect_equal(quantMSImageR:::.block_mode(FALSE), "off")
  expect_equal(quantMSImageR:::.block_mode(TRUE), "center")
  expect_equal(quantMSImageR:::.block_mode(NULL), "off")
  expect_equal(quantMSImageR:::.block_mode("centre"), "center")
  expect_equal(quantMSImageR:::.block_mode("baseline"), "baseline")
  expect_equal(quantMSImageR:::.block_mode("nonsense"), "nonsense")

  cfg <- blocked_cfg()
  cfg$parameters$block_centre <- FALSE
  expect_length(validateConfig(cfg)$errors, 0)
})

test_that("correcting within a block that is not there is pointless, and said so", {
  cfg <- blocked_cfg()
  cfg$parameters$block_centre <- "center"

  # No sample declares a block.
  bare <- cfg
  bare$samples <- lapply(bare$samples, function(s) { s$block <- NULL; s })
  expect_true(any(grepl("no sample has a block",
                        validateConfig(bare)$warnings)))

  # One block holding everything: correcting within it removes the comparison.
  one <- cfg
  one$samples <- lapply(one$samples, function(s) { s$block <- "P1"; s })
  expect_true(any(grepl("same block", validateConfig(one)$warnings)))

  # The ordinary case says nothing.
  expect_false(any(grepl("block_centre", validateConfig(cfg)$warnings)))
})

test_that("baseline centring names the blocks that lack the reference arm", {
  cfg <- blocked_cfg()
  cfg$parameters$block_centre <- "baseline"
  expect_false(any(grepl("own mean instead", validateConfig(cfg)$warnings)))

  # P2 loses its control, so it would fall back rather than be centred on one.
  cfg$samples[[3]]$group <- "Trt"
  w <- validateConfig(cfg)$warnings
  expect_true(any(grepl("P2", w) & grepl("own mean instead", w)))
})

test_that("stitch_overlap takes error, first or last", {
  cfg <- blocked_cfg()
  for (v in c("error", "first", "last")) {
    cfg$parameters$stitch_overlap <- v
    expect_length(validateConfig(cfg)$errors, 0)
  }
  cfg$parameters$stitch_overlap <- "average"
  expect_true(any(grepl("stitch_overlap", validateConfig(cfg)$errors)))
})
