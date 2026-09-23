good_cfg <- function(out = tempfile()) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  list(
    study = "v",
    paths = list(
      data_path    = system.file("extdata", package = "quantMSImageR"),
      out_path     = out,
      image_dir    = file.path(out, "images"),
      lib_ion_path = system.file("extdata", "example_ion_library.csv",
                                 package = "quantMSImageR")),
    samples = list(list(neg = "example", sections = list("section01"),
                        run_id = "S1", label = "A")),
    parameters = list(snr_thresh = 0),
    output = list(render_report = FALSE, output_txt = FALSE)
  )
}

test_that("a valid config passes with no errors", {
  v <- validateConfig(good_cfg())
  expect_s3_class(v, "quantValidation")
  expect_length(v$errors, 0)
  expect_true(nrow(v$summary) > 5)
})

test_that("the shipped template is structurally valid", {
  v <- validateConfig(system.file("config_template.yaml",
                                   package = "quantMSImageR"),
                       check_paths = FALSE)
  expect_length(v$errors, 0)
})

test_that("every problem is reported, not just the first", {
  cfg <- good_cfg()
  cfg$parameters$snr_thresh     <- -1
  cfg$parameters$average_method <- "geometric"
  cfg$parameters$is_mode        <- "rolling"
  cfg$colours                   <- list(heatmap = "not_a_palette")

  v <- validateConfig(cfg)
  expect_gte(length(v$errors), 4)
  expect_true(any(grepl("snr_thresh", v$errors)))
  expect_true(any(grepl("average_method", v$errors)))
  expect_true(any(grepl("is_mode", v$errors)))
  expect_true(any(grepl("palette", v$errors)))
})

test_that("is_norm_header is only demanded when the panel has several standards", {
  lib_path <- system.file("extdata", "example_ion_library.csv",
                          package = "quantMSImageR")
  lib <- read.csv(lib_path, check.names = FALSE)

  with_lib <- function(l, pars) {
    f <- tempfile(fileext = ".csv"); write.csv(l, f, row.names = FALSE)
    cfg <- good_cfg()
    cfg$paths$lib_ion_path <- f
    cfg$parameters <- utils::modifyList(cfg$parameters, pars)
    validateConfig(cfg)
  }

  # One standard in the library: None is right, and must not be an error.
  expect_length(with_lib(lib, list(is_name = "IS",
                                   is_norm_header = "None"))$errors, 0)

  # Two standards: None can no longer resolve which one each analyte uses.
  lib2 <- lib
  lib2$Type[lib2$transition_id == "12-HHTrE"] <- "IS"
  v <- with_lib(lib2, list(is_name = "IS", is_norm_header = "None"))
  expect_true(any(grepl("is_norm_header cannot be None", v$errors)))

  # Named but absent from the library.
  lib3 <- lib2; lib3$IS_norm <- NULL
  v3 <- with_lib(lib3, list(is_name = "IS", is_norm_header = "IS_norm"))
  expect_true(any(grepl("no 'IS_norm' column", v3$errors)))

  # Present but pointing at something that is not a standard.
  lib4 <- lib2
  lib4$IS_norm[lib4$Type != "IS"] <- "9-HOTE"
  v4 <- with_lib(lib4, list(is_name = "IS", is_norm_header = "IS_norm"))
  expect_true(any(grepl("is not a feature typed", v4$errors)))
})

test_that("cell_border is checked as a colour, not as a palette", {
  cfg <- good_cfg()

  # A palette name is exactly what cell_border must NOT be checked against.
  cfg$colours <- list(cell_border = "white")
  expect_length(validateConfig(cfg)$errors, 0)
  cfg$colours <- list(cell_border = "#1A2B3C")
  expect_length(validateConfig(cfg)$errors, 0)
  cfg$colours <- list(cell_border = "none")
  expect_length(validateConfig(cfg)$errors, 0)

  cfg$colours <- list(cell_border = "chartroose")
  expect_true(any(grepl("cell_border", validateConfig(cfg)$errors)))
})

test_that("is_name is checked against the ion library's type column", {
  cfg <- good_cfg()
  cfg$parameters$is_name <- "PGE2-d4"          # a transition name, not a type
  v <- validateConfig(cfg)
  expect_true(any(grepl("type label, not a transition name", v$errors)))

  cfg$parameters$is_name <- "IS"               # the correct type value
  expect_length(validateConfig(cfg)$errors, 0)
})

test_that("calibration requires cal_type and its files", {
  cfg <- good_cfg()
  cfg$calibration <- list(enabled = TRUE, cal_acquisition = "cal_run")
  v <- validateConfig(cfg)
  expect_true(any(grepl("cal_type", v$errors)))
  expect_true(any(grepl("cal_roi_csv", v$errors)))
  expect_true(any(grepl("cal_metadata", v$errors)))
})

test_that("missing blocks and duplicate run_ids are caught", {
  v <- validateConfig(list(study = "x"))
  expect_true(any(grepl("paths", v$errors)))
  expect_true(any(grepl("samples", v$errors)))

  cfg <- good_cfg()
  cfg$samples <- list(
    list(neg = "example", sections = list("section01"), run_id = "S1", label = "A"),
    list(neg = "example", sections = list("section02"), run_id = "S1", label = "B"))
  expect_true(any(grepl("duplicate run_id", validateConfig(cfg)$errors)))
})

test_that("runStudy refuses to start on an invalid config", {
  cfg <- good_cfg()
  cfg$parameters$average_method <- "geometric"
  f <- file.path(cfg$paths$out_path, "bad.yaml")
  yaml::write_yaml(cfg, f)

  expect_error(suppressWarnings(runStudy(f)),
               "Config validation failed")
})

test_that("print method summarises without erroring", {
  v <- validateConfig(good_cfg())
  expect_output(print(v), "quantMSImageR validation")
  expect_output(print(v), "All checks passed")
})

test_that("output$colocalisation names a level the report can draw", {
  cfg <- good_cfg()

  # One level only: off, by sample, or by region. True means "sample", so a
  # config written before the level existed keeps the section it had.
  for (v in list(FALSE, TRUE, "auto", "sample", "roi", "off")) {
    cfg$output$colocalisation <- v
    expect_length(validateConfig(cfg)$errors, 0)
  }

  # Anything else would otherwise be read as off, without a word.
  cfg$output$colocalisation <- "pixels"
  expect_true(any(grepl("colocalisation", validateConfig(cfg)$errors)))
})

test_that("a sample in pieces is checked for shape", {
  piece_sample <- function(...) c(list(run_id = "S2", label = "A"), list(...))
  two <- list(list(neg = "top", pos = list("top_p1")), list(neg = "btm"))

  cfg <- good_cfg()
  cfg$samples[[2]] <- piece_sample(pieces = two)
  expect_length(validateConfig(cfg)$errors, 0)

  # One piece is not a tissue in pieces.
  cfg$samples[[2]] <- piece_sample(pieces = two[1])
  expect_true(any(grepl("two or more", validateConfig(cfg)$errors)))

  # A piece with no acquisitions of its own.
  cfg$samples[[2]] <- piece_sample(pieces = list(two[[1]], list(label = "x")))
  expect_true(any(grepl("its own pos", validateConfig(cfg)$errors)))

  # Acquisitions beside pieces: which area do they cover?
  cfg$samples[[2]] <- piece_sample(pieces = two, neg = "other")
  expect_true(any(grepl("under each piece", validateConfig(cfg)$errors)))
})

test_that("the roi block takes only documented values", {
  cfg <- good_cfg()
  cfg$roi <- list(enabled = "auto", include_unassigned = FALSE,
                  compare = "between_groups", unit = "sample")
  expect_length(validateConfig(cfg)$errors, 0)

  for (cmp in c("between_rois", "both")) {
    cfg$roi$compare <- cmp
    expect_length(validateConfig(cfg)$errors, 0)
  }
  cfg$roi$compare <- "between_groups"
  cfg$roi$unit    <- "roi"
  expect_length(validateConfig(cfg)$errors, 0)

  # unit takes one of the two, never both: a comparison has a single unit.
  bad <- list(enabled = "yes please", include_unassigned = "no",
              compare = "groups", unit = "both")
  for (k in names(bad)) {
    cfg2 <- cfg
    cfg2$roi[[k]] <- bad[[k]]
    expect_true(any(grepl(paste0("roi$", k), validateConfig(cfg2)$errors,
                          fixed = TRUE)),
                info = k)
  }
})
