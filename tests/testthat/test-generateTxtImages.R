require(testthat)
require(quantMSImageR)


# Paths to the bundled pos04 test acquisition
extdata      <- system.file("extdata", package = "quantMSImageR")
pos04_folder <- extdata                          # data_path: contains pos04_test.raw/
pos04_name   <- "pos04_test"                     # acquisition name (folder = pos04_test.raw/)
lib_ion_path <- file.path(extdata, "ion_library_pos04.csv")

test_that("generateTxtImages returns the expected named list", {
  result <- generateTxtImages(
    fns          = pos04_name,
    data_path    = pos04_folder,
    image_dir    = tempdir(),
    lib_ion_path = lib_ion_path,
    output_txt   = FALSE,
    average_method = "median"
  )

  expect_true(all(c("combined", "combined_snr", "combined_snr_list",
                    "combined_NAbackground") %in% names(result)))
  expect_s4_class(result$combined,              "quant_MSImagingExperiment")
  expect_s4_class(result$combined_snr,          "quant_MSImagingExperiment")
  expect_s4_class(result$combined_NAbackground, "quant_MSImagingExperiment")
  expect_false("combined_FC" %in% names(result))   # tissue fold-change removed
})

test_that("output_txt = FALSE writes no files", {
  tmp <- file.path(tempdir(), paste0("qmsi_test_", as.integer(Sys.time())))
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  generateTxtImages(
    fns          = pos04_name,
    data_path    = pos04_folder,
    image_dir    = tmp,
    lib_ion_path = lib_ion_path,
    output_txt   = FALSE
  )

  expect_equal(length(list.files(tmp, recursive = TRUE)), 0L)
})

test_that("output_txt = TRUE creates expected subdirectories", {
  tmp <- file.path(tempdir(), paste0("qmsi_txt_", as.integer(Sys.time())))
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  result <- generateTxtImages(
    fns          = pos04_name,
    data_path    = pos04_folder,
    image_dir    = tmp,
    lib_ion_path = lib_ion_path,
    snr_thresh   = 1,
    output_txt   = TRUE
  )

  acq_dir <- file.path(tmp, pos04_name)
  expect_true(dir.exists(acq_dir))
  subdirs <- list.dirs(acq_dir, recursive = FALSE, full.names = FALSE)
  expect_true(any(grepl("intensity_SNRfiltered", subdirs)))
  expect_true(any(grepl("intensity_raw",         subdirs)))
})

test_that("combined_snr has fewer non-NA pixels than combined (SNR filter applied)", {
  result <- generateTxtImages(
    fns          = pos04_name,
    data_path    = pos04_folder,
    image_dir    = tempdir(),
    lib_ion_path = lib_ion_path,
    snr_thresh   = 3,
    output_txt   = FALSE
  )

  n_non_na_raw <- sum(!is.na(spectraData(result$combined)[["intensity"]]))
  n_non_na_snr <- sum(!is.na(spectraData(result$combined_snr)[["intensity"]]))

  expect_lte(n_non_na_snr, n_non_na_raw)
})

test_that("invalid average_method is rejected early", {
  expect_error(
    generateTxtImages(
      fns = "dummy", data_path = ".", image_dir = tempdir(),
      lib_ion_path = "dummy.csv", average_method = "geometric"
    ),
    regexp = "'arg' should be one of"
  )
})

test_that("regions of interest are off unless asked for", {
  result <- suppressMessages(generateTxtImages(
    fns = pos04_name, data_path = pos04_folder, image_dir = tempdir(),
    lib_ion_path = lib_ion_path, snr_thresh = 1, output_txt = FALSE))
  expect_false(result$rois)
  expect_false("roi_label" %in% names(pData(result$combined)))
})

test_that("roi_labels.csv is attached by pixel position", {
  # A copy of the bundled acquisition with ten tissue pixels labelled airway_01.
  tmp <- tempfile("qmsi_roi_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  file.copy(file.path(extdata, "pos04_test.raw"), tmp, recursive = TRUE)
  raw  <- file.path(tmp, "pos04_test.raw")
  obj  <- readRDS(file.path(raw, "MSImagingExperiment.rds"))
  mask <- read.csv(file.path(raw, "tissue_pixels.csv"))   # legacy: row order
  tiss <- which(as.logical(mask$tissue_pixels))[1:10]
  write.csv(data.frame(x = pData(obj)$x[tiss], y = pData(obj)$y[tiss],
                       roi_label = "airway", roi_id = "airway_01"),
            file.path(raw, "roi_labels.csv"), row.names = FALSE)

  result <- suppressMessages(generateTxtImages(
    fns = pos04_name, data_path = tmp, image_dir = tmp,
    lib_ion_path = lib_ion_path, snr_thresh = 1, output_txt = FALSE,
    rois = "auto"))
  pd <- pData(result$combined)
  bg <- as.character(pd$sample_name) == "background_pixels"

  expect_true(result$rois)
  expect_equal(sum(pd$roi_id == "airway_01", na.rm = TRUE), 10L)
  expect_true(all(is.na(pd$roi_label[bg])))              # background: no region
  expect_setequal(unique(pd$roi_label[!bg]), c("airway", "unassigned"))
  # The SNR-filtered objects, which the report reads, carry them too.
  expect_true("roi_label" %in% names(pData(result$combined_snr)))
})
