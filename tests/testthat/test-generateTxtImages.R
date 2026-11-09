require(testthat)
require(quantMSImageR)


# Paths to the bundled pos04 test acquisition
extdata      <- system.file("extdata", package = "quantMSImageR")
pos04_folder <- extdata                          # data_path: contains pos04_test.raw/
pos04_name   <- "pos04_test"                     # acquisition name (folder = pos04_test.raw/)
lib_ion_path <- file.path(extdata, "ion_library_pos04.csv")

# pos04_test has one feature whose background is mostly zero, so every SNR run
# on it carries the study-level imputed-noise warning. Tests about something
# else mute exactly that warning; it has its own test below.
quiet_imputed <- function(expr)
  withCallingHandlers(expr, quantMSImageR_imputed_noise_study = function(w)
    invokeRestart("muffleWarning"))

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

  result <- quiet_imputed(generateTxtImages(
    fns          = pos04_name,
    data_path    = pos04_folder,
    image_dir    = tmp,
    lib_ion_path = lib_ion_path,
    snr_thresh   = 1,
    output_txt   = TRUE
  ))

  acq_dir <- file.path(tmp, pos04_name)
  expect_true(dir.exists(acq_dir))
  subdirs <- list.dirs(acq_dir, recursive = FALSE, full.names = FALSE)
  expect_true(any(grepl("intensity_SNRfiltered", subdirs)))
  expect_true(any(grepl("intensity_raw",         subdirs)))
})

test_that("combined_snr has fewer non-NA pixels than combined (SNR filter applied)", {
  result <- quiet_imputed(generateTxtImages(
    fns          = pos04_name,
    data_path    = pos04_folder,
    image_dir    = tempdir(),
    lib_ion_path = lib_ion_path,
    snr_thresh   = 3,
    output_txt   = FALSE
  ))

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

test_that("thresh, perc and rot_clockwise take one valid value", {
  # Checked before anything is read, so dummy paths are enough.
  call_with <- function(...)
    generateTxtImages(fns = "dummy", data_path = ".", image_dir = tempdir(),
                      lib_ion_path = "dummy.csv", ...)

  # 180 is degrees; rot90() would read it as 180 %% 4 = 0, no rotation.
  expect_error(call_with(rot_clockwise = 180), "quarter turns")
  expect_error(call_with(rot_clockwise = 1.5), "rot_clockwise")
  expect_error(call_with(thresh = c(20, 25)), "one number")
  expect_error(call_with(perc = 101), "one number")
  expect_error(call_with(perc = "high"), "perc")
})

test_that("the imputed-noise warning is reported once for the study", {
  # int2SNR() warns once per SNR threshold; the study-level summary replaces
  # those with a single warning naming the sample by its label.
  ws <- list()
  withCallingHandlers(
    suppressMessages(generateTxtImages(
      fns = pos04_name, data_path = pos04_folder, image_dir = tempdir(),
      lib_ion_path = lib_ion_path, snr_thresh = c(1, 3), output_txt = FALSE)),
    warning = function(w) {
      ws[[length(ws) + 1L]] <<- w
      invokeRestart("muffleWarning")
    })

  study <- Filter(function(w) inherits(w, "quantMSImageR_imputed_noise_study"), ws)
  expect_length(study, 1L)
  expect_match(conditionMessage(study[[1]]), "pos04_test (1 of 7 features)",
               fixed = TRUE)
  expect_equal(study[[1]]$samples, "pos04_test")
  # The per-acquisition warnings it replaces do not leak out.
  expect_false(any(vapply(ws, inherits, logical(1),
                          what = "quantMSImageR_imputed_noise")))
})

test_that("regions of interest are off unless asked for", {
  result <- quiet_imputed(suppressMessages(generateTxtImages(
    fns = pos04_name, data_path = pos04_folder, image_dir = tempdir(),
    lib_ion_path = lib_ion_path, snr_thresh = 1, output_txt = FALSE)))
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

  result <- quiet_imputed(suppressMessages(generateTxtImages(
    fns = pos04_name, data_path = tmp, image_dir = tmp,
    lib_ion_path = lib_ion_path, snr_thresh = 1, output_txt = FALSE,
    rois = "auto")))
  pd <- pData(result$combined)
  bg <- as.character(pd$sample_name) == "background_pixels"

  expect_true(result$rois)
  expect_equal(sum(pd$roi_id == "airway_01", na.rm = TRUE), 10L)
  expect_true(all(is.na(pd$roi_label[bg])))              # background: no region
  expect_setequal(unique(pd$roi_label[!bg]), c("airway", "unassigned"))
  # The SNR-filtered objects, which the report reads, carry them too.
  expect_true("roi_label" %in% names(pData(result$combined_snr)))
})
