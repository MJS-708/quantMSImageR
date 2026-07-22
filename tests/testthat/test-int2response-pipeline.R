# The YAML `is_name:` option must actually reach int2response() and switch the
# whole pipeline onto the `response` layer -- previously int2response() was
# documented as part of the workflow but unreachable from a config.

is_fns <- list(list(neg = "example", section = "section01", label = "A"))

run_with <- function(...) {
  generate_txt_images(
    fns          = is_fns,
    data_path    = system.file("extdata", package = "quantMSImageR"),
    image_dir    = tempfile(),
    lib_ion_path = system.file("extdata", "example_ion_library.csv",
                               package = "quantMSImageR"),
    snr_thresh   = 0,
    output_txt   = FALSE,
    ...
  )
}

test_that("without is_name the pipeline stays on intensity", {
  res <- suppressMessages(run_with())
  expect_false("response" %in% names(spectraData(res$combined)))
})

test_that("is_name normalises and switches the pipeline to the response layer", {
  # is_name is the ion library's Type value (fData()$analyte), not a
  # transition name -- passing a name silently normalises each feature to
  # itself instead, which is why this is asserted explicitly.
  res <- suppressMessages(run_with(is_name = "IS", remove_IS = TRUE))

  expect_true("response" %in% names(spectraData(res$combined)))
  # The internal standard itself is dropped when remove_IS = TRUE
  expect_false("IS" %in% as.character(fData(res$combined)$analyte))
  # Normalised values are finite and not simply a copy of intensity
  r <- as.numeric(spectra(res$combined, "response"))
  expect_true(any(is.finite(r)))
  expect_false(isTRUE(all.equal(r, as.numeric(spectra(res$combined, "intensity")))))
})

test_that("an unknown internal standard fails loudly", {
  # And the message must say which column is consulted, since passing a
  # transition name here is the easy mistake.
  expect_error(suppressMessages(run_with(is_name = "not_a_feature")),
               "no feature with analyte")
  expect_error(suppressMessages(run_with(is_name = "not_a_feature")),
               "Type column")
})

test_that("the config template documents the IS options", {
  tpl <- yaml::read_yaml(system.file("config_template.yaml",
                                     package = "quantMSImageR"))
  expect_true(all(c("is_name", "is_mode", "remove_IS") %in%
                    names(tpl$parameters)))
  expect_equal(tpl$parameters$is_name, "None")   # off by default
})
