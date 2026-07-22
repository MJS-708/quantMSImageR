# The .txt export is the whole point of output_txt = TRUE, and nothing tested
# that files actually appear. make_txt_mat() wraps imageR() in tryCatch and
# returns an empty matrix on failure, which write_if_nonempty() then skips --
# so a broken imageR() call produced a successful run that wrote nothing.
# (It did: the call still passed the removed `overlay` argument.)

test_that("output_txt = TRUE actually writes per-feature matrices", {
  out <- tempfile("qmsi_txt_")
  dir.create(out, recursive = TRUE)

  res <- expect_no_warning(
    generate_txt_images(
      fns          = list(list(neg = "example", section = "section01",
                               label = "A")),
      data_path    = system.file("extdata", package = "quantMSImageR"),
      image_dir    = out,
      lib_ion_path = system.file("extdata", "example_ion_library.csv",
                                 package = "quantMSImageR"),
      snr_thresh   = 0,
      output_txt   = TRUE
    )
  )

  txt <- list.files(out, pattern = "\\.txt$", recursive = TRUE,
                    full.names = TRUE)
  expect_gt(length(txt), 0)

  # A written matrix must have the pixel grid, not be an empty placeholder.
  m <- as.matrix(read.table(txt[1], sep = "\t"))
  expect_gt(nrow(m), 1)
  expect_gt(ncol(m), 1)
  expect_true(any(is.finite(m)))

  expect_true("combined" %in% names(res))
})
