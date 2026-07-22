# Internal-standard selection. A named standard that matches nothing used to
# fall back to within-feature normalisation with only a message, and a name that
# matched several standards produced a matrix denominator that recycled.

make_obj <- function(types, names = NULL, map = NULL) {
  n <- length(types)
  if (is.null(names)) names <- paste0("f", seq_len(n))
  cols <- list(mz = seq_len(n), name = names, feature_type = types)
  if (!is.null(map)) cols$IS_norm <- map
  fdata <- do.call(MassDataFrame, cols)
  pdata <- PositionDataFrame(run   = rep("s1", 4),
                             coord = expand.grid(x = 1:2, y = 1:2))
  ints  <- matrix(seq_len(n * 4) * 10, nrow = n)
  as(MSImagingExperiment(spectraData = ints, featureData = fdata,
                          pixelData = pdata),
     "quant_MSImagingExperiment")
}

test_that("a standard that matches no feature is an error", {
  # Silently switching to within-feature normalisation changes the analytical
  # method, and every downstream value would carry that change unannounced.
  expect_error(
    int2response(make_obj(c("Analyte", "Analyte")), IS_name = "IS"),
    regexp = "no feature is typed 'IS'"
  )
})

test_that("several matching standards require an explicit mapping", {
  obj <- make_obj(c("IS", "IS", "Analyte", "Analyte"))
  expect_error(int2response(obj, IS_name = "IS"),
               regexp = "each analyte must name the standard it uses")

  # "None" is a legitimate setting, so the message must ask for a header
  # rather than for a column called "None".
  expect_error(int2response(obj, IS_name = "IS", is_norm_header = "None"),
               regexp = "Set is_norm_header to an ion-library column")
})

test_that("is_norm_header = None is correct for a single shared standard", {
  # One standard normalises everything, so there is nothing to map and the
  # column should not be demanded.
  obj <- make_obj(types = c("IS", "Analyte", "Analyte"),
                  names = c("IS_a", "an_a", "an_b"))

  v <- spectraData(obj)[["intensity"]]
  for (h in list("None", NULL, NA, "")) {
    out <- int2response(obj, IS_name = "IS", is_norm_header = h,
                        mode = "pixel", remove_IS = FALSE)
    expect_equal(as.numeric(spectraData(out)[["response"]][2, ]),
                 as.numeric(v[2, ] / v[1, ]))
  }
})

test_that("each analyte is divided by the standard it is mapped to", {
  obj <- make_obj(types = c("IS", "IS", "Analyte", "Analyte"),
                  names = c("IS_a", "IS_b", "an_a", "an_b"),
                  map   = c(NA, NA, "IS_a", "IS_b"))

  out <- int2response(obj, IS_name = "IS", mode = "pixel", remove_IS = FALSE)
  v   <- spectraData(obj)[["intensity"]]
  r   <- spectraData(out)[["response"]]

  expect_equal(as.numeric(r[3, ]), as.numeric(v[3, ] / v[1, ]))
  expect_equal(as.numeric(r[4, ]), as.numeric(v[4, ] / v[2, ]))
})

test_that("a mapping that points nowhere useful is rejected", {
  expect_error(
    int2response(make_obj(c("IS", "IS", "Analyte", "Analyte"),
                          names = c("IS_a", "IS_b", "an_a", "an_b"),
                          map   = c(NA, NA, "IS_a", "no_such_feature")),
                 IS_name = "IS"),
    regexp = "does not name a feature"
  )

  expect_error(
    int2response(make_obj(c("IS", "IS", "Analyte", "Analyte"),
                          names = c("IS_a", "IS_b", "an_a", "an_b"),
                          map   = c(NA, NA, "IS_a", "an_a")),
                 IS_name = "IS"),
    regexp = "not typed 'IS'"
  )
})

test_that("all standards are dropped when remove_IS is TRUE", {
  obj <- make_obj(types = c("IS", "IS", "Analyte", "Analyte"),
                  names = c("IS_a", "IS_b", "an_a", "an_b"),
                  map   = c(NA, NA, "IS_a", "IS_b"))
  out <- int2response(obj, IS_name = "IS", remove_IS = TRUE)
  expect_equal(as.character(fData(out)$name), c("an_a", "an_b"))
})

test_that("within_feature is stated rather than spelled 'None'", {
  obj <- make_obj(c("Analyte", "Analyte"))

  new <- int2response(obj, normalisation = "within_feature", mode = "pixel")
  old <- int2response(obj, IS_name = "None", mode = "pixel")   # deprecated form

  expect_equal(spectraData(new)[["response"]], spectraData(old)[["response"]])
  # Dividing a feature by itself per pixel is 1 everywhere.
  expect_true(all(spectraData(new)[["response"]] == 1))

  expect_error(int2response(obj, normalisation = "internal_standard"),
               regexp = "needs `IS_name`")
})

test_that("the shipped ion library defines IS_norm for every analyte", {
  # The mapping lives in the ion library, not in a calibration-only file,
  # because IS normalisation runs whether or not calibration is enabled.
  lib <- read.csv(system.file("extdata", "example_ion_library.csv",
                              package = "quantMSImageR"),
                  check.names = FALSE)

  expect_true("IS_norm" %in% names(lib))

  is_rows <- lib$Type == "IS"
  # Every analyte names a standard, and that standard is a real IS row.
  expect_true(all(nzchar(lib$IS_norm[!is_rows])))
  expect_true(all(lib$IS_norm[!is_rows] %in% lib$transition_id[is_rows]))
  # A standard has no standard of its own.
  expect_true(all(!nzchar(lib$IS_norm[is_rows])))
})

test_that("IS_norm reaches the analyte-to-standard lookup from the library", {
  # The library column is what a real study supplies, so drive the mapping
  # from the shipped CSV rather than from a hand-built fData.
  lib <- read.csv(system.file("extdata", "example_ion_library.csv",
                              package = "quantMSImageR"),
                  check.names = FALSE)

  obj <- make_obj(types = lib$Type, names = lib$transition_id,
                  map   = lib$IS_norm)

  out <- int2response(obj, IS_name = "IS", mode = "pixel", remove_IS = TRUE)
  # The one standard is consumed; the seven analytes survive.
  expect_equal(as.character(fData(out)$name),
               lib$transition_id[lib$Type != "IS"])

  # An optional mapping must not become a requirement: a library without the
  # column still works while there is only one standard.
  no_map <- make_obj(types = lib$Type, names = lib$transition_id)
  expect_s4_class(int2response(no_map, IS_name = "IS", mode = "pixel"),
                  "quant_MSImagingExperiment")
})

test_that("the legacy analyte column is still read as the feature type", {
  # Objects built by earlier versions store the type in fData()$analyte.
  fdata <- MassDataFrame(mz = 1:2, name = c("IS_a", "an_a"),
                         analyte = c("IS", "Analyte"))
  pdata <- PositionDataFrame(run = rep("s1", 4),
                             coord = expand.grid(x = 1:2, y = 1:2))
  obj <- as(MSImagingExperiment(spectraData = matrix(1:8 * 10, nrow = 2),
                                 featureData = fdata, pixelData = pdata),
            "quant_MSImagingExperiment")

  out <- int2response(obj, IS_name = "IS", mode = "pixel", remove_IS = TRUE)
  expect_equal(as.character(fData(out)$name), "an_a")
})
