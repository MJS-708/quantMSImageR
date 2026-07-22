# =============================================================================
# generate_example_data.R
# Run once (from the package root) to create the synthetic RDS sections
# stored in inst/extdata/example.raw/.
#
#   source("inst/generate_example_data.R")
#
# Requires: Cardinal (>= 3.6)
# =============================================================================

library(Cardinal)

set.seed(42)

# -----------------------------------------------------------------------------
# Grid geometry
# -----------------------------------------------------------------------------
nx <- 20L; ny <- 20L
x  <- rep(seq_len(nx), each = ny)
y  <- rep(seq_len(ny), times = nx)
n_pix <- nx * ny

# Tissue mask helpers — circle or square. Each section is shaped according to
# its `shape` field so users can visually distinguish SampleA (circle) from
# SampleB (square) in the rendered ion images.
cx <- (nx + 1) / 2; cy <- (ny + 1) / 2

make_mask <- function(shape) {
  switch(shape,
    circle = ((x - cx)^2 + (y - cy)^2) <= 7^2,
    square = abs(x - cx) <= 6 & abs(y - cy) <= 6,
    stop("Unknown shape: ", shape)
  )
}

dist_from_centre <- sqrt((x - cx)^2 + (y - cy)^2)
max_r            <- 7

# -----------------------------------------------------------------------------
# Feature definitions  (match example_ion_library.csv)
# -----------------------------------------------------------------------------
features <- data.frame(
  name         = c("12-HHTrE",
                   "12-HOPE",
                   "9-HOTE",
                   "13-HOTE",
                   "15-HOTE",
                   "8-isoPGE2 || PGE2 || PGD2",
                   "8-iso-PGF-2a || 11b-PGF-2a || PGF-2a",
                   "PGE2-d4"),
  precursor_mz = c(279L, 289L, 291L, 291L, 291L, 351L, 353L, 355L),
  product_mz   = c(179L, 135L, 169L, 193L, 221L, 271L, 193L, 275L),
  feature_type = c(rep("Analyte", 7L), "IS"),
  # Analyte -> the standard that normalises it. One standard here, so
  # int2response() never has to consult it; present so the bundled objects
  # match the shipped ion library column for column.
  IS_norm      = c(rep("PGE2-d4", 7L), ""),
  stringsAsFactors = FALSE
)
n_feat <- nrow(features)

# -----------------------------------------------------------------------------
# Helper: build one section RDS
# -----------------------------------------------------------------------------
make_section <- function(section_name, seed_offset = 0L,
                         tissue_mean_scale = 1, shape = "circle") {
  set.seed(42L + seed_offset)

  tissue_mask    <- make_mask(shape)
  spatial_weight <- pmax(0, 1 - dist_from_centre / (max_r + 2))

  imat <- matrix(NA_real_, nrow = n_feat, ncol = n_pix)

  for (f in seq_len(n_feat)) {
    if (features$feature_type[f] == "IS") {
      # IS: broadly uniform in tissue, low background
      imat[f, tissue_mask]  <- rnorm(sum(tissue_mask),  mean = 5000, sd = 400)
      imat[f, !tissue_mask] <- rnorm(sum(!tissue_mask), mean = 250,  sd = 80)
    } else {
      base <- runif(1L, 1000, 8000) * tissue_mean_scale
      # Analytes: concentrated in tissue with spatial gradient
      mu <- base * spatial_weight * ifelse(tissue_mask, 4, 0.4)
      imat[f, ] <- pmax(0, rnorm(n_pix, mean = mu, sd = base * 0.15))
    }
  }

  sample_name <- factor(
    ifelse(tissue_mask, "tissue_pixels", "background_pixels"),
    levels = c("tissue_pixels", "background_pixels")
  )

  pdata <- PositionDataFrame(
    run         = factor(rep(section_name, n_pix)),
    coord       = data.frame(x = x, y = y),
    sample_name = sample_name
  )

  fdata <- MassDataFrame(
    mz           = seq_len(n_feat),
    feature_type = features$feature_type,
    precursor_mz = as.character(features$precursor_mz),
    product_mz   = as.character(features$product_mz),
    name         = features$name,
    IS_norm      = features$IS_norm
  )

  # Acquisition metadata: pixelSize (micrometres) is required by int2conc() to
  # convert pg/pixel into pg/mm2, so the sections carry it like real data.
  exp_meta <- CardinalIO::ImzMeta()
  exp_meta$pixelSize <- 100

  obj <- MSImagingExperiment(
    spectraData    = SimpleList(intensity = imat),
    featureData    = fdata,
    pixelData      = pdata,
    experimentData = exp_meta
  )
  featureNames(obj) <- features$name
  obj
}

# -----------------------------------------------------------------------------
# Generate sections: 2 samples x 3 sections
#   SampleA (sections 01-03) — higher oxylipin levels
#   SampleB (sections 04-06) — ~60 % of SampleA levels
# -----------------------------------------------------------------------------
out_dir <- file.path("inst", "extdata", "example.raw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sections <- list(
  list(name = "section01", offset =  10L, scale = 1.0, shape = "circle"),
  list(name = "section02", offset =  20L, scale = 1.0, shape = "circle"),
  list(name = "section03", offset =  30L, scale = 1.0, shape = "circle"),
  list(name = "section04", offset = 110L, scale = 0.6, shape = "square"),
  list(name = "section05", offset = 120L, scale = 0.6, shape = "square"),
  list(name = "section06", offset = 130L, scale = 0.6, shape = "square")
)

for (s in sections) {
  obj <- make_section(s$name, seed_offset = s$offset,
                      tissue_mean_scale = s$scale, shape = s$shape)
  path <- file.path(out_dir, paste0(s$name, ".RDS"))
  saveRDS(obj, path)
  message("Saved ", path, " (", s$shape, ")")
}

message("Done — example data written to ", normalizePath(out_dir))
