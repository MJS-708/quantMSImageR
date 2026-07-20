# =============================================================================
# generate_cal_data.R
# Run once (from the package root) to create the synthetic calibration data
# bundled in inst/extdata/cal_example.raw/.
#
#   source("inst/generate_cal_data.R")
#
# Produces:
#   inst/extdata/cal_example.raw/cal_MSI.RDS            (MSImagingExperiment)
#   inst/extdata/cal_example.raw/calibration_metadata.csv
#
# The object carries the vocabulary the calibration functions require:
#   pData$sample_type in {"Cal","Tissue","Noise"} and pData$identifier per
#   calibration spot, fData$name matching cal_metadata$lipid, an intensity
#   slot scaling with amount, and experimentData$pixelSize. It powers the
#   runnable examples/tests for summarise_cal_levels -> create_cal_curve ->
#   int2conc and the vignette's quantification section.
# =============================================================================

library(Cardinal)

set.seed(7)

# ----- grid ------------------------------------------------------------------
nx <- 20L; ny <- 20L
x  <- rep(seq_len(nx), each = ny)
y  <- rep(seq_len(ny), times = nx)
n_pix <- nx * ny

# ----- features (3 lipid standards) ------------------------------------------
lipids <- c("SM 16:0", "PC 34:1", "LPC 16:0")
n_feat <- length(lipids)
slope     <- c(300, 500, 200)   # intensity per pg/pixel, per lipid
intercept <- c(50, 80, 30)

# ----- calibration design: 5 levels x 3 reps, 2x2 spots ----------------------
levels_amt <- c(L1 = 10, L2 = 20, L3 = 40, L4 = 80, L5 = 160)  # pg per spot
reps <- 1:3
pixels_per_spot <- 4L

sample_type <- rep("Noise", n_pix)
identifier  <- rep(NA_character_, n_pix)

for (i in seq_along(levels_amt)) {
  for (j in reps) {
    x0 <- 2 + (i - 1) * 3       # spot columns 2, 5, 8, 11, 14
    y0 <- 2 + (j - 1) * 3       # spot rows 2, 5, 8
    inds <- which(x %in% c(x0, x0 + 1) & y %in% c(y0, y0 + 1))
    sample_type[inds] <- "Cal"
    identifier[inds]  <- sprintf("%s_r%d", names(levels_amt)[i], j)
  }
}

# tissue block (x 15-20, y 12-17)
tinds <- which(x >= 15 & x <= 20 & y >= 12 & y <= 17)
sample_type[tinds] <- "Tissue"
ninds <- which(sample_type == "Noise")

# ----- intensity matrix ------------------------------------------------------
imat <- matrix(0, nrow = n_feat, ncol = n_pix)
for (f in seq_len(n_feat)) {
  for (i in seq_along(levels_amt)) {
    for (j in reps) {
      id   <- sprintf("%s_r%d", names(levels_amt)[i], j)
      inds <- which(identifier == id)
      pgpp <- levels_amt[i] / pixels_per_spot
      mu   <- slope[f] * pgpp + intercept[f]
      imat[f, inds] <- pmax(0, rnorm(length(inds), mu, mu * 0.05))
    }
  }
  # tissue pixels: a mid-range signal; noise pixels: low
  imat[f, tinds] <- pmax(0, rnorm(length(tinds), slope[f] * 5 + intercept[f], 50))
  imat[f, ninds] <- pmax(0, rnorm(length(ninds), 20, 10))
}

# ----- assemble MSImagingExperiment -----------------------------------------
pdata <- PositionDataFrame(
  run         = factor(rep("cal_slide", n_pix)),
  coord       = data.frame(x = x, y = y),
  sample_type = factor(sample_type, levels = c("Cal", "Tissue", "Noise")),
  identifier  = identifier
)

fdata <- MassDataFrame(
  mz           = seq_len(n_feat),
  name         = lipids,
  analyte      = rep("Analyte", n_feat),
  precursor_mz = as.character(c(703, 760, 496)),
  product_mz   = as.character(c(184, 184, 184))
)

exp_meta <- CardinalIO::ImzMeta()
exp_meta$pixelSize <- 100    # micrometres

obj <- MSImagingExperiment(
  spectraData    = SimpleList(intensity = imat),
  featureData    = fdata,
  pixelData      = pdata,
  experimentData = exp_meta
)
featureNames(obj) <- lipids

# ----- write outputs ---------------------------------------------------------
out_dir <- file.path("inst", "extdata", "cal_example.raw")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(obj, file.path(out_dir, "cal_MSI.RDS"))

# one row per (identifier, lipid): summarise_cal_levels joins on identifier and
# expands to per-lipid rows.
meta <- do.call(rbind, lapply(seq_along(levels_amt), function(i)
  do.call(rbind, lapply(reps, function(j)
    data.frame(
      identifier = sprintf("%s_r%d", names(levels_amt)[i], j),
      lipid      = lipids,
      amount_pg  = unname(levels_amt[i]),
      level      = names(levels_amt)[i],
      stringsAsFactors = FALSE
    )))))
write.csv(meta, file.path(out_dir, "calibration_metadata.csv"), row.names = FALSE)

message("Wrote ", file.path(out_dir, "cal_MSI.RDS"),
        " and calibration_metadata.csv")
