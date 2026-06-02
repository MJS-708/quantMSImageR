#' Read a Waters DESI-MRM acquisition into an MSImagingExperiment
#'
#' Reads one or more analyte text files from `<.raw>/imaging/`, matches
#' transitions to the ion library on rounded (precursor, product) m/z, and
#' constructs a [Cardinal::MSImagingExperiment]. Handles both the
#' single-file case (one combined `Analyte 1.txt`, the modern QuanOptimise
#' default) and the multi-file case (older exports that split transitions
#' across several `Analyte N.txt` files) transparently.
#'
#' Results are cached to `<.raw>/MSImagingExperiment.rds`; pass
#' `overwrite = TRUE` (the default) to re-parse from the raw text files.
#'
#' @import Cardinal
#' @import dplyr
#' @include setClasses.R
#'
#' @param name Acquisition name (without the `.raw` suffix).
#' @param folder Path to the directory that contains the `<name>.raw/` folder.
#' @param lib_ion_path Full path to the ion-library CSV. Must contain columns
#'   `transition_id`, `precursor_mz`, `product_mz`, `collision_eV`, `cone_V`,
#'   `Polarity`, `Type` (`"Analyte"` for in-tissue analytes, `"IS"` for
#'   internal standards). Additional metadata columns are preserved by
#'   downstream code via [build_feature_meta()].
#' @param overwrite Logical. When `TRUE` (default) the raw text files are
#'   re-parsed; when `FALSE` and a cached `MSImagingExperiment.rds` exists
#'   inside the `.raw` folder, it is returned instead.
#'
#' @return An `MSImagingExperiment` with feature metadata joined to the ion
#'   library, ready to pass into [generate_txt_images()] or
#'   [select_tissue_pixels()].
#'
#' @export
read_mrm <- function(name, folder, lib_ion_path, overwrite = TRUE) {

  # set Imaging folder
  imaging_folder <- sprintf("%s/%s.raw/imaging", folder, name)

  # Set rds filename
  rds_fn <- sprintf("%s/%s.raw/MSImagingExperiment.rds", folder, name)

  if (file.exists(rds_fn) && !overwrite) {
    return(readRDS(rds_fn))
  }

  # Read in ion library
  ion_lib <- read.csv(file = lib_ion_path, header = TRUE)

  # Find experimental parameters
  inf_file <- sprintf("%s/%s.raw/_extern.inf", folder, name)
  lines    <- readLines(inf_file)
  ystep    <- strsplit(lines[grep("DesiYStep", lines)], "\t")[[1]]
  ystep    <- as.numeric(ystep[length(ystep)]) * 1000
  polarity <- strsplit(lines[grep("Polarity", lines)], "\t")[[1]][2]
  polarity <- ifelse(polarity == "ES-", "Negative", "Positive")

  # Experimental metadata
  exp_metadata           <- CardinalIO::ImzMeta()
  exp_metadata$pixelSize <- ystep

  # ---- Read MRM imaging data --------------------------------------------
  # imaging/ may contain one analyte text file (modern QuanOptimise) or
  # several (older split exports). Iterate over all of them and combine.
  analyte_fns <- list.files(imaging_folder, full.names = TRUE)
  if (length(analyte_fns) == 0)
    stop("read_mrm: no analyte text files found in ", imaging_folder)

  result <- purrr::map(analyte_fns, function(analyte_fn) {
    temp_analyte <- read.table(analyte_fn, fill = TRUE, sep = "\t",
                                header = FALSE, blank.lines.skip = TRUE)[-1, ]

    # Extract transitions (first 3 rows: transition_id, precursor, product)
    temp_transitions <- t(temp_analyte[1:3, ]) |>
      `colnames<-`(c("transition_id", "precursor_mz", "product_mz")) |>
      na.omit() |>
      as.data.frame()

    col_heads <- c(
      "ind", "x_loci", "y_loci",
      sprintf("transition_%s", temp_transitions$transition_id),
      "sample", "pixel"
    )

    temp_analyte <- temp_analyte |>
      `colnames<-`(col_heads) |>
      dplyr::filter(!dplyr::row_number() %in% 1:3) |>
      dplyr::mutate(x = NA, y = NA, fn = basename(analyte_fn))

    # Map x_loci values to sequential integer x indices
    for (i in seq_along(sort(unique(temp_analyte$x_loci)))) {
      x_val <- unique(temp_analyte$x_loci)[i]
      temp_analyte$x[which(temp_analyte$x_loci == x_val)] <- i
    }

    list(temp_analyte = temp_analyte, temp_transitions = temp_transitions)
  })

  analyte_df <- dplyr::bind_rows(purrr::map(result, "temp_analyte")) |>
    dplyr::arrange(x_loci) |>
    dplyr::arrange(y_loci)

  transitions <- dplyr::bind_rows(purrr::map(result, "temp_transitions")) |>
    dplyr::distinct()

  # Guard against mismatched transition ordering across split files
  if (nrow(transitions) != max(transitions$transition_id))
    stop("read_mrm: merged transitions from split analyte files are inconsistent ",
         "(likely a different transition order across .txt files in ", imaging_folder, ").")

  # Map y_loci values to sequential integer y indices
  for (i in seq_along(sort(unique(analyte_df$y_loci)))) {
    y_val <- unique(analyte_df$y_loci)[i]
    analyte_df$y[which(analyte_df$y_loci == y_val)] <- i
  }

  # pixel metadata
  coord <- analyte_df |> dplyr::select(x, y)
  run   <- factor(rep(name, nrow(coord)))
  pdata <- PositionDataFrame(run = run, coord = coord)

  # Round precursors and products to 0 dp for the join (matches m/z key behaviour
  # used throughout the package — see build_feature_meta()).
  transitions <- transitions |>
    dplyr::mutate(precursor_mz = round(precursor_mz, digits = 0),
                  product_mz   = round(product_mz,   digits = 0))
  ion_lib <- ion_lib |>
    dplyr::mutate(precursor_mz = round(precursor_mz, digits = 0),
                  product_mz   = round(product_mz,   digits = 0))

  # Gather info of MRM transitions from the ion library
  ion_lib <- ion_lib |>
    subset(Polarity == polarity) |>
    dplyr::right_join(y = transitions, by = c("precursor_mz", "product_mz"),
                      suffix = c("_name", "_int")) |>
    dplyr::mutate(
      transition_id_name = ifelse(is.na(transition_id_name),
                                   transition_id_int, transition_id_name),
      Polarity = ifelse(is.na(Polarity), polarity, Polarity),
      Type     = ifelse(is.na(Type), "Unknown", Type)
    ) |>
    dplyr::arrange(transition_id_int) |>
    dplyr::group_by(precursor_mz, product_mz) |>
    dplyr::mutate(
      transition_id_name = paste0(transition_id_name, collapse = " || "),
      precursor_mz       = paste0(precursor_mz,       collapse = " || "),
      product_mz         = paste0(product_mz,         collapse = " || "),
      collision_eV       = paste0(collision_eV,       collapse = " || "),
      cone_V             = paste0(cone_V,             collapse = " || ")
    ) |>
    dplyr::distinct(transition_id_name, transition_id_int, .keep_all = TRUE) |>
    dplyr::mutate(
      duplicate_mrms     = dplyr::row_number(),
      transition_id_name = ifelse(duplicate_mrms > 1,
                                  sprintf("%s:- %s", duplicate_mrms, transition_id_name),
                                  transition_id_name)
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(precursor_mz, product_mz, transition_id_int) |>
    dplyr::mutate(new_transition_int = dplyr::row_number())

  # intensity data — pull each transition's column from the wide analyte_df
  idata <- t(sapply(ion_lib$transition_id_int, FUN = function(x) {
    analyte_df[, which(colnames(analyte_df) == sprintf("transition_%s", x))]
  }))

  # feature metadata
  fdata <- MassDataFrame(
    mz           = ion_lib$new_transition_int,
    analyte      = ion_lib$Type,
    precursor_mz = ion_lib$precursor_mz,
    product_mz   = ion_lib$product_mz,
    name         = ion_lib$transition_id_name
  )

  out <- MSImagingExperiment(
    spectraData    = idata,
    featureData    = fdata,
    pixelData      = pdata,
    experimentData = exp_metadata
  )

  featureNames(out) <- fData(out)$name

  # cache
  saveRDS(out, file = rds_fn)

  out
}
