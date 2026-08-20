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
#' @include setClasses.R
#'
#' @param name Acquisition name (without the `.raw` suffix).
#' @param folder Path to the directory that contains the `<name>.raw/` folder.
#' @param lib_ion_path Full path to the ion-library CSV. Must contain columns
#'   `transition_id`, `precursor_mz`, `product_mz`, `collision_eV`, `cone_V`,
#'   `Polarity`, `Type` (`"Analyte"` for in-tissue analytes, `"IS"` for
#'   internal standards). Additional metadata columns are preserved by
#'   downstream code via [buildFeatureMeta()].
#'
#'   Measured transitions are matched to it on precursor and product m/z within
#'   `mz_tolerance`; a transition the library does not describe is kept, named
#'   by its instrument index and typed `"Unknown"`.
#' @param type_header Character. Name of the ion-library column holding the
#'   feature type -- the values that mark internal standards versus analytes
#'   (default `"Type"`). Its contents become `fData()$feature_type`, which is what
#'   [int2response()] matches `IS_name` against.
#' @param is_norm_header Character. Name of the optional ion-library column that
#'   maps each analyte to the internal standard normalising it (default
#'   `"IS_norm"`). Carried through to `fData()` when present, so that
#'   [int2response()] can use it when a panel has more than one standard.
#'   `"None"` (or `NULL`), or simply omitting the column, is correct for a
#'   single-standard panel.
#' @param mz_tolerance Numeric. Half-width, in Da, within which a measured
#'   precursor and product are taken to be the library's. **Both** must agree;
#'   a row matching on precursor alone is not a match.
#'
#'   The default `0.4` is nominal-mass matching, because that is what the
#'   instrument does: a triple quadrupole running MRM selects Q1 and Q3 at unit
#'   resolution, so a product ion written as `308.1` and one written as `308.3`
#'   are the same measurement and no acquisition can separate them. Matching
#'   more tightly than the instrument resolves invents a distinction that is not
#'   in the data, and makes annotation depend on how many decimal places
#'   somebody typed into the library.
#'
#'   The consequence is that genuinely isobaric compounds now collide -- which
#'   is correct, and is what `ambiguity = "combine"` is for. Tighten it only for
#'   a high-resolution method where Q3 really is selective.
#' @param ambiguity Character. What to do when a measured transition matches
#'   more than one library entry within `mz_tolerance`, which m/z alone cannot
#'   resolve. Isomers sharing a nominal precursor and product are the usual
#'   cause, and at unit resolution they are indistinguishable by definition.
#'
#'   `"combine"` (default) names the feature for every entry it matched, joined
#'   with `" || "` -- `"LTC4 || 14_15-LTC4"`. The joined name says the
#'   measurement is one of these, which is what the data supports; naming it for
#'   one of them asserts more than was measured. Annotation columns come from
#'   the closest entry. `"error"` refuses, which is right when a collision means
#'   the library is wrong rather than the chemistry ambiguous; `"warn"` takes
#'   the closest and says so; `"nearest"` takes the closest quietly.
#'
#'   Two library entries for the *same* compound are redundancy, not ambiguity.
#'   `"combine"` reports those in its message so they can be removed from the
#'   library rather than carried as a joined name forever.
#' @param overwrite Logical. When `TRUE` (default) the raw text files are
#'   re-parsed; when `FALSE` and a cached `MSImagingExperiment.rds` exists
#'   inside the `.raw` folder, it is returned instead.
#'
#' @return An `MSImagingExperiment` with feature metadata joined to the ion
#'   library, ready to pass into [generateTxtImages()] or
#'   [selectTissuePixels()].
#'
#' @examples
#' # Return a previously parsed acquisition from its cached .rds
#' folder <- system.file("extdata", package = "quantMSImageR")
#' lib <- system.file("extdata", "ion_library_pos04.csv",
#'                    package = "quantMSImageR")
#' obj <- readMRM("pos04_test", folder = folder, lib_ion_path = lib,
#'                 overwrite = FALSE)
#'
#' @family acquisition
#' @export
readMRM <- function(name, folder, lib_ion_path, overwrite = TRUE,
                     type_header = "Type", is_norm_header = "IS_norm",
                     mz_tolerance = 0.4,
                     ambiguity = c("combine", "error", "warn", "nearest")) {

  ambiguity <- match.arg(ambiguity)

  # set Imaging folder
  imaging_folder <- sprintf("%s/%s.raw/imaging", folder, name)

  # Set rds filename
  rds_fn <- sprintf("%s/%s.raw/MSImagingExperiment.rds", folder, name)

  if (file.exists(rds_fn) && !overwrite) {
    return(readRDS(rds_fn))
  }

  # Read in ion library
  ion_lib <- read.csv(file = lib_ion_path, header = TRUE)

  # Find experimental parameters. Waters writes _extern.inf in Windows-1252
  # (non-ASCII bytes like 0xB0 for the degree symbol), so on a UTF-8 locale
  # readLines() emits "invalid UTF-8" warnings on the temperature lines --
  # harmless for the keyword scans below, but noisy. Decode as latin1.
  inf_file <- sprintf("%s/%s.raw/_extern.inf", folder, name)
  .con     <- file(inf_file, open = "r", encoding = "latin1")
  on.exit(close(.con), add = TRUE)
  lines    <- readLines(.con, warn = FALSE)
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
  if (length(analyte_fns) == 0) {
    # No raw text to (re)parse. Fall back to the cached object if present so
    # bundled acquisitions that ship only the .rds (no imaging/ folder) still
    # load, rather than erroring.
    if (file.exists(rds_fn)) {
      message("readMRM: no analyte text files in ", imaging_folder,
              "; returning cached MSImagingExperiment.rds.")
      return(readRDS(rds_fn))
    }
    stop("readMRM: no analyte text files found in ", imaging_folder)
  }

  result <- purrr::map(analyte_fns, function(analyte_fn) {
    temp_analyte <- read.table(analyte_fn, fill = TRUE, sep = "\t",
                                header = FALSE, blank.lines.skip = TRUE)[-1, ]

    # Extract transitions (first 3 rows: transition_id, precursor, product)
    temp_transitions <- t(temp_analyte[seq_len(3), ]) |>
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
      dplyr::filter(!dplyr::row_number() %in% seq_len(3)) |>
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
    stop("readMRM: merged transitions from split analyte files are inconsistent ",
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

  # Annotate the measured transitions from the ion library.
  ion_lib <- .join_ion_library(transitions, ion_lib, polarity = polarity,
                               type_header = type_header,
                               mz_tolerance = mz_tolerance,
                               ambiguity = ambiguity)

  # intensity data -- pull each transition's column from the wide analyte_df
  trans_cols <- match(sprintf("transition_%s", ion_lib$transition_id_int),
                      colnames(analyte_df))
  idata <- t(as.matrix(analyte_df[, trans_cols, drop = FALSE]))

  # feature metadata
  .fcols <- list(
    mz           = ion_lib$new_transition_int,
    feature_type = ion_lib[[type_header]],
    precursor_mz = ion_lib$precursor_mz,
    product_mz   = ion_lib$product_mz,
    name         = ion_lib$transition_id_name
  )
  # The analyte-to-standard map only matters for panels with several internal
  # standards, so it is carried when the library defines it and omitted
  # otherwise rather than being required of every library.
  .map_col <- .none(is_norm_header)
  if (!is.null(.map_col) && .map_col %in% names(ion_lib))
    .fcols[[.map_col]] <- as.character(ion_lib[[.map_col]])

  fdata <- do.call(MassDataFrame, .fcols)

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
