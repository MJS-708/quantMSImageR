#' Build per-feature metadata by joining an ion library on m/z
#'
#' Joins each row of `fData(combined)` to the ion library on the rounded
#' (`precursor_mz`, `product_mz`) pair. Robust to feature renames in the YAML
#' (`features.rename`) because m/z is stable across name changes.
#'
#' Used by [`run_study`] (via `inst/run_study.R`) and [`run_example`] to build
#' a per-feature metadata frame the heatmap report can consume — e.g. for
#' splitting / colour-bar annotation by an ion-library column such as `Met-1`.
#'
#' @param combined An MSI object whose `fData` has columns `name`,
#'   `precursor_mz` and `product_mz`. Multi-transition features stored with
#'   `" || "`-joined m/z strings are matched on the first token.
#' @param ion_lib_meta A `data.frame` of ion library metadata (must contain
#'   `precursor_mz` and `product_mz`). Pass `NULL` to short-circuit and return
#'   `NULL`.
#' @param verbose Logical. Emit a `message()` reporting how many features
#'   matched. Default `TRUE`.
#'
#' @return A `data.frame` of metadata, one row per feature in `fData(combined)`,
#'   with row names set to the feature name and a `name` column appended.
#'   Unmatched features become a row of NAs. Returns `NULL` when
#'   `ion_lib_meta` is `NULL`.
#'
#' @export
build_feature_meta <- function(combined, ion_lib_meta, verbose = TRUE) {
  if (is.null(ion_lib_meta)) return(NULL)

  # Cardinal's `fData` returns a MassDataFrame; its as.data.frame method
  # warns on extra args, so we don't pass any.
  fpd <- as.data.frame(Cardinal::fData(combined))

  first_num <- function(x) {
    suppressWarnings(as.numeric(
      sapply(strsplit(as.character(x), " || ", fixed = TRUE), `[`, 1)
    ))
  }

  key_fd  <- paste(first_num(fpd$precursor_mz),
                   first_num(fpd$product_mz), sep = "_")
  key_lib <- paste(round(suppressWarnings(as.numeric(ion_lib_meta$precursor_mz)), 0),
                   round(suppressWarnings(as.numeric(ion_lib_meta$product_mz)), 0),
                   sep = "_")

  midx <- match(key_fd, key_lib)
  feature_meta <- ion_lib_meta[midx, , drop = FALSE]
  feature_meta$name <- fpd$name
  rownames(feature_meta) <- fpd$name

  # Coerce all character columns to valid UTF-8. CSVs saved by Excel on
  # Windows are usually Windows-1252; characters like en-dash (0x96) and
  # smart quotes break trimws()/sub() downstream. iconv with sub="" drops
  # any byte that can't be represented in UTF-8.
  .chr_cols <- vapply(feature_meta, is.character, logical(1))
  feature_meta[.chr_cols] <- lapply(feature_meta[.chr_cols], function(x)
    iconv(x, from = "", to = "UTF-8", sub = ""))

  if (verbose) {
    message(sprintf(
      "Feature metadata: %d / %d features matched to ion library by m/z.",
      sum(!is.na(midx)), nrow(fpd)
    ))
  }

  feature_meta
}
