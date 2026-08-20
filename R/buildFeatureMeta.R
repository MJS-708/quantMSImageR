#' Build per-feature metadata by joining an ion library on m/z
#'
#' Joins each row of `fData(combined)` to the ion library on the
#' (`precursor_mz`, `product_mz`) pair, within `mz_tolerance`. Robust to feature
#' renames in the YAML (`features.rename`) because m/z is stable across name
#' changes.
#'
#' Matching uses the same rule as [readMRM()], so the two cannot disagree about
#' which library entry a feature is. They previously built their m/z keys
#' independently -- one rounding to whole numbers, the other not -- and a
#' feature matching several library entries silently took the first one's
#' annotation.
#'
#' Used by `runStudy.R` (via `inst/runStudy.R`) and [runExample()] to build
#' a per-feature metadata frame the heatmap report can consume -- e.g. for
#' splitting / colour-bar annotation by an ion-library column such as `Met-1`.
#'
#' @param combined An MSI object whose `fData` has columns `name`,
#'   `precursor_mz` and `product_mz`. Multi-transition features stored with
#'   `" || "`-joined m/z strings are matched on the first token.
#' @param ion_lib_meta A `data.frame` of ion library metadata (must contain
#'   `precursor_mz` and `product_mz`). Pass `NULL` to short-circuit and return
#'   `NULL`.
#' @param mz_tolerance Numeric. Half-width in Da within which a feature's
#'   precursor and product -- both, not either -- are taken to be the library's.
#'   Default `0.4`, nominal-mass matching, because MRM selects Q1 and Q3 at unit
#'   resolution. See [readMRM()].
#' @param ambiguity Character. What to do when a feature matches more than one
#'   library entry: `"combine"` (default, name it for all of them joined with
#'   `" || "`), `"error"`, `"warn"` or `"nearest"`. See [readMRM()].
#' @param verbose Logical. Emit a `message()` reporting how many features
#'   matched. Default `TRUE`.
#'
#' @return A `data.frame` of metadata, one row per feature in `fData(combined)`,
#'   with row names set to the feature name and a `name` column appended.
#'   Unmatched features become a row of NAs. Returns `NULL` when
#'   `ion_lib_meta` is `NULL`.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- readRDS(p)
#' lib <- read.csv(system.file("extdata", "example_ion_library.csv",
#'                             package = "quantMSImageR"), check.names = FALSE)
#' fm <- buildFeatureMeta(obj, lib)
#'
#' @seealso [quantileHm()], [generateTxtImages()]
#' @export
buildFeatureMeta <- function(combined, ion_lib_meta,
                               mz_tolerance = 0.4,
                               ambiguity = c("combine", "error", "warn",
                                             "nearest"),
                               verbose = TRUE) {
  ambiguity <- match.arg(ambiguity)
  if (is.null(ion_lib_meta)) return(NULL)

  # Without these, .match_transitions receives numeric(0), every which() is
  # empty, and the function returns a full-height data frame of NA rows. The
  # only signal is a verbose message reporting "0 / N features matched", which
  # verbose = FALSE suppresses entirely -- a silent wrong answer.
  if (!is.data.frame(ion_lib_meta))
    stop("buildFeatureMeta: ion_lib_meta must be a data.frame; got ",
         class(ion_lib_meta)[1L], ".", call. = FALSE)
  .need <- c("precursor_mz", "product_mz")
  if (!all(.need %in% names(ion_lib_meta)))
    stop("buildFeatureMeta: ion_lib_meta must carry ",
         paste(.need, collapse = " and "), "; missing ",
         paste(setdiff(.need, names(ion_lib_meta)), collapse = ", "),
         ". Features are matched on transition m/z.", call. = FALSE)

  # Cardinal's `fData` returns a MassDataFrame; its as.data.frame method
  # warns on extra args, so we don't pass any.
  fpd <- as.data.frame(Cardinal::fData(combined))

  first_num <- function(x) {
    suppressWarnings(as.numeric(
      vapply(strsplit(as.character(x), " || ", fixed = TRUE), `[`, character(1), 1)
    ))
  }

  midx <- .match_transitions(
    first_num(fpd$precursor_mz), first_num(fpd$product_mz),
    ion_lib_meta$precursor_mz, ion_lib_meta$product_mz,
    tolerance  = mz_tolerance,
    ambiguity  = ambiguity,
    labels     = as.character(fpd$name),
    ref_labels = as.character(ion_lib_meta$transition_id),
    context    = "buildFeatureMeta")
  feature_meta <- ion_lib_meta[midx, , drop = FALSE]

  # Same rule as readMRM(): a feature matching several library entries is
  # named for all of them. Kept in step here so a report cannot show the joined
  # name in one column and one of its halves in another.
  amb_sets <- attr(midx, "matches")
  if (ambiguity == "combine" && !is.null(amb_sets) &&
      "transition_id" %in% names(feature_meta)) {
    for (i in seq_along(amb_sets)) {
      h <- amb_sets[[i]]
      if (is.null(h)) next
      feature_meta$transition_id[i] <-
        paste(unique(as.character(ion_lib_meta$transition_id[h])),
              collapse = " || ")
    }
  }

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
