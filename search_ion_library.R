#!/usr/bin/env Rscript
# =============================================================================
# search_ion_library.R
#
# Search the combined ion library by name and/or m/z and report hits.
#
# -- Interactive use (RStudio, R console) --------------------------------
#   source("D:/__STUDIES/search_ion_library.R")
#   search_library(name = "PGE2")
#   search_library(precursor = 351.2, tol = 0.01)
#   search_library(precursor = 351.2, product = 271.1, tol = 0.05)
#   search_library(name = "HETE", polarity = "Negative")
#   search_library(lipid = "Oxylipin", class = "HETE")
#
# -- Command-line use ----------------------------------------------------
#   Rscript search_ion_library.R --name PGE2
#   Rscript search_ion_library.R --precursor 351.2 --tol 0.01
#   Rscript search_ion_library.R --precursor 351.2 --product 271.1 --tol 0.05
#   Rscript search_ion_library.R --name HETE --polarity Negative
#   Rscript search_ion_library.R --lipid Oxylipin --class HETE
#
# Result is printed to the console; pass --csv path/to/out.csv to also save.
# =============================================================================

suppressMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

LIB_PATH <- "D:/__STUDIES/ion_library_combined.csv"

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Search the combined ion library
#'
#' Any subset of arguments can be combined; criteria are AND-ed together.
#'
#' @param name      Substring (case-insensitive) matched against `transition_id`.
#'                  Use `regex = TRUE` to interpret as a regular expression.
#' @param precursor Numeric precursor m/z. Matches within `tol` (Daltons).
#' @param product   Numeric product   m/z. Matches within `tol`.
#' @param tol       m/z tolerance in Daltons (default 0.01).
#' @param polarity  "Positive" or "Negative" (case-insensitive).
#' @param lipid     Coarse Lipid category (e.g. "Oxylipin", "Phospholipid").
#' @param class     Specific Class (e.g. "HETE", "PC", "Cer").
#' @param flag      Restrict to "ok" / "conflict" / "blank" / "conflict_blank".
#' @param regex     Treat `name` as a regex (default FALSE = substring match).
#' @param lib_path  Path to the combined ion library CSV (default `LIB_PATH`).
#' @param verbose   Print a one-line summary to the console (default TRUE).
#'
#' @return A tibble of matching rows. The hit count is printed when
#'   `verbose = TRUE`.
search_library <- function(name = NULL,
                           precursor = NULL,
                           product = NULL,
                           tol = 0.01,
                           polarity = NULL,
                           lipid = NULL,
                           class = NULL,
                           flag = NULL,
                           regex = FALSE,
                           lib_path = LIB_PATH,
                           verbose = TRUE) {

  if (!file.exists(lib_path))
    stop("Library not found: ", lib_path,
         "\nRun build_ion_library.R first.")

  lib <- suppressMessages(read_csv(lib_path, show_col_types = FALSE))

  hits <- lib

  # ---- Name match (substring or regex on transition_id) ------------------
  if (!is.null(name) && nzchar(name)) {
    if (regex) {
      hits <- hits |> filter(str_detect(transition_id, regex(name, ignore_case = TRUE)))
    } else {
      hits <- hits |> filter(str_detect(transition_id, fixed(name, ignore_case = TRUE)))
    }
  }

  # ---- m/z windows -------------------------------------------------------
  if (!is.null(precursor)) {
    precursor <- as.numeric(precursor)
    hits <- hits |> filter(!is.na(precursor_mz),
                            abs(precursor_mz - precursor) <= tol)
  }
  if (!is.null(product)) {
    product <- as.numeric(product)
    hits <- hits |> filter(!is.na(product_mz),
                            abs(product_mz - product) <= tol)
  }

  # ---- Categorical filters ----------------------------------------------
  if (!is.null(polarity) && nzchar(polarity)) {
    hits <- hits |> filter(toupper(Polarity) == toupper(polarity))
  }
  if (!is.null(lipid) && nzchar(lipid)) {
    hits <- hits |> filter(!is.na(Lipid),
                            tolower(Lipid) == tolower(lipid))
  }
  if (!is.null(class) && nzchar(class)) {
    hits <- hits |> filter(!is.na(Class),
                            str_detect(Class, fixed(class, ignore_case = TRUE)))
  }
  if (!is.null(flag) && nzchar(flag)) {
    hits <- hits |> filter(.data$flag == .env$flag)
  }

  if (verbose) {
    crit <- character()
    if (!is.null(name))      crit <- c(crit, sprintf("name~'%s'", name))
    if (!is.null(precursor)) crit <- c(crit, sprintf("precursor=%.4f±%.3f", precursor, tol))
    if (!is.null(product))   crit <- c(crit, sprintf("product=%.4f±%.3f", product, tol))
    if (!is.null(polarity))  crit <- c(crit, sprintf("polarity=%s", polarity))
    if (!is.null(lipid))     crit <- c(crit, sprintf("lipid=%s", lipid))
    if (!is.null(class))     crit <- c(crit, sprintf("class~'%s'", class))
    if (!is.null(flag))      crit <- c(crit, sprintf("flag=%s", flag))
    if (length(crit) == 0)   crit <- "(no criteria — full library)"
    cat(sprintf("[%d hit%s] %s\n",
                nrow(hits), if (nrow(hits) == 1) "" else "s",
                paste(crit, collapse = " & ")))
  }

  hits
}

#' Batch-search a table of input transitions against the library
#'
#' For each input row, find library rows where precursor and product m/z fall
#' within `tol`. Reports per-row hit count, whether any of the hits also
#' share a name token (split on `||`, `/`, or `,`), and the matching
#' library transition_ids and sources.
#'
#' @param queries  data.frame / tibble with at minimum `precursor_mz` and
#'                 `product_mz`. An optional `transition_id` column enables
#'                 the name-match check.
#' @param tol      m/z tolerance in Da (default 0.01).
#' @param lib_path Path to combined library CSV.
#'
#' @return The input table augmented with:
#'   * `n_hits`         number of library rows within `tol`
#'   * `n_name_match`   how many of those also share a name token (NA if no
#'                      `transition_id` column was provided)
#'   * `hit_ids`        `;`-separated transition_ids of all hits
#'   * `hit_sources`    `;`-separated unique source panels of all hits
#'   * `match_status`   "exact" (name+mz), "mz_only", or "miss"
batch_search <- function(queries,
                         tol = 0.01,
                         lib_path = LIB_PATH) {

  if (!file.exists(lib_path))
    stop("Library not found: ", lib_path,
         "\nRun build_ion_library.R first.")

  needed <- c("precursor_mz", "product_mz")
  miss <- setdiff(needed, names(queries))
  if (length(miss) > 0)
    stop("Input table is missing required columns: ",
         paste(miss, collapse = ", "))

  lib <- suppressMessages(read_csv(lib_path, show_col_types = FALSE))

  # Tokenise a name on common separators, lowercase, trim.
  tokens <- function(x) {
    if (is.na(x) || !nzchar(x)) return(character(0))
    toks <- unlist(strsplit(x, "\\s*(\\|\\||/|,)\\s*"))
    tolower(trimws(toks[nzchar(toks)]))
  }

  has_name_col <- "transition_id" %in% names(queries)

  q_pre <- as.numeric(queries$precursor_mz)
  q_pro <- as.numeric(queries$product_mz)
  q_nam <- if (has_name_col) as.character(queries$transition_id) else rep(NA_character_, nrow(queries))

  out <- queries
  out$n_hits       <- 0L
  out$n_name_match <- if (has_name_col) 0L else NA_integer_
  out$hit_ids      <- ""
  out$hit_sources  <- ""

  for (i in seq_len(nrow(queries))) {
    if (is.na(q_pre[i]) || is.na(q_pro[i])) next

    hits <- lib |>
      filter(!is.na(precursor_mz), !is.na(product_mz),
             abs(precursor_mz - q_pre[i]) <= tol,
             abs(product_mz   - q_pro[i]) <= tol)

    out$n_hits[i]      <- nrow(hits)
    out$hit_ids[i]     <- paste(unique(hits$transition_id), collapse = "; ")
    out$hit_sources[i] <- paste(unique(unlist(strsplit(hits$sources, "\\s*\\|\\s*"))),
                                collapse = "; ")

    if (has_name_col && nrow(hits) > 0 && !is.na(q_nam[i])) {
      qtok  <- tokens(q_nam[i])
      n_ok  <- sum(vapply(hits$transition_id, function(x) {
        any(tokens(x) %in% qtok)
      }, logical(1)))
      out$n_name_match[i] <- n_ok
    }
  }

  out$match_status <- with(out, dplyr::case_when(
    n_hits == 0                                     ~ "miss",
    has_name_col & (n_name_match %||% 0L) > 0       ~ "exact",
    TRUE                                            ~ "mz_only"
  ))

  message(sprintf(
    "[batch] %d queries: %d exact, %d mz-only, %d miss (tol=%.3f Da)",
    nrow(out),
    sum(out$match_status == "exact"),
    sum(out$match_status == "mz_only"),
    sum(out$match_status == "miss"),
    tol
  ))

  out
}

# =============================================================================
# CLI entry point
# =============================================================================
parse_cli <- function(args) {
  if (length(args) == 0) return(NULL)
  out <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (!startsWith(a, "--")) {
      stop("Unexpected positional arg '", a, "'. Use --flag value.")
    }
    key <- sub("^--", "", a)
    if (key %in% c("regex", "help", "h")) {
      out[[key]] <- TRUE; i <- i + 1L
    } else {
      if (i + 1L > length(args))
        stop("Missing value for --", key)
      out[[key]] <- args[i + 1L]; i <- i + 2L
    }
  }
  out
}

print_help <- function() {
  cat("Usage: Rscript search_ion_library.R [options]\n\n",
      "Single-query mode:\n",
      "  --name      <str>      substring match against transition_id\n",
      "  --regex                treat --name as regex\n",
      "  --precursor <mz>       precursor m/z (Da)\n",
      "  --product   <mz>       product   m/z (Da)\n",
      "  --polarity  <Pos|Neg>  filter by polarity\n",
      "  --lipid     <cat>      coarse Lipid category (Polar, Oxylipin, ...)\n",
      "  --class     <cls>      specific Class (HETE, PC, Cer, ...)\n",
      "  --flag      <flag>     ok | conflict | blank | conflict_blank\n\n",
      "Batch mode (one query per row):\n",
      "  --batch     <path>     TSV/CSV with cols: transition_id (opt),\n",
      "                         precursor_mz, product_mz\n\n",
      "Common:\n",
      "  --tol       <Da>       m/z tolerance (default 0.01)\n",
      "  --csv       <path>     also write results to CSV\n",
      "  --lib       <path>     override library CSV path\n",
      "  --help                 print this help and exit\n",
      sep = "")
}

if (!interactive() && sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- tryCatch(parse_cli(args), error = function(e) {
    cat("Error:", e$message, "\n\n"); print_help(); quit(status = 1)
  })

  if (isTRUE(opts$help) || isTRUE(opts$h) || length(args) == 0) {
    print_help(); quit(status = 0)
  }

  # ---- Batch mode --------------------------------------------------------
  if (!is.null(opts$batch) && nzchar(opts$batch)) {
    if (!file.exists(opts$batch))
      stop("Batch input not found: ", opts$batch)

    queries <- if (grepl("\\.tsv$|\\.txt$", opts$batch, ignore.case = TRUE))
                 suppressMessages(read_tsv(opts$batch, show_col_types = FALSE))
               else
                 suppressMessages(read_csv(opts$batch, show_col_types = FALSE))

    res <- batch_search(
      queries  = queries,
      tol      = if (!is.null(opts$tol)) as.numeric(opts$tol) else 0.01,
      lib_path = if (!is.null(opts$lib)) opts$lib else LIB_PATH
    )

    print(res, n = nrow(res))

    if (!is.null(opts$csv) && nzchar(opts$csv)) {
      write_csv(res, opts$csv)
      cat("Wrote", nrow(res), "rows to", opts$csv, "\n")
    }
    quit(status = 0)
  }

  # ---- Single-query mode -------------------------------------------------
  hits <- search_library(
    name      = opts$name,
    precursor = if (!is.null(opts$precursor)) as.numeric(opts$precursor) else NULL,
    product   = if (!is.null(opts$product))   as.numeric(opts$product)   else NULL,
    tol       = if (!is.null(opts$tol))       as.numeric(opts$tol)       else 0.01,
    polarity  = opts$polarity,
    lipid     = opts$lipid,
    class     = opts$class,
    flag      = opts$flag,
    regex     = isTRUE(opts$regex),
    lib_path  = opts$lib %||% LIB_PATH
  )

  if (nrow(hits) > 0) {
    print(hits, n = 50)
  }

  if (!is.null(opts$csv) && nzchar(opts$csv)) {
    write_csv(hits, opts$csv)
    cat("Wrote", nrow(hits), "hits to", opts$csv, "\n")
  }
}
