setGeneric("int2response", function(MSIobject, ...) standardGeneric("int2response"))

#' Normalise pixel intensities to internal-standard response
#'
#' Divides each feature's intensity by an internal standard measured in the
#' same pixel, line or sample, which suppresses drift and much of the local
#' variation in ionisation efficiency.
#'
#' @details
#' `normalisation` states what is being done, so that no argument value has to
#' be read as an instruction:
#' \describe{
#'   \item{`"internal_standard"`}{Divide by the standard named in `IS_name`.
#'     This is the quantitative path.}
#'   \item{`"within_feature"`}{Divide each feature by its own summary at the
#'     chosen `mode`. This removes between-line or between-sample scale
#'     differences within a feature, but it is *not* internal-standard
#'     normalisation and does not make features comparable to one another.}
#' }
#' To leave values untouched, do not call this function. There is deliberately
#' no "none" value here: skipping a step is the caller's decision, not a mode of
#' the step. The YAML workflow exposes `normalisation: none` for that, and
#' simply does not call this function.
#'
#' @section Multiple internal standards:
#' When `IS_name` matches exactly one feature, every analyte is normalised to
#' it. When it matches several -- a panel carrying one class-specific standard
#' per lipid class, for instance -- an explicit analyte-to-standard mapping is
#' required, because dividing by several standards at once has no defined
#' meaning. Supply it as an ion-library column (named by `is_norm_header`,
#' default `"IS_norm"`) whose value on each analyte row is the transition name
#' of the standard that normalises it. Any analyte left unmapped is an error
#' rather than a silent fallback.
#'
#' With a single standard the column is never read, so `is_norm_header = "None"`
#' is the right setting for a panel where every analyte shares one standard.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @param MSIobject A `quant_MSImagingExperiment` object.
#' @param val_slot Character. Spectra slot to normalise (default
#'   `"intensity"`).
#' @param normalisation Character. `"internal_standard"` (default) or
#'   `"within_feature"`. See Details.
#' @param IS_name Character. Value identifying the internal standard in the
#'   feature-type column of `fData(MSIobject)` -- the ion library's `Type`
#'   column, typically `"IS"`. This is a type label, not a transition name:
#'   passing a feature name will not match. Required when
#'   `normalisation = "internal_standard"`; a value that matches no feature is
#'   an error, since falling back to a different normalisation would silently
#'   change the analytical method.
#' @param is_norm_header Character. Feature-metadata column holding the
#'   analyte-to-standard mapping (default `"IS_norm"`). It is consulted **only**
#'   when `IS_name` matches more than one feature, so `"None"` (or `NULL`) is
#'   correct for the common case of a single standard shared by every analyte.
#' @param mode Character. Level at which the standard is summarised:
#'   \describe{
#'     \item{`"line"`}{Median across a whole acquisition line (the default). A
#'       *line* is one horizontal raster row -- constant `y`, varying `x` --
#'       which is the order DESI acquires in, so it is also the axis along which
#'       source drift accumulates.}
#'     \item{`"sample"`}{Median across the whole acquisition.}
#'     \item{`"pixel"`}{The standard in that pixel alone -- responsive to local
#'       suppression, but carries the standard's own shot noise.}
#'     \item{`"window"`}{Rolling median over `window` consecutive pixels along
#'       the same horizontal row, ordered by `x`. A compromise between `"pixel"`
#'       and `"line"`: it smooths the standard's own shot noise while still
#'       tracking drift across the row. The window is truncated at the ends of a
#'       row rather than wrapping, so pixels from different rows are never
#'       mixed.}
#'   }
#' @param window Integer. Number of consecutive pixels averaged when
#'   `mode = "window"` (default `15`). Ignored for the other modes.
#' @param remove_IS Logical. Drop the internal-standard features from the
#'   returned object (default `TRUE`).
#' @param ... Additional arguments (currently unused).
#' @return The input object with a `response` spectra slot holding the
#'   normalised values.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' # No internal standard in the example panel: normalise each feature to
#' # itself, per acquisition line.
#' obj <- int2response(obj, val_slot = "intensity",
#'                     normalisation = "within_feature")
#'
#' @family filtering
#' @aliases int2response
#' @export
setMethod("int2response", "quant_MSImagingExperiment",
          function(MSIobject, val_slot = "intensity",
                   normalisation = c("internal_standard", "within_feature"),
                   IS_name = NULL, is_norm_header = "IS_norm",
                   mode = c("line", "sample", "pixel", "window"),
                   window = 15, remove_IS = TRUE, ...){

            mode <- match.arg(mode)
            # "None" is how a YAML config spells "not set"; a single-standard
            # panel never needs a mapping column, so this is the normal case.
            is_norm_header <- .none(is_norm_header)

            # Back-compatibility: IS_name = "None" used to select within-feature
            # normalisation. That overloaded a standard's name with a mode, and
            # meant the same token skipped normalisation in the YAML while
            # performing it here.
            if(!is.null(IS_name) && identical(as.character(IS_name), "None")){
              if(missing(normalisation)) normalisation <- "within_feature"
              IS_name <- NULL
            }
            normalisation <- match.arg(normalisation)

            n_feat  <- nrow(fData(MSIobject))
            type_col <- .feature_type(MSIobject)

            # is_row[f] is the feature index of the standard that normalises
            # feature f, or NA when f is normalised to itself.
            is_row  <- rep(NA_integer_, n_feat)
            IS_ind  <- integer(0)

            if(normalisation == "internal_standard"){

              if(is.null(IS_name) || !nzchar(as.character(IS_name)))
                stop("int2response: normalisation = \"internal_standard\" needs ",
                     "`IS_name`. Use normalisation = \"within_feature\" to ",
                     "normalise each feature to itself.", call. = FALSE)

              if(is.null(type_col))
                stop("int2response: fData(MSIobject) has no feature-type column, ",
                     "so IS_name = '", IS_name, "' cannot be matched. This ",
                     "column comes from the ion library's `Type` column via ",
                     "read_mrm(type_header = ).", call. = FALSE)

              IS_ind <- which(type_col == IS_name)

              # A missing standard is an error: silently normalising to
              # something else changes the analytical method behind the user's
              # back, and every downstream value would carry that change.
              if(length(IS_ind) == 0L)
                stop("int2response: no feature is typed '", IS_name, "'. ",
                     "Types present: ",
                     paste(sort(unique(stats::na.omit(type_col))), collapse = ", "),
                     ".", call. = FALSE)

              if(length(IS_ind) == 1L){
                is_row[] <- IS_ind
              } else {
                # Several standards: an explicit mapping is the only defined
                # interpretation.
                map <- .feature_col(MSIobject, is_norm_header)
                if(is.null(map))
                  stop("int2response: ", length(IS_ind), " features are typed '",
                       IS_name, "' (", paste(fData(MSIobject)$name[IS_ind],
                                              collapse = ", "),
                       "), so each analyte must name the standard it uses. ",
                       if(is.null(is_norm_header))
                         "Set is_norm_header to an ion-library column giving "
                       else
                         paste0("Add an ion-library column '", is_norm_header,
                                "' giving "),
                       "that standard's transition name for every analyte.",
                       call. = FALSE)

                nm  <- as.character(fData(MSIobject)$name)
                tgt <- match(trimws(as.character(map)), nm)
                # Standards themselves need no mapping.
                tgt[IS_ind] <- IS_ind

                bad <- setdiff(which(is.na(tgt)), IS_ind)
                if(length(bad))
                  stop("int2response: '", is_norm_header, "' does not name a ",
                       "feature in this object for: ",
                       paste(nm[bad], collapse = ", "), ".", call. = FALSE)

                not_is <- setdiff(unique(tgt[-IS_ind]), IS_ind)
                if(length(not_is))
                  stop("int2response: '", is_norm_header, "' points at features ",
                       "that are not typed '", IS_name, "': ",
                       paste(nm[not_is], collapse = ", "), ".", call. = FALSE)

                is_row <- tgt
              }
            }

            spectra(MSIobject, "response") = matrix(nrow = nrow(MSIobject), ncol = ncol(MSIobject))

            # Iterate over samples in study
            for(sample in unique(pData(MSIobject)$run)){

              # Select IS mz and pixels for specific sample
              sample_pixels = which(pData(MSIobject)$run == sample)

              tempMSIobject = MSIobject[, sample_pixels]

              vals = spectraData(tempMSIobject)[[val_slot]]

              for(mz_ind in seq_len(nrow(fData(tempMSIobject)))){

                ints = vals[mz_ind, ]

                # Denominator vector: this feature's own standard where one was
                # mapped, otherwise the feature itself (within-feature
                # normalisation).
                denom_src = if(is.na(is_row[mz_ind])) ints else vals[is_row[mz_ind], ]

                # Build `response` by index rather than by concatenation. The
                # per-line and per-window modes group pixels, and groups are not
                # guaranteed to be contiguous or in pixel order, so appending
                # results group-by-group would misalign them on assignment.
                response = rep(NA_real_, length(ints))

                if(mode == "pixel"){
                  response = ints / denom_src

                } else if(mode == "sample"){
                  response = ints / median(denom_src, na.rm = TRUE)

                } else if(mode == "line"){
                  for(line in unique(pData(tempMSIobject)$y)){
                    lp = which(pData(tempMSIobject)$y == line)
                    response[lp] = ints[lp] / median(denom_src[lp], na.rm = TRUE)
                  }

                } else if(mode == "window"){
                  # Rolling median of the standard over `window` consecutive
                  # pixels along each acquisition line, ordered by x. The window
                  # is truncated at the ends of a line rather than wrapping, so
                  # pixels from different lines are never mixed.
                  for(line in unique(pData(tempMSIobject)$y)){
                    lp = which(pData(tempMSIobject)$y == line)
                    lp = lp[order(pData(tempMSIobject)$x[lp])]
                    d  = .roll_median(denom_src[lp], window)
                    response[lp] = ints[lp] / d
                  }
                }

                spectra(MSIobject, "response")[mz_ind, sample_pixels] = response

              }
            }

            # Remove the standard features once every analyte has been divided
            # by them.
            if(remove_IS == TRUE && length(IS_ind) > 0L){
              MSIobject = MSIobject[-IS_ind, ]
            }

            return(MSIobject)
          })
