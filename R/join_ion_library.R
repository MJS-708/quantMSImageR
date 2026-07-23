# Annotate measured MRM transitions from an ion library.
#
# Split out of read_mrm() so it can be exercised without a .raw folder: the
# join is the part with the analytical risk in it, and the surrounding function
# needs instrument text files that no test can supply.
#
# One row out per measured transition, in the order the features will appear.
# A transition the library does not describe is kept, named by its instrument
# index and typed "Unknown", because dropping a measured channel silently would
# be worse than carrying it unannotated.
#
# `transitions` needs transition_id, precursor_mz, product_mz (as read from the
# instrument). `ion_lib` is the library CSV. Returns `ion_lib`'s columns for the
# matched row, plus:
#   transition_id_int   instrument transition index (selects the data column)
#   transition_id_name  library name, or the instrument index when unmatched
#   precursor_mz/product_mz  as measured
#   new_transition_int  1..n, the mz key
.join_ion_library <- function(transitions, ion_lib, polarity,
                              type_header = "Type",
                              mz_tolerance = 0.05,
                              ambiguity = c("error", "warn", "nearest")) {

  ambiguity <- match.arg(ambiguity)

  lib <- ion_lib[as.character(ion_lib$Polarity) == polarity, , drop = FALSE]
  if (nrow(lib) == 0)
    warning("read_mrm: the ion library has no ", polarity,
            " entries, so no transition can be annotated.", call. = FALSE)

  idx <- .match_transitions(
    transitions$precursor_mz, transitions$product_mz,
    lib$precursor_mz, lib$product_mz,
    tolerance  = mz_tolerance,
    ambiguity  = ambiguity,
    labels     = paste0("transition ", transitions$transition_id),
    ref_labels = as.character(lib$transition_id),
    context    = "read_mrm")

  # Library columns for the matched row; NA throughout where nothing matched.
  out <- lib[idx, , drop = FALSE]
  rownames(out) <- NULL

  out$transition_id_int  <- transitions$transition_id
  out$precursor_mz       <- as.character(transitions$precursor_mz)
  out$product_mz         <- as.character(transitions$product_mz)
  out$transition_id_name <- ifelse(is.na(idx),
                                   as.character(transitions$transition_id),
                                   as.character(lib$transition_id[idx]))
  out$Polarity <- polarity
  if (type_header %in% names(out))
    out[[type_header]] <- ifelse(is.na(out[[type_header]]), "Unknown",
                                 as.character(out[[type_header]]))

  # Two measured channels can legitimately point at one library entry (the same
  # transition acquired twice). Names index features downstream, so they are
  # made unique the way they always were.
  dup <- stats::ave(seq_len(nrow(out)), out$transition_id_name, FUN = seq_along)
  out$transition_id_name <- ifelse(dup > 1,
                                   sprintf("%s:- %s", dup, out$transition_id_name),
                                   out$transition_id_name)

  ord <- order(suppressWarnings(as.numeric(out$precursor_mz)),
               suppressWarnings(as.numeric(out$product_mz)),
               out$transition_id_int)
  out <- out[ord, , drop = FALSE]
  rownames(out) <- NULL
  out$new_transition_int <- seq_len(nrow(out))
  out
}
