#' Colour palettes used by quantMSImageR
#'
#' The palettes available to [imageR()], [quantile_hm()] and the metadata
#' colour bars, returned as named character vectors of hex colours.
#'
#' @details
#' Four of these -- `heatmap0`, `heatmap2`, `hat` and `reading` -- are taken
#' from the `ltc` package by Loukas Theodosiou
#' (<https://github.com/loukesio/ltc-color-palettes>), which is MIT licensed.
#' `heatmap2` and `hat` are re-ordered relative to ltc -- see below -- and
#' `hat` substitutes a grey for ltc's black; the remaining colours are
#' unchanged.
#' They are reproduced here rather than depended on: `ltc` brings in 34
#' recursive dependencies, which is a large amount of installation surface for
#' four colour vectors in a package that is otherwise light.
#'
#' Which palette suits what:
#' \describe{
#'   \item{`viridis`}{Sequential, perceptually uniform and colour-blind safe.
#'     The default for ion images, where the quantity is one-directional.}
#'   \item{`heatmap0`}{Nine-colour sequential ramp running dark blue to teal to
#'     sand to red. An alternative for ion images.}
#'   \item{`heatmap2`}{Five-colour diverging ramp, blue to white to red. Suits
#'     z-scores, where the midpoint is meaningful.}
#'   \item{`hat`, `reading`}{Qualitative, for categorical metadata such as study
#'     group or pathway class. `hat` carries 10 well-separated hues, re-sequenced
#'     from ltc's ordering so that the first few are maximally distinct, with
#'     ltc's black replaced by a neutral grey placed last; `reading` is a softer
#'     8-colour set.}
#' }
#'
#' @param name Optional palette name. When `NULL` (default) the whole list is
#'   returned.
#' @param n Optional number of colours. Qualitative palettes are truncated,
#'   continuous ones are interpolated. `NULL` (default) returns the palette at
#'   its native length.
#'
#' @return A named list of character vectors, or a single character vector of
#'   hex colours when `name` is given.
#'
#' @examples
#' names(quant_palettes())
#' quant_palettes("heatmap2")
#' quant_palettes("hat", n = 3)
#'
#' @family visualisation
#' @export
quant_palettes <- function(name = NULL, n = NULL) {

  pals <- list(
    # Sequential ramps
    heatmap0 = c("#001219", "#005F73", "#0A9396", "#94D2BD", "#E9D8A6",
                 "#EE9B00", "#CA6702", "#AE2012", "#9B2226"),
    # Diverging ramp, blue (low) -> white -> red (high). Reversed relative to
    # ltc's ordering so that low values are blue, as z-score maps expect.
    heatmap2 = rev(c("#ca0020", "#f4a582", "#f7f7f7", "#92c5de", "#0571b0")),
    # Qualitative sets for categorical metadata. `hat` is re-sequenced from
    # ltc's ordering, which runs around the colour wheel and so gives adjacent
    # yellows and oranges to the first few classes of a pathway annotation.
    # This ordering interleaves hue and lightness, so the two to five classes
    # a study actually has are distinguishable at a glance.
    #
    # ltc's black is replaced by a neutral grey and moved to the end: black
    # reads as a border or as text rather than as a category, and a grey is
    # the natural colour for whatever class ends up last (often "Other").
    hat      = c("#4e54ac", "#e8351e", "#17a769", "#efb306", "#852f88",
                 "#0f8096", "#cd023d", "#7db954", "#eb990c", "#7f7f7f"),
    reading  = c("#EFBC68", "#919F89", "#EDBDAE", "#57717C", "#5F97A4",
                 "#CAEAC8", "#95A1AE", "#C8CFD6")
  )

  if (is.null(name)) return(pals)

  name <- match.arg(name, names(pals))
  out  <- pals[[name]]

  if (!is.null(n)) {
    if (name %in% c("hat", "reading")) {
      # Qualitative: recycle rather than interpolate, so hues stay distinct.
      out <- rep_len(out, n)
    } else {
      out <- grDevices::colorRampPalette(out)(n)
    }
  }
  out
}

# Colour vector for a categorical annotation, named by level. Falls back to
# hcl.colors when a palette name is not one of the vendored sets.
.anno_cols <- function(levels, palette = "hat") {
  n <- max(length(levels), 2)
  cols <- if (palette %in% names(quant_palettes()))
            quant_palettes(palette, n = n)
          else
            grDevices::hcl.colors(n, palette = palette)
  stats::setNames(cols[seq_along(levels)], levels)
}
