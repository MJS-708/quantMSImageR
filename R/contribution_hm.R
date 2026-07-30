#' Group-mean heatmap with per-sample contribution as opacity
#'
#' The same heatmap as [quantile_hm()] -- same orientation, same group and
#' feature colour bars, same square cells -- with only the cell fill changed.
#' Each cell encodes two things at once:
#'
#' \describe{
#'   \item{Hue}{the sample's **group** mean for that feature -- what the group
#'     did.}
#'   \item{Opacity}{that **sample's** contribution to the group mean -- who
#'     made it do that.}
#' }
#'
#' A plain group-mean heatmap cannot show that a group effect rests on one
#' replicate. Here it is visible directly: a uniformly solid block is a real
#' group effect, whereas one opaque tile among faded ones means the mean is
#' carried by a single sample.
#'
#' @section How the two channels are computed:
#'
#' Values are the per-sample quantile of pixel intensities, as in
#' [quantile_hm()], then z-scored **per feature across all samples**, not
#' within group -- so a cell reads "high or low for this feature" rather than
#' "abundant or not". Features with zero variance, and missing values, stay
#' `NA` and are drawn in `na_col`; they are never imputed to the row mean.
#'
#' The group mean of those z-scores is the hue. The contribution is
#' `z * sign(group mean)`, so positive means the sample moved with its group
#' and zero or negative means it did not; only the positive part drives
#' opacity, rescaled onto `[alpha_floor, 1]`. `alpha_floor` exists because a
#' nearly transparent tile reads as missing data rather than as weak agreement.
#'
#' @section Why there are two colour limits:
#'
#' The hue is a group **mean**, and averaging \eqn{n} replicates shrinks it by
#' roughly \eqn{\sqrt{n}}. Reusing the per-sample z limit for the group mean
#' therefore produces a washed-out panel that no amount of alpha tuning fixes.
#' Both limits are computed from the quantity actually drawn, each as the
#' `1 - saturate` quantile of its own absolute values, so a similar small
#' fraction of cells saturates on each channel and the bulk of the data gets
#' the full ramp. Expect the two to differ by a factor of two or more. Note
#' this differs from [quantile_hm()], which clips to a fixed \[-1, 1\].
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @inheritParams quantile_hm
#' @param quant_val Numeric in (0, 1). Quantile of pixel intensities summarised
#'   per feature per sample (default `0.5`, the median).
#' @param alpha_floor Numeric in \[0, 1). Opacity given to a sample that did not
#'   contribute (default `0.6`). Below roughly this a non-contributing sample
#'   starts to look like missing data.
#' @param saturate Numeric in (0, 0.5). Fraction of cells allowed to saturate
#'   at the colour limit (default `0.04`, i.e. 4%).
#' @param na_col Colour for features that cannot be scored -- no variance
#'   across samples, or missing values (default `"grey88"`).
#'
#' @return A `ComplexHeatmap::Heatmap` (rows = samples, columns = features),
#'   carrying the opacity key in `attr(x, "contribution_legend")`. Pass that to
#'   `ComplexHeatmap::draw(annotation_legend_list = )` so the opacity channel
#'   is labelled; the report does this.
#'
#' @examples
#' p <- system.file("extdata", "example.raw", "section01.RDS",
#'                  package = "quantMSImageR")
#' obj <- as(readRDS(p), "quant_MSImagingExperiment")
#' hm <- contribution_hm(obj, quant_val = 0.5,
#'                       heatmap_order = "section01", heatmap_labs = "A")
#'
#' @seealso [quantile_hm()] for the per-sample view of the same matrix.
#' @family visualisation
#' @export
contribution_hm <- function(MSIobject, quant_val = 0.5,
                            heatmap_order = NA, heatmap_labs = NA,
                            feature_split = NULL,
                            feature_split_name = "Pathway",
                            group_split_name = "Group",
                            cell_size = 8, max_aspect = 1.5,
                            cell_border = "white", fontsize = 8,
                            palette = c("heatmap2", "heatmap0"),
                            group_palette = "hat",
                            feature_palette = "reading",
                            alpha_floor = 0.6, saturate = 0.04,
                            na_col = "grey88") {

  palette <- match.arg(palette)
  if (!is.numeric(alpha_floor) || length(alpha_floor) != 1L ||
      is.na(alpha_floor) || alpha_floor < 0 || alpha_floor >= 1)
    stop("contribution_hm: alpha_floor must be a single number in [0, 1).",
         call. = FALSE)
  if (!is.numeric(saturate) || length(saturate) != 1L || is.na(saturate) ||
      saturate <= 0 || saturate >= 0.5)
    stop("contribution_hm: saturate must be a single number in (0, 0.5).",
         call. = FALSE)

  mat <- .quantile_matrix(MSIobject, quant_val, heatmap_order)
  samples <- colnames(mat)

  grp <- if (all(is.na(heatmap_labs))) rep("All", length(samples))
         else as.character(heatmap_labs)
  if (length(grp) != length(samples))
    stop("contribution_hm: heatmap_labs has ", length(grp),
         " entries but there are ", length(samples), " samples.", call. = FALSE)

  v <- .contribution_values(mat, grp, alpha_floor = alpha_floor,
                            saturate = saturate)

  ## ---- same orientation and layout as quantile_hm --------------------------
  fill_t  <- t(v$fill)
  alpha_t <- t(v$alpha)

  gs <- if (!all(is.na(heatmap_labs)))
          factor(heatmap_labs, levels = unique(heatmap_labs)) else NULL
  fs <- if (!is.null(feature_split))
          factor(feature_split, levels = unique(feature_split)) else NULL

  layout <- .hm_layout(n_row = nrow(fill_t), n_col = ncol(fill_t),
                       gs = gs, fs = fs,
                       group_split_name = group_split_name,
                       feature_split_name = feature_split_name,
                       group_palette = group_palette,
                       feature_palette = feature_palette,
                       cell_size = cell_size, max_aspect = max_aspect,
                       cell_border = cell_border, fontsize = fontsize)

  cols   <- quant_palettes(palette)
  col_fn <- circlize::colorRamp2(
    seq(-v$cap_fill, v$cap_fill, length.out = length(cols)), cols)

  # Every cell is drawn by hand, because the opacity varies per cell and a
  # colour-mapping function sees only the value. rect_gp type = "none"
  # suppresses ComplexHeatmap's own rectangles so they are not drawn twice.
  .border <- layout$rect_gp$col
  cell_fun <- function(j, i, x, y, width, height, fill) {
    a <- alpha_t[i, j]
    if (!is.finite(a)) a <- 1
    grid::grid.rect(x, y, width, height,
                    gp = grid::gpar(
                      fill = grDevices::adjustcolor(fill, alpha.f = a),
                      col = .border, lwd = 0.5))
  }

  layout$rect_gp <- grid::gpar(type = "none")

  hm <- do.call(ComplexHeatmap::Heatmap, c(list(
    fill_t, name = "Mean\nz-score", col = col_fn, na_col = na_col,
    cell_fun = cell_fun), layout))

  # The opacity channel is an encoding, so it needs a key. ComplexHeatmap
  # legends for a heatmap come from its colour mapping, which knows nothing
  # about alpha, so this one is built here and attached for draw() to place.
  ref <- cols[length(cols)]
  attr(hm, "contribution_legend") <- ComplexHeatmap::Legend(
    title = "Sample\ncontribution",
    at = c(1, 2), labels = c("not carrying it", "carrying it"),
    legend_gp = grid::gpar(
      fill = c(grDevices::adjustcolor(ref, alpha.f = alpha_floor), ref)),
    title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = fontsize))
  hm
}


# Per-feature, per-sample summary of pixel intensities.
#
# Shared with quantile_hm(), so the two views of one study cannot disagree
# about the underlying numbers.
#
# Tolerates heatmap_order = NA (discovery order) and a run named in
# heatmap_order that has no pixels (left NA rather than erroring); the previous
# inline version summarised every run and then indexed by `heatmap_order`, so
# the documented NA default subscripted by NA.
.quantile_matrix <- function(MSIobject, quant_val, heatmap_order = NA) {
  runs <- as.character(pData(MSIobject)$run)
  sample_names <- if (all(is.na(heatmap_order))) unique(runs)
                  else as.character(heatmap_order)
  features <- as.character(fData(MSIobject)$name)

  out <- matrix(NA_real_, nrow = length(features), ncol = length(sample_names),
                dimnames = list(features, sample_names))
  for (j in seq_along(sample_names)) {
    idx <- which(runs == sample_names[j])
    if (!length(idx)) next
    inten <- spectraData(MSIobject[, idx])[["intensity"]]
    out[, j] <- apply(inten, 1L, stats::quantile, probs = quant_val,
                      na.rm = TRUE)
  }
  out
}


# Colour limit that lets a small fraction of cells saturate.
#
# A limit set to the maximum spends the ramp on a handful of outliers and
# leaves the bulk of the data in a narrow band of near-identical colour.
.saturating_cap <- function(v, saturate = 0.04) {
  v <- abs(as.vector(v))
  v <- v[is.finite(v)]
  if (!length(v)) return(1)
  cap <- unname(stats::quantile(v, 1 - saturate, na.rm = TRUE))
  if (!is.finite(cap) || cap <= 0) cap <- max(v)
  if (!is.finite(cap) || cap <= 0) cap <- 1
  cap
}


# The two encoded channels, separated from the drawing so they can be tested
# without rendering a heatmap.
#
# Returns the per-cell group mean (hue), the per-cell opacity, and the two
# colour limits -- which are computed from different distributions on purpose;
# see the function docs.
.contribution_values <- function(mat, grp, alpha_floor = 0.6,
                                 saturate = 0.04) {
  # z-score per feature across ALL samples, not within group: within-group
  # scaling would centre every group on itself and erase the difference being
  # drawn.
  z <- t(apply(mat, 1L, function(x) {
    s <- stats::sd(x, na.rm = TRUE)
    # A feature measured identically everywhere carries no information here.
    # NA says so; 0 would draw it mid-ramp, indistinguishable from a feature
    # that genuinely sits at the mean.
    if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
    (x - mean(x, na.rm = TRUE)) / s
  }))
  if (ncol(mat) == 1L) z <- matrix(z, nrow = nrow(mat))
  dimnames(z) <- dimnames(mat)

  grp_levels <- unique(grp)
  gmean <- vapply(grp_levels, function(g)
    rowMeans(z[, which(grp == g), drop = FALSE], na.rm = TRUE),
    numeric(nrow(z)))
  gmean <- matrix(gmean, nrow = nrow(z))
  # rowMeans over an all-NA row returns NaN; keep it NA so it draws as missing.
  gmean[!is.finite(gmean)] <- NA_real_

  fill <- gmean[, match(grp, grp_levels), drop = FALSE]
  dimnames(fill) <- dimnames(z)

  # Signed towards the group mean, so "high" always means "moved with the
  # group" whichever way the group went. sign(0) is 0, zeroing the
  # contribution where the group mean is exactly 0 -- correct, since there is
  # no direction to agree with.
  contrib <- z * sign(fill)
  contrib[!is.finite(contrib)] <- NA_real_
  contrib_pos <- pmax(contrib, 0)

  cap_fill  <- .saturating_cap(gmean,       saturate)
  cap_alpha <- .saturating_cap(contrib_pos, saturate)

  alpha <- alpha_floor + (1 - alpha_floor) * pmin(contrib_pos / cap_alpha, 1)
  # A cell drawn in na_col must be opaque, or a missing measurement fades and
  # reads as a weak one.
  alpha[is.na(alpha)] <- 1
  dimnames(alpha) <- dimnames(z)

  list(z = z, gmean = gmean, fill = fill, alpha = alpha,
       cap_fill = cap_fill, cap_alpha = cap_alpha)
}


# Which of the two heatmaps section 2 draws.
#
# The contribution encoding answers "did the whole group do this, or one
# sample?" -- a question that only arises once a group has enough replicates
# for a mean to hide behind. At two per group the per-sample heatmap already
# shows both values directly and puts no summarising step between the reader
# and the data; at three or more, a group mean can rest on one replicate
# without saying so, which is exactly what the opacity channel exposes.
#
# Judged on the LARGEST group: if any group can hide a single-sample effect,
# the panel that reveals it is the more informative one for the whole study.
#
# Deliberately a named function rather than an inline condition: it is a
# judgement about study design, and it should be one line to change.
.heatmap_style <- function(style = "auto", heatmap_labs = NULL) {
  style <- match.arg(as.character(style),
                     c("auto", "per_sample", "contribution"))
  if (style != "auto") return(style)
  labs <- stats::na.omit(as.character(heatmap_labs))
  if (!length(labs)) return("per_sample")
  if (max(table(labs)) >= 3L) "contribution" else "per_sample"
}
