#' Group-mean heatmap with per-sample contribution as opacity
#'
#' The same heatmap as [quantileHm()] -- same orientation, same group and
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
#' [quantileHm()], then z-scored **per feature across all samples**, not
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
#' this differs from [quantileHm()], which clips to a fixed \[-1, 1\].
#'
#' Both limits above are single numbers covering the whole panel, which is what
#' lets one column be compared with another. `scale = "feature"` replaces them
#' with a limit per feature, so each column is normalised by its own extreme
#' and fills the ramp whatever its actual spread. That is the same reading
#' [quantileHm()] gives, but keeping the group-mean hue and the contribution
#' opacity: use it to look along a single feature, never across features.
#'
#' @section Comparing features with each other:
#'
#' One cap or many. That is the whole of it.
#'
#' `scale = "shared"` computes a single limit for the panel and draws every
#' column on it, so a feature whose group mean reaches the cap saturates and one
#' at a fifth of it stays pale. The pale column is pale because it moved less,
#' and that difference is the information: features can be ranked against each
#' other by how far their groups separate, in units of each feature's own
#' spread.
#'
#' `scale = "feature"` gives every column its own limit, its own extreme, so
#' each fills the ramp whatever it actually did. Nothing is pale, so nothing can
#' be ranked. A feature with no group effect at all looks exactly as vivid as
#' the strongest one in the study -- its largest group mean might be 0.05, and
#' it is still drawn at full intensity. Use it to read one feature's profile
#' across samples; never to argue that a feature matters.
#'
#' What is being compared under `"shared"` is **effect size, not fold change**.
#' The hue is a group separation measured in standard deviations of that
#' feature, so a feature that doubles but varies wildly between replicates can
#' sit below one that moves 20% in every sample. For magnitude, read the
#' fold-change table instead. The two orderings are both correct and answer
#' different questions.
#'
#' For a balanced two-group design the group-mean z-score is bounded at
#' \eqn{\pm 1}: perfectly separated groups with no within-group scatter put
#' every sample at \eqn{z = \pm 1}, so the group means land on the bound, and
#' identical groups put them on 0. A cap of 0.5 therefore means the strongest
#' feature in the panel reaches about half the separation the design can show.
#'
#' Two limits on the comparison, both worth remembering. The standard deviation
#' is estimated from the samples in hand, so with a handful per group those
#' estimates are noisy and the ranking is coarse -- read it as "which features
#' stand out", not as a precise order. And the deviation is pooled across
#' groups, so a large separation inflates its own denominator; the scale
#' compresses at the top and is not linear in effect size there.
#'
#' @import Cardinal
#' @include setClasses.R
#'
#' @inheritParams quantileHm
#' @param scale Character. How the colour scale is set.
#'   `"shared"` (default) puts every feature on ONE scale, capped at a
#'   quantile of all the group means, so a strongly separating feature is
#'   vivid and a flat one is pale -- colours mean the same thing in every
#'   column and features can be compared with each other.
#'   `"feature"` normalises each feature by its own extreme, so every column
#'   fills the ramp regardless of how far it actually moved. That shows the
#'   pattern within a feature and makes comparison BETWEEN features
#'   meaningless -- a feature varying by a few percent looks exactly like one
#'   that doubles. Use `"shared"` to ask which features changed most, and
#'   `"feature"` to read each feature's own profile.
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
#' hm <- contributionHm(obj, quant_val = 0.5,
#'                       heatmap_order = "section01", heatmap_labs = "A")
#'
#' @seealso [quantileHm()] for the per-sample view of the same matrix.
#' @family visualisation
#' @export
contributionHm <- function(MSIobject, quant_val = 0.5,
                            scale = c("shared", "feature"),
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
    stop("contributionHm: alpha_floor must be a single number in [0, 1).",
         call. = FALSE)
  if (!is.numeric(saturate) || length(saturate) != 1L || is.na(saturate) ||
      saturate <= 0 || saturate >= 0.5)
    stop("contributionHm: saturate must be a single number in (0, 0.5).",
         call. = FALSE)

  mat <- .quantile_matrix(MSIobject, quant_val, heatmap_order,
                          context = "contributionHm")
  samples <- colnames(mat)

  grp <- if (all(is.na(heatmap_labs))) rep("All", length(samples))
         else as.character(heatmap_labs)
  if (length(grp) != length(samples))
    stop("contributionHm: heatmap_labs has ", length(grp),
         " entries but there are ", length(samples), " samples.", call. = FALSE)

  scale <- match.arg(scale)
  v <- .contribution_values(mat, grp, alpha_floor = alpha_floor,
                            saturate = saturate, scale = scale)

  ## ---- same orientation and layout as quantileHm --------------------------
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

  cols   <- quantPalettes(palette)
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
    fill_t,
        # The two scales are different readings; the key must say which.
        name = if (identical(scale, "feature"))
                 "Group mean\n(per feature)" else "Group mean\nz-score",
        col = col_fn, na_col = na_col,
    cell_fun = cell_fun), layout))

  # The opacity channel is an encoding, so it needs a key. ComplexHeatmap
  # legends for a heatmap come from its colour mapping, which knows nothing
  # about alpha, so this one is built here and attached for draw() to place.
  # Neutral grey, NOT the palette's top colour. The fill scale already uses red
  # for a high z-score, so drawing the opacity key in red made the same hue mean
  # two different things -- a reader sees a red swatch and reads "high value"
  # when it means "this sample agreed with its group". Grey carries the
  # light-to-dark reading without claiming a direction.
  ref <- "#3F3F3F"
  attr(hm, "contribution_legend") <- ComplexHeatmap::Legend(
    title = "Contribution to\ngroup mean",
    at = c(1, 2), labels = c("low", "high"),
    legend_gp = grid::gpar(
      fill = c(grDevices::adjustcolor(ref, alpha.f = alpha_floor), ref)),
    title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = fontsize))
  hm
}


# Per-feature, per-sample summary of pixel intensities.
#
# Shared with quantileHm(), so the two views of one study cannot disagree
# about the underlying numbers.
#
# Tolerates heatmap_order = NA (discovery order) and a run named in
# heatmap_order that has no pixels (left NA rather than erroring); the previous
# inline version summarised every run and then indexed by `heatmap_order`, so
# the documented NA default subscripted by NA.
.quantile_matrix <- function(MSIobject, quant_val, heatmap_order = NA,
                             context = ".quantile_matrix") {
  # missing() first, and before any use of quant_val: quantileHm() declares it
  # with no default, so evaluating it when the caller omitted it raises R's own
  # "argument is missing" from inside stats::quantile(), naming neither the
  # heatmap nor the argument the user actually has to supply.
  if (missing(quant_val))
    stop(context, ": quant_val is required -- the quantile of pixel ",
         "intensities to summarise each sample by, e.g. 0.5 for the median.",
         call. = FALSE)
  if (!is(MSIobject, "MSImagingExperiment"))
    stop(context, ": MSIobject must be an MSImagingExperiment ",
         "(quant_MSImagingExperiment included); got ",
         class(MSIobject)[1L], ".", call. = FALSE)
  if (!is.numeric(quant_val) || length(quant_val) != 1L ||
      is.na(quant_val) || quant_val < 0 || quant_val > 1)
    stop(context, ": quant_val must be a single number in [0, 1].",
         call. = FALSE)

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
                                 saturate = 0.04,
                                 scale = c("shared", "feature")) {
  scale <- match.arg(scale)
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

  if (identical(scale, "feature")) {
    # Every feature normalised by its OWN extreme, so each column uses the full
    # ramp whatever its actual spread. That is the point: it shows the pattern
    # within a feature and makes comparison BETWEEN features meaningless, which
    # is the opposite of the shared cap below. A feature that varies by a few
    # percent then looks exactly like one that doubles.
    .rowcap <- function(m) {
      v <- apply(abs(m), 1L, function(r) {
        r <- r[is.finite(r)]
        if (!length(r)) return(NA_real_)
        max(r)
      })
      v[!is.finite(v) | v <= 0] <- NA_real_
      v
    }
    # m / v with length(v) == nrow(m) divides each ROW by its own cap.
    .cf <- .rowcap(gmean);       fill        <- fill        / ifelse(is.na(.cf), 1, .cf)
    .ca <- .rowcap(contrib_pos); contrib_pos <- contrib_pos / ifelse(is.na(.ca), 1, .ca)
    cap_fill  <- 1
    cap_alpha <- 1
  } else {
    cap_fill  <- .saturating_cap(gmean,       saturate)
    cap_alpha <- .saturating_cap(contrib_pos, saturate)
  }

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
