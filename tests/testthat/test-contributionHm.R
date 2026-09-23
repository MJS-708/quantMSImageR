# The point of this plot is that a group mean never hides which sample made
# it. The maths is tested through .contribution_values(), which is separated
# from the drawing precisely so it can be checked without rendering.

mk <- function(runs = c("s1", "s2", "s3", "s4"), px = 3, values = NULL) {
  n_feat <- if (is.null(values)) 4L else nrow(values)
  m <- if (is.null(values)) matrix(1, n_feat, length(runs) * px) else values
  obj <- MSImagingExperiment(
    spectraData = m,
    featureData = MassDataFrame(mz = seq_len(n_feat),
                                name = paste0("f", seq_len(n_feat))),
    pixelData = PositionDataFrame(
      coord = expand.grid(x = seq_len(px), y = seq_along(runs)),
      run = factor(rep(runs, each = px), levels = runs)))
  as(obj, "quant_MSImagingExperiment")
}

# f1: s4 alone is high -- the group mean rests on one sample.
# f2: s3 and s4 both high  -- a real group effect.
# f3: flat everywhere      -- no variance.
vals <- function(px = 3) {
  rbind(
    f1 = c(rep(1, px), rep(1, px), rep(1, px), rep(9, px)),
    f2 = c(rep(1, px), rep(1, px), rep(9, px), rep(9, px)),
    f3 = rep(5, 4 * px),
    f4 = c(rep(2, px), rep(3, px), rep(4, px), rep(5, px)))
}

ord <- c("s1", "s2", "s3", "s4")
lab <- c("ctrl", "ctrl", "trt", "trt")

vv <- function(...) {
  m <- .quantile_matrix(mk(values = vals()), 0.5, ord)
  .contribution_values(m, lab, ...)
}

test_that("hue is the group mean and is constant within a group", {
  v <- vv()
  # Both treated samples carry the same value: their group's mean.
  expect_equal(v$fill["f1", "s3"], v$fill["f1", "s4"])
  expect_equal(v$fill["f2", "s3"], v$fill["f2", "s4"])
  # And it is genuinely the mean of their z-scores.
  expect_equal(v$fill["f1", "s3"], mean(v$z["f1", c("s3", "s4")]))
})

test_that("opacity separates a one-sample effect from a whole-group effect", {
  a <- vv()$alpha
  # f1: only s4 is high, so within the treated group s4 carries the mean and
  # s3 does not.
  expect_gt(a["f1", "s4"], a["f1", "s3"])
  # f2: both treated samples are high, so neither is singled out.
  expect_equal(a["f2", "s3"], a["f2", "s4"])
  # ... and the pair agrees more than the split pair does.
  expect_gt(a["f2", "s3"], a["f1", "s3"])
})

test_that("a zero-variance feature is NA, not zero", {
  # Drawing it as 0 would put it mid-ramp, indistinguishable from a feature
  # that genuinely sits at the mean.
  v <- vv()
  expect_true(all(is.na(v$fill["f3", ])))
  # NA cells are drawn opaque, or the missing-value colour fades into nothing.
  expect_true(all(v$alpha["f3", ] == 1))
})

test_that("alpha never falls below the floor", {
  for (fl in c(0, 0.6, 0.9)) {
    a <- vv(alpha_floor = fl)$alpha
    expect_gte(min(a), fl)
    expect_lte(max(a), 1)
  }
})

test_that("the two colour limits are computed separately", {
  # Averaging replicates shrinks the values by roughly sqrt(n), so a limit
  # inherited from the per-sample z-scores would leave the panel washed out.
  # The group-mean limit must therefore be the smaller of the two.
  v <- vv()
  expect_lt(v$cap_fill, max(abs(v$z), na.rm = TRUE))
  expect_gt(v$cap_fill, 0)
  expect_gt(v$cap_alpha, 0)
})

test_that("bad arguments are refused rather than silently coerced", {
  obj <- mk(values = vals())
  expect_error(contributionHm(obj, heatmap_order = ord,
                               heatmap_labs = c("a", "b")),
               regexp = "heatmap_labs has 2 entries")
  expect_error(contributionHm(obj, heatmap_order = ord, alpha_floor = 1),
               regexp = "alpha_floor")
  expect_error(contributionHm(obj, heatmap_order = ord, saturate = 0),
               regexp = "saturate")
})

test_that("it is the same panel as quantileHm, differing only in the cells", {
  # This is the whole design constraint: a reader must be able to put the two
  # side by side and trust that a colour bar means the same thing in both.
  obj <- mk(values = vals())
  fs  <- c("A", "A", "B", "B")
  q <- quantileHm(obj, quant_val = 0.5, heatmap_order = ord,
                   heatmap_labs = lab, feature_split = fs)
  cn <- contributionHm(obj, quant_val = 0.5, heatmap_order = ord,
                        heatmap_labs = lab, feature_split = fs)

  expect_s4_class(cn, "Heatmap")
  expect_equal(dim(cn@matrix), dim(q@matrix))
  expect_equal(rownames(cn@matrix), rownames(q@matrix))
  expect_equal(colnames(cn@matrix), colnames(q@matrix))
  # Same splits, so the slices line up.
  expect_equal(cn@row_order_list, q@row_order_list)
  expect_equal(cn@column_order_list, q@column_order_list)
})

test_that("the opacity key travels with the heatmap", {
  # Opacity is an encoded channel; ComplexHeatmap builds legends from the
  # colour mapping, which knows nothing about alpha, so the key is attached.
  cn <- contributionHm(mk(values = vals()), heatmap_order = ord,
                        heatmap_labs = lab)
  lgd <- attr(cn, "contribution_legend")
  expect_false(is.null(lgd))
  expect_s4_class(lgd, "Legends")
})

test_that("a single ungrouped sample yields no scorable feature", {
  # One acquisition cannot be z-scored, so every cell is missing rather than
  # silently drawn at the ramp midpoint. The report checks for this and says
  # so instead of drawing an empty panel.
  m <- .quantile_matrix(mk(runs = "s1",
                           values = matrix(c(1, 2, 3, 4), 4, 3)), 0.5, "s1")
  v <- .contribution_values(m, "All")
  expect_true(all(is.na(v$fill)))
})

test_that("auto picks the heatmap from group size, not group count", {
  # The contribution encoding matters once a mean has enough replicates to
  # hide behind. Two groups of three is exactly the case that motivated it,
  # and an earlier version of this rule read it as "three groups" and drew
  # the wrong panel for the example study.
  expect_equal(.heatmap_style("auto", c("A", "A", "A", "B", "B", "B")),
               "contribution")
  expect_equal(.heatmap_style("auto", c("A", "A", "B", "B")), "per_sample")
  # Judged on the largest group: if any group can hide it, show it.
  expect_equal(.heatmap_style("auto", c("A", "A", "B", "B", "B")),
               "contribution")
  # Many groups of one is not replication.
  expect_equal(.heatmap_style("auto", c("A", "B", "C", "D")), "per_sample")

  # An explicit choice always wins, and nothing else is accepted.
  expect_equal(.heatmap_style("per_sample", c("A", "A", "A")), "per_sample")
  expect_equal(.heatmap_style("contribution", c("A", "B")), "contribution")
  expect_error(.heatmap_style("sideways", c("A", "B")))
})

# ---- Blocks (regions of interest) ------------------------------------------

test_that("a block is scored against itself, not against the panel", {
  # Two blocks whose levels differ wildly: scored across the panel, block b
  # would be uniformly high and say nothing about which sample within it is.
  m <- .quantile_matrix(mk(values = vals()), 0.5, ord)
  m[, c(3, 4)] <- m[, c(3, 4)] * 100          # block b on a different scale

  blocked <- .contribution_values(m, grp = c("ctrl", "trt", "ctrl", "trt"),
                                  block = c("a", "a", "b", "b"))
  # Within each block of two, the two z-scores are mirror images: each block
  # was centred on its own mean rather than on the panel's. f4 is the feature
  # that varies inside both blocks.
  expect_equal(unname(blocked$z["f4", 1] + blocked$z["f4", 2]), 0,
               tolerance = 1e-8)
  expect_equal(unname(blocked$z["f4", 3] + blocked$z["f4", 4]), 0,
               tolerance = 1e-8)
  # f1 varies across the panel (s4 alone is high) but is flat inside block a,
  # so it cannot be scored there -- scoring within a block means exactly this.
  expect_true(all(is.na(blocked$z["f1", 1:2])))
  expect_false(anyNA(blocked$z["f1", 3:4]))
  # And multiplying one block by 100 cannot reach the other block's scores.
  plain <- .contribution_values(m, grp = c("ctrl", "trt", "ctrl", "trt"),
                                block = c("a", "a", "a", "a"))
  expect_false(isTRUE(all.equal(blocked$z[, 1], plain$z[, 1])))
})

test_that("a group mean belongs to its group within its block", {
  m <- .quantile_matrix(mk(values = vals()), 0.5, ord)
  v <- .contribution_values(m, grp = c("ctrl", "trt", "ctrl", "trt"),
                            block = c("a", "a", "b", "b"))
  # Four block-group combinations, so the fill of column 1 (block a, ctrl) is
  # that cell's own z -- a mean over one sample -- not a mean over both ctrls.
  expect_equal(unname(v$fill[, 1]), unname(v$z[, 1]))
})

test_that("contributionHm draws the region bar beside the group bar", {
  obj <- mk(values = vals())
  hm <- contributionHm(obj, quant_val = 0.5, heatmap_order = ord,
                       heatmap_labs = c("ctrl", "trt", "ctrl", "trt"),
                       sample_block = c("a", "a", "b", "b"))
  expect_s4_class(hm, "Heatmap")
  expect_equal(names(hm@left_annotation@anno_list), c("Region", "Group"))
  expect_equal(as.character(hm@matrix_param$row_split[[1]]),
               c("a", "a", "b", "b"))
})

test_that("sample_block must have one entry per sample", {
  expect_error(contributionHm(mk(values = vals()), quant_val = 0.5,
                              heatmap_order = ord, heatmap_labs = lab,
                              sample_block = c("a", "b")),
               regexp = "sample_block has")
})
