# quantMSImageR 0.99.7

* `contributionHm()` gains `scale`. `"shared"` (default, unchanged behaviour)
  keeps every feature on one colour scale so features can be compared with each
  other; `"feature"` normalises each feature by its own extreme, so every
  column fills the ramp whatever its actual spread. The second is the reading
  [quantileHm()] gives, but keeping the group-mean hue and the contribution
  opacity -- use it to look along one feature, never across features. The
  legend states which scale is in force. Exposed in the study YAML as
  `parameters$heatmap_scale`, and checked by `validateConfig()`.
* Documented the colour system on both heatmaps. `?quantileHm` previously
  mentioned only that z-scores are "clipped to [-1, 1]" and never said the
  z-scoring is per feature, so nothing warned a reader that its colours are not
  comparable between features.
* `?quantileHm` now records that a feature with no variance across samples is
  drawn at the bottom of the ramp rather than as missing, so an entirely
  single-coloured column should be checked against the values. (The behaviour
  itself is unchanged in this release.)

# quantMSImageR 0.99.6

Changes in response to the Bioconductor package review (Contributions #115).

* **Function names are now camelCase throughout**, following the Bioconductor
  style guidance. The convention applied is: acronyms upper case (`MRM`, `MSI`,
  `SNR`, `NA`), other abbreviations camel (`Hm`, `Cal`). Renamed:
  `align_features()` to `alignFeatures()`, `bind_panels()` to `bindPanels()`,
  `build_feature_meta()` to `buildFeatureMeta()`, `combine_MSIs()` to
  `combineMSIs()`, `contribution_hm()` to `contributionHm()`,
  `create_cal_curve()` to `createCalCurve()`, `createMSIDatamatrix()` to
  `createMSIDataMatrix()`, `generate_txt_images()` to `generateTxtImages()`,
  `int2snr()` to `int2SNR()`, `plot_cal_coverage()` to `plotCalCoverage()`,
  `quant_palettes()` to `quantPalettes()`, `quantile_hm()` to `quantileHm()`,
  `read_mrm()` to `readMRM()`, `remove_blank_mzs()` to `removeBlankMzs()`,
  `run_example()` to `runExample()`, `run_study()` to `runStudy()`,
  `select_tissue_pixels()` to `selectTissuePixels()`, `stitch_acquisitions()`
  to `stitchAcquisitions()`, `summarise_cal_levels()` to
  `summariseCalLevels()`, `trim_MSI()` to `trimMSI()`, `validate_config()` to
  `validateConfig()`, and `zero2na()` to `zero2NA()`. The S3 class
  `quant_validation` is now `quantValidation`.

  **No deprecated aliases are provided.** The package has never appeared in a
  Bioconductor release, so there are no downstream users to protect and
  Bioconductor's deprecation cycle does not apply. Existing scripts must be
  updated. The S4 class `quant_MSImagingExperiment` is unchanged.

* **Argument checking made consistent.** The exported plain functions that did
  not validate their inputs now do: `alignFeatures()`, `bindPanels()`,
  `trimMSI()`, `buildFeatureMeta()`, and the summarising helper shared by
  `quantileHm()` and `contributionHm()`. Checks test for
  `MSImagingExperiment`, which `quant_MSImagingExperiment` extends, so both are
  accepted. `buildFeatureMeta()` previously returned a full-height frame of
  `NA` rows when the ion library lacked `precursor_mz` / `product_mz`; it now
  errors.

* `?alignFeatures` documented `obj1` and `obj2` as `quant_MSImagingExperiment`
  when the code has always accepted any `MSImagingExperiment`. The
  documentation was wrong and has been corrected.

* **Fixed two warnings emitted during report generation.**
  `out-of-range values treated as 0 in coercion to raw` came from rmarkdown's
  resource discovery: passing `intermediates_dir=` to `rmarkdown::render()`
  makes it render a preview copy of the template with a `.md` extension, which
  skips knitr, so an unexpanded `sprintf` placeholder in an HTML attribute
  reached `utils::URLdecode()`. The tag is now built with `paste0()`. The
  pandoc `Div ... unclosed, closing implicitly` warning came from HTML block
  fences emitted without surrounding blank lines, whose parse depended on the
  `markdown_in_html_blocks` extension; the fences are now padded, and the
  statements between them are guarded so an error cannot skip a closing tag.

* The sample-contribution key on `contributionHm()` is drawn in neutral grey
  and labelled low / high. It previously used the top colour of the fill
  palette, so red meant both "high z-score" and "contributed to the group
  mean".

* Long feature-class and sample-group legend labels are wrapped, and the report
  heatmap device is no longer capped at a width that silently cropped the
  right-hand legends off wide studies.

# quantMSImageR 0.99.5

* Added `contributionHm()`, a heatmap that encodes the group mean as hue and
  each sample's agreement with it as opacity, so a group mean resting on one
  replicate is visible rather than hidden.
* Added `stitchAcquisitions()` for assembling several acquisitions of one
  section into a single image.
* Fixed the `combine` check in `validateConfig()`, which referenced the wrong
  object and sat in the ion-library block rather than the samples block.

# quantMSImageR 0.99.4

* Transitions are now matched to the ion library at unit resolution
  (`mz_tolerance` default `0.4`, was `0.05`). A triple quadrupole running MRM
  selects Q1 and Q3 at unit resolution, so a product ion written as `308.1` and
  one written as `308.3` are the same measurement; matching more tightly than
  the instrument resolves made annotation depend on how many decimal places
  were typed into the library. Applies to `readMRM()`, `buildFeatureMeta()`,
  `alignFeatures()` and `bindPanels()`. Precursor **and** product must still
  both agree.
* `ambiguity` gains a `"combine"` mode, now the default. A transition matching
  several library entries is named for all of them, joined with `" || "` --
  e.g. `"LTC4 || 14_15-LTC4"`, isomers 0.03 Da apart on the product ion that no
  acquisition can separate. Annotation columns come from the closest entry.
  `"error"`, `"warn"` and `"nearest"` are unchanged.
* Combining across entries of differing `Type` is refused: one feature cannot
  be both analyte and internal standard.
* Fixed the `show` methods export, which was declared in NAMESPACE without the
  generic being imported from `methods`.

# quantMSImageR 0.99.3

* Fixed `runExample()` / `generateTxtImages()`: trimming empty background
  borders no longer removes the last background pixels when the tissue fills a
  perfect rectangle, which had left `int2SNR()` with nothing to reference.
* Fixed the HTML report for single-acquisition studies: the run-to-group map is
  now built in an always-run chunk, so single-sample reports render.

# quantMSImageR 0.99.2

* Report heatmap now shows the sample-group colour legend (the group colour bar
  was previously drawn without a key). Legends are merged into one column, which
  also avoids a ComplexHeatmap drawing error seen with multiple annotations.

# quantMSImageR 0.99.1

* Added `show()` methods for the `calibrationInfo`, `tissueInfo` and
  `quant_MSImagingExperiment` classes.
* Enabled R-universe build checking.
* Vignette, help-page and package-title clarifications; documentation tidy-up.

# quantMSImageR 0.99.0

* First Bioconductor submission.
* Added per-analyte SNR overrides (`snr_overrides`) and per-report `SNR_used`
  reporting.
* Added `bindPanels()` for coordinate-matched merging of cross-panel /
  cross-polarity acquisitions of the same tissue.
* Added `buildFeatureMeta()` to join ion-library metadata to features by m/z.
* Coordinate-aware tissue masks (`selectTissuePixels()` now writes `x`/`y`),
  portable across MRM panels of the same sample.
* Consolidated `readMRM()` to handle single- and multi-analyte imaging folders.
* Calibration (`summariseCalLevels()`, `createCalCurve()`, `int2conc()`)
  can be driven from the study YAML via a `calibration:` block.
