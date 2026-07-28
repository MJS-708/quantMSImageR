# quantMSImageR 0.99.4

* Transitions are now matched to the ion library at unit resolution
  (`mz_tolerance` default `0.4`, was `0.05`). A triple quadrupole running MRM
  selects Q1 and Q3 at unit resolution, so a product ion written as `308.1` and
  one written as `308.3` are the same measurement; matching more tightly than
  the instrument resolves made annotation depend on how many decimal places
  were typed into the library. Applies to `read_mrm()`, `build_feature_meta()`,
  `align_features()` and `bind_panels()`. Precursor **and** product must still
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

* Fixed `run_example()` / `generate_txt_images()`: trimming empty background
  borders no longer removes the last background pixels when the tissue fills a
  perfect rectangle, which had left `int2snr()` with nothing to reference.
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
* Added `bind_panels()` for coordinate-matched merging of cross-panel /
  cross-polarity acquisitions of the same tissue; `bind_polarities()` is now a
  deprecated alias.
* Added `build_feature_meta()` to join ion-library metadata to features by m/z.
* Coordinate-aware tissue masks (`select_tissue_pixels()` now writes `x`/`y`),
  portable across MRM panels of the same sample.
* Consolidated `read_mrm()` to handle single- and multi-analyte imaging folders.
* Calibration (`summarise_cal_levels()`, `create_cal_curve()`, `int2conc()`)
  can be driven from the study YAML via a `calibration:` block.
