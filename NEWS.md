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
