# quantMSImageR 0.99.11

First Bioconductor release (Bioconductor 3.24). quantMSImageR processes and
quantifies targeted DESI-MRM mass spectrometry imaging data, extending
Cardinal's `MSImagingExperiment` with a `quant_MSImagingExperiment` class that
carries calibration and tissue metadata.

## Reading and assembling acquisitions

* `readMRM()` reads Waters DESI-MRM acquisitions and names each transition
  from an ion-library CSV, matching precursor and product at unit resolution.
  Isomers the instrument cannot separate are named for all of them
  (`"A || B"`). Each pixel keeps its stage position.
* `selectTissuePixels()` draws the tissue/background mask, and `labelROIs()`
  draws regions of interest within the tissue (`airway_01`, `airway_02`, ...).
  Both are saved beside the acquisition and matched by pixel position, so one
  mask serves every MRM panel of the same area.
* `combineMSIs()` holds several sections in one object, each its own `run`;
  `bindPanels()` merges panels or polarities acquired over the same area; and
  `stitchAcquisitions()` joins pieces of one tissue by stage position.
  `alignFeatures()`, `trimMSI()` and `removeBlankMzs()` prepare them.

## Normalisation and filtering

* `int2response()` normalises to an internal standard per acquisition line,
  rolling window, sample or pixel, with an ion-library column mapping each
  analyte to its standard when there are several.
* `int2SNR()` and `applySNR()` filter each feature against the acquisition's
  own background pixels (median or mean), with optional per-feature thresholds.
  `int2SNR()` warns when a mostly-zero background leaves the median noise level
  set by the zero imputation rather than by the background.
  `back2NA()` and `zero2NA()` mark background and empty values as missing.
* `blockCentre()` corrects for a blocking factor -- patient, slide, sex -- by
  dividing each block by its own geometric mean or by its baseline arm.

## Quantification

* `summariseCalLevels()`, `createCalCurve()` and `int2conc()` convert response
  to amount (pg per pixel, and pg/mm^2 where the pixel size is known) against
  on-slide (`"cal"`) or on-tissue standard-addition (`"std_addition"`)
  standards, with 1/x, 1/x^2 or no weighting.
* `calibrationDiagnostics()`, `calibrationR2()` and `plotCalCoverage()` report
  the fit and whether the standards bracket the tissue signal. `int2conc()`
  stops when more than 10% of a feature's pixels would be extrapolated, unless
  told otherwise.

## Visualisation and summaries

* `imageR()` draws ion images with capped, histogram-equalised or square-root
  colour scales, in the palettes of `quantPalettes()`.
* `quantileHm()` draws per-sample quantile heatmaps. `contributionHm()` colours
  each cell by the group mean and shades it by the sample's contribution to it,
  so a group effect resting on one replicate is visible.
* `createMSIDataMatrix()`, `tissueData()` and `tissueMatrix()` return tables
  for downstream statistics.

## Whole-study workflow

* `runStudy()` runs a study from one YAML file: loading, normalisation, a sweep
  of SNR thresholds, optional calibration and ion ratios, per-feature `.txt`
  images for ImageJ and image registration, and one HTML report per threshold.
  The report holds heatmaps, group comparisons, colocalisation, ion images
  with an animated overview (also written as a `.gif`), regions of interest and
  an embedded tables workbook.
* `validateConfig()` checks a config -- paths, ion-library columns, enumerated
  values, calibration and internal-standard settings -- and reports every
  problem at once. `runStudy()` will not start on an invalid config.
* `generateTxtImages()` is the engine behind `runStudy()`. `runExample()` runs
  the workflow end to end on the bundled synthetic data: six sections in two
  groups, plus a calibration acquisition.

Changes made during development and review (0.99.0 to 0.99.10) are listed in
full at <https://github.com/MJS-708/quantMSImageR/blob/d1987f2/NEWS.md>.
Scripts written against a pre-release version from GitHub need updating for
one rename: function names became camelCase in 0.99.6 (`run_study()` is now
`runStudy()`, and so on). Study configs need no changes: the old YAML keys
`label:` and `baseline_label:` are still read.
