# Provenance of the bundled example data

Everything under `inst/extdata` is **synthetic**. No patient, animal or
unpublished experimental data is distributed with this package. The two scripts
in this directory generate it, and are the authoritative record of how it was
made.

Run either from the package root:

```r
source("inst/scripts/generate_example_data.R")
source("inst/scripts/generate_cal_data.R")
```

## `generate_example_data.R` → `inst/extdata/example.raw/section0[1-6].RDS`

Six synthetic DESI-MRM sections on a 20 x 20 pixel grid, used by
`run_example()`, the vignettes and most of the runnable examples.

- Two groups of three sections: **A** (`section01`-`03`, circular tissue) and
  **B** (`section04`-`06`, square tissue at ~60% of A's levels). The distinct
  shapes make the groups immediately distinguishable in rendered ion images.
- Eight negative-mode oxylipin features: seven analytes plus one internal
  standard, matching `inst/extdata/example_ion_library.csv` by precursor and
  product m/z.
- Pixels are labelled `tissue_pixels` / `background_pixels` in
  `pData()$sample_name`, which is what `int2snr()` and `back2NA()` reference.
- `experimentData()$pixelSize` is set to 100 micrometres, without which
  `int2conc()` cannot produce the `pg_mm2` layer.
- Seeded (`set.seed(42)`, plus a per-section offset) so the data is
  reproducible.

Note that the `.raw` folder contains **only** RDS files -- there are no Waters
acquisition files, so any code path calling `read_mrm()` against it will fail.
The section-aware path in `generate_txt_images()` is what makes the example
work.

## `generate_cal_data.R` → `inst/extdata/cal_example.raw/`

A synthetic calibration standards acquisition (`cal_MSI.RDS`) plus its design
table (`calibration_metadata.csv`), used by the quantification vignette and the
calibration examples and tests.

- Five calibration levels (10-160 pg per spot) x three replicates, each spot
  covering 4 pixels, on a 20 x 20 grid.
- Pixels labelled `Cal`, `Tissue` or `Background` in `pData()$sample_type`, with
  a unique `identifier` per standard spot.
- Three analytes -- `12-HHTrE`, `12-HOPE`, `9-HOTE` -- deliberately named to
  match features in `example.raw`, so the fitted models can be carried across
  and applied to a study section.
- Response is linear in amount with known slope and intercept per analyte plus
  5% noise, so the recovered calibration models are checkable against the
  values in the script.
- Seeded (`set.seed(7)`).

If either script is changed, regenerate **both** sets of files and re-run
`devtools::test()`: several tests derive their expectations from the object
itself, but the analyte names are shared between the two datasets.
