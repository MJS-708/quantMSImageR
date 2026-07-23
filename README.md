[![](https://badgen.net/static/Publication/10.1021.acs.analchem.4c02350/green?.svg)](https://doi.org/10.1021/acs.analchem.4c02350) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10807654.svg)](https://doi.org/10.5281/zenodo.10807654)

# quantMSImageR <img src="man/figures/logo.png" align="right" height="139" alt="quantMSImageR logo" />

Software tools for processing and quantifying targeted multiple reaction monitoring (MRM) mass spectrometry imaging (MSI) data. Largely an extension to [Cardinal](https://bioconductor.org/packages/Cardinal/) for targeted data acquired by DESI-MRM, with additional functionality for quantification and outputs that facilitate multimodal imaging.

## Publication

**Development of a Desorption Electrospray Ionization–Multiple-Reaction-Monitoring Mass Spectrometry (DESI-MRM) Workflow for Spatially Mapping Oxylipins in Pulmonary Tissue**

Matthew J. Smith, Mu Nie, Mikael Adner, Jesper Säfholm, Craig E. Wheelock

*Analytical Chemistry* (2024). DOI: [10.1021/acs.analchem.4c02350](https://doi.org/10.1021/acs.analchem.4c02350)

## Requirements

- R ≥ 4.5.0
- [Cardinal](https://bioconductor.org/packages/Cardinal/) ≥ 3.6.2 (Bioconductor)
- [ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap/) (Bioconductor)

Install Bioconductor dependencies first:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Cardinal", "ComplexHeatmap"))
```

### macOS notes

- Install ImageMagick before R packages so the `magick` dependency
  compiles cleanly:

  ```bash
  brew install imagemagick
  ```

- In `select_tissue_pixels()`, finish ROI point selection with **Esc**
  (not right-click — that is the Windows convention). Make sure the plot
  opens in an external Quartz window, not the RStudio Plots pane, or
  clicks will not register.

## Installation

Install the latest release from GitHub:

```r
# install.packages("remotes")
remotes::install_github("MJS-708/quantMSImageR", ref = "main")
```

For the development version:

```r
remotes::install_github("MJS-708/quantMSImageR", ref = "dev")
```

Local development install:

```r
devtools::install("path/to/quantMSImageR")
```

## Quick start

Run the built-in example to verify your installation. The example uses a
synthetic dataset of two samples × three 20×20-pixel sections, with seven
oxylipins (negative-ion mode) and one internal standard:

```r
quantMSImageR::run_example()
```

This renders an HTML report and opens it directly in the RStudio viewer (or
your default browser). No file paths to configure.

## Usage

### Step 1 — Select tissue pixels (once per acquisition)

Before running the YAML workflow, each acquisition needs a `tissue_pixels.csv`
saved inside its `.raw` folder.  This is generated interactively by
`select_tissue_pixels()`, which displays all ion images so you can choose the
feature that best separates tissue from background, then opens an ROI-selection
window:

```r
select_tissue_pixels(
  name         = "sample_1",
  data_path    = "path/to/raw",
  lib_ion_path = "path/to/ion_library.csv"
)
```

This displays all features first — choose the one with the clearest
tissue/background contrast — then draws a selection window on that feature.
The resulting mask is saved automatically to
`{data_path}/sample_1.raw/tissue_pixels.csv`.

Repeat for each `.raw` acquisition. For acquisitions loaded from pre-saved
per-section RDS files (`sections:` in the YAML), the tissue mask is embedded
in the RDS and this step is not required.

### Step 2 — Run the workflow via YAML

Analysis is configured via a YAML file and run with a single command:

```r
library(quantMSImageR)

# Check the configuration before anything is read or written
validate_config("path/to/config.yaml")

run_study("path/to/config.yaml")
```

or from a shell:

```bash
Rscript -e 'quantMSImageR::run_study("path/to/config.yaml")'
```

This will:
1. Load and process all acquisitions defined in the YAML
2. Apply SNR filtering, IS normalisation and other user specified functionality (modular to easy to add other functionality)
3. Export per-feature ion image `.txt` files for each sample for instance for multimodal integration workflows
4. Render an HTML report with heatmaps, comparison plots, and ion images

A template YAML config is provided:

```r
file.copy(system.file("config_template.yaml", package = "quantMSImageR"), "config.yaml")
```

### YAML structure

```yaml
study: "study_name"

paths:
  data_path:    "path/to/raw/data"
  out_path:     "path/to/output"
  image_dir:    "path/to/ion/images"
  lib_ion_path: "path/to/ion_library.csv"
  markdown_rmd: "path/to/quantMSImageR___general_heatmap.Rmd"

samples:
  - neg: "acquisition_filename"   # without .raw extension
    label: "GroupA"
  - neg: "acquisition_filename2"
    label: "GroupB"
  - neg: "acquisition_filename3"
    label: "GroupA"               # duplicate labels = biological replicates

parameters:
  snr_thresh:  3
  baseline_label: "GroupA"

output:
  render_report: true
  report_fn: "my_study_SNR3"
```

### Multi-section acquisitions

When a single `.raw` folder contains multiple pre-saved per-section RDS files
(e.g. several tissue sections on one slide), list them under `sections:`. Each
section becomes its own run; sections sharing a `label` are treated as
biological replicates and pooled in comparison plots.

```yaml
samples:
  - neg: "SlideA"
    sections: ["section-01", "section-02", "section-03"]
    label: "GroupA"
```

Per-section labels are also supported (each section becomes its own group):

```yaml
samples:
  - neg: "SlideB"
    sections:
      - { file: "top", label: "Section01" }
      - { file: "mid", label: "Section02" }
      - { file: "btm", label: "Section03" }
```

### Dual-polarity acquisitions

A sample scanned in both polarities is given `pos:` and `neg:` together. The two
panels are matched on pixel coordinates and combined into one object, so
positive- and negative-mode features sit side by side in the same heatmaps and
reports:

```yaml
samples:
  - pos: "sample_1_pos"
    neg: "sample_1_neg"
    run_id: "S1"
    label:  "Ctrl"
```

Either polarity may also be a list, when a panel was split across several
acquisitions:

```yaml
samples:
  - pos:
      - "sample_2_pos_panelA"
      - "sample_2_pos_panelB"
    neg: "sample_2_neg"
    run_id: "S2"
    label:  "Treatment"
```

## Contact

For questions or comments please contact Matthew J. Smith (matthew.smith@ki.se || mattyjsmith123@gmail.com).
