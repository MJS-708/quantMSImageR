#!/usr/bin/env Rscript
# =============================================================================
# run_study.R  —  quantMSImageR DESI-MRM study runner (compatibility shim)
# =============================================================================
#
# The study runner is now the exported function `quantMSImageR::run_study()`.
# Prefer calling it directly:
#
#   library(quantMSImageR)
#   run_study("path/to/config.yaml")
#
# This script is kept so existing workflows keep working:
#
#   Rscript run_study.R path/to/config.yaml
#
#   CONFIG_FILE <- "path/to/config.yaml"
#   source(system.file("run_study.R", package = "quantMSImageR"))
#
# The YAML must follow the structure in inst/config_template.yaml.
# =============================================================================

library(quantMSImageR)

if (!exists("CONFIG_FILE")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0)
    stop("Provide a YAML config path:\n  Rscript run_study.R config.yaml")
  CONFIG_FILE <- args[1]
}

run_study(CONFIG_FILE)
