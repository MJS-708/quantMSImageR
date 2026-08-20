#!/usr/bin/env Rscript
# =============================================================================
# runStudy.R  —  quantMSImageR DESI-MRM study runner (compatibility shim)
# =============================================================================
#
# The study runner is now the exported function `quantMSImageR::runStudy()`.
# Prefer calling it directly:
#
#   library(quantMSImageR)
#   runStudy("path/to/config.yaml")
#
# This script is kept so existing workflows keep working:
#
#   Rscript runStudy.R path/to/config.yaml
#
#   CONFIG_FILE <- "path/to/config.yaml"
#   source(system.file("runStudy.R", package = "quantMSImageR"))
#
# The YAML must follow the structure in inst/config_template.yaml.
# =============================================================================

library(quantMSImageR)

if (!exists("CONFIG_FILE")) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0)
    stop("Provide a YAML config path:\n  Rscript runStudy.R config.yaml")
  CONFIG_FILE <- args[1]
}

runStudy(CONFIG_FILE)
