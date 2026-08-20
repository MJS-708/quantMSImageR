remove.packages("quantMSImageR")
unlink("C:/R/R-4.5.3/library/quantMSImageR", recursive = TRUE, force = TRUE)
.rs.restartR()

"quantMSImageR" %in% loadedNamespaces()   # must now be FALSE

devtools::document()
devtools::build(vignettes = TRUE)     # tarball with vignettes built from .Rmd
devtools::install(build_vignettes = TRUE, upgrade = F)

devtools::check(vignettes = TRUE)
BiocCheck::BiocCheck()

vignette("quantMSImageR", package = "quantMSImageR")   # the installed one
pkgdown::build_site()                                  # writes docs/, local only

quantMSImageR::runExample()

#remotes::install_github("MJS-708/quantMSImageR", ref = "main")

library(quantMSImageR)

quantMSImageR::runExample()
