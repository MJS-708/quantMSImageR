remove.packages("quantMSImageR")
unlink("C:/R/R-4.5.3/library/quantMSImageR", recursive = TRUE, force = TRUE)
.rs.restartR()

"quantMSImageR" %in% loadedNamespaces()   # must now be FALSE

devtools::document(
  "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/Bioinformatics/quantMSImageR"
)

devtools::test(
  "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/Bioinformatics/quantMSImageR"
)

devtools::install(
  "C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/Bioinformatics/quantMSImageR",
  upgrade = F,
  build = TRUE,
  force = TRUE
)

#remotes::install_github("MJS-708/quantMSImageR", ref = "main")

library(quantMSImageR)

quantMSImageR::run_example()
