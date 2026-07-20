source("C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/Bioinformatics/quantMSImageR/build_ion_library.R")
source("C:/Users/matsmi/OneDrive - Karolinska Institutet/Dokument/Bioinformatics/quantMSImageR/search_ion_library.R")

search_library(name = "PGE2")
search_library(precursor = 351.2, product = 271.1, tol = 0.05)
search_library(lipid = "Oxylipin", polarity = "Negative")

# TEST the fn
fn = "D:/__STUDIES/021_HDM_mice/ion_library.csv"

q <- readr::read_csv(fn)
a = batch_search(q, tol = 0.05)
