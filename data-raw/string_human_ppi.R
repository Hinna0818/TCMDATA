# Build package data object for the default STRING human PPI background.
#
# The processed RDS is created by dev/string/prepare_string_human_ppi.R from:
# https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz
# https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz

path <- "dev/string/string_human_ppi_v12.0_score700_lcc_igraph.rds"
if (!file.exists(path)) {
  stop("Processed STRING PPI RDS not found: ", path)
}

string_human_ppi <- readRDS(path)
usethis::use_data(string_human_ppi, overwrite = TRUE, compress = "xz")
