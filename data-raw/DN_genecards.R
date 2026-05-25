## load DN targets from GeneCards
dir <- "/Users/hinna/Desktop/yulab/casestudy/DN_GeneCards.csv"
raw <- read.csv(dir, check.names = FALSE, stringsAsFactors = FALSE)

dn_gcds_tbl <- data.frame(
  symbol = raw[["Gene Symbol"]],
  score = suppressWarnings(as.numeric(raw[["Relevance score"]])),
  source = "GeneCards",
  description = raw[["Description"]],
  category = raw[["Category"]],
  uniprot_id = raw[["Uniprot ID"]],
  genecards_id = raw[["GC Id"]],
  stringsAsFactors = FALSE
)

dn_gcds_tbl <- dn_gcds_tbl[!is.na(dn_gcds_tbl$symbol) &
                             nzchar(dn_gcds_tbl$symbol) &
                             is.finite(dn_gcds_tbl$score), , drop = FALSE]
dn_gcds_tbl <- dn_gcds_tbl[order(-dn_gcds_tbl$score, dn_gcds_tbl$symbol), , drop = FALSE]
dn_gcds_tbl <- dn_gcds_tbl[!duplicated(dn_gcds_tbl$symbol), , drop = FALSE]
rownames(dn_gcds_tbl) <- NULL

## Backward-compatible gene vector used by existing tutorials.
dn_gcds <- dn_gcds_tbl$symbol

usethis::use_data(dn_gcds, overwrite = TRUE)
usethis::use_data(dn_gcds_tbl, overwrite = TRUE)
