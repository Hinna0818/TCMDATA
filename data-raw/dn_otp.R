#' download target genes associated with DN from Open Targets Platform
dir <- "/Users/hinna/Desktop/yulab/casestudy/OT-EFO_0000401-associated-targets-2_21_2026-v25_12.tsv"
raw <- read.delim(dir, check.names = FALSE, stringsAsFactors = FALSE)

dn_otp_tbl <- data.frame(
  symbol = raw[["symbol"]],
  score = suppressWarnings(as.numeric(raw[["globalScore"]])),
  source = "OpenTargets",
  stringsAsFactors = FALSE
)

channel_cols <- setdiff(names(raw), c("symbol", "globalScore"))
for (col in channel_cols) {
  dn_otp_tbl[[col]] <- suppressWarnings(as.numeric(raw[[col]]))
}

dn_otp_tbl <- dn_otp_tbl[!is.na(dn_otp_tbl$symbol) &
                           nzchar(dn_otp_tbl$symbol) &
                           is.finite(dn_otp_tbl$score), , drop = FALSE]
dn_otp_tbl <- dn_otp_tbl[order(-dn_otp_tbl$score, dn_otp_tbl$symbol), , drop = FALSE]
dn_otp_tbl <- dn_otp_tbl[!duplicated(dn_otp_tbl$symbol), , drop = FALSE]
rownames(dn_otp_tbl) <- NULL

## Backward-compatible gene vector used by existing tutorials.
dn_otp <- dn_otp_tbl$symbol

usethis::use_data(dn_otp, overwrite = TRUE)
usethis::use_data(dn_otp_tbl, overwrite = TRUE)
