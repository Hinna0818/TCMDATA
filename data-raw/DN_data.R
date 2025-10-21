## prepare GSEA example data

DN_genes <- openxlsx::read.xlsx("./GSE142025_advanced_control.xlsx")
DN_genes <- DN_genes[, c(1, 2, 6, 8)]

# ---- full gene vector for GSEA ----
DN_genes <- DN_genes[!duplicated(DN_genes$Symbol), ]
DN_data <- DN_genes$log2FoldChange
names(DN_data) <- DN_genes$Symbol
DN_data <- sort(DN_data, decreasing = TRUE)

# ---- filter significant DEGs ----
DN_deg <- subset(DN_genes, padj < 0.05 & abs(log2FoldChange) > 1.5)

usethis::use_data(DN_data, overwrite = TRUE, compress = "xz")
usethis::use_data(DN_deg, overwrite = TRUE, compress = "xz")
