#' @title Example data: ranked gene expression changes between early DN and controls
#'
#' @description This example data represents a ranked vector of gene expression 
#' changes (log2 fold change) between early diabetic nephropathy (DN) patients 
#' and healthy controls. It is based on DEG data curated from a published study
#' (doi: 10.2337/db19-0204).
#'
#' @name DN_data
#' @aliases DN_data
#' @format A vector with 748 DEGs and their log2fc
#' @docType data
#' @keywords data GSEA DEG diabetic_nephropathy
#' @returns \url{https://diabetesjournals.org
#' /diabetes/article/68/12/2301/39855
#' /Comparison-of-Kidney-Transcriptomic-Profiles-of}
#' @examples
#' data(DN_data)

DN_degs <- openxlsx::read.xlsx("DN_control.xlsx")
DN_degs <- DN_degs[, c("Symbol", "Log2Rat", "p.adj")]
DN_data <- DN_degs$Log2Rat
names(DN_data) <- DN_degs$Symbol
DN_data <- sort(DN_data, decreasing = TRUE)
save(DN_data, file = "DN_data.rda")


