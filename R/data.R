#' Differential expression results for diabetic nephropathy (DN)
#'
#' A named numeric vector of log2 fold changes derived from GSE142025,
#' comparing kidney transcriptomes between advanced diabetic nephropathy
#' (DN) patients and control individuals.  
#'  
#' The differential expression analysis was performed using GEO2R
#' based on the raw count matrix from GSE142025, representing the
#' **full gene set** detected in human kidney tissue.
#'
#' @format A named numeric vector:
#' \describe{
#'   \item{name}{gene symbol}
#'   \item{value}{log2 fold change (advanced DN vs control)}
#' }
#'
#' @source GEO dataset [GSE142025](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142025),
#' analyzed using GEO2R (R 4.2.2, DESeq2 pipeline)
#'
#' @examples
#' data(DN_data)
#' head(DN_data)
"DN_data"


#' Differentially expressed genes (DEGs) in diabetic nephropathy
#'
#' A data frame of significantly differentially expressed genes from
#' GSE142025, comparing advanced diabetic nephropathy (DN) with control samples.
#'
#' @format A data frame with 4 columns:
#' \describe{
#'   \item{GeneID}{Ensembl gene identifier}
#'   \item{padj}{adjusted p-value (FDR)}
#'   \item{log2FoldChange}{log2 fold change (advanced DN vs control)}
#'   \item{Symbol}{gene symbol}
#' }
#'
#' @details DEGs were filtered using padj < 0.05 and |log2FC| > 1.5.
#'
#' @source GEO dataset [GSE142025](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142025)
#'
#' @examples
#' data(DN_deg)
#' head(DN_deg)
"DN_deg"


#' Docking example dataset
#'
#' A dataset containing randomly selected docking scores and predicted affinities for compound–target pairs.
#'
#' @format A data frame with N rows and M columns:
#' \describe{
#'   \item{compound}{Compound name or ID.}
#'   \item{target}{Protein target name.}
#'   \item{binding_energy}{Docking score (kcal/mol).}
#' }
#' @source Molecular docking results from TCMDATA internal workflow.
"dock_data"

#' Venn diagram demo dataset
#'
#' A small and random dataset demonstrating Venn analysis across gene sets.
#'
#' @format A named list of four character vectors:
#' \describe{
#'   \item{TCM}{Character vector of 1,000 genes related to traditional Chinese medicine ingredients.}
#'   \item{DEGs}{Character vector of 800 differentially expressed genes identified from transcriptomic analysis.}
#'   \item{Markers}{Character vector of 500 known biomarker genes.}
#'   \item{DB}{Character vector of 700 targets collected from curated databases.}
#' }
#' @source Generated for TCMDATA package demonstration.
"venn_demo_data"
