#' Example dataset: DN_data
#'
#' @description
#' A ranked numeric vector representing differential gene expression (log₂ fold change)
#' between early diabetic nephropathy (DN) patients and healthy controls.
#' The gene names correspond to human gene symbols, and the numeric values indicate
#' the direction and magnitude of expression changes.
#'
#'
#' @format A named numeric vector of length *N*:
#' \describe{
#'   \item{names}{Character vector of gene symbols (e.g., `"TUSC5"`, `"ADIPOQ"`, `"CIDEC"`).}
#'   \item{values}{Numeric values indicating log₂ fold change between DN and control groups.}
#' }
#'
#' @details
#' Higher positive values indicate genes upregulated in DN patients, while
#' negative values (if present) would indicate downregulation.  
#' The dataset is curated from publicly available transcriptomic data 
#' (see: \doi{10.2337/db19-0204}) and standardized for demonstration use.
#'
#' @source Curated and processed within the TCMDATA analysis workflow.
"DN_data"


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
