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

#' Combined Human and Mouse TF-Target Interactions
#'
#' A comprehensive dataset containing high-confidence (levels A, B, and C) transcription factor (TF)
#' and target gene interactions for both Human and Mouse, derived from the DoRothEA database.
#'
#' @format A tibble with `r nrow(tf_targets)` rows and 5 columns:
#' \describe{
#'   \item{Species}{Character. The organism ("Human" or "Mouse").}
#'   \item{TF}{Character. Gene symbol of the transcription factor.}
#'   \item{Target}{Character. Gene symbol of the target gene.}
#'   \item{Confidence}{Character. Confidence level of the interaction source:
#'     \itemize{
#'       \item \strong{A}: High confidence (Literature-curated).
#'       \item \strong{B}: Moderate confidence (ChIP-seq evidence).
#'       \item \strong{C}: Medium confidence (TFBS predictions + ChIP-seq).
#'     }
#'   }
#'   \item{Mode_of_Regulation}{Numeric. Indicates if the TF activates or inhibits the target:
#'     \itemize{
#'       \item \strong{1}: Activation
#'       \item \strong{-1}: Inhibition
#'     }
#'   }
#' }
#' @source \url{https://saezlab.github.io/dorothea/}
#' @references
#' Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J.
#' Benchmark and integration of resources for the estimation of human transcription factor activities.
#' Genome Research. 2019. DOI: 10.1101/gr.240663.118.
"tf_targets"

#' Gut Microbiota-Metabolite-Target Axis Data (GutMGene v2.0)
#'
#' @description
#' A comprehensive list containing datasets that describe the relationships between gut microbiota,
#' their metabolites, and host target genes. This data allows for the exploration of the "Gut-Target" axis
#' by linking bacteria to metabolites and metabolites to host gene targets.
#'
#' The data is sourced from the GutMGene v2.0 database and processed into a "wide" format
#' using full outer joins.
#'
#' @details
#' This object is a named list containing two data frames:
#' \describe{
#'   \item{\code{human}}{Data derived from human samples.}
#'   \item{\code{mouse}}{Data derived from mouse samples.}
#' }
#'
#' You can access the specific dataset using \code{gutMGene$human} or \code{gutMGene$mouse}.
#'
#' Missing values (\code{NA}) in the \code{Bacteria} column indicate that the metabolite-target
#' relationship is known, but the specific producing bacteria species is not recorded in the
#' current database context. Similarly, \code{NA} in the \code{Target} column indicates a
#' known bacteria-metabolite pair with no currently recorded target gene.
#'
#' @format A list with 2 elements, where each element is a data frame with the following columns:
#' \describe{
#'   \item{Bacteria}{Name of the gut bacterium (e.g., "Akkermansia muciniphila").}
#'   \item{Bacteria_ID}{NCBI Taxonomy ID of the bacterium.}
#'   \item{Metabolite}{Name of the microbial metabolite (e.g., "Acetate").}
#'   \item{Metabolite_ID}{PubChem CID of the metabolite.}
#'   \item{Target}{Symbol of the host target gene (e.g., "FFAR2").}
#'   \item{Target_ID}{Entrez Gene ID of the target gene.}
#'   \item{Interaction}{Type of interaction (e.g., "activation", "inhibition").}
#'   \item{PMID_Bac_Met}{PubMed ID for the evidence linking Bacteria to Metabolite.}
#'   \item{PMID_Met_Target}{PubMed ID for the evidence linking Metabolite to Target.}
#' }
#'
#' @source \url{https://bio-computing.hrbmu.edu.cn/gutmgene/}
#'
#' @references
#' Qi, C., He, G., Qian, K., Guan, S., Li, Z., Liang, S., Liu, J., Ke, X., Zhang, S., Lu, M., Cheng, L., & Zhang, X. (2025).
#' GutMGene v2.0: an updated comprehensive database for target genes of gut microbes and microbial metabolites.
#' \emph{Nucleic Acids Research}, 53(D1), D783–D788. \doi{10.1093/nar/gkae1002}
#'
#' @note
#' This data is distributed under the CC BY-NC 4.0 license.
#' Please cite the original GutMGene publication when using this data.
#'
#' @usage data(gutMGene)
"gutMGene"

#'
#' A demo subset of a protein–protein interaction (PPI) network
#' derived from a diabetic kidney disease (DKD) dataset.
#'
#' @details
#' This object is an \code{igraph} graph containing a small
#' protein–protein interaction network, useful for demonstrating
#' downstream network analysis workflows.
#'
#' @format An \code{igraph} object with:
#' \describe{
#'   \item{58 vertices}{Each vertex corresponds to a gene symbol.}
#'   \item{Edges}{Represent protein–protein interactions.}
#'   \item{Vertex attributes}{Gene symbol and optional annotations.}
#'   \item{Edge attributes}{Interaction score or weight, if available.}
#' }
#'
#' @usage data(demo_ppi)
#'
#' @source Internal demo dataset.
"demo_ppi"
