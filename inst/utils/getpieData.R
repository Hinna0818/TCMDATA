#' Generate Pie Chart Data for Network Visualization from Enrichment Results
#'
#' This function extracts gene-term relationships from a clusterProfiler enrichment result
#' and formats them into a matrix suitable for scatterpie-based node pie chart visualization.
#'
#' @param enrich_obj An object of class `enrichResult`, typically from `enrichGO()`, `enrichKEGG()`, or `enricher()`.
#' @param ppi_genes A character vector of gene symbols or IDs present in the PPI network (i.e., `V(ppi)$name`).
#' @param top_n Integer. The number of top enrichment terms (e.g., GO terms) to include. Default is 5.
#' @param use_weight Logical. If `TRUE`, uses enrichment significance as weights (e.g., -log10(p.adjust)); otherwise binary (0/1). Default is `FALSE`.
#' @param weight_scale Character. One of `"logp"` (default) or `"invp"`. Defines the weighting method if `use_weight = TRUE`:
#'  - `"logp"`: use `-log10(p.adjust)`
#'  - `"invp"`: use `1 / p.adjust`
#'
#' @return A data frame with genes in rows and selected enrichment terms in columns. Values represent either binary membership or weighted scores.
#'

getPieData <- function(
    enrich_obj, 
    ppi_genes, 
    top_n = 5, 
    use_weight = FALSE, 
    weight_scale = c("logp", "invp")) {
  
  stopifnot(inherits(enrich_obj, "enrichResult"))
  
  # 1. extract top n pathways
  enrich_df <- enrich_obj@result %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n)
  
  # 2. extract the name of term and gene id
  enrich_terms <- enrich_df$Description
  gene_lists <- strsplit(enrich_df$geneID, "/")
  names(gene_lists) <- enrich_terms
  
  all_genes <- unique(unlist(gene_lists))
  
  # 3. set gene-term df
  score_df <- data.frame(name = all_genes, stringsAsFactors = FALSE)
  
  for (i in seq_along(enrich_terms)) {
    term <- enrich_terms[i]
    genes_in_term <- gene_lists[[i]]
    
    if (use_weight) {
      if (weight_scale[1] == "logp") {
        score_df[[term]] <- ifelse(score_df$name %in% genes_in_term,
                                   -log10(enrich_df$p.adjust[i] + 1e-10), 0)
      } else if (weight_scale[1] == "invp") {
        score_df[[term]] <- ifelse(score_df$name %in% genes_in_term,
                                   1 / (enrich_df$p.adjust[i] + 1e-10), 0)
      } else {
        stop("Unknown weight_scale. Choose 'logp' or 'invp'")
      }
    } else {
      score_df[[term]] <- as.integer(score_df$name %in% genes_in_term)
    }
  }
  
  # 4. filter the nodes in ppi
  pie_data <- score_df %>%
    dplyr::filter(name %in% ppi_genes)
  
  return(pie_data)
}
