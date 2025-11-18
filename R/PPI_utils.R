#' Subset PPI network by edge score and top-n node degree
#'
#' @param ppi_obj An igraph object with edge attribute 'score'
#' @param n Number of top-degree nodes to keep. Default is NULL.
#' @param score_cutoff Minimum edge score to keep. Default is 0.7.
#'
#' @return A subgraph of the original PPI network
#' @importFrom igraph degree E subgraph_from_edges vcount induced_subgraph
#' @export

ppi_subset <- function(ppi_obj, n = NULL, score_cutoff = 0.7) {

  stopifnot(inherits(ppi_obj, "igraph"))
  score <- E(ppi_obj)$score

  # edge score filter
  if (is.null(score)) {
    stop("Edges must have a 'score' attribute.")
  }

  ppi_filtered <- igraph::subgraph_from_edges(ppi_obj, eids = E(ppi_obj)[score >= score_cutoff], delete.vertices = TRUE)

  if (vcount(ppi_filtered) == 0) {
    warning("No nodes left after edge score filtering.")
    return(ppi_filtered)
  }

  # degree filter
  if (!is.null(n)){
    deg <- igraph::degree(ppi_filtered, mode = "all")
    deg <- deg[!is.na(deg)]
    top_nodes <- names(sort(deg, decreasing = TRUE))[1:min(n, length(deg))]
    ppi_filtered <- igraph::induced_subgraph(ppi_filtered, vids = top_nodes)
  }

  return(ppi_filtered)
}


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
#' @importFrom dplyr arrange slice_head filter
#' @export

getPieData <- function(
    enrich_obj,
    ppi_genes,
    top_n = 5,
    use_weight = FALSE,
    weight_scale = c("logp", "invp")) {

  stopifnot(inherits(enrich_obj, "enrichResult"))
  weight_scale <- match.arg(weight_scale)

  # 1. extract top n pathways
  enrich_df <- enrich_obj@result %>%
    dplyr::arrange(.data$p.adjust) %>%
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
    dplyr::filter(.data$name %in% ppi_genes)

  if (nrow(pie_data) == 0) {
    warning("No overlapping genes found between enrichment results and PPI network.")
  }

  return(pie_data)
}


#' compute MCC (maximum clique centrality) of PPI network.
#' @param graph A igraph object from PPI network.
#' @importFrom igraph max_cliques vcount V V<-
#'
#' @return An igraph object containing MCC.
#' @export

compute_MCC <- function(graph) {

  max_cliques_list <- igraph::max_cliques(graph)
  mcc <- setNames(rep(0, igraph::vcount(graph)), igraph::V(graph)$name)

  for (clq in max_cliques_list) {
    size <- length(clq)
    weight <- factorial(size - 1)
    mcc[igraph::V(graph)[clq]$name] <- mcc[igraph::V(graph)[clq]$name] + weight
  }

  igraph::V(graph)$MCC <- mcc

  return(graph)
}


#' Compute some useful metrics for PPI nodes
#' @param g An igraph object containing PPI information.
#' @param weight_attr Character; The attribute of weight in this PPI object.
#' @importFrom igraph E V edge_attr_names edge_attr
#' @importFrom igraph degree strength betweenness closeness
#' @importFrom igraph eigen_centrality page_rank coreness
#' @importFrom igraph transitivity eccentricity articulation_points
#' @return The input \code{igraph} object with additional vertex
#' @export

compute_nodeinfo <- function(g, weight_attr = "score") {
  stopifnot(inherits(g, "igraph"))

  ## PPI score
  w <- NULL
  if (!is.null(weight_attr) && weight_attr %in% igraph::edge_attr_names(g)) {
    w <- igraph::edge_attr(g, "score")
  }

  w_dist <- NULL
  if (!is.null(w)) {
    w_dist <- 1 / pmax(w, .Machine$double.eps)
  }

  ## compute node info
  message("Calculating degree / strength ...")
  V(g)$degree <- igraph::degree(g, mode = "all")
  V(g)$strength <- if (!is.null(w)) igraph::strength(g, mode = "all", weights = w) else NA_real_

  message("Calculating betweenness (unweighted) ...")
  V(g)$betweenness <- igraph::betweenness(g, directed = FALSE, normalized = TRUE)

  if (!is.null(w_dist)) {
    message("Calculating betweenness (weighted) ...")
    V(g)$betweenness_w <- igraph::betweenness(g, directed = FALSE,
                                              weights = w_dist,
                                              normalized = TRUE)
  } else {
    V(g)$betweenness_w <- NA_real_
  }

  message("Calculating closeness (unweighted) ...")
  V(g)$closeness <- igraph::closeness(g, normalized = TRUE)

  if (!is.null(w_dist)) {
    message("Calculating closeness (weighted) ...")
    V(g)$closeness_w <- igraph::closeness(g, normalized = TRUE,
                                          weights = w_dist)
  } else {
    V(g)$closeness_w <- NA_real_
  }

  message("Calculating eigenvector centrality ...")
  V(g)$eigen_centrality <- igraph::eigen_centrality(g, weights = w)$vector

  message("Calculating PageRank ...")
  V(g)$pagerank <- igraph::page_rank(g, weights = w)$vector

  message("Calculating coreness (k-core) ...")
  V(g)$coreness <- igraph::coreness(g, mode = "all")

  message("Calculating local clustering coefficient ...")
  V(g)$clustering_coef <- igraph::transitivity(g, type = "local", isolates = "zero")

  message("Calculating eccentricity ...")
  V(g)$eccentricity <- igraph::eccentricity(g)

  message("Detecting articulation points ...")
  art <- igraph::articulation_points(g)
  V(g)$is_articulation <- FALSE
  V(g)$is_articulation[art] <- TRUE

  ## compute MCC
  message("Calculating MCC ...")
  g <- compute_MCC(g)

  return(g)
}
