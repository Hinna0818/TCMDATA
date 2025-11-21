#' Compute MCL Clustering for PPI Network
#'
#' This function performs Markov Cluster Algorithm (MCL) on an igraph object
#' and assigns the cluster labels to the nodes.
#'
#' @param g An \code{igraph} object.
#' @param inflation Numeric. The inflation parameter controls the granularity of clusters.
#'   Common values for PPI: 2.0 - 3.0.
#'   Larger values (e.g., 3.0) produce smaller, tighter clusters.
#'   Smaller values (e.g., 1.5) produce larger, coarser clusters.
#'   Default is 2.5.
#' @param addLoops Logical; Self-loops with weight 1 are added to each vertex of g when TRUE; Default is TRUE (necessary).
#' @param ... Additional parameters in function `mcl()`.
#' @importFrom igraph as_adjacency_matrix V
#' @importFrom MCL mcl
#'
#' @return The input \code{igraph} object with a new vertex attribute \code{mcl_cluster}.
#' @export
run_MCL <- function(g, inflation = 2.5, addLoops = TRUE, ...) {
  
  stopifnot(inherits(g, "igraph"))
  message(sprintf("Running MCL with inflation = %.1f ...", inflation))
  
  adj_mat <- igraph::as_adjacency_matrix(g, sparse = TRUE)
  
  ## MCL running
  mcl_res <- MCL::mcl(x = adj_mat, addLoops = addLoops, inflation = inflation, ...)
  
  igraph::V(g)$mcl_cluster <- as.factor(mcl_res$Cluster)
  message(sprintf("Done! Identified %d modules (Iterations: %d).", 
                  mcl_res$K, mcl_res$n.iterations))
  
  return(g)
}


#' Compute Louvain Clustering for PPI Network
#'
#' This function performs Louvain clustering (multi-level modularity optimization)
#' on an igraph object. It is suitable for detecting larger, macro-scale functional
#' modules in PPI networks.
#'
#' @param g An \code{igraph} object.
#' @param resolution Numeric. The resolution parameter controls the granularity of clusters. Default is 1.0 (standard modularity).
#' @param weights Numeric vector or NULL. Edge weights to use for clustering.
#'   If NULL (default), the function attempts to use the 'weight' or 'score' edge attribute.
#'   Set to NA to perform unweighted clustering.
#'
#' @return The input \code{igraph} object with a new vertex attribute \code{louvain_cluster}.
#' @importFrom igraph cluster_louvain V edge_attr_names E membership
#' @export
run_louvain <- function(g, resolution = 1.0, weights = NULL) {
  
  stopifnot(inherits(g, "igraph"))
  
  # PPI networks often have 'score' or 'weight'. Using them improves accuracy.
  if (is.null(weights)) {
    edge_attrs <- igraph::edge_attr_names(g)
    if ("weight" %in% edge_attrs) {
      weights <- igraph::E(g)$weight
      message("Using edge attribute 'weight' for clustering.")
    } else if ("score" %in% edge_attrs) {
      weights <- igraph::E(g)$score
      message("Using edge attribute 'score' for clustering.")
    } else {
      weights <- NULL # Unweighted
      message("No weights found. Running unweighted clustering.")
    }
  } else if (identical(weights, NA)) {
    weights <- NULL # Explicitly unweighted
  }
  
  message(sprintf("Running Louvain (Resolution: %.1f)...", resolution))
  
  # Run Louvain Algorithm
  louvain_res <- igraph::cluster_louvain(g, weights = weights, resolution = resolution)
  
  # Assign Cluster Labels
  igraph::V(g)$louvain_cluster <- as.factor(igraph::membership(louvain_res))
  num_clusters <- length(unique(igraph::membership(louvain_res)))
  modularity_score <- max(louvain_res$modularity, na.rm = TRUE)
  
  message(sprintf("Done! Identified %d modules (Modularity: %.3f).", 
                  num_clusters, modularity_score))
  
  return(g)
}


