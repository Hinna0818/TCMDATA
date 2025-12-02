#' Molecular Complex Detection (MCODE) method for identifying PPI clusters.
#'
#' @param g An igraph object (PPI network).
#' @param vwp Numeric. Node Score Cutoff (Vertex Weight Percentage). Default 0.2.
#' @param degree_cutoff Numeric. Nodes with degree < cutoff will get a score of 0. Default 2.
#' @param k_core_threshold Numeric. Filters out clusters that do not contain a k-core of at least this level. Default 2.
#' @param haircut Logical. Whether to prune singly-connected nodes. Default TRUE.
#' @param fluff Logical. Whether to expand cluster by adding dense neighbors. Default FALSE.
#' @param fdt Numeric. Fluff Node Density Cutoff. Used if fluff=TRUE. Default 0.1.
#' @param loops Logical. Whether to include self-loops in scoring. Default FALSE.
#' @param max_depth Numeric. Maximum recursion depth for cluster finding (to prevent stack overflow on huge networks). Default 100.
#'
#' @importFrom igraph is_igraph simplify as_undirected induced_subgraph edge_density vcount
#' @return A list containing:
#' \item{complexes}{List of gene symbols for each module.}
#' \item{scores}{Numeric vector of node scores.}
#' \item{module_scores}{Score (Density * Size) for each module.}
#' @export
runMCODE <- function(g, 
                  vwp = 0.2, 
                  degree_cutoff = 2, 
                  k_core_threshold = 2, 
                  haircut = TRUE, 
                  fluff = FALSE, 
                  fdt = 0.1, 
                  loops = FALSE,
                  max_depth = 100) {
  if (!is_igraph(g)) {
    stop("Input 'g' must be an igraph object.")
  }
  
  if (!loops) {
    g <- simplify(as_undirected(g), remove.multiple = TRUE, remove.loops = TRUE)
  } else {
    g <- simplify(as_undirected(g), remove.multiple = TRUE, remove.loops = FALSE)
  }
  
  if (vcount(g) < 2) {
    warning("Graph is too small.")
    return(list(complexes = list(), scores = numeric(0), module_scores = numeric(0)))
  }
  
  #  Stage 1: Vertex Weighting (Scoring) 
  message("Stage 1: Vertex Weighting (k-core * density)...")
  node_scores <- .mcode_score_nodes(g, degree_cutoff)
  
  # Stage 2: Cluster Prediction (Seeding)
  message("Stage 2: Molecular Complex Prediction...")
  raw_clusters <- .mcode_find_clusters(g, node_scores, vwp, max_depth)
  
  # Stage 3: Post-processing & Filtering
  message("Stage 3: Post-processing (Haircut, Fluff, Filter)...")
  final_clusters <- .mcode_post_process(g, raw_clusters, node_scores, 
                                        haircut, fluff, fdt, k_core_threshold)
  
  # combine results
  message(paste("Identified", length(final_clusters), "complexes."))
  
  module_scores <- numeric(0)
  if (length(final_clusters) > 0) {
    module_scores <- sapply(final_clusters, function(nodes) {
      if (length(nodes) < 2) return(0)
      sub <- induced_subgraph(g, nodes)
      edge_density(sub) * length(nodes)
    })
    
    # order
    ord <- order(module_scores, decreasing = TRUE)
    final_clusters <- final_clusters[ord]
    module_scores <- module_scores[ord]
    names(final_clusters) <- paste0("Module_", seq_along(final_clusters))
  }
  
  res <- list(
    complexes = final_clusters,
    scores = node_scores,
    module_scores = module_scores)
  
  return(res)
}


#' Vertix weighting (step 1 for MCODE method)
#' @param g An igraph object (PPi networks).
#' @param degree_cutoff Numeric. Nodes with degree < cutoff will get a score of 0. Default 2.
#' @importFrom igraph V degree neighbors induced_subgraph coreness edge_density
#' @return A numeric vectors containing each vertex's score.
#' @keywords internal
.mcode_score_nodes <- function(g, degree_cutoff) {
  nodes <- V(g)$name
  scores <- numeric(length(nodes))
  names(scores) <- nodes
  degrees <- degree(g) 
  
  for (i in seq_along(nodes)) {
    v <- nodes[i]

    # Degree Cutoff
    if (degrees[v] < degree_cutoff) {
      scores[i] <- 0
      next
    }
    
    # get neighbors
    neis <- names(neighbors(g, v))
    sub_nodes <- c(v, neis)
    sub_g <- induced_subgraph(g, sub_nodes)
    
    # k-core calculation
    k_cores <- coreness(sub_g)
    k_max <- max(k_cores)
    
    if (k_max < 1) {
      scores[i] <- 0; next
    }
  
    core_nodes <- names(k_cores[k_cores == k_max])
    if (length(core_nodes) < 2) { 
      scores[i] <- 0; next 
    }
    k_core_subg <- induced_subgraph(sub_g, core_nodes)
    
    d <- edge_density(k_core_subg)
    scores[i] <- k_max * d
  }
  return(scores)
}


#' Molecular complex prediction (step 2 for MCODE method)
#' @param g An igraph object (PPI networks).
#' @param scores A numeric vector containing each vertex's score after running `.mcode_score_nodes()`.
#' @param vwp Numeric. Node Score Cutoff (Vertex Weight Percentage). Default 0.2.
#' @param max_depth Numeric. Maximum recursion depth for cluster finding (to prevent stack overflow on huge networks). Default 100.
#' @importFrom igraph neighbors
#' @return A list containing clusters.
#' @keywords internal
.mcode_find_clusters <- function(g, scores, vwp, max_depth) {
  sorted_nodes <- names(sort(scores, decreasing = TRUE))
  seen <- rep(FALSE, length(sorted_nodes))
  names(seen) <- names(scores)
  
  clusters <- list()
  
  for (seed in sorted_nodes) {
    if (seen[seed]) next
    if (scores[seed] == 0) next
    
    # BFS initialization
    complex_nodes <- c(seed)
    seen[seed] <- TRUE
    queue <- list(list(node = seed, depth = 0)) 
    
    threshold <- scores[seed] * (1 - vwp)
    
    while(length(queue) > 0) {
      curr_item <- queue[[1]]
      queue <- queue[-1]
      
      curr_node <- curr_item$node
      curr_depth <- curr_item$depth
      
      # Max Depth Protection
      if (curr_depth >= max_depth) next
      
      curr_neis <- names(neighbors(g, curr_node))
      
      for (v in curr_neis) {
        if (!seen[v]) {
          # Node Score Cutoff (VWP)
          if (scores[v] > threshold) {
            complex_nodes <- c(complex_nodes, v)
            seen[v] <- TRUE
            queue[[length(queue) + 1]] <- list(node = v, depth = curr_depth + 1)
          }
        }
      }
    }
    
    if (length(complex_nodes) >= 2) {
      clusters[[seed]] <- complex_nodes
    }
  }
  return(clusters)
}

#' Post-processing clusters (step 3 for MCODE method)
#' @param g An igraph object (PPI networks)
#' @param clusters A list containing cluster results after running `.mcode_find_clusters()`.
#' @param scores A numeric vector containing each vertex's score after running `.mcode_score_nodes()`.
#' @param haircut Logical. Whether to prune singly-connected nodes. Default TRUE.
#' @param fluff Logical. Whether to expand cluster by adding dense neighbors. Default FALSE.
#' @param fdt Numeric. Fluff Node Density Cutoff. Used if fluff=TRUE. Default 0.1.
#' @param k_core_threshold Numeric. Filters out clusters that do not contain a k-core of at least this level. Default 2.
#' @importFrom igraph neighbors induced_subgraph degree V delete_vertices coreness
#' @return A list containing final clusters after post-processing.
#' @keywords internal
.mcode_post_process <- function(g, clusters, scores, haircut, fluff, fdt, k_core_threshold) {
  final_clusters <- list()
  
  for (seed in names(clusters)) {
    c_nodes <- clusters[[seed]]
    
    # step1. Fluff
    if (fluff) {
      potential_nodes <- c()
      for (u in c_nodes) {
        u_neis <- names(neighbors(g, u))
        candidates <- setdiff(u_neis, c_nodes)
        # Fluff Node Density Cutoff
        valid_candidates <- candidates[scores[candidates] > fdt]
        potential_nodes <- c(potential_nodes, valid_candidates)
      }
      c_nodes <- unique(c(c_nodes, potential_nodes))
    }
    
    # step2. Haircut
    if (haircut && length(c_nodes) > 2) {
      sub_c <- induced_subgraph(g, c_nodes)
      while(TRUE) {
        deg <- degree(sub_c)
        leaf_nodes <- names(deg[deg < 2])
        if (length(leaf_nodes) == 0) break
        if (length(V(sub_c)) <= length(leaf_nodes)) {
          c_nodes <- character(0); break
        }
        sub_c <- delete_vertices(sub_c, leaf_nodes)
      }
      if (length(c_nodes) > 0) c_nodes <- V(sub_c)$name
    }
    
    # Filter (K-Core Threshold)
    if (length(c_nodes) >= 2) {
      sub_final <- induced_subgraph(g, c_nodes)
      if (vcount(sub_final) > 0) {
        max_coreness <- max(coreness(sub_final))
        if (max_coreness >= k_core_threshold) {
          final_clusters[[seed]] <- c_nodes
        }
      }
    }
  }
  
  return(final_clusters)
}

