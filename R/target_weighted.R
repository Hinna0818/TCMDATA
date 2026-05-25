# Target weighting and PPI proximity ranking

#' Prepare weighted disease targets from one or more evidence sources
#'
#' Combines disease-target evidence from scored sources such as Open Targets
#' Platform and GeneCards. By default, every disease-target source must provide
#' an explicit score column. This avoids replacing missing disease evidence with
#' artificial binary scores.
#'
#' @param ... Disease-target inputs. Each input must be a scored data frame by
#'   default. Data frames may contain columns such as \code{symbol},
#'   \code{gene_symbol}, \code{target}, \code{gene}, or \code{names}; scored
#'   inputs may contain \code{score},
#'   \code{overallAssociationScore}, \code{Relevance.Score}, or
#'   \code{disease_weight}.
#' @param source_weights Named numeric vector giving reliability weights for
#'   evidence sources. Missing source names default to the \code{user} weight,
#'   or 1 when \code{user} is absent.
#' @param gene_col Optional column name to use as the gene symbol column for
#'   all data-frame inputs.
#' @param score_col Optional column name to use as the raw score column for
#'   all data-frame inputs.
#' @param source_col Column containing source labels when present.
#' @param deg_score Logical. If \code{TRUE}, data frames with
#'   \code{log2FoldChange} and \code{padj} but no explicit score column are
#'   scored with an optional DEG heuristic. Defaults to \code{FALSE} because
#'   this is not part of the TCMNet network-proximity method.
#' @param deg_score_method DEG heuristic used when \code{deg_score = TRUE}.
#' @param require_score Logical. If \code{TRUE}, unscored disease-target inputs
#'   are rejected. Set to \code{FALSE} only for exploratory analyses where a
#'   source-level binary evidence fallback is intentional.
#' @param method How to combine multiple evidence rows per gene.
#'
#' @return A data frame with one row per gene and columns \code{symbol},
#'   \code{disease_weight}, \code{sources}, \code{n_sources},
#'   \code{n_evidence}, and \code{source_detail}.
#'
#' @importFrom stats quantile
#' @importFrom utils data
#'
#' @examples
#' otp_targets <- data.frame(gene_symbol = c("IL6", "VEGFA"),
#'                           score = c(0.91, 0.72))
#' gcds_targets <- data.frame(symbol = c("IL6", "TNF"),
#'                            score = c(67.2, 43.1))
#' prepare_disease_weights(OpenTargets = otp_targets, GeneCards = gcds_targets)
#'
#' @export
prepare_disease_weights <- function(...,
                                    source_weights = c(
                                      OpenTargets = 1.0,
                                      OTP = 1.0,
                                      GeneCards = 1.0,
                                      DOSE = 1.0,
                                      DisGeNET = 1.0,
                                      DEG = 1.0,
                                      user = 1.0
                                    ),
                                    gene_col = NULL,
                                    score_col = NULL,
                                    source_col = "source",
                                    deg_score = FALSE,
                                    deg_score_method = c("effect_significance", "wald_stat"),
                                    require_score = TRUE,
                                    method = c("noisy_or", "weighted_mean", "max")) {
  method <- match.arg(method)
  deg_score_method <- match.arg(deg_score_method)
  inputs <- list(...)
  input_names <- names(inputs)
  if (is.null(input_names)) input_names <- rep("", length(inputs))

  if (length(inputs) == 0L) {
    stop("At least one disease-target input is required.", call. = FALSE)
  }

  evidence <- list()
  for (i in seq_along(inputs)) {
    source_name <- input_names[[i]]
    if (is.null(source_name) || !nzchar(source_name)) source_name <- "user"
    evidence[[i]] <- .tw_standardise_disease_input(
      inputs[[i]],
      source_name = source_name,
      gene_col = gene_col,
      score_col = score_col,
      source_col = source_col,
      deg_score = deg_score,
      deg_score_method = deg_score_method,
      require_score = require_score
    )
  }

  evidence <- do.call(rbind, evidence)
  evidence$symbol <- .tw_clean_symbol(evidence$symbol)
  evidence <- evidence[!is.na(evidence$symbol) & nzchar(evidence$symbol), , drop = FALSE]
  if (nrow(evidence) == 0L) {
    stop("No valid disease-target symbols were found.", call. = FALSE)
  }

  evidence$source_weight <- .tw_source_weight(evidence$source, source_weights)
  evidence$evidence_weight <- pmin(1, pmax(0, evidence$score_norm * evidence$source_weight))

  split_ev <- split(evidence, evidence$symbol)
  out <- lapply(split_ev, function(x) {
    sources <- sort(unique(x$source))
    disease_weight <- switch(
      method,
      noisy_or = 1 - prod(1 - x$evidence_weight),
      weighted_mean = {
        denom <- sum(x$source_weight)
        if (denom == 0) 0 else sum(x$score_norm * x$source_weight) / denom
      },
      max = max(x$evidence_weight)
    )

    data.frame(
      symbol = x$symbol[[1]],
      disease_weight = pmin(1, pmax(0, disease_weight)),
      sources = paste(sources, collapse = ";"),
      n_sources = length(sources),
      n_evidence = nrow(x),
      source_detail = paste(
        sprintf("%s:%.3f", x$source, x$evidence_weight),
        collapse = ";"
      ),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out <- out[order(-out$disease_weight, out$symbol), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Prepare target-side weights from herb-compound-target records
#'
#' @param herb_df A data frame containing herb, molecule, and target columns,
#'   typically returned by \code{\link{search_herb}}.
#' @param target_col Column containing target gene symbols.
#' @param molecule_col Column containing compound names.
#' @param herb_col Column containing herb names.
#' @param method Weighting method: \code{"compound_count"} counts unique
#'   compounds per target, \code{"record_count"} counts rows, and
#'   \code{"binary"} gives all targets weight 1. Default is \code{"binary"}.
#'
#' @return A data frame with \code{target}, \code{herb_target_weight},
#'   \code{n_compounds}, \code{n_herbs}, and \code{n_records}.
#'
#' @importFrom stats na.omit
#'
#' @export
prepare_herb_target_weights <- function(herb_df,
                                        target_col = "target",
                                        molecule_col = "molecule",
                                        herb_col = "herb",
                                        method = c("binary", "record_count", "compound_count")) {
  method <- match.arg(method)

  if (is.character(herb_df) && !is.data.frame(herb_df)) {
    targets <- .tw_clean_symbol(herb_df)
    targets <- unique(targets[!is.na(targets) & nzchar(targets)])
    return(data.frame(
      target = targets,
      herb_target_weight = 1,
      n_compounds = NA_integer_,
      n_herbs = NA_integer_,
      n_records = 1L,
      stringsAsFactors = FALSE
    ))
  }

  if (!is.data.frame(herb_df)) {
    stop("herb_df must be a data frame or character vector.", call. = FALSE)
  }
  if (!target_col %in% names(herb_df)) {
    stop("Target column not found: ", target_col, call. = FALSE)
  }

  df <- herb_df
  df[[target_col]] <- .tw_clean_symbol(df[[target_col]])
  df <- df[!is.na(df[[target_col]]) & nzchar(df[[target_col]]), , drop = FALSE]
  if (nrow(df) == 0L) {
    stop("No valid herb target symbols were found.", call. = FALSE)
  }

  split_df <- split(df, df[[target_col]])
  out <- lapply(split_df, function(x) {
    n_compounds <- if (molecule_col %in% names(x)) length(unique(na.omit(x[[molecule_col]]))) else NA_integer_
    n_herbs <- if (herb_col %in% names(x)) length(unique(na.omit(x[[herb_col]]))) else NA_integer_
    n_records <- nrow(x)
    herb_target_weight <- switch(
      method,
      compound_count = if (is.na(n_compounds)) n_records else n_compounds,
      record_count = n_records,
      binary = 1
    )

    data.frame(
      target = x[[target_col]][[1]],
      herb_target_weight = herb_target_weight,
      n_compounds = n_compounds,
      n_herbs = n_herbs,
      n_records = n_records,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  out <- out[order(-out$herb_target_weight, out$target), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Rank TCM targets by weighted PPI proximity to a disease module
#'
#' @param herb_targets A herb target data frame, the output of
#'   \code{\link{prepare_herb_target_weights}}, or a character vector of target
#'   genes.
#' @param disease_targets A scored disease target data frame or the output of
#'   \code{\link{prepare_disease_weights}}.
#' @param ppi An undirected \code{igraph} PPI background network. If
#'   \code{NULL}, the function uses the package's built-in
#'   \code{string_human_ppi} background, or a path supplied through
#'   \code{ppi_path} / \code{getOption("TCMDATA.string_ppi")}.
#' @param ppi_path Optional path to a saved igraph RDS used when \code{ppi} is
#'   \code{NULL}.
#' @param target_col Target column in \code{herb_targets}.
#' @param disease_gene_col Gene column in \code{disease_targets}.
#' @param disease_weight_col Disease weight column in \code{disease_targets}.
#' @param weight_attr Edge confidence attribute. Used to derive distances when
#'   \code{distance_attr} is absent.
#' @param distance_attr Edge distance attribute. Defaults to \code{"distance"}.
#' @param n_perm Number of degree-matched random draws per candidate target.
#' @param degree_bins Number of degree quantile bins for background matching.
#' @param seed Random seed.
#' @param eps Deprecated numeric kept for backward compatibility. The current
#'   proximity kernel is bounded as \code{1 / (1 + distance)}, so direct overlap
#'   no longer produces unbounded scores.
#' @param warn_coverage Numeric coverage threshold below which a warning is
#'   issued for herb or disease targets present in the PPI network.
#'
#' @return A \code{tcm_ppi_rank} object containing \code{result},
#'   \code{disease_module}, \code{diagnostics}, \code{ppi_summary}, and
#'   \code{params}.
#'
#' @importFrom igraph V distances degree vcount ecount components edge_attr_names edge_attr
#' @importFrom stats p.adjust sd setNames
#'
#' @export
rank_tcm_targets_by_ppi <- function(herb_targets,
                                    disease_targets,
                                    ppi = NULL,
                                    ppi_path = NULL,
                                    target_col = "target",
                                    disease_gene_col = "symbol",
                                    disease_weight_col = "disease_weight",
                                    weight_attr = "score",
                                    distance_attr = "distance",
                                    n_perm = 1000L,
                                    degree_bins = 20L,
                                    seed = 42L,
                                    eps = 1e-8,
                                    warn_coverage = 0.5) {
  if (is.null(ppi)) {
    ppi <- .tw_load_default_ppi(ppi_path)
  }
  if (!inherits(ppi, "igraph")) {
    stop("ppi must be an igraph object.", call. = FALSE)
  }
  n_perm <- as.integer(n_perm)
  if (n_perm < 1L) stop("n_perm must be at least 1.", call. = FALSE)

  herb_w <- .tw_as_herb_weights(herb_targets, target_col = target_col)
  disease_w <- .tw_as_disease_weights(
    disease_targets,
    disease_gene_col = disease_gene_col,
    disease_weight_col = disease_weight_col
  )

  vertex_names <- V(ppi)$name
  if (is.null(vertex_names)) {
    stop("ppi must have vertex names containing gene symbols.", call. = FALSE)
  }

  herb_w$target_query <- herb_w$target
  disease_w$symbol_query <- disease_w$symbol

  herb_w$target <- .tw_match_graph_symbols(herb_w$target, vertex_names)
  disease_w$symbol <- .tw_match_graph_symbols(disease_w$symbol, vertex_names)

  dropped_herb <- herb_w$target_query[is.na(herb_w$target)]
  dropped_disease <- disease_w$symbol_query[is.na(disease_w$symbol)]
  herb_w <- herb_w[!is.na(herb_w$target), , drop = FALSE]
  disease_w <- disease_w[!is.na(disease_w$symbol), , drop = FALSE]

  if (nrow(herb_w) == 0L) {
    stop("None of the herb targets are present in the PPI network.", call. = FALSE)
  }
  if (nrow(disease_w) == 0L) {
    stop("None of the disease targets are present in the PPI network.", call. = FALSE)
  }

  diagnostics <- list(
    n_herb_targets = length(unique(herb_w$target_query)),
    n_disease_targets = length(unique(disease_w$symbol_query)),
    n_herb_in_ppi = length(unique(herb_w$target)),
    n_disease_in_ppi = length(unique(disease_w$symbol)),
    dropped_herb_targets = unique(dropped_herb),
    dropped_disease_targets = unique(dropped_disease)
  )
  diagnostics$herb_ppi_coverage <- diagnostics$n_herb_in_ppi / diagnostics$n_herb_targets
  diagnostics$disease_ppi_coverage <- diagnostics$n_disease_in_ppi / diagnostics$n_disease_targets
  .tw_warn_coverage(diagnostics, warn_coverage = warn_coverage)
  if (n_perm < 100L) {
    warning("n_perm < 100; empirical P-values will be coarse.", call. = FALSE)
  }

  disease_w <- .tw_collapse_disease_weights(disease_w)
  herb_w <- .tw_collapse_herb_weights(herb_w)

  distance_weights <- .tw_edge_distances(ppi, weight_attr, distance_attr)
  disease_nodes <- unique(disease_w$symbol)
  disease_weights <- setNames(disease_w$disease_weight, disease_w$symbol)

  message("Computing weighted distances from PPI background to disease module ...")
  dist_mat <- distances(
    ppi,
    v = V(ppi),
    to = disease_nodes,
    weights = distance_weights,
    mode = "all"
  )
  all_proximity <- .tw_weighted_proximity_matrix(dist_mat, disease_weights[colnames(dist_mat)])
  all_proximity_no_self <- .tw_weighted_proximity_matrix_no_self(
    dist_mat,
    disease_weights[colnames(dist_mat)]
  )
  names(all_proximity) <- vertex_names
  names(all_proximity_no_self) <- vertex_names

  candidates <- herb_w$target
  observed <- all_proximity[candidates]
  observed_no_self <- all_proximity_no_self[candidates]
  candidate_dist <- dist_mat[candidates, disease_nodes, drop = FALSE]
  nearest_idx <- max.col(-candidate_dist, ties.method = "first")
  nearest_distance <- candidate_dist[cbind(seq_along(candidates), nearest_idx)]
  nearest_disease_target <- colnames(candidate_dist)[nearest_idx]
  nearest_disease_target[!is.finite(nearest_distance)] <- NA_character_
  nearest_distance[!is.finite(nearest_distance)] <- NA_real_

  degree <- degree(ppi, mode = "all")
  bins <- .tw_degree_bins(degree, degree_bins)
  names(bins) <- names(degree)

  set.seed(seed)
  random_stats <- lapply(candidates, function(target) {
    pool <- names(bins)[bins == bins[[target]] & is.finite(all_proximity[names(bins)])]
    pool <- setdiff(pool, target)
    matched_pool_size <- length(pool)
    if (length(pool) < 10L) {
      pool <- names(all_proximity)[is.finite(all_proximity)]
      pool <- setdiff(pool, target)
    }
    draws <- sample(pool, n_perm, replace = length(pool) < n_perm)
    random <- all_proximity[draws]
    random <- random[is.finite(random)]
    if (length(random) == 0L) random <- 0
    random_no_self <- all_proximity_no_self[draws]
    random_no_self <- random_no_self[is.finite(random_no_self)]
    if (length(random_no_self) == 0L) random_no_self <- 0
    c(
      random_mean = mean(random),
      random_sd = sd(random),
      p_empirical = (sum(random >= observed[[target]]) + 1) / (length(random) + 1),
      random_mean_no_self = mean(random_no_self),
      random_sd_no_self = sd(random_no_self),
      p_empirical_no_self = (sum(random_no_self >= observed_no_self[[target]]) + 1) / (length(random_no_self) + 1),
      random_pool_size = length(pool),
      degree_matched_pool_size = matched_pool_size
    )
  })
  random_stats <- as.data.frame(do.call(rbind, random_stats), stringsAsFactors = FALSE)

  random_sd_safe <- ifelse(is.na(random_stats$random_sd) | random_stats$random_sd == 0, 1, random_stats$random_sd)
  z_proximity <- (observed - random_stats$random_mean) / random_sd_safe
  random_sd_no_self_safe <- ifelse(
    is.na(random_stats$random_sd_no_self) | random_stats$random_sd_no_self == 0,
    1,
    random_stats$random_sd_no_self
  )
  z_proximity_no_self <- (observed_no_self - random_stats$random_mean_no_self) / random_sd_no_self_safe

  direct_overlap <- candidates %in% disease_nodes
  disease_overlap_weight <- disease_weights[candidates]
  disease_overlap_weight[is.na(disease_overlap_weight)] <- 0

  result <- data.frame(
    target = candidates,
    direct_overlap = direct_overlap,
    degree = as.numeric(degree[candidates]),
    herb_target_weight = herb_w$herb_target_weight,
    disease_weight_if_overlap = as.numeric(disease_overlap_weight),
    nearest_disease_target = nearest_disease_target,
    nearest_distance = as.numeric(nearest_distance),
    proximity_raw = as.numeric(observed),
    proximity_no_self = as.numeric(observed_no_self),
    random_mean = random_stats$random_mean,
    random_sd = random_stats$random_sd,
    z_proximity = as.numeric(z_proximity),
    p_empirical = random_stats$p_empirical,
    p_adjust = p.adjust(random_stats$p_empirical, method = "BH"),
    random_mean_no_self = random_stats$random_mean_no_self,
    random_sd_no_self = random_stats$random_sd_no_self,
    z_proximity_no_self = as.numeric(z_proximity_no_self),
    p_empirical_no_self = random_stats$p_empirical_no_self,
    p_adjust_no_self = p.adjust(random_stats$p_empirical_no_self, method = "BH"),
    random_pool_size = as.integer(random_stats$random_pool_size),
    degree_matched_pool_size = as.integer(random_stats$degree_matched_pool_size),
    stringsAsFactors = FALSE
  )
  for (extra_col in intersect(c("n_compounds", "n_herbs", "n_records"), names(herb_w))) {
    result[[extra_col]] <- herb_w[[extra_col]]
  }

  result$proximity_no_self_norm <- .tw_norm01(result$z_proximity_no_self)
  result$Score_final <- result$z_proximity_no_self
  result$Rank_final <- rank(-result$Score_final, ties.method = "min")
  result <- result[order(result$Rank_final, result$target), , drop = FALSE]
  rownames(result) <- NULL

  ppi_summary <- list(
    nodes = vcount(ppi),
    edges = ecount(ppi),
    components = components(ppi)$no,
    weight_attr = weight_attr,
    distance_attr = if (distance_attr %in% edge_attr_names(ppi)) distance_attr else NA_character_
  )

  out <- list(
    result = result,
    diagnostics = diagnostics,
    ppi_summary = ppi_summary,
    params = list(
      n_perm = n_perm,
      degree_bins = degree_bins,
      seed = seed,
      eps = eps,
      warn_coverage = warn_coverage
    )
  )
  out$disease_module <- disease_w[order(-disease_w$disease_weight, disease_w$symbol), , drop = FALSE]
  class(out) <- c("tcm_ppi_rank", "list")
  out
}

#' Select targets from a TCM PPI ranking result
#'
#' @param x A \code{tcm_ppi_rank} object returned by
#'   \code{\link{rank_tcm_targets_by_ppi}}, or a ranking data frame.
#' @param top_n Maximum number of targets to keep. If \code{NULL}, no top-n
#'   limit is applied.
#' @param min_score Minimum \code{Score_final}.
#' @param min_z Minimum \code{z_proximity}.
#' @param max_p Maximum empirical P-value.
#' @param max_p_adjust Maximum BH-adjusted empirical P-value.
#' @param direct_overlap Optional logical filter for direct disease-target
#'   overlap.
#' @param return Character. Return the selected table or only target symbols.
#'
#' @return A data frame or character vector.
#'
#' @importFrom utils head
#'
#' @export
select_tcm_targets <- function(x,
                               top_n = 50L,
                               min_score = NULL,
                               min_z = NULL,
                               max_p = NULL,
                               max_p_adjust = NULL,
                               direct_overlap = NULL,
                               return = c("table", "targets")) {
  return <- match.arg(return)
  df <- if (inherits(x, "tcm_ppi_rank")) x$result else x
  if (!is.data.frame(df)) {
    stop("x must be a tcm_ppi_rank object or a ranking data frame.", call. = FALSE)
  }
  required <- c("target", "Score_final", "Rank_final")
  if (!all(required %in% names(df))) {
    stop("Ranking table must contain target, Score_final, and Rank_final.", call. = FALSE)
  }

  df <- df[order(df$Rank_final, -df$Score_final, df$target), , drop = FALSE]
  keep <- rep(TRUE, nrow(df))
  if (!is.null(min_score)) keep <- keep & df$Score_final >= min_score
  if (!is.null(min_z) && "z_proximity" %in% names(df)) keep <- keep & df$z_proximity >= min_z
  if (!is.null(max_p) && "p_empirical" %in% names(df)) keep <- keep & df$p_empirical <= max_p
  if (!is.null(max_p_adjust) && "p_adjust" %in% names(df)) keep <- keep & df$p_adjust <= max_p_adjust
  if (!is.null(direct_overlap) && "direct_overlap" %in% names(df)) {
    keep <- keep & df$direct_overlap %in% direct_overlap
  }
  df <- df[keep, , drop = FALSE]
  if (!is.null(top_n)) {
    top_n <- max(0L, as.integer(top_n))
    df <- head(df, top_n)
  }
  rownames(df) <- NULL

  if (return == "targets") df$target else df
}

#' @export
print.tcm_ppi_rank <- function(x, ...) {
  cat("TCM PPI target ranking\n")
  cat("Targets ranked:", nrow(x$result), "\n")
  cat("PPI background:", x$ppi_summary$nodes, "nodes,", x$ppi_summary$edges, "edges\n")
  cat("Herb targets in PPI:", x$diagnostics$n_herb_in_ppi, "/", x$diagnostics$n_herb_targets, "\n")
  cat("Disease targets in PPI:", x$diagnostics$n_disease_in_ppi, "/", x$diagnostics$n_disease_targets, "\n")
  if (!is.null(x$disease_module)) {
    cat("Disease module:", nrow(x$disease_module), "weighted targets\n")
  }
  print(head(x$result, 10), row.names = FALSE)
  invisible(x)
}

#' Evaluate target-set proximity to a weighted disease module on a PPI network
#'
#' This set-level analysis is closer to the original TCMNet question than
#' single-target ranking: does a herb, formula, compound group, or any target
#' set lie closer to the disease module than expected from degree-matched random
#' target sets?
#'
#' @param target_sets A character vector, a target data frame, or a named list
#'   of character vectors/data frames. When a data frame contains multiple
#'   groups, use \code{set_col} to split it into target sets.
#' @param disease_targets A scored disease target data frame or the output of
#'   \code{\link{prepare_disease_weights}}.
#' @param ppi An undirected \code{igraph} PPI background network. If
#'   \code{NULL}, the package's built-in STRING background is used.
#' @param ppi_path Optional path to a saved igraph RDS used when \code{ppi} is
#'   \code{NULL}.
#' @param set_col Optional grouping column in \code{target_sets}.
#' @param target_col Target column in data-frame target-set inputs.
#' @param disease_gene_col Gene column in \code{disease_targets}.
#' @param disease_weight_col Disease weight column in \code{disease_targets}.
#' @param weight_attr Edge confidence attribute. Used to derive distances when
#'   \code{distance_attr} is absent.
#' @param distance_attr Edge distance attribute. Defaults to \code{"distance"}.
#' @param n_perm Number of degree-matched random target sets per observed set.
#' @param degree_bins Number of degree quantile bins for background matching.
#' @param seed Random seed.
#' @param warn_coverage Numeric coverage threshold below which a warning is
#'   issued for target-set or disease targets present in the PPI network.
#'
#' @return A \code{tcm_ppi_set_eval} object containing \code{result},
#'   \code{set_details}, \code{disease_module}, \code{diagnostics},
#'   \code{ppi_summary}, and \code{params}.
#'
#' @importFrom igraph V distances degree vcount ecount components edge_attr_names edge_attr
#' @importFrom stats p.adjust sd setNames
#'
#' @export
evaluate_tcm_target_sets_by_ppi <- function(target_sets,
                                            disease_targets,
                                            ppi = NULL,
                                            ppi_path = NULL,
                                            set_col = NULL,
                                            target_col = "target",
                                            disease_gene_col = "symbol",
                                            disease_weight_col = "disease_weight",
                                            weight_attr = "score",
                                            distance_attr = "distance",
                                            n_perm = 1000L,
                                            degree_bins = 20L,
                                            seed = 42L,
                                            warn_coverage = 0.5) {
  if (is.null(ppi)) {
    ppi <- .tw_load_default_ppi(ppi_path)
  }
  if (!inherits(ppi, "igraph")) {
    stop("ppi must be an igraph object.", call. = FALSE)
  }
  n_perm <- as.integer(n_perm)
  if (n_perm < 1L) stop("n_perm must be at least 1.", call. = FALSE)
  if (n_perm < 100L) {
    warning("n_perm < 100; empirical P-values will be coarse.", call. = FALSE)
  }

  set_list <- .tw_as_target_set_list(target_sets, set_col = set_col, target_col = target_col)
  disease_w <- .tw_as_disease_weights(
    disease_targets,
    disease_gene_col = disease_gene_col,
    disease_weight_col = disease_weight_col
  )

  vertex_names <- V(ppi)$name
  if (is.null(vertex_names)) {
    stop("ppi must have vertex names containing gene symbols.", call. = FALSE)
  }

  disease_w$symbol_query <- disease_w$symbol
  disease_w$symbol <- .tw_match_graph_symbols(disease_w$symbol, vertex_names)
  dropped_disease <- disease_w$symbol_query[is.na(disease_w$symbol)]
  disease_w <- disease_w[!is.na(disease_w$symbol), , drop = FALSE]
  if (nrow(disease_w) == 0L) {
    stop("None of the disease targets are present in the PPI network.", call. = FALSE)
  }
  disease_w <- .tw_collapse_disease_weights(disease_w)
  disease_nodes <- unique(disease_w$symbol)
  disease_weights <- setNames(disease_w$disease_weight, disease_w$symbol)

  mapped_sets <- lapply(set_list, function(x) {
    x$target_query <- x$target
    x$target <- .tw_match_graph_symbols(x$target, vertex_names)
    dropped <- x$target_query[is.na(x$target)]
    x <- x[!is.na(x$target), , drop = FALSE]
    if (nrow(x) > 0L) x <- .tw_collapse_herb_weights(x)
    attr(x, "dropped_targets") <- unique(dropped)
    x
  })

  if (any(vapply(mapped_sets, nrow, integer(1)) == 0L)) {
    empty <- names(mapped_sets)[vapply(mapped_sets, nrow, integer(1)) == 0L]
    stop(
      "The following target sets have no targets in the PPI network: ",
      paste(empty, collapse = ", "),
      call. = FALSE
    )
  }

  diagnostics <- list(
    n_sets = length(mapped_sets),
    n_disease_targets = length(unique(disease_w$symbol_query)),
    n_disease_in_ppi = length(unique(disease_w$symbol)),
    dropped_disease_targets = unique(dropped_disease),
    set_coverage = data.frame(
      set = names(mapped_sets),
      n_targets = vapply(set_list, function(x) length(unique(x$target)), integer(1)),
      n_targets_in_ppi = vapply(mapped_sets, function(x) length(unique(x$target)), integer(1)),
      stringsAsFactors = FALSE
    )
  )
  diagnostics$disease_ppi_coverage <- diagnostics$n_disease_in_ppi / diagnostics$n_disease_targets
  diagnostics$set_coverage$target_ppi_coverage <-
    diagnostics$set_coverage$n_targets_in_ppi / diagnostics$set_coverage$n_targets
  diagnostics$dropped_set_targets <- lapply(mapped_sets, attr, which = "dropped_targets")
  .tw_warn_set_coverage(diagnostics, warn_coverage = warn_coverage)

  distance_weights <- .tw_edge_distances(ppi, weight_attr, distance_attr)
  message("Computing weighted distances from PPI background to disease module ...")
  dist_mat <- distances(
    ppi,
    v = V(ppi),
    to = disease_nodes,
    weights = distance_weights,
    mode = "all"
  )
  all_proximity <- .tw_weighted_proximity_matrix(dist_mat, disease_weights[colnames(dist_mat)])
  names(all_proximity) <- vertex_names

  degree <- degree(ppi, mode = "all")
  bins <- .tw_degree_bins(degree, degree_bins)
  names(bins) <- names(degree)

  set.seed(seed)
  result <- lapply(names(mapped_sets), function(set_name) {
    set_w <- mapped_sets[[set_name]]
    set_nodes <- set_w$target
    set_weights <- set_w$herb_target_weight
    set_weights[!is.finite(set_weights) | set_weights < 0] <- 0
    if (sum(set_weights) == 0) set_weights[] <- 1

    observed <- weighted.mean(all_proximity[set_nodes], set_weights, na.rm = TRUE)
    random <- replicate(
      n_perm,
      {
        random_nodes <- .tw_sample_degree_matched_set(
          set_nodes = set_nodes,
          bins = bins,
          all_proximity = all_proximity
        )
        weighted.mean(all_proximity[random_nodes], set_weights, na.rm = TRUE)
      }
    )
    random <- random[is.finite(random)]
    if (length(random) == 0L) random <- 0
    random_mean <- mean(random)
    random_sd <- sd(random)
    random_sd_safe <- ifelse(is.na(random_sd) || random_sd == 0, 1, random_sd)
    z_proximity <- (observed - random_mean) / random_sd_safe
    p_empirical <- (sum(random >= observed) + 1) / (length(random) + 1)

    overlap_nodes <- intersect(set_nodes, disease_nodes)
    union_nodes <- union(set_nodes, disease_nodes)
    direct <- .tw_directlink_stats(ppi, set_nodes = set_nodes, disease_nodes = disease_nodes)
    target_coverage <- diagnostics$set_coverage$target_ppi_coverage[
      match(set_name, diagnostics$set_coverage$set)
    ]

    data.frame(
      set = set_name,
      n_targets = diagnostics$set_coverage$n_targets[match(set_name, diagnostics$set_coverage$set)],
      n_targets_in_ppi = length(unique(set_nodes)),
      target_ppi_coverage = target_coverage,
      n_disease_targets = diagnostics$n_disease_targets,
      n_disease_in_ppi = length(disease_nodes),
      overlap_n = length(overlap_nodes),
      overlap_targets = paste(sort(overlap_nodes), collapse = ";"),
      overlap_target_fraction = length(overlap_nodes) / length(unique(set_nodes)),
      overlap_disease_fraction = length(overlap_nodes) / length(disease_nodes),
      jaccard = length(overlap_nodes) / length(union_nodes),
      directlink_edges = direct$directlink_edges,
      directlink_target_n = direct$directlink_target_n,
      directlink_disease_n = direct$directlink_disease_n,
      first_order_target_n = length(unique(c(overlap_nodes, direct$directlink_targets))),
      directlink_target_fraction = direct$directlink_target_n / length(unique(set_nodes)),
      first_order_target_coverage = length(unique(c(overlap_nodes, direct$directlink_targets))) / length(unique(set_nodes)),
      weighted_proximity = as.numeric(observed),
      random_mean = random_mean,
      random_sd = random_sd,
      z_proximity = as.numeric(z_proximity),
      p_empirical = p_empirical,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, result)
  result$p_adjust <- p.adjust(result$p_empirical, method = "BH")
  result$Score_set <- result$z_proximity
  result$Rank_set <- rank(-result$Score_set, ties.method = "min")
  result <- result[order(result$Rank_set, result$set), , drop = FALSE]
  rownames(result) <- NULL

  ppi_summary <- list(
    nodes = vcount(ppi),
    edges = ecount(ppi),
    components = components(ppi)$no,
    weight_attr = weight_attr,
    distance_attr = if (distance_attr %in% edge_attr_names(ppi)) distance_attr else NA_character_
  )

  out <- list(
    result = result,
    set_details = mapped_sets,
    disease_module = disease_w[order(-disease_w$disease_weight, disease_w$symbol), , drop = FALSE],
    diagnostics = diagnostics,
    ppi_summary = ppi_summary,
    params = list(
      n_perm = n_perm,
      degree_bins = degree_bins,
      seed = seed,
      warn_coverage = warn_coverage
    )
  )
  class(out) <- c("tcm_ppi_set_eval", "list")
  out
}

#' @export
print.tcm_ppi_set_eval <- function(x, ...) {
  cat("TCM PPI target-set evaluation\n")
  cat("Target sets:", nrow(x$result), "\n")
  cat("PPI background:", x$ppi_summary$nodes, "nodes,", x$ppi_summary$edges, "edges\n")
  cat("Disease targets in PPI:", x$diagnostics$n_disease_in_ppi, "/", x$diagnostics$n_disease_targets, "\n")
  print(head(x$result, 10), row.names = FALSE)
  invisible(x)
}

.tw_standardise_disease_input <- function(x,
                                          source_name,
                                          gene_col = NULL,
                                          score_col = NULL,
                                          source_col = "source",
                                          deg_score = FALSE,
                                          deg_score_method = "effect_significance",
                                          require_score = TRUE) {
  if (is.null(x)) {
    return(data.frame(symbol = character(), source = character(), score_raw = numeric(),
                      score_norm = numeric(), has_score = logical(), stringsAsFactors = FALSE))
  }

  if (is.character(x) && !is.data.frame(x)) {
    if (isTRUE(require_score)) {
      stop(
        "Disease-target source '", source_name, "' has no target-level score. ",
        "Pass a scored data frame, for example Open Targets or GeneCards with a score column, ",
        "or set require_score = FALSE for exploratory binary evidence.",
        call. = FALSE
      )
    }
    return(data.frame(
      symbol = x,
      source = source_name,
      score_raw = 1,
      score_norm = 1,
      has_score = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  if (!is.data.frame(x)) {
    stop("Disease-target inputs must be data frames or character vectors.", call. = FALSE)
  }

  gcol <- gene_col
  if (is.null(gcol)) {
    gcol <- .tw_first_existing(names(x), c("symbol", "gene_symbol", "approvedSymbol", "target", "gene", "names", "Gene.Symbol", "Gene"))
  }
  if (is.na(gcol)) {
    stop("Cannot find a gene symbol column in disease-target input: ", source_name, call. = FALSE)
  }

  src <- rep(source_name, nrow(x))
  if (!is.null(source_col) && source_col %in% names(x)) {
    src <- as.character(x[[source_col]])
    src[is.na(src) | !nzchar(src)] <- source_name
  }

  score <- NULL
  has_score <- TRUE
  scol <- score_col
  if (is.null(scol)) {
    scol <- .tw_first_existing(names(x), c("score", "overallAssociationScore", "overall_association_score", "Relevance.Score", "relevance_score", "disease_weight"))
  }
  if (!is.na(scol)) {
    score <- suppressWarnings(as.numeric(x[[scol]]))
  } else if (isTRUE(deg_score) && all(c("log2FoldChange", "padj") %in% names(x))) {
    if (identical(deg_score_method, "wald_stat") && "stat" %in% names(x)) {
      score <- abs(suppressWarnings(as.numeric(x$stat)))
    } else {
      score <- abs(suppressWarnings(as.numeric(x$log2FoldChange))) *
        -log10(pmax(suppressWarnings(as.numeric(x$padj)), 1e-300))
    }
    src <- ifelse(src == source_name & source_name == "user", "DEG", src)
  } else {
    if (isTRUE(require_score)) {
      stop(
        "Disease-target source '", source_name, "' has no score column. ",
        "Expected one of: score, overallAssociationScore, overall_association_score, ",
        "Relevance.Score, relevance_score, disease_weight. ",
        "Set require_score = FALSE only if a binary fallback is intentional.",
        call. = FALSE
      )
    }
    score <- rep(1, nrow(x))
    has_score <- FALSE
  }

  score[!is.finite(score)] <- NA_real_
  score_norm <- .tw_score_norm(score, has_score = has_score)

  data.frame(
    symbol = x[[gcol]],
    source = src,
    score_raw = score,
    score_norm = score_norm,
    has_score = has_score,
    stringsAsFactors = FALSE
  )
}

.tw_load_default_ppi <- function(ppi_path = NULL) {
  candidates <- unlist(c(
    ppi_path,
    getOption("TCMDATA.string_ppi")
  ), use.names = FALSE)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  existing <- if (length(candidates) > 0L) candidates[file.exists(candidates)] else character()
  if (length(existing) > 0L) {
    return(readRDS(existing[[1]]))
  }

  ns <- asNamespace("TCMDATA")
  if (exists("string_human_ppi", envir = ns, inherits = FALSE)) {
    return(get("string_human_ppi", envir = ns, inherits = FALSE))
  }

  data_env <- new.env(parent = emptyenv())
  loaded <- try(data("string_human_ppi", package = "TCMDATA", envir = data_env), silent = TRUE)
  if (!inherits(loaded, "try-error") && exists("string_human_ppi", envir = data_env, inherits = FALSE)) {
    return(get("string_human_ppi", envir = data_env, inherits = FALSE))
  }

  dev_path <- file.path(getwd(), "dev/string/string_human_ppi_v12.0_score700_lcc_igraph.rds")
  if (file.exists(dev_path)) {
    return(readRDS(dev_path))
  }

  stop(
    "No PPI background was supplied and no built-in/default STRING background could be loaded. ",
    "Pass an igraph object with 'ppi', set options(TCMDATA.string_ppi = ...), ",
    "or provide 'ppi_path'.",
    call. = FALSE
  )
}

.tw_clean_symbol <- function(x) {
  x <- trimws(as.character(x))
  toupper(x)
}

.tw_first_existing <- function(existing, candidates) {
  hit <- candidates[candidates %in% existing]
  if (length(hit) == 0L) NA_character_ else hit[[1]]
}

.tw_score_norm <- function(score, has_score) {
  if (!has_score) return(rep(1, length(score)))
  finite <- is.finite(score)
  out <- rep(NA_real_, length(score))
  if (!any(finite)) {
    out[] <- 1
    return(out)
  }
  s <- score[finite]
  if (all(s >= 0 & s <= 1, na.rm = TRUE)) {
    out[finite] <- s
  } else {
    rng <- range(s, na.rm = TRUE)
    if (diff(rng) == 0) {
      out[finite] <- 1
    } else {
      out[finite] <- (s - rng[[1]]) / diff(rng)
    }
  }
  out[is.na(out)] <- 0
  pmin(1, pmax(0, out))
}

.tw_source_weight <- function(source, source_weights) {
  if (is.null(source_weights) || length(source_weights) == 0L) {
    return(rep(1, length(source)))
  }
  nms <- names(source_weights)
  fallback <- if ("user" %in% nms) source_weights[["user"]] else 1
  idx <- match(tolower(source), tolower(nms))
  out <- rep(fallback, length(source))
  matched <- !is.na(idx)
  out[matched] <- as.numeric(source_weights[idx[matched]])
  out[!is.finite(out)] <- fallback
  pmin(1, pmax(0, out))
}

.tw_as_herb_weights <- function(herb_targets, target_col = "target") {
  if (is.character(herb_targets) && !is.data.frame(herb_targets)) {
    return(prepare_herb_target_weights(herb_targets))
  }
  if (!is.data.frame(herb_targets)) {
    stop("herb_targets must be a data frame or character vector.", call. = FALSE)
  }
  if (!"herb_target_weight" %in% names(herb_targets)) {
    return(prepare_herb_target_weights(herb_targets, target_col = target_col))
  }
  if (!target_col %in% names(herb_targets)) {
    stop("Target column not found in herb_targets: ", target_col, call. = FALSE)
  }
  out <- herb_targets
  names(out)[names(out) == target_col] <- "target"
  out$target <- .tw_clean_symbol(out$target)
  out
}

.tw_as_disease_weights <- function(disease_targets, disease_gene_col = "symbol", disease_weight_col = "disease_weight") {
  if (is.character(disease_targets) && !is.data.frame(disease_targets)) {
    return(prepare_disease_weights(user = disease_targets))
  }
  if (!is.data.frame(disease_targets)) {
    stop("disease_targets must be a data frame or character vector.", call. = FALSE)
  }
  if (!disease_gene_col %in% names(disease_targets)) {
    alt <- .tw_first_existing(names(disease_targets), c("symbol", "gene_symbol", "target", "gene", "names"))
    if (is.na(alt)) stop("Cannot find disease gene column.", call. = FALSE)
    disease_gene_col <- alt
  }
  if (!disease_weight_col %in% names(disease_targets)) {
    return(prepare_disease_weights(user = disease_targets, gene_col = disease_gene_col))
  }
  out <- disease_targets
  names(out)[names(out) == disease_gene_col] <- "symbol"
  names(out)[names(out) == disease_weight_col] <- "disease_weight"
  out$symbol <- .tw_clean_symbol(out$symbol)
  out$disease_weight <- suppressWarnings(as.numeric(out$disease_weight))
  out$disease_weight[!is.finite(out$disease_weight)] <- 0
  out$disease_weight <- pmin(1, pmax(0, out$disease_weight))
  out
}

.tw_as_target_set_list <- function(target_sets, set_col = NULL, target_col = "target") {
  if (is.list(target_sets) && !is.data.frame(target_sets)) {
    if (length(target_sets) == 0L) {
      stop("target_sets must contain at least one target set.", call. = FALSE)
    }
    set_names <- names(target_sets)
    if (is.null(set_names)) set_names <- rep("", length(target_sets))
    missing_names <- is.na(set_names) | !nzchar(set_names)
    set_names[missing_names] <- paste0("set", which(missing_names))
    out <- lapply(target_sets, .tw_as_herb_weights, target_col = target_col)
    names(out) <- make.unique(set_names)
    return(out)
  }

  if (is.character(target_sets) && !is.data.frame(target_sets)) {
    return(list(set1 = prepare_herb_target_weights(target_sets)))
  }

  if (!is.data.frame(target_sets)) {
    stop("target_sets must be a character vector, data frame, or named list.", call. = FALSE)
  }

  if (!is.null(set_col)) {
    if (!set_col %in% names(target_sets)) {
      stop("set_col not found in target_sets: ", set_col, call. = FALSE)
    }
    split_df <- split(target_sets, target_sets[[set_col]])
    out <- lapply(split_df, .tw_as_herb_weights, target_col = target_col)
    names(out) <- make.unique(names(split_df))
    return(out)
  }

  list(set1 = .tw_as_herb_weights(target_sets, target_col = target_col))
}

.tw_match_graph_symbols <- function(symbols, vertex_names) {
  symbols_clean <- .tw_clean_symbol(symbols)
  vertex_clean <- .tw_clean_symbol(vertex_names)
  idx <- match(symbols_clean, vertex_clean)
  out <- vertex_names[idx]
  out[is.na(idx)] <- NA_character_
  out
}

.tw_collapse_disease_weights <- function(x) {
  split_x <- split(x, x$symbol)
  out <- lapply(split_x, function(z) {
    row <- data.frame(
      symbol = z$symbol[[1]],
      symbol_query = paste(unique(z$symbol_query), collapse = ";"),
      disease_weight = max(z$disease_weight, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if ("sources" %in% names(z)) {
      row$sources <- paste(sort(unique(unlist(strsplit(paste(z$sources, collapse = ";"), ";", fixed = TRUE)))), collapse = ";")
    }
    if ("n_sources" %in% names(z)) row$n_sources <- max(z$n_sources, na.rm = TRUE)
    if ("n_evidence" %in% names(z)) row$n_evidence <- sum(z$n_evidence, na.rm = TRUE)
    if ("source_detail" %in% names(z)) row$source_detail <- paste(unique(z$source_detail), collapse = ";")
    row
  })
  do.call(rbind, out)
}

.tw_collapse_herb_weights <- function(x) {
  split_x <- split(x, x$target)
  out <- lapply(split_x, function(z) {
    row <- data.frame(
      target = z$target[[1]],
      target_query = paste(unique(z$target_query), collapse = ";"),
      herb_target_weight = sum(z$herb_target_weight, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    if ("n_compounds" %in% names(z)) row$n_compounds <- .tw_sum_optional(z$n_compounds)
    if ("n_herbs" %in% names(z)) row$n_herbs <- .tw_sum_optional(z$n_herbs)
    if ("n_records" %in% names(z)) row$n_records <- .tw_sum_optional(z$n_records)
    row
  })
  do.call(rbind, out)
}

.tw_sum_optional <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
}

.tw_edge_distances <- function(ppi, weight_attr, distance_attr) {
  edge_attrs <- edge_attr_names(ppi)
  if (!is.null(distance_attr) && distance_attr %in% edge_attrs) {
    w <- edge_attr(ppi, distance_attr)
    w[!is.finite(w) | w <= 0] <- 1
    return(w)
  }
  if (!weight_attr %in% edge_attrs) {
    return(NULL)
  }
  score <- suppressWarnings(as.numeric(edge_attr(ppi, weight_attr)))
  if (max(score, na.rm = TRUE) > 1) score <- score / 1000
  score[!is.finite(score) | score <= 0] <- .Machine$double.eps
  1 / score
}

.tw_weighted_proximity_matrix <- function(dist_mat, disease_weights) {
  disease_weights <- as.numeric(disease_weights)
  disease_weights[!is.finite(disease_weights)] <- 0
  if (sum(disease_weights) == 0) disease_weights[] <- 1
  inv <- 1 / (1 + dist_mat)
  inv[!is.finite(inv)] <- 0
  as.numeric(inv %*% disease_weights / sum(disease_weights))
}

.tw_weighted_proximity_matrix_no_self <- function(dist_mat, disease_weights) {
  disease_weights <- as.numeric(disease_weights)
  disease_weights[!is.finite(disease_weights)] <- 0
  if (sum(disease_weights) == 0) disease_weights[] <- 1

  inv <- 1 / (1 + dist_mat)
  inv[!is.finite(inv)] <- 0

  denom <- sum(disease_weights)
  numerator <- as.numeric(inv %*% disease_weights)

  row_ids <- rownames(dist_mat)
  col_ids <- colnames(dist_mat)
  self_idx <- match(row_ids, col_ids)
  has_self <- !is.na(self_idx)
  if (any(has_self)) {
    self_weight <- disease_weights[self_idx[has_self]]
    numerator[has_self] <- numerator[has_self] - self_weight
    denom_self <- denom - self_weight
    out <- numerator / denom
    out[has_self] <- ifelse(denom_self > 0, numerator[has_self] / denom_self, NA_real_)
  } else {
    out <- numerator / denom
  }
  out[!is.finite(out)] <- NA_real_
  out
}

.tw_degree_bins <- function(degree, degree_bins) {
  degree_bins <- max(1L, as.integer(degree_bins))
  if (length(unique(degree)) <= degree_bins) {
    return(as.integer(factor(degree, levels = sort(unique(degree)))))
  }
  probs <- seq(0, 1, length.out = degree_bins + 1L)
  brks <- unique(quantile(degree, probs = probs, na.rm = TRUE, type = 8))
  if (length(brks) <= 2L) {
    return(rep(1L, length(degree)))
  }
  as.integer(cut(degree, breaks = brks, include.lowest = TRUE, labels = FALSE))
}

.tw_sample_degree_matched_set <- function(set_nodes, bins, all_proximity) {
  vnames <- names(all_proximity)
  random_nodes <- vapply(set_nodes, function(target) {
    pool <- names(bins)[bins == bins[[target]] & is.finite(all_proximity[names(bins)])]
    if (length(pool) < 10L) {
      pool <- vnames[is.finite(all_proximity)]
    }
    sample(pool, 1L)
  }, character(1))
  random_nodes
}

.tw_directlink_stats <- function(ppi, set_nodes, disease_nodes) {
  set_nodes <- setdiff(set_nodes, disease_nodes)
  if (length(set_nodes) == 0L) {
    return(list(
      directlink_edges = 0L,
      directlink_target_n = 0L,
      directlink_disease_n = 0L,
      directlink_targets = character(),
      directlink_disease = character()
    ))
  }
  unweighted_dist <- distances(
    ppi,
    v = set_nodes,
    to = disease_nodes,
    weights = NA,
    mode = "all"
  )
  direct_mat <- unweighted_dist == 1
  direct_mat[is.na(direct_mat)] <- FALSE
  directlink_targets <- rownames(direct_mat)[rowSums(direct_mat) > 0]
  directlink_disease <- colnames(direct_mat)[colSums(direct_mat) > 0]
  list(
    directlink_edges = sum(direct_mat),
    directlink_target_n = length(directlink_targets),
    directlink_disease_n = length(directlink_disease),
    directlink_targets = directlink_targets,
    directlink_disease = directlink_disease
  )
}

.tw_warn_coverage <- function(diagnostics, warn_coverage = 0.5) {
  if (!is.null(warn_coverage) && is.finite(warn_coverage)) {
    if (diagnostics$herb_ppi_coverage < warn_coverage) {
      warning(
        sprintf(
          "Only %.1f%% of herb targets are present in the PPI background.",
          100 * diagnostics$herb_ppi_coverage
        ),
        call. = FALSE
      )
    }
    if (diagnostics$disease_ppi_coverage < warn_coverage) {
      warning(
        sprintf(
          "Only %.1f%% of disease targets are present in the PPI background.",
          100 * diagnostics$disease_ppi_coverage
        ),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

.tw_warn_set_coverage <- function(diagnostics, warn_coverage = 0.5) {
  if (!is.null(warn_coverage) && is.finite(warn_coverage)) {
    low_sets <- diagnostics$set_coverage$set[
      diagnostics$set_coverage$target_ppi_coverage < warn_coverage
    ]
    if (length(low_sets) > 0L) {
      warning(
        "Some target sets have low PPI coverage: ",
        paste(low_sets, collapse = ", "),
        call. = FALSE
      )
    }
    if (diagnostics$disease_ppi_coverage < warn_coverage) {
      warning(
        sprintf(
          "Only %.1f%% of disease targets are present in the PPI background.",
          100 * diagnostics$disease_ppi_coverage
        ),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

.tw_norm01 <- function(x) {
  x <- as.numeric(x)
  finite <- is.finite(x)
  if (!any(finite)) return(rep(0.5, length(x)))
  rng <- range(x[finite], na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0.5, length(x)))
  out <- (x - rng[[1]]) / diff(rng)
  out[!is.finite(out)] <- 0.5
  out
}
