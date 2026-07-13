test_that("compute_nodeinfo stores reciprocal eccentricity centrality", {
  g <- igraph::graph_from_edgelist(
    matrix(c(1, 2, 2, 3, 3, 4, 4, 5), ncol = 2, byrow = TRUE),
    directed = FALSE
  )
  igraph::V(g)$name <- LETTERS[1:5]
  igraph::E(g)$score <- 1

  annotated <- suppressMessages(compute_nodeinfo(g, seed = 1))
  ecc <- igraph::vertex_attr(annotated, "eccentricity")
  ecc_centrality <- igraph::vertex_attr(annotated, "eccentricity_centrality")

  expect_equal(ecc_centrality, 1 / ecc)
  expect_gt(ecc_centrality[[3]], ecc_centrality[[1]])
})

test_that("calculate_metric_weights shares weight across positive correlation groups", {
  g <- igraph::make_empty_graph(8)
  igraph::V(g)$m1 <- 1:8
  igraph::V(g)$m2 <- (1:8)^2
  igraph::V(g)$m3 <- -(1:8)
  igraph::V(g)$m4 <- c(1, 8, 2, 7, 3, 6, 4, 5)

  result <- calculate_metric_weights(
    g,
    metrics = c("m1", "m2", "m3", "m4"),
    use_weight = FALSE,
    correlation_threshold = 0.8
  )

  expect_s3_class(result, "tcm_metric_weights")
  expect_equal(sum(result$weights), 1)
  expect_equal(
    result$weights,
    c(m1 = 1 / 6, m2 = 1 / 6, m3 = 1 / 3, m4 = 1 / 3)
  )
  groups <- stats::setNames(
    result$metric_groups$cluster,
    result$metric_groups$metric
  )
  expect_equal(groups[["m1"]], groups[["m2"]])
  expect_false(groups[["m1"]] == groups[["m3"]])
  expect_false(groups[["m1"]] == groups[["m4"]])
  expect_equal(
    names(result$metric_groups),
    c("metric", "cluster", "cluster_size", "weight")
  )
  expect_equal(result$correlation[["m1", "m3"]], -1)
  expect_equal(result$params$correlation_method, "spearman")
  expect_equal(result$params$clustering_method, "complete")
  expect_equal(dim(result$correlation), c(4, 4))
})

test_that("calculate_metric_weights requires three paired observations", {
  g <- igraph::make_empty_graph(6)
  igraph::V(g)$m1 <- c(1, 2, 3, NA, NA, NA)
  igraph::V(g)$m2 <- c(1, 2, NA, 3, NA, NA)

  expect_warning(
    result <- calculate_metric_weights(
      g,
      metrics = c("m1", "m2"),
      use_weight = FALSE
    ),
    "fewer than 3 paired observations"
  )

  expect_true(is.na(result$correlation[["m1", "m2"]]))
  expect_equal(result$pairwise_n[["m1", "m2"]], 2)
  expect_false(
    result$metric_groups$cluster[[1]] ==
      result$metric_groups$cluster[[2]]
  )
})

test_that("complete linkage enforces the positive correlation threshold", {
  correlation <- matrix(
    c(
      1.00, 0.91, 0.70,
      0.91, 1.00, 0.91,
      0.70, 0.91, 1.00
    ),
    nrow = 3,
    dimnames = list(paste0("m", 1:3), paste0("m", 1:3))
  )
  basis <- qr.Q(qr(scale(cbind(
    seq_len(1000),
    sin(seq_len(1000)),
    cos(seq_len(1000) / 7)
  ))))
  metric_values <- basis %*% chol(correlation)
  colnames(metric_values) <- colnames(correlation)
  g <- igraph::make_empty_graph(nrow(metric_values))
  for (metric in colnames(metric_values)) {
    g <- igraph::set_vertex_attr(g, metric, value = metric_values[, metric])
  }

  result <- calculate_metric_weights(
    g,
    metrics = colnames(metric_values),
    use_weight = FALSE,
    correlation_method = "pearson",
    correlation_threshold = 0.8
  )

  groups <- stats::setNames(
    result$metric_groups$cluster,
    result$metric_groups$metric
  )
  expect_false(
    groups[["m1"]] == groups[["m2"]] &&
      groups[["m2"]] == groups[["m3"]]
  )
})

test_that("calculate_metric_weights excludes infinite values from correlation", {
  g <- igraph::make_empty_graph(5)
  igraph::V(g)$m1 <- 1:5
  igraph::V(g)$m2 <- c(2, 4, 6, 8, Inf)

  result <- calculate_metric_weights(
    g,
    metrics = c("m1", "m2"),
    use_weight = FALSE,
    correlation_method = "pearson"
  )

  expect_equal(result$correlation[["m1", "m2"]], 1)
  expect_equal(
    result$metric_groups$cluster[[1]],
    result$metric_groups$cluster[[2]]
  )
})

test_that("calculate_metric_weights resolves metrics for ranking", {
  g <- igraph::make_empty_graph(6)
  igraph::V(g)$name <- LETTERS[1:6]
  igraph::V(g)$degree <- c(1, 4, 2, 6, 3, 5)
  igraph::V(g)$betweenness <- c(3, 1, 5, 2, 6, 4)
  igraph::V(g)$betweenness_w <- c(1, 3, 6, 2, 5, 4)
  igraph::V(g)$closeness <- c(6, 4, 1, 5, 2, 3)
  igraph::V(g)$closeness_w <- c(2, 6, 3, 1, 5, 4)
  igraph::V(g)$eccentricity <- c(4, 2, 5, 1, 6, 3)
  igraph::V(g)$eccentricity_centrality <-
    1 / igraph::V(g)$eccentricity

  result <- calculate_metric_weights(
    g,
    metrics = c("degree", "betweenness", "closeness", "eccentricity"),
    use_weight = TRUE
  )

  expect_equal(
    names(result$weights),
    c("degree", "betweenness_w", "closeness_w", "eccentricity_centrality")
  )
  ranked <- suppressMessages(rank_ppi_nodes(
    g,
    metrics = names(result$weights),
    weights = result$weights,
    use_weight = FALSE
  ))
  expect_equal(ranked$scoring$weights, result$weights)
})

test_that("calculate_metric_weights drops unavailable and constant metrics", {
  g <- igraph::make_empty_graph(5)
  igraph::V(g)$varying <- 1:5
  igraph::V(g)$constant <- rep(1, 5)

  expect_warning(
    result <- calculate_metric_weights(
      g,
      metrics = c("varying", "constant", "missing"),
      use_weight = FALSE
    ),
    "dropped"
  )

  expect_equal(result$weights, c(varying = 1))
  expect_null(result$clustering)
  expect_equal(
    result$dropped_metrics$reason,
    c("non_informative", "unavailable")
  )
})

test_that("calculate_metric_weights validates its inputs", {
  g <- igraph::make_empty_graph(5)
  igraph::V(g)$metric <- 1:5

  expect_error(
    calculate_metric_weights(g, "metric", correlation_threshold = 1.1),
    "correlation_threshold"
  )
  expect_error(
    calculate_metric_weights(g, character()),
    "metrics"
  )
  expect_error(
    calculate_metric_weights(data.frame(x = 1:5), "x"),
    "igraph"
  )
})

test_that("rank_ppi_nodes preserves requested metrics", {
  g <- igraph::make_star(5, mode = "undirected")
  igraph::V(g)$name <- LETTERS[1:5]
  igraph::V(g)$degree <- igraph::degree(g)
  igraph::V(g)$MCC <- c(1, 4, 3, 2, 1)

  ranked <- suppressMessages(rank_ppi_nodes(
    g,
    metrics = c("degree", "MCC"),
    use_weight = TRUE
  ))

  expect_equal(ranked$scoring$metrics, c("degree", "MCC"))
  expect_true(all(c("degree_norm", "MCC_norm") %in% names(ranked$table)))
  expect_false("eccentricity_centrality_norm" %in% names(ranked$table))
  expect_equal(length(ranked), 3L)
})

test_that("rank_ppi_nodes maps named weights after metric resolution", {
  g <- igraph::graph_from_edgelist(
    matrix(c(1, 2, 2, 3, 3, 4), ncol = 2, byrow = TRUE),
    directed = FALSE
  )
  igraph::V(g)$name <- LETTERS[1:4]
  igraph::V(g)$degree <- igraph::degree(g)
  ecc <- igraph::eccentricity(g)
  igraph::V(g)$eccentricity_centrality <- 1 / ecc

  ranked <- suppressMessages(rank_ppi_nodes(
    g,
    metrics = c("degree", "eccentricity_centrality"),
    weights = c(eccentricity_centrality = 3, degree = 1),
    use_weight = FALSE
  ))

  expect_equal(
    ranked$scoring$weights,
    c(degree = 0.25, eccentricity_centrality = 0.75)
  )
  expected <- 0.25 * ranked$table$degree_norm +
    0.75 * ranked$table$eccentricity_centrality_norm
  expect_equal(ranked$table$Score_network, expected)
})

test_that("rank_ppi_nodes validates weights", {
  g <- igraph::make_ring(4)
  igraph::V(g)$degree <- igraph::degree(g)
  igraph::V(g)$MCC <- seq_len(4)

  expect_error(
    rank_ppi_nodes(g, c("degree", "MCC"), weights = c(1, -1), use_weight = FALSE),
    "non-negative"
  )
  expect_error(
    rank_ppi_nodes(g, c("degree", "MCC"), weights = c(0, 0), use_weight = FALSE),
    "positive sum"
  )
  expect_error(
    rank_ppi_nodes(g, c("degree", "MCC"), weights = 1, use_weight = FALSE),
    "match"
  )
})

test_that("rank_ppi_nodes warns and falls back when weighted metrics are absent", {
  g <- igraph::graph_from_edgelist(
    matrix(c(1, 2, 2, 3, 3, 4), ncol = 2, byrow = TRUE),
    directed = FALSE
  )
  igraph::V(g)$betweenness <- igraph::betweenness(g)
  igraph::V(g)$closeness <- igraph::closeness(g)

  expect_warning(
    ranked <- suppressMessages(rank_ppi_nodes(
      g,
      metrics = c("betweenness", "closeness"),
      use_weight = TRUE
    )),
    "falling back"
  )
  expect_equal(ranked$scoring$metrics, c("betweenness", "closeness"))
})
