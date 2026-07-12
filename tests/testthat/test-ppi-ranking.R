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
