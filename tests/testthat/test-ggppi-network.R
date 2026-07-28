test_that("ggppi_network maps degree to node colour and size", {
  graph <- igraph::make_ring(6)
  plot <- ggppi_network(graph, metric = "degree", show_text = FALSE)

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$data[[".tcm_colour_metric"]], rep(2, 6))
  expect_equal(plot$data[[".tcm_size_metric"]], rep(2, 6))
  expect_length(plot$layers, 2)
})

test_that("ggppi_network accepts numeric vertex attributes", {
  graph <- igraph::make_star(6, mode = "undirected")
  graph <- igraph::set_vertex_attr(graph, "hub_score", value = 1:6)
  plot <- ggppi_network(graph, metric = "hub_score", label_top = 2)

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$data[[".tcm_colour_metric"]], 1:6)
  expect_equal(plot$data[[".tcm_size_metric"]], 1:6)
  expect_length(plot$layers, 3)
  expect_equal(nrow(plot$layers[[3]]$data), igraph::vcount(graph))
  expect_equal(
    sum(nzchar(plot$layers[[3]]$data[[".tcm_label"]])),
    2
  )
})

test_that("ggppi_network maps a separate metric to node size", {
  graph <- igraph::make_star(6, mode = "undirected")
  plot <- ggppi_network(
    graph,
    metric = "degree",
    size_metric = "betweenness",
    show_text = FALSE
  )

  expect_equal(plot$data[[".tcm_colour_metric"]], igraph::degree(graph))
  expect_equal(
    plot$data[[".tcm_size_metric"]],
    igraph::betweenness(graph, normalized = TRUE)
  )
  expect_s3_class(plot$guides$guides$size, "GuideLegend")
})

test_that("ggppi_network validates graph and metric inputs", {
  graph <- igraph::make_ring(5)
  expect_error(ggppi_network(data.frame()), "igraph object")
  expect_error(ggppi_network(graph, metric = "not_a_metric"), "was not found")
  expect_error(ggppi_network(graph, size_metric = "not_a_metric"), "was not found")

  graph <- igraph::set_vertex_attr(graph, "group", value = letters[1:5])
  expect_error(ggppi_network(graph, metric = "group"), "must be numeric")
  expect_error(ggppi_network(graph, show_text = NA), "TRUE or FALSE")
  expect_error(ggppi_network(graph, edge_width = -1), "non-negative")
})

test_that("ggppi_network exposes edge colour and width controls", {
  graph <- igraph::make_ring(5)
  plot <- ggppi_network(
    graph,
    show_text = FALSE,
    edge_color = "#333333",
    edge_width = 0.8
  )

  expect_equal(plot$layers[[1]]$aes_params$colour, "#333333")
  expect_equal(plot$layers[[1]]$aes_params$linewidth, 0.8)
})

test_that("ggppi_network adds optional node and edge mappings on request", {
  graph <- igraph::make_ring(6)
  graph <- igraph::set_vertex_attr(
    graph,
    "cluster",
    value = rep(c("Module 1", "Module 2"), 3)
  )
  graph <- igraph::set_vertex_attr(
    graph,
    "candidate",
    value = rep(c(TRUE, FALSE), 3)
  )
  graph <- igraph::set_edge_attr(
    graph,
    "score",
    value = seq(0.4, 0.9, length.out = igraph::ecount(graph))
  )

  plot <- ggppi_network(
    graph,
    cluster = "cluster",
    candidate = "candidate",
    score = "score",
    show_text = FALSE
  )

  expect_equal(
    plot$data[[".tcm_cluster"]],
    rep(c("Module 1", "Module 2"), 3)
  )
  expect_equal(
    plot$data[[".tcm_candidate"]],
    rep(c("TRUE", "FALSE"), 3)
  )
  expect_false(is.null(plot$scales$get_scales("colour")))
  expect_false(is.null(plot$scales$get_scales("shape")))
  expect_false(is.null(plot$scales$get_scales("linewidth")))
  expect_silent(ggplot2::ggplot_build(plot))
})

test_that("ggppi_network keeps optional mappings disabled by default", {
  expect_null(formals(ggppi_network)$cluster)
  expect_null(formals(ggppi_network)$candidate)
  expect_null(formals(ggppi_network)$score)

  plot <- ggppi_network(igraph::make_ring(5), show_text = FALSE)
  expect_null(plot$scales$get_scales("colour"))
  expect_null(plot$scales$get_scales("shape"))
  expect_null(plot$scales$get_scales("linewidth"))
})

test_that("ggppi_network validates optional mapping attributes", {
  graph <- igraph::make_ring(5)

  expect_error(
    ggppi_network(graph, cluster = "cluster"),
    "requested by 'cluster'"
  )
  expect_error(
    ggppi_network(graph, candidate = "candidate"),
    "requested by 'candidate'"
  )
  expect_error(
    ggppi_network(graph, score = "score"),
    "Edge attribute 'score' was not found"
  )

  graph <- igraph::set_edge_attr(graph, "score", value = letters[1:5])
  expect_error(
    ggppi_network(graph, score = "score"),
    "must be numeric"
  )
})
