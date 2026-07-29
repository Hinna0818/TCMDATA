.edge_score_by_pair <- function(graph) {
  edge_ends <- igraph::ends(graph, igraph::E(graph), names = TRUE)
  edge_keys <- apply(
    edge_ends,
    1L,
    function(x) paste(sort(x), collapse = "--")
  )
  scores <- stats::setNames(igraph::edge_attr(graph, "score"), edge_keys)
  scores[order(names(scores))]
}

test_that("PPI clustering outputs retain score and are plot-ready", {
  data(demo_ppi)
  runners <- list(
    mcode = function(graph) run_mcode(graph),
    louvain = function(graph) run_louvain(graph),
    fastgreedy = function(graph) run_fastgreedy(graph)
  )
  if (requireNamespace("Matrix", quietly = TRUE)) {
    runners$mcl <- function(graph) run_MCL(graph)
  }

  expected_scores <- .edge_score_by_pair(demo_ppi)
  results <- lapply(runners, function(run_cluster) {
    suppressMessages(run_cluster(demo_ppi))
  })

  for (result in results) {
    expect_true("score" %in% igraph::edge_attr_names(result))
    expect_equal(.edge_score_by_pair(result), expected_scores)
    expect_true("cluster" %in% igraph::vertex_attr_names(result))
    expect_equal(
      length(igraph::vertex_attr(result, "cluster")),
      igraph::vcount(result)
    )

    plot <- ggppi_network(
      result,
      size_metric = "betweenness",
      cluster = "cluster",
      candidate = if ("candidate" %in% igraph::vertex_attr_names(result)) {
        "candidate"
      } else {
        NULL
      },
      score = "score",
      show_text = FALSE
    )
    expect_silent(ggplot2::ggplot_build(plot))
    expect_false(is.null(plot$scales$get_scales("colour")))
    expect_false(is.null(plot$scales$get_scales("linewidth")))
  }

  expect_true("candidate" %in% igraph::vertex_attr_names(results$mcode))
  expect_equal(
    igraph::vertex_attr(results$mcode, "candidate"),
    igraph::vertex_attr(results$mcode, "mcode_is_seed")
  )
  expect_false(is.null(
    ggppi_network(
      results$mcode,
      cluster = "cluster",
      candidate = "candidate",
      score = "score",
      show_text = FALSE
    )$scales$get_scales("shape")
  ))
})

test_that("computed node metrics flow directly into ggppi_network", {
  graph <- igraph::make_ring(6)
  graph <- igraph::set_edge_attr(
    graph,
    "score",
    value = seq(0.4, 0.9, length.out = igraph::ecount(graph))
  )
  graph <- suppressMessages(run_louvain(graph))
  invisible(capture.output(
    graph <- suppressMessages(compute_nodeinfo(graph, seed = 1))
  ))

  expect_true(all(c("cluster", "score", "MCC", "betweenness") %in% c(
    igraph::vertex_attr_names(graph),
    igraph::edge_attr_names(graph)
  )))

  plot <- ggppi_network(
    graph,
    metric = "MCC",
    size_metric = "betweenness",
    cluster = "cluster",
    score = "score",
    show_text = FALSE
  )
  expect_silent(ggplot2::ggplot_build(plot))
})

test_that("run_mcode keeps an igraph return type for small graphs", {
  graph <- igraph::make_empty_graph(1)
  result <- NULL
  expect_warning(result <- run_mcode(graph), "too small")

  expect_s3_class(result, "igraph")
  expect_true(all(c("cluster", "candidate") %in% igraph::vertex_attr_names(result)))
})
