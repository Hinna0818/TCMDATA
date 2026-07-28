test_that("ggtcm_network plots a layered tripartite data frame", {
  data <- data.frame(
    herb = c("H1", "H1", "H2", "H2"),
    molecule = c("C1", "C2", "C1", "C3"),
    target = c("T1", "T2", "T1", "T2")
  )

  p <- ggtcm_network(data, layout = "layered", show_text = FALSE)

  expect_s3_class(p, "ggplot")
  expect_setequal(unique(p$data$.tcm_type), c("Herb", "Molecule", "Target"))
  expect_equal(
    sort(unique(p$data$x)),
    c(0, 1, 2)
  )
})

test_that("ggtcm_network uses concentric type layers and size hierarchy", {
  data <- data.frame(
    herb = c("H1", "H1", "H2", "H2"),
    molecule = c("C1", "C2", "C1", "C3"),
    target = c("T1", "T2", "T1", "T3")
  )

  p <- ggtcm_network(data, show_text = FALSE)
  radius <- sqrt(p$data$x^2 + p$data$y^2)

  expect_equal(
    unique(round(radius[p$data$.tcm_type == "Herb"], 2)),
    0.22
  )
  expect_equal(
    unique(round(radius[p$data$.tcm_type == "Molecule"], 2)),
    1.05
  )
  expect_equal(
    unique(round(radius[p$data$.tcm_type == "Target"], 2)),
    2
  )
  expect_gt(
    min(p$data$.tcm_plot_size[p$data$.tcm_type == "Herb"]),
    max(p$data$.tcm_plot_size[p$data$.tcm_type == "Molecule"])
  )
  expect_gt(
    min(p$data$.tcm_plot_size[p$data$.tcm_type == "Molecule"]),
    max(p$data$.tcm_plot_size[p$data$.tcm_type == "Target"])
  )
  expect_equal(
    unname(p$scales$get_scales("shape")$palette(3)),
    rep(21, 3)
  )
})

test_that("a single herb is placed at the concentric origin", {
  data <- data.frame(
    herb = c("H1", "H1"),
    molecule = c("C1", "C2"),
    target = c("T1", "T2")
  )

  p <- ggtcm_network(data, show_text = FALSE)
  herb <- p$data[p$data$.tcm_type == "Herb", , drop = FALSE]

  expect_equal(herb$x, 0)
  expect_equal(herb$y, 0)
})

test_that("saved concentric coordinates retain type-specific node sizes", {
  data <- data.frame(
    herb = c("H1", "H1", "H2", "H2"),
    molecule = c("C1", "C2", "C1", "C3"),
    target = c("T1", "T2", "T1", "T3")
  )
  graph <- .ggtcm_as_graph(
    data,
    herb_col = "herb",
    compound_col = "molecule",
    target_col = "target"
  )
  graph <- .ggtcm_prepare_attributes(graph, "degree")
  xy <- .ggtcm_concentric_layout(
    graph,
    c(Herb = 0.22, Molecule = 1.05, Target = 2)
  )
  coordinates <- data.frame(
    name = rownames(xy),
    x = xy[, 1] * 240,
    y = xy[, 2] * 240
  )
  class(coordinates) <- c("tcm_network_layout", "data.frame")
  attr(coordinates, "layout_type") <- "concentric"

  p <- ggtcm_network(data, layout = coordinates, show_text = FALSE)

  expect_true(".tcm_plot_size" %in% names(p$data))
  expect_gt(
    min(p$data$.tcm_plot_size[p$data$.tcm_type == "Herb"]),
    max(p$data$.tcm_plot_size[p$data$.tcm_type == "Molecule"])
  )
})

test_that("ggtcm_network keeps identical labels in different node types", {
  data <- data.frame(
    herb = "shared",
    molecule = "shared",
    target = "shared"
  )

  p <- ggtcm_network(data, show_text = FALSE)

  expect_equal(nrow(p$data), 3L)
  expect_equal(length(unique(p$data$name)), 3L)
  expect_equal(unique(p$data$.tcm_label), "shared")
})

test_that("ggtcm_network accepts prepared igraph objects", {
  data <- data.frame(
    herb = c("H1", "H1"),
    molecule = c("C1", "C2"),
    target = c("T1", "T2")
  )
  graph <- prepare_herb_graph(data)

  p <- ggtcm_network(graph, node_size = "betweenness", show_text = FALSE)

  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), igraph::vcount(graph))
})
