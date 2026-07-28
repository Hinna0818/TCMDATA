#' Plot a herb-compound-target network
#'
#' @param x A data frame with herb, compound and target columns, or an
#'   `igraph` object with a `type` vertex attribute.
#' @param herb_col,compound_col,target_col Column names used when `x` is a
#'   data frame.
#' @param layout `"concentric"`, `"layered"`, an igraph layout, or saved node
#'   coordinates.
#' @param node_size Numeric vertex attribute or centrality name mapped to node
#'   size within each node type; use `NULL` for type-specific constant sizes.
#' @param show_text Whether to label selected nodes.
#' @param label_top Number of labels by node type.
#' @param highlight Optional node labels to highlight with an accent border.
#' @param colors,shapes Named values for herbs, compounds and targets.
#' @param type_sizes Node sizes for herbs, compounds and targets in the
#'   concentric layout.
#' @param ring_radii Radii of the three concentric node layers. A single herb
#'   is placed exactly at the origin.
#' @param show_rings Whether to draw subtle guides for the compound and target
#'   rings.
#' @param size_range Node-size range for non-concentric layouts.
#' @param edge_color,edge_alpha,edge_width Edge appearance.
#' @param seed Random seed used by layouts and labels.
#' @param base_size,base_family Base text size and family.
#' @param title Optional title.
#' @param ... Arguments passed to the selected layout.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' hct <- search_herb(c("lingzhi", "huangqi"), "Herb_pinyin_name")
#' ggtcm_network(hct, label_top = c(Herb = Inf, Molecule = 8, Target = 10))
#' }
#'
#' @importFrom ggplot2 aes annotate coord_equal element_text expansion geom_path geom_point ggplot guide_legend guides labs margin scale_colour_manual scale_fill_manual scale_shape_manual scale_size_continuous scale_size_identity scale_x_continuous scale_y_continuous theme
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggtangle geom_edge
#' @importFrom igraph graph_from_data_frame is_igraph layout_with_sugiyama simplify set_vertex_attr vcount vertex_attr vertex_attr_names
#' @importFrom rlang .data
#' @importFrom stats complete.cases setNames
#' @importFrom utils head
#' @export
ggtcm_network <- function(
    x,
    herb_col = "herb",
    compound_col = "molecule",
    target_col = "target",
    layout = "concentric",
    node_size = "degree",
    show_text = TRUE,
    label_top = c(Herb = Inf, Molecule = 12, Target = 12),
    highlight = NULL,
    colors = c(
      Herb = "#D94B2B",
      Molecule = "#F07A3E",
      Target = "#F5B06A"
    ),
    shapes = c(Herb = 21, Molecule = 21, Target = 21),
    type_sizes = c(Herb = 11, Molecule = 6, Target = 4),
    ring_radii = c(Herb = 0.22, Molecule = 1.05, Target = 2),
    show_rings = TRUE,
    size_range = c(3, 8),
    edge_color = "#C8C8C8",
    edge_alpha = 0.35,
    edge_width = 0.3,
    seed = 42L,
    base_size = 7,
    base_family = "Arial",
    title = NULL,
    ...) {
  base_family <- .resolve_tcm_font_family(base_family)
  graph <- .ggtcm_as_graph(
    x,
    herb_col = herb_col,
    compound_col = compound_col,
    target_col = target_col
  )
  graph <- .ggtcm_prepare_attributes(graph, node_size)

  type_levels <- c("Herb", "Molecule", "Target")
  colors <- .ggtcm_validate_named_values(colors, type_levels, "colors")
  shapes <- .ggtcm_validate_named_values(shapes, type_levels, "shapes")
  type_sizes <- .ggtcm_validate_named_values(
    type_sizes,
    type_levels,
    "type_sizes"
  )
  ring_radii <- .ggtcm_validate_named_values(
    ring_radii,
    type_levels,
    "ring_radii"
  )
  label_top <- .ggtcm_label_limits(label_top, type_levels)

  if (!is.logical(show_text) || length(show_text) != 1L || is.na(show_text)) {
    stop("'show_text' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(highlight) && !is.character(highlight)) {
    stop("'highlight' must be NULL or a character vector.", call. = FALSE)
  }
  if (!is.logical(show_rings) || length(show_rings) != 1L ||
      is.na(show_rings)) {
    stop("'show_rings' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(type_sizes) || any(!is.finite(type_sizes)) ||
      any(type_sizes <= 0)) {
    stop("'type_sizes' must contain positive numeric values.", call. = FALSE)
  }
  if (!is.numeric(ring_radii) || any(!is.finite(ring_radii)) ||
      any(ring_radii < 0) ||
      !(ring_radii[["Herb"]] < ring_radii[["Molecule"]] &&
        ring_radii[["Molecule"]] < ring_radii[["Target"]])) {
    stop("'ring_radii' must be non-negative and increase by node type.", call. = FALSE)
  }
  if (!is.numeric(size_range) || length(size_range) != 2L ||
      any(!is.finite(size_range)) || any(size_range <= 0) ||
      size_range[[1L]] > size_range[[2L]]) {
    stop("'size_range' must contain two positive, increasing values.", call. = FALSE)
  }
  if (!is.numeric(edge_alpha) || length(edge_alpha) != 1L ||
      !is.finite(edge_alpha) || edge_alpha < 0 || edge_alpha > 1) {
    stop("'edge_alpha' must be between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(edge_width) || length(edge_width) != 1L ||
      !is.finite(edge_width) || edge_width < 0) {
    stop("'edge_width' must be a non-negative number.", call. = FALSE)
  }

  is_concentric_name <- is.character(layout) && length(layout) == 1L &&
    identical(layout, "concentric")
  is_concentric_saved <- inherits(layout, "tcm_network_layout") &&
    identical(attr(layout, "layout_type"), "concentric")
  is_concentric <- is_concentric_name || is_concentric_saved
  is_layered <- is.character(layout) && length(layout) == 1L &&
    identical(layout, "layered")
  if (is_concentric) {
    graph <- set_vertex_attr(
      graph,
      ".tcm_plot_size",
      value = .ggtcm_concentric_sizes(graph, type_sizes, node_size)
    )
  }
  if (is_concentric_name) {
    layout <- .ggtcm_concentric_layout(graph, ring_radii)
  }
  if (is_layered) {
    layout <- .ggtcm_layered_layout(graph)
  }
  layout <- .ggnetwork_resolve_layout(graph, layout)

  node_labels <- vertex_attr(graph, ".tcm_label")
  highlight_values <- node_labels %in% highlight
  graph <- set_vertex_attr(
    graph,
    ".tcm_highlight",
    value = ifelse(highlight_values, "Highlighted", "Standard")
  )

  size_attribute <- if (is_concentric) ".tcm_plot_size" else ".tcm_size"
  node_layer <- if (any(highlight_values)) {
    geom_point(
      aes(
        fill = .data[[".tcm_type"]],
        shape = .data[[".tcm_type"]],
        size = .data[[size_attribute]],
        colour = .data[[".tcm_highlight"]]
      ),
      stroke = 0.65,
      alpha = 0.97
    )
  } else {
    geom_point(
      aes(
        fill = .data[[".tcm_type"]],
        shape = .data[[".tcm_type"]],
        size = .data[[size_attribute]]
      ),
      colour = "white",
      stroke = 0.35,
      alpha = 0.97
    )
  }

  p <- .ggnetwork_with_seed(seed, ggplot(graph, layout = layout, ...))
  if (is_concentric_name && isTRUE(show_rings)) {
    p <- p + .ggtcm_ring_layer(ring_radii)
  }
  p <- p +
    geom_edge(
      colour = edge_color,
      alpha = edge_alpha,
      linewidth = edge_width,
      lineend = "round"
    ) +
    node_layer +
    scale_fill_manual(
      values = colors,
      breaks = type_levels,
      labels = c(Herb = "Herb", Molecule = "Compound", Target = "Target"),
      name = "Node type"
    ) +
    scale_shape_manual(
      values = shapes,
      breaks = type_levels,
      labels = c(Herb = "Herb", Molecule = "Compound", Target = "Target"),
      name = "Node type"
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        override.aes = list(
          size = c(6, 4.5, 3.5),
          shape = unname(shapes)
        )
      ),
      shape = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.18))) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.16))) +
    coord_equal(clip = "off") +
    labs(title = title) +
    .theme_tcm_void(base_size = base_size, base_family = base_family) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "plain"),
      plot.margin = margin(7, 10, 7, 10)
    )

  if (is_concentric) {
    p <- p + scale_size_identity(guide = "none")
  } else {
    p <- p + scale_size_continuous(
      range = size_range,
      transform = "sqrt",
      name = if (is.null(node_size)) NULL else .ggnetwork_metric_title(node_size),
      guide = if (is.null(node_size)) "none" else guide_legend(order = 2)
    )
  }

  if (any(highlight_values)) {
    p <- p + scale_colour_manual(
      values = c(Standard = "#FFFFFF", Highlighted = "#B85C5C"),
      breaks = "Highlighted",
      name = NULL,
      guide = guide_legend(
        order = 3,
        override.aes = list(fill = "white", shape = 21, size = 4)
      )
    )
  }

  if (is_layered) {
    p <- p + annotate(
      "text",
      x = c(0, 1, 2),
      y = 1.13,
      label = c("Herb", "Compound", "Target"),
      family = base_family,
      fontface = "plain",
      size = base_size / 2.85,
      colour = "#3F3F3F"
    )
  }

  if (isTRUE(show_text)) {
    label_data <- .ggtcm_label_data(
      p$data,
      label_top,
      size_range,
      size_attribute = size_attribute,
      size_identity = is_concentric
    )
    if (nrow(label_data) > 0L) {
      p <- p + geom_text_repel(
        data = label_data,
        aes(
          label = .data[[".tcm_plot_label"]],
          point.size = .data[[".tcm_point_size"]]
        ),
        size = max(2, base_size / 3),
        family = base_family,
        fontface = "plain",
        colour = "#272727",
        box.padding = 0.35,
        point.padding = 0.2,
        min.segment.length = 0,
        segment.colour = "#A0A0A0",
        segment.size = 0.22,
        force = 1.6,
        max.time = 3,
        max.overlaps = Inf,
        seed = seed,
        show.legend = FALSE
      )
    }
  }

  p
}

.ggtcm_as_graph <- function(x, herb_col, compound_col, target_col) {
  if (is_igraph(x)) {
    if (vcount(x) == 0L) {
      stop("'x' must contain at least one node.", call. = FALSE)
    }
    return(simplify(
      x,
      remove.multiple = TRUE,
      remove.loops = TRUE,
      edge.attr.comb = "first"
    ))
  }
  if (!is.data.frame(x)) {
    stop("'x' must be a data frame or an igraph object.", call. = FALSE)
  }

  required <- c(herb_col, compound_col, target_col)
  if (!all(required %in% names(x))) {
    stop(
      "Input data must contain columns: ",
      paste(required, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  data <- x[, required, drop = FALSE]
  names(data) <- c("herb", "compound", "target")
  data[] <- lapply(data, function(value) trimws(as.character(value)))
  keep <- complete.cases(data) &
    data$herb != "" & data$compound != "" & data$target != ""
  data <- unique(data[keep, , drop = FALSE])
  if (nrow(data) == 0L) {
    stop("No complete herb-compound-target relationships were found.", call. = FALSE)
  }

  herb_id <- paste0("Herb::", data$herb)
  compound_id <- paste0("Molecule::", data$compound)
  target_id <- paste0("Target::", data$target)
  edges <- unique(rbind(
    data.frame(
      from = herb_id,
      to = compound_id,
      relation = "Herb-Compound",
      stringsAsFactors = FALSE
    ),
    data.frame(
      from = compound_id,
      to = target_id,
      relation = "Compound-Target",
      stringsAsFactors = FALSE
    )
  ))
  vertices <- unique(rbind(
    data.frame(
      name = herb_id,
      display_label = data$herb,
      type = "Herb",
      stringsAsFactors = FALSE
    ),
    data.frame(
      name = compound_id,
      display_label = data$compound,
      type = "Molecule",
      stringsAsFactors = FALSE
    ),
    data.frame(
      name = target_id,
      display_label = data$target,
      type = "Target",
      stringsAsFactors = FALSE
    )
  ))

  graph_from_data_frame(edges, directed = TRUE, vertices = vertices)
}

.ggtcm_prepare_attributes <- function(graph, node_size) {
  if (!"type" %in% vertex_attr_names(graph)) {
    stop("An igraph input must contain a 'type' vertex attribute.", call. = FALSE)
  }
  type <- as.character(vertex_attr(graph, "type"))
  type[tolower(type) %in% c("compound", "molecule")] <- "Molecule"
  type[tolower(type) == "herb"] <- "Herb"
  type[tolower(type) == "target"] <- "Target"
  invalid <- setdiff(unique(type), c("Herb", "Molecule", "Target"))
  if (length(invalid) > 0L) {
    stop(
      "Unsupported node types: ",
      paste(invalid, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  label <- if ("display_label" %in% vertex_attr_names(graph)) {
    as.character(vertex_attr(graph, "display_label"))
  } else if ("label" %in% vertex_attr_names(graph)) {
    as.character(vertex_attr(graph, "label"))
  } else {
    sub("^[^:]+::", "", as.character(vertex_attr(graph, "name")))
  }
  if ("label" %in% vertex_attr_names(graph)) {
    graph <- set_vertex_attr(
      graph,
      "label",
      value = as.character(vertex_attr(graph, "name"))
    )
  }
  size <- if (is.null(node_size)) {
    rep(1, vcount(graph))
  } else {
    if (!is.character(node_size) || length(node_size) != 1L ||
        is.na(node_size) || !nzchar(node_size)) {
      stop("'node_size' must be NULL or one metric name.", call. = FALSE)
    }
    .ggnetwork_clean_metric(.ggnetwork_metric(graph, node_size), node_size)
  }

  graph <- set_vertex_attr(graph, ".tcm_type", value = type)
  graph <- set_vertex_attr(graph, ".tcm_label", value = label)
  set_vertex_attr(graph, ".tcm_size", value = as.numeric(size))
}

.ggtcm_layered_layout <- function(graph) {
  type <- vertex_attr(graph, ".tcm_type")
  layer <- unname(c(Herb = 1L, Molecule = 2L, Target = 3L)[type])
  initial <- layout_with_sugiyama(graph, layers = layer)$layout[, 1L]
  label <- vertex_attr(graph, ".tcm_label")
  coordinates <- matrix(0, nrow = vcount(graph), ncol = 2L)
  coordinates[, 1L] <- unname(c(Herb = 0, Molecule = 1, Target = 2)[type])

  for (node_type in c("Herb", "Molecule", "Target")) {
    index <- which(type == node_type)
    index <- index[order(initial[index], label[index])]
    coordinates[index, 2L] <- if (length(index) == 1L) {
      0
    } else {
      seq(1, -1, length.out = length(index))
    }
  }
  rownames(coordinates) <- vertex_attr(graph, "name")
  colnames(coordinates) <- c("x", "y")
  coordinates
}

.ggtcm_concentric_layout <- function(graph, ring_radii) {
  type <- vertex_attr(graph, ".tcm_type")
  layer <- unname(c(Herb = 1L, Molecule = 2L, Target = 3L)[type])
  initial <- layout_with_sugiyama(graph, layers = layer)$layout[, 1L]
  label <- vertex_attr(graph, ".tcm_label")
  coordinates <- matrix(0, nrow = vcount(graph), ncol = 2L)

  for (node_type in c("Herb", "Molecule", "Target")) {
    index <- which(type == node_type)
    if (length(index) == 0L) {
      next
    }
    index <- index[order(initial[index], label[index])]
    radius <- ring_radii[[node_type]]
    if (node_type == "Herb" && length(index) == 1L) {
      radius <- 0
    }
    angles <- pi / 2 - 2 * pi * (seq_along(index) - 1) / length(index)
    coordinates[index, 1L] <- radius * cos(angles)
    coordinates[index, 2L] <- radius * sin(angles)
  }

  rownames(coordinates) <- vertex_attr(graph, "name")
  colnames(coordinates) <- c("x", "y")
  coordinates
}

.ggtcm_concentric_sizes <- function(graph, type_sizes, node_size) {
  type <- vertex_attr(graph, ".tcm_type")
  metric <- vertex_attr(graph, ".tcm_size")
  plot_size <- numeric(vcount(graph))

  for (node_type in c("Herb", "Molecule", "Target")) {
    index <- which(type == node_type)
    if (length(index) == 0L) {
      next
    }
    base <- type_sizes[[node_type]]
    values <- metric[index]
    if (is.null(node_size) || diff(range(values, finite = TRUE)) == 0) {
      plot_size[index] <- base
    } else {
      relative <- (values - min(values)) / diff(range(values))
      plot_size[index] <- base * (0.85 + 0.3 * relative)
    }
  }

  plot_size
}

.ggtcm_ring_layer <- function(ring_radii) {
  angle <- seq(0, 2 * pi, length.out = 361L)
  radii <- ring_radii[c("Molecule", "Target")]
  ring_data <- do.call(
    rbind,
    lapply(names(radii), function(node_type) {
      data.frame(
        x = radii[[node_type]] * cos(angle),
        y = radii[[node_type]] * sin(angle),
        ring = node_type
      )
    })
  )

  geom_path(
    data = ring_data,
    aes(x = .data[["x"]], y = .data[["y"]], group = .data[["ring"]]),
    inherit.aes = FALSE,
    colour = "#E6E6E6",
    linewidth = 0.25,
    linetype = "dashed"
  )
}

.ggtcm_validate_named_values <- function(values, required, argument) {
  if (is.null(names(values)) && length(values) == length(required)) {
    names(values) <- required
  }
  names(values)[names(values) == "Compound"] <- "Molecule"
  missing <- setdiff(required, names(values))
  if (length(missing) > 0L) {
    stop(
      sprintf("'%s' must provide values for: %s.", argument, paste(required, collapse = ", ")),
      call. = FALSE
    )
  }
  values[required]
}

.ggtcm_label_limits <- function(label_top, types) {
  if (!is.numeric(label_top) || length(label_top) == 0L ||
      any(is.na(label_top)) || any(label_top < 0)) {
    stop("'label_top' must contain non-negative numbers.", call. = FALSE)
  }
  if (length(label_top) == 1L) {
    return(setNames(rep(label_top, length(types)), types))
  }
  names(label_top)[names(label_top) == "Compound"] <- "Molecule"
  if (is.null(names(label_top)) || !all(types %in% names(label_top))) {
    stop(
      "'label_top' must be one number or values named Herb, Molecule and Target.",
      call. = FALSE
    )
  }
  label_top[types]
}

.ggtcm_label_data <- function(plot_data,
                              label_top,
                              size_range,
                              size_attribute = ".tcm_size",
                              size_identity = FALSE) {
  selected <- integer(0)
  for (node_type in names(label_top)) {
    index <- which(plot_data[[".tcm_type"]] == node_type)
    limit <- label_top[[node_type]]
    if (!is.infinite(limit)) {
      index <- head(
        index[order(
          plot_data[[".tcm_size"]][index],
          decreasing = TRUE,
          na.last = NA
        )],
        as.integer(limit)
      )
    }
    selected <- c(selected, index)
  }
  data <- plot_data[selected, , drop = FALSE]
  data[[".tcm_plot_label"]] <- data[[".tcm_label"]]
  data[[".tcm_point_size"]] <- if (isTRUE(size_identity)) {
    data[[size_attribute]]
  } else {
    .ggnetwork_rescale_size(data[[".tcm_size"]], size_range)
  }
  data
}
