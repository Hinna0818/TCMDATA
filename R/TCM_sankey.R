#' Create a Sankey Diagram for Hierarchical TCM or Omics Data
#'
#' @param data A `data.frame` for sankey_plot.
#' @param text.size Numeric. Size of node label text. Default is 3.
#' @param text.position Numeric. Horizontal justification of node labels; typically 0 (left),
#'                      0.5 (center), or 1 (right). Default is 1.
#' @param x.axis.text.size Numeric. Font size of x-axis labels. Default is 12.
#' @param flow.alpha Numeric (0–1). Transparency of the Sankey flows. Default is 0.5.
#' @param palette Character. Name of the color palette to use from `cols4all::c4a_pals()`.
#'                Default is `"rainbow_wh_rd"`.
#'
#' @return A `ggplot` object.
#'
#' @import ggplot2
#' @importFrom cols4all c4a
#' @importFrom ggsankey make_long geom_sankey geom_sankey_text theme_sankey
#'
#' @export
TCM_sankey <- function(data,
                       text.size = 3,
                       text.position = 1,
                       x.axis.text.size = 12,
                       flow.alpha = 0.5,
                       palette = "rainbow_wh_rd",
                       ...) {
  
  
  stopifnot(is.data.frame(data), ncol(data) >= 2)
  
  layer_cols <- colnames(data)
  
  df_long <- ggsankey::make_long(data, !!!rlang::syms(layer_cols))
  
  node_colors <- cols4all::c4a(palette, length(unique(df_long$node)))
  names(node_colors) <- unique(df_long$node)
  
  p <- ggplot(df_long,
              aes(x = .data[["x"]],
                  next_x = .data[["next_x"]],
                  node = .data[["node"]],
                  next_node = .data[["next_node"]],
                  fill = .data[["node"]])) +
    ggsankey::geom_sankey(flow.alpha = flow.alpha, node.color = NA, show.legend = FALSE) +
    ggsankey::geom_sankey_text(
                     aes(label = .data[["node"]]),
                     size = text.size,
                     hjust = text.position,
                     color = "black",
                     position = position_nudge(x = -0.05)) +
    scale_fill_manual(values = node_colors) +
    ggsankey::theme_sankey(base_size = 16) +
    theme(axis.title = element_blank(),
          axis.text.x = element_text(size = x.axis.text.size,
                                     color = "black",
                                     hjust = 0.5))
  
  return(p)
  
}
