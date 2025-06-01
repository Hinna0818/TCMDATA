#' Draw an Alluvial Diagram for Hierarchical Data
#'
#' @param data A `data.frame` with at least two columns representing hierarchical levels.
#' @param text.size Numeric. Font size of the node labels. Default is 3.
#' @param text.position Numeric. Horizontal justification of the node labels. Typically 0 (left), 0.5 (center), or 1 (right). Default is 1.
#' @param x.axis.text.size Numeric. Font size of the x-axis text labels. Default is 12.
#' @param flow.alpha Numeric (between 0 and 1). Transparency of the flow paths. Default is 0.5.
#' @param palette Character. Name of the color palette to use from `cols4all::c4a`. Default is `"rainbow_wh_rd"`.
#'
#' @return A `ggplot` object representing the alluvial diagram.
#'
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom cols4all c4a
#' @importFrom ggsankey make_long geom_alluvial geom_alluvial_text theme_sankey
#'
#' @export
TCM_alluvial <- function(data,
                       text.size = 3,
                       text.position = 1,
                       x.axis.text.size = 12,
                       flow.alpha = 0.5,
                       palette = "rainbow_wh_rd",
                       ...) {
  
  library(ggplot2)
  library(cols4all)
  library(dplyr)
  
  stopifnot(is.data.frame(data), ncol(data) >= 2)
  
  layer_cols <- colnames(data)
  
  df_long <- ggsankey::make_long(data, !!!rlang::syms(layer_cols))
  
  node_colors <- cols4all::c4a(palette, length(unique(df_long$node)))
  names(node_colors) <- unique(df_long$node)
  
  p <- ggplot(df_long,
              aes(x = x,
                  next_x = next_x,
                  node = node,
                  next_node = next_node,
                  fill = node,
                  label = node)) +
    ggsankey::geom_alluvial(flow.alpha = flow.alpha, node.color = NA, show.legend = FALSE) +
    ggsankey::geom_alluvial_text(size = text.size, 
                                 color = "black", 
                                 hjust = text.position,
                                 position = position_nudge(x = -0.05)) +
    ggsankey::theme_sankey(base_size = 16) +
    scale_fill_manual(values = node_colors) +
    theme(axis.title = element_blank()) +
    theme(axis.text.x = element_text(
      size = x.axis.text.size,
      hjust = 0.5, vjust = 10,
      colour = "black"
    ))
  
  return(p)
}
