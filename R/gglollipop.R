#' Lollipop Plot for Enrichment Results
#'
#' @param enrich_obj An enrichment result object from clusterProfiler.
#' @param x Character. The variable used for the x-axis. Default is `"RichFactor"`.
#' @param top_n Integer. Number of top enriched terms to display. Default is 10.
#' @param orderBy Character. Variable used to order the y-axis terms. Default is `"x"`.
#' @param text.col Character. Colors for text. Default is `black`.
#' @param text.size Numeric. Base font size for plot text. Default is 7.
#' @param font.family Character. Font family used throughout the plot.
#'   Default is `"Arial"`.
#' @param text.width Numeric. Font width for axis text and title. Default is 35.
#' @param palette Character. Color palette name from `RColorBrewer` to use for dot color. Default is `"RdBu"`.
#' @param line.col Character. Color of the segment lines. Default is `"grey60"`.
#' @param line.type Character. Line type for segments. Default is `"solid"`.
#' @param line.size Numeric. Line width for segments. Default is 0.9.
#' @param plot_title Character. Optional plot title. Default is `NULL`.
#' @param show_count Logical. Whether to display the count value as a text label next to each dot. Default is TRUE
#' @param ... Additional arguments passed to internal helper functions.

#' @return A `ggplot` object showing a lollipop-style enrichment plot.
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom enrichplot dotplot
#' @importFrom rlang .data
#'
#' @export

gglollipop <- function(enrich_obj,
                               x = "RichFactor",
                               top_n = 10,
                               orderBy = NULL,
                               text.col = "black",
                               text.size = 7,
                               font.family = "Arial",
                               text.width = 35,
                               palette = "RdBu",
                               line.col = "grey60",
                               line.type = "solid",
                               line.size = 0.9,
                               plot_title = NULL,
                               show_count = TRUE,
                               ...) {

  font.family <- .resolve_tcm_font_family(font.family)
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required for gglollipop(). Please install it with: install.packages('RColorBrewer')")
  }

  if (is.null(orderBy)) {
    orderBy <- "x"
  }

  p <- enrichplot::dotplot(enrich_obj,
                           showCategory = top_n,
                           x = x,
                           orderBy = orderBy,
                           ...)

  df_plot <- p$data
  n_color <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]

  ## remove former color_scale
  p$scales$scales <- Filter(function(s) {
    !"ScaleContinuous" %in% class(s) || s$aesthetics != "fill"
  }, p$scales$scales)

  p <- p +
    geom_segment(data = df_plot,
                 aes(x = 0, xend = .data[[x]],
                     y = .data[["Description"]], yend = .data[["Description"]]),
                 color = line.col,
                 linetype = line.type,
                 linewidth = line.size,
                 inherit.aes = FALSE) +
    scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(n_color, palette),
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE)
    )

  p$layers <- rev(p$layers)

  if (show_count){
    p <- p + geom_text_repel(data = df_plot,
              aes(x = .data[[x]],
                  y = .data[["Description"]],
              label = .data[["Count"]]),
              size = text.size / 2,
              color = text.col,
              family = font.family,
              fontface = "plain",
              hjust = 0,
              direction = "x",
              nudge_x = 0.025 * max(df_plot[[x]], na.rm = TRUE),
              segment.color = NA,
              inherit.aes = FALSE)
  }

  p <- p +
    .theme_tcm_pub(
      base_size = text.size,
      base_family = font.family,
      grid = "none"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "plain")
    ) +
    ggtitle(plot_title)

  return(p)
}
