#' @param enrich_obj An enrichment result object of class `enrichResult` from clusterProfiler 
#' or a `data.frame` containing at least the columns: `Description`, `RichFactor`, `p.adjust`, and `Count`.
#' @param top_n Integer. Number of top enriched terms to display. Default is 10.
#' @param text.size Numeric. Font size of axis text and title. Default is 10.
#' @param text.width Integer. Line wrapping width for y-axis labels. Default is 35 characters.
#' @param color Character. Color palette name from RColorBrewer to use for bubble color. Default is "RdBu".
#' @param plot_title Character. Optional title for the plot. Default is NULL (no title).
#'
#' @return A `ggplot` object representing the lollipop plot.
#'
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stringr str_wrap
#' @importFrom dplyr arrange slice_head %>% 
#' @importFrom tidyr drop_na
#' 
#' @export
#' 
gglollipop <- function(
    enrich_obj,
    top_n = 10,
    text.size = 10,
    text.width = 35,
    color = "RdBu",
    plot_title = NULL) 
  {
  
  library(ggplot2)
  library(RColorBrewer)
  library(stringr)
  library(dplyr)
  
  if (inherits(enrich_obj, "enrichResult")) {
    df <- enrich_obj@result %>% tidyr::drop_na()
  } else if (is.data.frame(enrich_obj)) {
    df <- enrich_obj %>% tidyr::drop_na()
  } else {
    stop("Input must be an enrichResult object or a data.frame.")
  }
  
  required_cols <- c("Description", "RichFactor", "p.adjust", "Count")
  if (!all(required_cols %in% colnames(df))) {
    stop(paste("Data must contain the following columns:", paste(required_cols, collapse = ", ")))
  }
  
  df <- df %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n)
  
  df$RichFactor <- df$RichFactor %>% round(digits = 2)
  
  label_wrap <- function(labels) str_wrap(labels, width = text.width)
  
  p <- ggplot(
    df,
    aes(RichFactor, stats::reorder(Description, RichFactor))
  ) +
    geom_segment(aes(xend = 0, yend = Description)) +
    geom_point(aes(color = p.adjust, size = Count)) +
    geom_text(
      aes(label = Count, x = RichFactor + max(RichFactor) * 0.03),  
      size = text.size / 2,
      color = "black"
    ) + 
    scale_color_gradientn(
      colours = RColorBrewer::brewer.pal(8, color),
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    scale_size_continuous(range = c(2, 10)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = text.size, colour = "black", vjust = 1),
      axis.text.y = element_text(size = text.size, colour = "black", hjust = 1),
      axis.title = element_text(margin = margin(10, 5, 0, 0), color = "black", size = text.size),
      plot.title = element_text(size = text.size + 2, face = "bold", hjust = 0.5)
    ) +
    xlab("Rich Factor") +
    ylab(NULL) +
    ggtitle(plot_title) +
    scale_y_discrete(labels = label_wrap)
  
  return(p)
}
