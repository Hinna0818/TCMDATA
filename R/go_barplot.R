#' Bar Plot for GO Enrichment Results
#'
#' @param enrich_obj An enrichment result object from clusterProfiler::enrichGO(),
#'   or a data.frame containing at least three columns: ONTOLOGY, Description, and Count.
#' @param top_n Integer. Number of top enriched terms to display per ontology. Default is 10.
#' @param colors Named character vector. Colors for BP, CC, MF. Default uses "#E64B35", "#4DBBD5", and "#00A087" colors.
#' @param bar.width Numeric. Width of bars. Default is 0.7.
#' @param text.size Numeric. Font size for axis text. Default is 9.
#' @param label.size Numeric. Font size for count labels on bars. Default is 3.2.
#' @param show_count Logical. Whether to display the count value as a text label above each bar. Default is TRUE.
#' @param x.angle Numeric. Angle of x-axis text. Default is 60.
#' @param colored.label Logical. Whether to color x-axis labels by ontology.
#'   Default is FALSE. Note: When TRUE, this uses vectorized input to element_text()
#'   which may produce a warning in ggplot2 but works correctly.
#' @param legend.position Character. Position of legend. Default is "top".
#' @param plot_title Character. Optional plot title. Default is NULL.
#' @param y_expand Numeric. Expansion factor for y-axis upper limit. Default is 1.1.
#' @param ... Additional arguments (currently unused).
#'
#' @return A `ggplot` object showing a bar-style GO enrichment plot.
#'
#' @import ggplot2
#' @importFrom dplyr select mutate group_by slice_max ungroup arrange desc all_of %>%
#' @importFrom rlang .data
#' @importFrom scales pretty_breaks
#'
#' @examples
#' \dontrun{
#' ## demo data
#' herbs <- c("灵芝")
#' lz <- search_herb(herb = herbs, type = "Herb_cn_name")
#' set.seed(2025)
#' g <- sample(lz$target, 200)
#'
#' ## GO enrichment analysis
#' library(clusterProfiler)
#' x <- enrichGO(g, ont="all", OrgDb='org.Hs.eg.db', keyType="SYMBOL")
#' p1 <- go_barplot(x)
#' print(p1)
#'
#' ## Or use a custom data.frame
#' df <- data.frame(
#'   ONTOLOGY = c("BP", "BP", "CC", "MF"),
#'   Description = c("term1", "term2", "term3", "term4"),
#'   Count = c(10, 8, 15, 12)
#' )
#' p2 <- go_barplot(df)
#' print(p2)
#' }
#'
#' @export

go_barplot <- function(enrich_obj,
                       top_n = 10,
                       colors = NULL,
                       bar.width = 0.7,
                       text.size = 9,
                       label.size = 3.2,
                       show_count = TRUE,
                       x.angle = 60,
                       colored.label = FALSE,
                       legend.position = "top",
                       plot_title = NULL,
                       y_expand = 1.1,
                       ...) {

  # Check input data and extract result
  required_cols <- c("ONTOLOGY", "Description", "Count")

  if (is.data.frame(enrich_obj)) {
    # User provided a data.frame
    missing_cols <- setdiff(required_cols, colnames(enrich_obj))
    if (length(missing_cols) > 0) {
      stop("Input data.frame is missing required columns: ",
           paste(missing_cols, collapse = ", "),
           "\nRequired columns are: ONTOLOGY, Description, Count",
           call. = FALSE)
    }
    go_res <- enrich_obj
  } else if (inherits(enrich_obj, "enrichResult")) {
    # clusterProfiler enrichResult object
    go_res <- enrich_obj@result
  } else {
    stop("enrich_obj must be either a clusterProfiler enrichResult object ",
         "or a data.frame with columns: ONTOLOGY, Description, Count",
         call. = FALSE)
  }

  # Check for empty data
  if (nrow(go_res) == 0) {
    stop("Input data contains no rows.", call. = FALSE)
  }

  # Use default colors if necessary
  if (is.null(colors)) {
    colors <- c(
      "BP" = "#E64B35",
      "CC" = "#4DBBD5",
      "MF" = "#00A087"
    )
  }

  # Label mapping for legend
  label_mapping <- c(
    "BP" = "Biological Process",
    "CC" = "Cellular Component",
    "MF" = "Molecular Function"
  )

  # Process data
  plot_data <- go_res %>%
    select(all_of(required_cols)) %>%
    mutate(Count = as.numeric(.data$Count)) %>%
    group_by(.data$ONTOLOGY) %>%
    slice_max(order_by = .data$Count, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(.data$ONTOLOGY, desc(.data$Count)) %>%
    mutate(Description = factor(.data$Description, levels = .data$Description))

  # Check for empty plot data after filtering
  if (nrow(plot_data) == 0) {
    stop("No data remaining after filtering.", call. = FALSE)
  }

  # Build plot
  p <- ggplot(data = plot_data, aes(x = .data$Description, y = .data$Count, fill = .data$ONTOLOGY)) +
    geom_bar(stat = "identity", position = position_dodge(width = bar.width), width = bar.width)

  # Add count labels if requested
  if (show_count) {
    p <- p + geom_text(
      aes(label = .data$Count),
      position = position_dodge(width = bar.width),
      vjust = -0.5,
      size = label.size,
      fontface = "bold"
    )
  }

  # Determine x-axis label color
  if (colored.label) {
    x_label_color <- colors[plot_data$ONTOLOGY]
  } else {
    x_label_color <- "black"
  }

  # Apply scales and theme
  p <- p +
    scale_fill_manual(values = colors, labels = label_mapping) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(
        angle = x.angle,
        hjust = 1,
        vjust = 1,
        size = text.size,
        face = "bold"
      ),
      axis.ticks.x = element_line(color = "black"),
      axis.title.y = element_text(
        margin = margin(r = 10),
        face = "bold",
        size = text.size + 2
      ),
      axis.line.y = element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = element_line(color = "black"),
      axis.text.y = element_text(margin = margin(r = 5), size = text.size + 1),
      panel.grid = element_blank(),
      legend.position = legend.position,
      legend.title = element_blank(),
      legend.key.size = unit(5, "mm"),
      legend.text = element_text(size = text.size + 1),
      legend.justification = "center",
      legend.spacing.x = unit(8, "mm"),
      plot.margin = margin(10, 15, 10, 10, unit = "mm"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = text.size + 3, face = "bold", hjust = 0.5)
    ) +
    labs(x = "", y = "Number of Genes") +
    scale_y_continuous(
      limits = c(0, max(plot_data$Count) * y_expand),
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0, 0))
    )

  # Apply colored x-axis labels using ggplot_build (CRAN-compliant approach)
  if (colored.label) {
    p <- p + theme(axis.text.x = element_text(
      angle = x.angle,
      hjust = 1,
      vjust = 1,
      size = text.size,
      face = "bold",
      color = x_label_color
    ))
  }

  # Add title if provided
  if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title)
  }

  return(p)
}
