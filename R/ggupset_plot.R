#' Plot an UpSet diagram directly from getvenndata() output with colored matrix points
#'
#' @param venn_df A logical matrix (output of getvenndata), first column is element ID, others are sets
#' @param sets Character vector: set names (columns). Default: all except first column.
#' @param palette RColorBrewer palette name (e.g., "Set2", "Dark2")
#' @param text.size Font size scaling
#' @param order.by Order type ("freq", "degree")
#' @param mb.ratio Top bar vs matrix height
#' @param keep.order Whether to preserve input set order
#'
#' @return A static UpSet plot
#' @export
#'
ggupset_plot <- function(venn_df,
                         sets = NULL,
                         palette = "Set2",
                         text.size = 3.5,
                         order.by = "freq",
                         mb.ratio = c(0.6, 0.4),
                         keep.order = FALSE) {
  
  if (is.null(sets)) {
    sets <- colnames(venn_df)[-1]
  }
  
  # convert logical matrix into 0/1 matrix for upset input
  matrix_data <- as.data.frame(lapply(venn_df[, sets, drop = FALSE], as.integer))
  
  # get unique intersection groups
  intersection_df <- unique(matrix_data)
  intersection_df <- intersection_df[rowSums(intersection_df) > 0, , drop = FALSE]
  n_combinations <- nrow(intersection_df)
  
  # color prepare
  color_vec <- RColorBrewer::brewer.pal(n = min(8, max(n_combinations, 3)), palette)
  if (n_combinations > length(color_vec)) {
    color_vec <- grDevices::colorRampPalette(color_vec)(n_combinations)
  }
  
  queries_list <- lapply(seq_len(n_combinations), function(i) {
    set_names <- names(intersection_df)[which(intersection_df[i, ] == 1)]
    list(
      query = UpSetR::intersects,
      params = list(set_names),
      color = color_vec[i],
      active = TRUE
    )
  })
  
  sets.bar.color <- color_vec[seq_len(length(sets))]
  
  p <- UpSetR::upset(
    matrix_data,
    nsets = length(sets),
    sets = if (keep.order) sets else NULL,
    sets.bar.color = sets.bar.color,
    text.scale = text.size / 2.5,
    mb.ratio = mb.ratio,
    order.by = order.by,
    keep.order = keep.order,
    decreasing = TRUE,
    queries = queries_list,
    shade.color = NA
  )
  
  return(p)
}
