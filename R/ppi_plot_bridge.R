.ppi_add_plot_attributes <- function(g,
                                     cluster_attr,
                                     candidate_attr = NULL,
                                     unassigned = NULL) {
  if (!cluster_attr %in% vertex_attr_names(g)) {
    stop(
      sprintf("Internal cluster attribute '%s' was not found.", cluster_attr),
      call. = FALSE
    )
  }

  cluster <- vertex_attr(g, cluster_attr)
  cluster <- as.character(cluster)
  if (!is.null(unassigned)) {
    cluster[cluster %in% as.character(unassigned)] <- NA_character_
  }
  cluster[!nzchar(cluster)] <- NA_character_
  g <- set_vertex_attr(g, "cluster", value = cluster)

  if (!is.null(candidate_attr)) {
    if (!candidate_attr %in% vertex_attr_names(g)) {
      stop(
        sprintf(
          "Internal candidate attribute '%s' was not found.",
          candidate_attr
        ),
        call. = FALSE
      )
    }
    g <- set_vertex_attr(
      g,
      "candidate",
      value = vertex_attr(g, candidate_attr)
    )
  }

  g
}
