#' Herb-target ORA (single herb) with TCM targets
#'
#' @param genes Character vector of query genes.
#' @param herb_targets Character vector of targets for a single herb.
#' @param universe Optional. Either an integer background size or a character vector
#'   of background genes. If `NULL`, it defaults to the union of all targets in `tcmdata$target`.
#'
#' @importFrom stats phyper
#' @return List with fields: ngene, setsize, indrg, pvalue, overlap
herb_ora_one <- function(genes, herb_targets, universe = NULL) {
  
  res <- list(ngene = 0L, setsize = 0L, indrg = 0L, pvalue = NA_real_, overlap = NA_character_)
  
  # sanitize inputs
  genes <- unique(na.omit(as.character(genes)))
  herb_targets <- unique(na.omit(as.character(herb_targets)))
  
  # default universe: all TCM targets
  if (is.null(universe)) {
    universe <- unique(na.omit(as.character(tcm_data$target)))
  }
  
  # universe handling: integer size or explicit background vector
  if (is.numeric(universe) && length(universe) == 1L) {
    U <- as.integer(universe)
    genes_in <- genes
    set_in   <- herb_targets
  } else {
    bg <- unique(na.omit(as.character(universe)))
    genes_in <- intersect(genes, bg)
    set_in   <- intersect(herb_targets, bg)
    U <- length(bg)
  }
  
  K <- length(genes_in)                       # query size
  M <- length(set_in)                         # set size
  k <- length(intersect(genes_in, set_in))    # overlap
  
  res$ngene   <- K
  res$setsize <- M
  res$indrg   <- k
  
  if (U <= 0 || K <= 0 || M <= 0 || k <= 0) {
    return(res)
  }
  
  q <- k - 1
  m <- M
  n <- U - M
  res$pvalue <- stats::phyper(q, m, n, K, lower.tail = FALSE)
  res$overlap <- paste0(sort(intersect(genes_in, set_in)), collapse = "/")
  return(res)
}


#' Herb-target ORA (multiple herbs) with TCM targets
#'
#' @param genes Character vector of query genes.
#' @param herb_sets A vector of herbs' names.
#' @param type A vector of search types of herb.
#' @param universe Optional background size or background gene vector.
#'   If `NULL`, defaults to the union of `tcmdata$target`.
#' @param p_adjust_method Multiple-testing correction method (default "BH").
#' 
#' @importFrom stats p.adjust
#' @importFrom tibble as_tibble
#' @return A tibble with columns: herb, ngene, setsize, indrg, pvalue, padj, overlap
herb_ora <- function(genes, 
                     herb_sets = NULL, 
                     type = c("Herb_cn_name","Herb_pinyin_name","Herb_en_name"),
                     universe = NULL, 
                     p_adjust_method = "BH") {
  
  type <- match.arg(type)
  herb_sets <- get_herbsets(herb_name = herb_sets, type = type)
  
  # default universe: all TCM targets
  if (is.null(universe)) {
    universe <- unique(na.omit(as.character(tcm_data$target)))
  }
  
  out <- lapply(names(herb_sets), function(h) {
    r <- herb_ora_one(genes, herb_sets[[h]], universe = universe)
    data.frame(
      herb    = h,
      ngene   = r$ngene,
      setsize = r$setsize,
      indrg   = r$indrg,
      pvalue  = r$pvalue,
      overlap = r$overlap,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, out)
  
  # adjust p only where finite
  p_ok <- is.finite(out$pvalue)
  out$padj <- NA_real_
  if (any(p_ok)) {
    out$padj[p_ok] <- p.adjust(out$pvalue[p_ok], method = p_adjust_method)
  }
  
  # order by padj then pvalue
  ord <- order(ifelse(is.na(out$padj), Inf, out$padj),
               ifelse(is.na(out$pvalue), Inf, out$pvalue))
  out <- out[ord, , drop = FALSE]
  
  res <- tibble::as_tibble(out)
  return(res)
}

#' Build herb target sets as the input for herb_ora analysis
#'
#' @param herb_name Character vector of herb names to query. If NULL, use all herbs in tcm_data.
#' @param type One of "Herb_cn_name", "Herb_pinyin_name", or "Herb_en_name".
#' @return A named list: names = herb (pinyin), values = target vectors.
#' @keywords internal
get_herbsets <- function(herb_name = NULL, type = c("Herb_cn_name","Herb_pinyin_name","Herb_en_name")){
  type <- match.arg(type)

  if (is.null(herb_name)) {
    df <- tcm_data[, c("Herb_pinyin_name", "target")]
    names(df) <- c("herb", "target")
  } else {
    df <- search_herb(herb = herb_name, type = type)[, c("herb", "target")]
  }
  
  df <- df[!is.na(df$target) & !is.na(df$herb), , drop = FALSE]
  df <- unique(df)
  
  herb_sets <- split(df$target, df$herb)
  herb_sets <- lapply(herb_sets, unique)
  
  return(herb_sets)
}


#' Get top enriched herbs
#' @param ora_res A tibble/data.frame from herb_ora(), must contain columns: herb, pvalue, padj.
#' @param top Integer, number of top herbs to keep after sorting by FDR. Default 10.
#' @importFrom tibble as_tibble
#' @keywords internal
get_top_herb <- function(ora_res, top = 10){
  if (is.null(ora_res)) stop("`ora_res` cannot be NULL.")
  res <- tibble::as_tibble(ora_res)
  
  # compute -log10(FDR)
  if (!all(c("herb","pvalue","padj") %in% names(res))) {
    stop("`ora_res` must contain columns: herb, pvalue, padj.")
  }
  res$mlp <- -log10(res$padj)
  res <- res[is.finite(res$mlp), , drop = FALSE]
  if (nrow(res) == 0) return(tibble::tibble())
  
  # sort by FDR then p-value
  res <- res[order(res$padj, res$pvalue), , drop = FALSE]
  
  # keep top N
  if (!is.null(top) && top > 0L && nrow(res) > top) {
    res <- res[seq_len(top), , drop = FALSE]
  }
  
  return(res)
}


#' Bar plot for herb ORA results
#'
#' @param ora_res A tibble containing herb-ora results.
#' @param top Integer. Number of herbs to plot (after sorting by FDR). Default 10.
#' @param base_size Base font size.
#' @param show_overlap Logical, whether to annotate bars with overlap counts.
#' @param ... Additional arguments passed to internal helper functions.
#' 
#' @import ggplot2
#' @importFrom forcats fct_reorder
#' @importFrom rlang .data
#' @return A ggplot object.
#' @export
herb_ora_plot <- function(ora_res, 
                          top = 10,
                          base_size = 12, 
                          show_overlap = TRUE,
                          ...){
  
  ## extract top herbs
  df <- get_top_herb(ora_res = ora_res, top = top)
  if (nrow(df) == 0) {
    warning("No finite FDR values to plot.")
    return(invisible(NULL))
  }
  
  ## calculate mlp
  if (!"mlp" %in% names(df)) df$mlp <- -log10(df$padj)
  df <- df[is.finite(df$mlp), , drop = FALSE]
  df$herb <- forcats::fct_reorder(df$herb, df$mlp)
  
  p <- ggplot(data = df, aes(x = .data[["herb"]], y = .data[["mlp"]])) +
    geom_col() +
    coord_flip() +
    labs(x = NULL, y = expression(-log[10](FDR))) +
    theme_minimal(base_size = base_size)
  
  if (isTRUE(show_overlap) && all(c("indrg") %in% names(df))) {
    p <- p + geom_text(
      aes(label = .data[["indrg"]]),
      hjust = -0.1, size = 3
    ) +
      expand_limits(y = max(df$mlp, na.rm = TRUE) * 1.08)
  }
  
  return(p)
}