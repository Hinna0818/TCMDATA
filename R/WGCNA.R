#' Run weighted gene co-expression network analysis
#'
#' A lightweight wrapper around [WGCNA::pickSoftThreshold()] and
#' [WGCNA::blockwiseModules()]. Expression data use the same orientation as
#' other TCMDATA workflows: genes in rows and samples in columns.
#'
#' @param expr_mat Numeric expression matrix with genes in rows and samples in
#'   columns.
#' @param traits Optional sample-level vector, factor, matrix, or data frame.
#'   Rows should correspond to samples; row names are used for matching when
#'   available. Categorical variables are converted to indicator variables.
#' @param power Positive soft-thresholding power. If `NULL`, it is selected with
#'   [WGCNA::pickSoftThreshold()].
#' @param powers Candidate powers used when `power = NULL`.
#' @param scale_free_R2 Target scale-free topology fit used for power selection.
#' @param top_n Optional number of most variable genes to retain before network
#'   construction.
#' @param network_type Network type, either `"unsigned"` or `"signed"`.
#' @param tom_type Topological-overlap type. Defaults to `network_type`.
#' @param cor_type Correlation method, either `"pearson"` or `"bicor"`.
#' @param max_block_size Maximum number of genes per block.
#' @param min_module_size Minimum module size.
#' @param merge_cut_height Height used to merge similar modules.
#' @param deep_split Module-detection sensitivity passed to
#'   [WGCNA::blockwiseModules()].
#' @param reassign_threshold P-value threshold for reassigning genes between
#'   modules.
#' @param pam_respects_dendro Whether PAM assignments must respect the
#'   dendrogram.
#' @param seed Random seed used by blockwise module detection.
#' @param n_threads Number of WGCNA worker threads. Use `1` for deterministic
#'   serial execution.
#' @param verbose WGCNA verbosity level.
#'
#' @return A `tcm_wgcna` object containing the selected power, network,
#'   module assignments, eigengenes, gene-level module membership, module gene
#'   lists, and optional module-trait and gene-trait correlations.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' expr <- matrix(rnorm(80 * 20), nrow = 80)
#' rownames(expr) <- paste0("Gene", seq_len(nrow(expr)))
#' colnames(expr) <- paste0("Sample", seq_len(ncol(expr)))
#' group <- factor(rep(c("Control", "Disease"), each = 10))
#'
#' fit <- run_wgcna(
#'   expr,
#'   traits = group,
#'   power = 6,
#'   min_module_size = 10
#' )
#' fit$gene_info
#' }
#'
#' @export
run_wgcna <- function(
    expr_mat,
    traits = NULL,
    power = NULL,
    powers = c(1:10, seq(12, 20, by = 2)),
    scale_free_R2 = 0.85,
    top_n = NULL,
    network_type = c("unsigned", "signed"),
    tom_type = NULL,
    cor_type = c("pearson", "bicor"),
    max_block_size = 5000,
    min_module_size = 30,
    merge_cut_height = 0.25,
    deep_split = 2,
    reassign_threshold = 0,
    pam_respects_dendro = FALSE,
    seed = 2025,
    n_threads = 1,
    verbose = 1) {
  .check_wgcna_dependency()

  network_type <- match.arg(network_type)
  cor_type <- match.arg(cor_type)
  if (is.null(tom_type)) {
    tom_type <- network_type
  }
  tom_type <- match.arg(tom_type, c("unsigned", "signed"))
  .validate_wgcna_parameters(
    power = power,
    powers = powers,
    scale_free_R2 = scale_free_R2,
    top_n = top_n,
    max_block_size = max_block_size,
    min_module_size = min_module_size,
    merge_cut_height = merge_cut_height,
    deep_split = deep_split,
    n_threads = n_threads
  )

  prepared <- .prepare_wgcna_expression(expr_mat, top_n = top_n)
  dat_expr <- prepared$data
  trait_info <- .prepare_wgcna_traits(traits, rownames(dat_expr))

  # Several WGCNA internals resolve correlation functions in the caller.
  # Keep these local bindings so namespace-only use does not fall back to
  # stats::cor, which lacks WGCNA's weight and cosine arguments.
  cor <- WGCNA::cor
  bicor <- WGCNA::bicor

  power_fit <- .select_wgcna_power(
    dat_expr = dat_expr,
    power = power,
    powers = powers,
    scale_free_R2 = scale_free_R2,
    network_type = network_type,
    cor_type = cor_type,
    verbose = verbose
  )

  if (verbose > 0) {
    message(
      "Constructing ", network_type, " network with power ",
      power_fit$power, " (", nrow(dat_expr), " samples; ",
      ncol(dat_expr), " genes)."
    )
  }

  set.seed(seed)
  network <- WGCNA::blockwiseModules(
    datExpr = dat_expr,
    power = power_fit$power,
    maxBlockSize = as.integer(max_block_size),
    randomSeed = as.integer(seed),
    corType = cor_type,
    networkType = network_type,
    TOMType = tom_type,
    deepSplit = deep_split,
    minModuleSize = as.integer(min_module_size),
    reassignThreshold = reassign_threshold,
    mergeCutHeight = merge_cut_height,
    numericLabels = FALSE,
    pamRespectsDendro = pam_respects_dendro,
    saveTOMs = FALSE,
    nThreads = as.integer(n_threads),
    verbose = verbose
  )

  module_colors <- stats::setNames(
    as.character(network$colors),
    colnames(dat_expr)
  )
  eigengenes <- WGCNA::moduleEigengenes(
    dat_expr,
    colors = module_colors
  )$eigengenes
  eigengenes <- WGCNA::orderMEs(eigengenes)

  kme <- WGCNA::signedKME(
    dat_expr,
    eigengenes,
    outputColumnName = "MM.",
    corFnc = if (cor_type == "pearson") {
      "WGCNA::cor"
    } else {
      "WGCNA::bicor"
    },
    corOptions = "use = 'p'"
  )
  kme <- as.data.frame(kme, check.names = FALSE)

  own_kme_column <- paste0("MM.", module_colors)
  own_kme <- vapply(
    seq_along(module_colors),
    function(i) {
      column <- own_kme_column[[i]]
      if (column %in% colnames(kme)) kme[i, column] else NA_real_
    },
    numeric(1)
  )
  gene_info <- data.frame(
    gene = names(module_colors),
    module = unname(module_colors),
    module_membership = own_kme,
    stringsAsFactors = FALSE,
    row.names = names(module_colors)
  )
  module_genes <- split(gene_info$gene, gene_info$module)

  module_trait_cor <- NULL
  module_trait_p <- NULL
  gene_trait_cor <- NULL
  gene_trait_p <- NULL
  if (!is.null(trait_info$design)) {
    module_trait_cor <- .wgcna_cor(
      eigengenes,
      trait_info$design,
      cor_type = cor_type
    )
    module_trait_p <- WGCNA::corPvalueStudent(
      module_trait_cor,
      nSamples = nrow(dat_expr)
    )
    gene_trait_cor <- .wgcna_cor(
      dat_expr,
      trait_info$design,
      cor_type = cor_type
    )
    gene_trait_p <- WGCNA::corPvalueStudent(
      gene_trait_cor,
      nSamples = nrow(dat_expr)
    )
  }

  out <- list(
    power = power_fit$power,
    power_source = power_fit$source,
    power_diagnostics = power_fit$diagnostics,
    expression = dat_expr,
    traits = trait_info$design,
    trait_metadata = trait_info$metadata,
    network = network,
    module_colors = module_colors,
    module_eigengenes = eigengenes,
    module_membership = kme,
    gene_info = gene_info,
    module_genes = module_genes,
    module_trait_cor = module_trait_cor,
    module_trait_p = module_trait_p,
    gene_trait_cor = gene_trait_cor,
    gene_trait_p = gene_trait_p,
    qc = prepared$qc,
    params = list(
      network_type = network_type,
      tom_type = tom_type,
      cor_type = cor_type,
      max_block_size = max_block_size,
      min_module_size = min_module_size,
      merge_cut_height = merge_cut_height,
      deep_split = deep_split,
      reassign_threshold = reassign_threshold,
      pam_respects_dendro = pam_respects_dendro,
      seed = seed,
      n_threads = n_threads
    )
  )
  class(out) <- c("tcm_wgcna", "list")
  out
}

#' Extract genes from a WGCNA module
#'
#' @param x A `tcm_wgcna` object returned by `run_wgcna()`.
#' @param module Module colour or label.
#' @param min_module_membership Optional minimum absolute module membership.
#' @param trait Optional trait name from `colnames(x$gene_trait_cor)`.
#' @param min_gene_significance Optional minimum absolute gene-trait
#'   correlation. Used only when `trait` is supplied.
#'
#' @return A character vector of gene names.
#'
#' @export
get_wgcna_module_genes <- function(
    x,
    module,
    min_module_membership = NULL,
    trait = NULL,
    min_gene_significance = NULL) {
  if (!inherits(x, "tcm_wgcna")) {
    stop("'x' must be a tcm_wgcna object.", call. = FALSE)
  }
  if (!is.character(module) || length(module) != 1L || is.na(module)) {
    stop("'module' must be one module name.", call. = FALSE)
  }

  info <- x$gene_info[x$gene_info$module == module, , drop = FALSE]
  if (nrow(info) == 0L) {
    stop("Module '", module, "' was not found.", call. = FALSE)
  }
  if (!is.null(min_module_membership)) {
    .check_unit_threshold(
      min_module_membership,
      "min_module_membership"
    )
    info <- info[
      abs(info$module_membership) >= min_module_membership,
      ,
      drop = FALSE
    ]
  }

  if (!is.null(trait)) {
    if (is.null(x$gene_trait_cor) ||
        !trait %in% colnames(x$gene_trait_cor)) {
      stop("Trait '", trait, "' was not found.", call. = FALSE)
    }
    if (is.null(min_gene_significance)) {
      min_gene_significance <- 0
    }
    .check_unit_threshold(
      min_gene_significance,
      "min_gene_significance"
    )
    significance <- x$gene_trait_cor[info$gene, trait]
    info <- info[
      is.finite(significance) &
        abs(significance) >= min_gene_significance,
      ,
      drop = FALSE
    ]
  } else if (!is.null(min_gene_significance)) {
    stop(
      "'trait' is required when 'min_gene_significance' is supplied.",
      call. = FALSE
    )
  }

  unname(info$gene)
}

#' @export
print.tcm_wgcna <- function(x, ...) {
  module_sizes <- sort(table(x$module_colors), decreasing = TRUE)
  cat("<tcm_wgcna>\n")
  cat("  Samples:", nrow(x$expression), "\n")
  cat("  Genes:", ncol(x$expression), "\n")
  cat("  Power:", x$power, paste0(" (", x$power_source, ")"), "\n")
  cat("  Modules:", length(module_sizes), "\n")
  cat(
    "  Module sizes:",
    paste(names(module_sizes), module_sizes, sep = "=", collapse = ", "),
    "\n"
  )
  invisible(x)
}

#' Plot WGCNA modules
#'
#' @param x A `tcm_wgcna` object.
#' @param block Block number to plot.
#' @param linewidth Dendrogram line width.
#' @param base_size,base_family Text size and font family.
#'
#' @return A ggplot-compatible object.
#'
#' @importFrom aplot insert_top
#' @importFrom ggfun theme_noxaxis
#' @importFrom ggplot2 aes geom_tile ggplot labs margin scale_fill_identity theme theme_void
#' @importFrom ggtree ggtree
#' @export
plot_wgcna_modules <- function(
    x,
    block = 1L,
    linewidth = 0.25,
    base_size = 7,
    base_family = "Arial") {
  .check_wgcna_plot_input(x)
  block <- as.integer(block)
  if (block < 1L || block > length(x$network$dendrograms)) {
    stop("'block' is outside the available range.", call. = FALSE)
  }

  genes <- x$network$blockGenes[[block]]
  modules <- unname(x$module_colors[genes])
  strip_data <- data.frame(
    label = seq_along(genes),
    module = modules
  )
  family <- .resolve_tcm_font_family(base_family)

  tree <- suppressWarnings(
    ggtree(
      x$network$dendrograms[[block]],
      layout = "dendrogram",
      ladderize = FALSE,
      linewidth = linewidth
    ) +
      labs(x = "Height") +
      theme_noxaxis() +
      .theme_tcm_text(base_size, family)
  )

  strip <- ggplot(
    strip_data,
    aes(x = .data$label, y = "Module", fill = .data$module)
  ) +
    geom_tile() +
    scale_fill_identity() +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = base_size, base_family = family) +
    theme(plot.margin = margin(0, 5.5, 0, 5.5))

  insert_top(strip, tree, height = 8)
}

#' Plot WGCNA eigengene relationships
#'
#' @param x A `tcm_wgcna` object.
#' @param show_heatmap Show a correlation heatmap instead of the dendrogram.
#' @param linewidth Tree line width.
#' @param label_size Tip-label size.
#' @param base_size,base_family Text size and font family.
#'
#' @return A ggplot-compatible object.
#'
#' @importFrom ggfun theme_noxaxis
#' @importFrom ggplot2 aes coord_equal element_blank element_text geom_tile ggplot labs scale_fill_gradient2 theme theme_minimal
#' @importFrom ggtree geom_tiplab ggtree
#' @export
plot_wgcna_eigengenes <- function(
    x,
    show_heatmap = FALSE,
    linewidth = 0.4,
    label_size = 2.5,
    base_size = 7,
    base_family = "Arial") {
  .check_wgcna_plot_input(x)
  eigengenes <- x$module_eigengenes
  correlation <- stats::cor(eigengenes, use = "pairwise.complete.obs")
  clustering <- stats::hclust(
    stats::as.dist(1 - correlation),
    method = "average"
  )
  family <- .resolve_tcm_font_family(base_family)

  if (show_heatmap) {
    module_order <- clustering$labels[clustering$order]
    heat_data <- as.data.frame(as.table(correlation))
    names(heat_data) <- c("module_x", "module_y", "correlation")
    heat_data$module_x <- factor(heat_data$module_x, levels = module_order)
    heat_data$module_y <- factor(
      heat_data$module_y,
      levels = rev(module_order)
    )

    return(
      ggplot(
        heat_data,
        aes(
          x = .data$module_x,
          y = .data$module_y,
          fill = .data$correlation
        )
      ) +
        geom_tile(colour = "white", linewidth = 0.25) +
        scale_fill_gradient2(
          low = "#4E79A7",
          mid = "white",
          high = "#B85C5C",
          midpoint = 0,
          limits = c(-1, 1),
          name = "Correlation"
        ) +
        coord_equal() +
        labs(x = NULL, y = NULL) +
        theme_minimal(base_size = base_size, base_family = family) +
        .theme_tcm_text(base_size, family) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          )
        )
    )
  }

  suppressWarnings(
    ggtree(
      clustering,
      layout = "dendrogram",
      ladderize = FALSE,
      linewidth = linewidth
    ) +
      geom_tiplab(
        size = label_size,
        family = family,
        fontface = "plain"
      ) +
      labs(x = "Height") +
      theme_noxaxis() +
      .theme_tcm_text(base_size, family)
  )
}

#' Plot WGCNA module-trait associations
#'
#' @param x A `tcm_wgcna` object containing trait associations.
#' @param show_labels Show correlations and P values in cells.
#' @param digits Number of correlation digits.
#' @param label_size Label size.
#' @param base_size,base_family Text size and font family.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 aes element_blank element_text geom_text geom_tile ggplot labs scale_fill_gradient2 theme theme_minimal
#' @export
plot_wgcna_traits <- function(
    x,
    show_labels = TRUE,
    digits = 2,
    label_size = 2.4,
    base_size = 7,
    base_family = "Arial") {
  .check_wgcna_plot_input(x)
  if (is.null(x$module_trait_cor) || is.null(x$module_trait_p)) {
    stop("No trait associations are available in 'x'.", call. = FALSE)
  }

  correlation <- x$module_trait_cor
  p_value <- x$module_trait_p
  plot_data <- expand.grid(
    module = rownames(correlation),
    trait = colnames(correlation),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  plot_data$correlation <- as.vector(correlation)
  plot_data$p_value <- as.vector(p_value)
  plot_data$label <- paste0(
    formatC(
      plot_data$correlation,
      format = "f",
      digits = as.integer(digits)
    ),
    "\n(",
    format.pval(plot_data$p_value, digits = 2, eps = 0.001),
    ")"
  )
  plot_data$module <- factor(
    plot_data$module,
    levels = rev(rownames(correlation))
  )
  plot_data$trait <- factor(
    plot_data$trait,
    levels = colnames(correlation)
  )
  family <- .resolve_tcm_font_family(base_family)

  plot <- ggplot(
    plot_data,
    aes(
      x = .data$trait,
      y = .data$module,
      fill = .data$correlation
    )
  ) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "deepskyblue",
      mid = "white",
      high = "orangered",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Correlation"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = base_size, base_family = family) +
    .theme_tcm_text(base_size, family) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )
  )

  if (show_labels) {
    plot <- plot + geom_text(
      aes(label = .data$label),
      size = label_size,
      family = family,
      fontface = "plain"
    )
  }
  plot
}

.check_wgcna_plot_input <- function(x) {
  if (!inherits(x, "tcm_wgcna")) {
    stop("'x' must be a tcm_wgcna object.", call. = FALSE)
  }
}

.check_wgcna_dependency <- function() {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop(
      "Package 'WGCNA' is required. Install it with ",
      "install.packages('WGCNA').",
      call. = FALSE
    )
  }
}

.validate_wgcna_parameters <- function(
    power,
    powers,
    scale_free_R2,
    top_n,
    max_block_size,
    min_module_size,
    merge_cut_height,
    deep_split,
    n_threads) {
  if (!is.null(power) &&
      (!is.numeric(power) || length(power) != 1L ||
       !is.finite(power) || power <= 0)) {
    stop("'power' must be NULL or one positive number.", call. = FALSE)
  }
  if (!is.numeric(powers) || length(powers) == 0L ||
      any(!is.finite(powers)) || any(powers <= 0)) {
    stop("'powers' must contain positive finite numbers.", call. = FALSE)
  }
  .check_unit_threshold(scale_free_R2, "scale_free_R2")
  if (!is.null(top_n) &&
      (!is.numeric(top_n) || length(top_n) != 1L ||
       !is.finite(top_n) || top_n < 2)) {
    stop("'top_n' must be NULL or an integer of at least 2.", call. = FALSE)
  }
  positive_integer <- list(
    max_block_size = max_block_size,
    min_module_size = min_module_size,
    n_threads = n_threads
  )
  for (name in names(positive_integer)) {
    value <- positive_integer[[name]]
    if (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 1) {
      stop("'", name, "' must be a positive integer.", call. = FALSE)
    }
  }
  .check_unit_threshold(merge_cut_height, "merge_cut_height")
  if (!is.numeric(deep_split) || length(deep_split) != 1L ||
      !is.finite(deep_split) || deep_split < 0 || deep_split > 4) {
    stop("'deep_split' must be between 0 and 4.", call. = FALSE)
  }
}

.check_unit_threshold <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L ||
      !is.finite(value) || value < 0 || value > 1) {
    stop("'", name, "' must be between 0 and 1.", call. = FALSE)
  }
}

.prepare_wgcna_expression <- function(expr_mat, top_n = NULL) {
  expr_mat <- as.matrix(expr_mat)
  if (!is.numeric(expr_mat) || length(dim(expr_mat)) != 2L) {
    stop("'expr_mat' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(expr_mat) < 2L || ncol(expr_mat) < 4L) {
    stop(
      "'expr_mat' must contain at least 2 genes and 4 samples.",
      call. = FALSE
    )
  }
  if (any(is.infinite(expr_mat))) {
    stop("'expr_mat' contains infinite values.", call. = FALSE)
  }
  if (is.null(rownames(expr_mat))) {
    rownames(expr_mat) <- paste0("Gene", seq_len(nrow(expr_mat)))
    warning("Gene names were missing and have been generated.", call. = FALSE)
  }
  if (is.null(colnames(expr_mat))) {
    colnames(expr_mat) <- paste0("Sample", seq_len(ncol(expr_mat)))
    warning("Sample names were missing and have been generated.", call. = FALSE)
  }
  if (anyDuplicated(rownames(expr_mat))) {
    stop("Gene names in 'expr_mat' must be unique.", call. = FALSE)
  }
  if (anyDuplicated(colnames(expr_mat))) {
    stop("Sample names in 'expr_mat' must be unique.", call. = FALSE)
  }

  gene_variance <- apply(expr_mat, 1L, stats::var, na.rm = TRUE)
  keep_variance <- is.finite(gene_variance) & gene_variance > 0
  removed_zero_variance <- rownames(expr_mat)[!keep_variance]
  expr_mat <- expr_mat[keep_variance, , drop = FALSE]
  gene_variance <- gene_variance[keep_variance]

  selected_by_variance <- character()
  if (!is.null(top_n) && nrow(expr_mat) > top_n) {
    selected <- order(gene_variance, decreasing = TRUE)[
      seq_len(as.integer(top_n))
    ]
    selected_by_variance <- rownames(expr_mat)[selected]
    expr_mat <- expr_mat[selected, , drop = FALSE]
  }

  dat_expr <- as.data.frame(t(expr_mat), check.names = FALSE)
  qc <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
  removed_samples <- rownames(dat_expr)[!qc$goodSamples]
  removed_genes <- colnames(dat_expr)[!qc$goodGenes]
  dat_expr <- dat_expr[
    qc$goodSamples,
    qc$goodGenes,
    drop = FALSE
  ]

  if (nrow(dat_expr) < 4L || ncol(dat_expr) < 2L) {
    stop(
      "Too few samples or genes remain after WGCNA quality control.",
      call. = FALSE
    )
  }
  if (nrow(dat_expr) < 15L) {
    warning(
      "WGCNA is being run with fewer than 15 samples; module stability may be limited.",
      call. = FALSE
    )
  }

  list(
    data = dat_expr,
    qc = list(
      all_ok = qc$allOK &&
        length(removed_zero_variance) == 0L,
      removed_zero_variance_genes = removed_zero_variance,
      removed_samples = removed_samples,
      removed_genes = removed_genes,
      variance_selected_genes = selected_by_variance
    )
  )
}

.prepare_wgcna_traits <- function(traits, sample_names) {
  if (is.null(traits)) {
    return(list(design = NULL, metadata = NULL))
  }
  if (is.atomic(traits) && is.null(dim(traits))) {
    trait_name <- deparse(substitute(traits))
    if (!nzchar(trait_name) || trait_name == "traits") {
      trait_name <- "trait"
    }
    traits <- data.frame(
      stats::setNames(list(traits), trait_name),
      check.names = FALSE
    )
  } else {
    traits <- as.data.frame(traits, check.names = FALSE)
  }
  if (nrow(traits) != length(sample_names) &&
      (is.null(rownames(traits)) ||
       !all(sample_names %in% rownames(traits)))) {
    stop(
      "Rows of 'traits' must match the expression samples.",
      call. = FALSE
    )
  }
  if (!is.null(rownames(traits)) &&
      all(sample_names %in% rownames(traits))) {
    traits <- traits[sample_names, , drop = FALSE]
  } else {
    rownames(traits) <- sample_names
  }

  design_parts <- lapply(names(traits), function(name) {
    value <- traits[[name]]
    if (is.numeric(value) || is.integer(value) || is.logical(value)) {
      out <- matrix(as.numeric(value), ncol = 1L)
      colnames(out) <- name
      return(out)
    }
    value <- factor(value)
    out <- stats::model.matrix(~ 0 + value)
    colnames(out) <- paste0(
      make.names(name),
      "_",
      make.names(levels(value))
    )
    out
  })
  design <- do.call(cbind, design_parts)
  rownames(design) <- sample_names
  storage.mode(design) <- "double"

  trait_variance <- apply(design, 2L, stats::var, na.rm = TRUE)
  keep <- is.finite(trait_variance) & trait_variance > 0
  if (!all(keep)) {
    warning(
      "Constant or non-informative trait columns were removed: ",
      paste(colnames(design)[!keep], collapse = ", "),
      call. = FALSE
    )
    design <- design[, keep, drop = FALSE]
  }
  if (ncol(design) == 0L) {
    stop("No informative trait columns remain.", call. = FALSE)
  }

  list(design = design, metadata = traits)
}

.select_wgcna_power <- function(
    dat_expr,
    power,
    powers,
    scale_free_R2,
    network_type,
    cor_type,
    verbose) {
  if (!is.null(power)) {
    return(list(
      power = as.numeric(power),
      source = "user",
      diagnostics = NULL
    ))
  }

  pick_power <- function() {
    WGCNA::pickSoftThreshold(
      dat_expr,
      RsquaredCut = scale_free_R2,
      powerVector = sort(unique(as.numeric(powers))),
      corFnc = if (cor_type == "pearson") WGCNA::cor else WGCNA::bicor,
      corOptions = list(use = "p"),
      networkType = network_type,
      verbose = verbose
    )
  }
  if (verbose > 0) {
    fit <- pick_power()
  } else {
    fit <- NULL
    invisible(utils::capture.output(fit <- pick_power()))
  }
  diagnostics <- as.data.frame(fit$fitIndices, check.names = FALSE)
  diagnostics$signed_R2 <- -sign(diagnostics$slope) * diagnostics$SFT.R.sq

  selected_power <- fit$powerEstimate
  source <- "powerEstimate"
  if (length(selected_power) == 0L ||
      !is.finite(selected_power)) {
    acceptable <- which(
      is.finite(diagnostics$signed_R2) &
        diagnostics$signed_R2 >= scale_free_R2
    )
    if (length(acceptable) > 0L) {
      selected_power <- diagnostics$Power[acceptable[[1L]]]
      source <- "first_R2_match"
    } else {
      best <- which.max(replace(
        diagnostics$signed_R2,
        !is.finite(diagnostics$signed_R2),
        -Inf
      ))
      selected_power <- diagnostics$Power[[best]]
      source <- "best_available_R2"
      warning(
        "No candidate power reached scale_free_R2 = ",
        scale_free_R2,
        "; using power ",
        selected_power,
        " with the highest available signed R2.",
        call. = FALSE
      )
    }
  }

  list(
    power = as.numeric(selected_power),
    source = source,
    diagnostics = diagnostics
  )
}

.wgcna_cor <- function(x, y, cor_type) {
  if (cor_type == "bicor") {
    return(WGCNA::bicor(x, y, use = "pairwise.complete.obs"))
  }
  stats::cor(x, y, use = "pairwise.complete.obs", method = "pearson")
}
