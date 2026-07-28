.simulate_wgcna_expression <- function(seed = 11) {
  set.seed(seed)
  n_samples <- 24
  group <- factor(rep(c("Control", "Disease"), each = n_samples / 2))
  signal_a <- scale(as.numeric(group == "Disease"))[, 1] +
    stats::rnorm(n_samples, sd = 0.25)
  signal_b <- stats::rnorm(n_samples)
  signal_c <- stats::rnorm(n_samples)

  make_module <- function(signal, n_genes) {
    vapply(
      seq_len(n_genes),
      function(i) signal + stats::rnorm(n_samples, sd = 0.35),
      numeric(n_samples)
    )
  }
  sample_by_gene <- cbind(
    make_module(signal_a, 30),
    make_module(signal_b, 30),
    make_module(signal_c, 30)
  )
  expr <- t(sample_by_gene)
  rownames(expr) <- paste0("Gene", seq_len(nrow(expr)))
  colnames(expr) <- paste0("Sample", seq_len(ncol(expr)))

  traits <- data.frame(
    subtype = group,
    row.names = colnames(expr)
  )
  list(expr = expr, traits = traits)
}

test_that("run_wgcna returns modules and trait associations", {
  skip_if_not_installed("WGCNA")
  simulated <- .simulate_wgcna_expression()

  fit <- suppressMessages(run_wgcna(
    simulated$expr,
    traits = simulated$traits,
    power = 6,
    max_block_size = 200,
    min_module_size = 10,
    merge_cut_height = 0.15,
    seed = 7,
    n_threads = 1,
    verbose = 0
  ))

  expect_s3_class(fit, "tcm_wgcna")
  expect_identical(fit$power, 6)
  expect_identical(fit$power_source, "user")
  expect_equal(nrow(fit$expression), ncol(simulated$expr))
  expect_equal(length(fit$module_colors), ncol(fit$expression))
  expect_named(fit$module_colors, colnames(fit$expression))
  expect_true(all(c(
    "gene", "module", "module_membership"
  ) %in% colnames(fit$gene_info)))
  expect_true(length(fit$module_genes) >= 1L)
  expect_equal(
    sort(unlist(fit$module_genes, use.names = FALSE)),
    sort(fit$gene_info$gene)
  )
  expect_true(is.matrix(fit$module_trait_cor))
  expect_equal(dim(fit$module_trait_cor), dim(fit$module_trait_p))
  expect_true(is.matrix(fit$gene_trait_cor))
})

test_that("get_wgcna_module_genes filters a fitted module", {
  skip_if_not_installed("WGCNA")
  simulated <- .simulate_wgcna_expression(seed = 22)
  fit <- suppressMessages(run_wgcna(
    simulated$expr,
    traits = simulated$traits,
    power = 6,
    max_block_size = 200,
    min_module_size = 10,
    seed = 8,
    n_threads = 1,
    verbose = 0
  ))

  module <- names(sort(table(fit$module_colors), decreasing = TRUE))[[1L]]
  genes <- get_wgcna_module_genes(
    fit,
    module = module,
    min_module_membership = 0
  )

  expect_type(genes, "character")
  expect_setequal(genes, fit$module_genes[[module]])
  expect_error(
    get_wgcna_module_genes(fit, module = "not-a-module"),
    "was not found"
  )
})

test_that("run_wgcna can select a power from candidate values", {
  skip_if_not_installed("WGCNA")
  simulated <- .simulate_wgcna_expression(seed = 33)

  fit <- suppressWarnings(suppressMessages(run_wgcna(
    simulated$expr,
    power = NULL,
    powers = 1:4,
    scale_free_R2 = 0.5,
    max_block_size = 200,
    min_module_size = 10,
    seed = 9,
    n_threads = 1,
    verbose = 0
  )))

  expect_true(fit$power %in% 1:4)
  expect_true(is.data.frame(fit$power_diagnostics))
  expect_true("signed_R2" %in% colnames(fit$power_diagnostics))
})

test_that("WGCNA plotting helpers accept run_wgcna output", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("ggtree")
  skip_if_not_installed("aplot")
  simulated <- .simulate_wgcna_expression(seed = 44)
  fit <- suppressMessages(run_wgcna(
    simulated$expr,
    traits = simulated$traits,
    power = 6,
    max_block_size = 200,
    min_module_size = 10,
    seed = 10,
    n_threads = 1,
    verbose = 0
  ))

  module_plot <- plot_wgcna_modules(fit)
  eigengene_plot <- plot_wgcna_eigengenes(fit)
  eigengene_heatmap <- plot_wgcna_eigengenes(
    fit,
    show_heatmap = TRUE
  )
  trait_plot <- plot_wgcna_traits(fit)

  expect_true(inherits(module_plot, "ggplot"))
  expect_s3_class(eigengene_plot, "ggplot")
  expect_true(inherits(eigengene_heatmap, "ggplot"))
  expect_s3_class(trait_plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(eigengene_plot))
  expect_silent(ggplot2::ggplot_build(trait_plot))
})
