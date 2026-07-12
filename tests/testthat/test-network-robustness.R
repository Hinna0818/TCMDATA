test_that("ppi_knock_impact returns node-level perturbation table", {
  data(demo_ppi)
  target <- igraph::V(demo_ppi)$name[[1]]

  res <- suppressWarnings(ppi_knock_impact(
    demo_ppi,
    targets = target,
    n_perm = 100,
    seed = 1
  ))

  expect_s3_class(res, "tcm_ppi_knock_impact")
  expect_true(all(c("impact", "protein_summary", "params") %in% names(res)))
  expect_equal(res$params$targets, target)
  expect_true(all(c(
    "knocked_targets", "affected_protein", "metric",
    "before", "after", "delta", "relative_change",
    "impact_change", "formula", "metric_sign",
    "random_mean", "random_sd", "random_n",
    "z_delta", "Pvalue", "P_normal", "P_empirical", "P_adjust",
    "P_adjust_normal", "P_adjust_empirical",
    "distance_to_knockout", "is_direct_neighbor"
  ) %in% names(res$impact)))
  expect_false(target %in% res$impact$affected_protein)
  expect_equal(sort(unique(res$impact$metric)), c("AD", "ASPL", "CC", "DC"))
  expect_equal(nrow(res$protein_summary), igraph::vcount(demo_ppi) - 1L)
})

test_that("ppi_knock_impact supports multiple knockout targets", {
  data(demo_ppi)
  targets <- igraph::V(demo_ppi)$name[seq_len(2)]

  res <- suppressWarnings(ppi_knock_impact(
    demo_ppi,
    targets = targets,
    metrics = "AD",
    n_perm = 100,
    seed = 2
  ))

  expect_s3_class(res, "tcm_ppi_knock_impact")
  expect_equal(sort(res$params$targets), sort(targets))
  expect_false(any(targets %in% res$impact$affected_protein))
  expect_equal(nrow(res$protein_summary), igraph::vcount(demo_ppi) - length(targets))
})

test_that("ppi_knock supports exploratory whole-network metrics", {
  data(demo_ppi)
  target <- igraph::V(demo_ppi)$name[[1]]

  res <- suppressWarnings(ppi_knock(
    demo_ppi,
    targets = target,
    metrics = c("ASPL", "AD", "density", "components"),
    n_perm = 100,
    seed = 3
  ))

  expect_true(all(c("ASPL", "AD", "density", "components") %in% res$Summary$Metric))
  expect_true(is.na(res$Total_Score))
  expect_true(is.na(res$Total_Pvalue))
})

test_that("ppi_knock_impact supports extra PPI node metrics", {
  data(demo_ppi)
  target <- igraph::V(demo_ppi)$name[[1]]

  res <- suppressWarnings(ppi_knock_impact(
    demo_ppi,
    targets = target,
    metrics = c("ASPL", "pagerank", "MCC"),
    n_perm = 100,
    seed = 4
  ))

  expect_true(all(c("ASPL", "pagerank", "MCC") %in% res$impact$metric))
  expect_true(all(res$impact$formula[res$impact$metric == "ASPL"] == "(after - before) / abs(before)"))
  expect_true(all(res$impact$formula[res$impact$metric == "pagerank"] == "-(after - before) / abs(before)"))
})

test_that("ppi_knock_impact validates requested metrics", {
  data(demo_ppi)
  target <- igraph::V(demo_ppi)$name[[1]]

  expect_error(
    ppi_knock_impact(demo_ppi, target, metrics = "not_a_metric", n_perm = 100),
    "Unsupported node-level metric"
  )
})

test_that("network perturbation uses 1000 permutations by default", {
  expect_identical(formals(ppi_knock)$n_perm, 1000L)
  expect_identical(formals(ppi_knock_impact)$n_perm, 1000L)
})

test_that("network perturbation validates low permutation counts", {
  data(demo_ppi)
  target <- igraph::V(demo_ppi)$name[[1]]

  expect_error(ppi_knock(demo_ppi, target, n_perm = 99), "at least 100")
  expect_error(ppi_knock_impact(demo_ppi, target, n_perm = 99), "at least 100")
  expect_warning(
    ppi_knock(demo_ppi, target, n_perm = 100, seed = 10),
    "fewer than 1,000"
  )
})

test_that("ppi_knock reports empirical and selected p-values", {
  data(demo_ppi)
  target <- igraph::V(demo_ppi)$name[[1]]

  res <- suppressWarnings(ppi_knock(
    demo_ppi,
    target,
    n_perm = 100,
    seed = 11,
    p_method = "empirical"
  ))

  expect_true(all(c(
    "P_normal", "P_empirical", "P_adjust_empirical", "Pvalue"
  ) %in% names(res$Summary)))
  expect_equal(res$Summary$Pvalue, res$Summary$P_empirical)
  expect_equal(
    res$Summary$P_adjust_empirical,
    stats::p.adjust(res$Summary$P_empirical, method = "BH")
  )
  expect_true(all(res$Summary$P_empirical >= 1 / 101))
  expect_true(all(res$Summary$P_empirical <= 1))
  expect_equal(res$Total_Pvalue, res$Total_P_empirical)
  expect_equal(res$params$p_method, "empirical")
  expect_equal(res$params$n_perm, 100L)
})

test_that("ppi_knock_impact adjusts the selected p-value family", {
  data(demo_ppi)
  target <- igraph::V(demo_ppi)$name[[1]]

  empirical <- suppressWarnings(ppi_knock_impact(
    demo_ppi, target, metrics = "AD", n_perm = 100,
    seed = 12, p_method = "empirical"
  ))
  normal <- suppressWarnings(ppi_knock_impact(
    demo_ppi, target, metrics = "AD", n_perm = 100,
    seed = 12, p_method = "normal"
  ))

  expect_equal(empirical$impact$Pvalue, empirical$impact$P_empirical)
  expect_equal(empirical$impact$P_adjust, empirical$impact$P_adjust_empirical)
  expect_equal(normal$impact$Pvalue, normal$impact$P_normal)
  expect_equal(normal$impact$P_adjust, normal$impact$P_adjust_normal)
  expect_equal(empirical$impact$P_empirical, normal$impact$P_empirical)
})

test_that("zero-variance null distributions do not create artificial z-scores", {
  expect_warning(
    standardized <- .standardize_against_null(
      observed = c(2, 3),
      null_mean = c(1, 1),
      null_sd = c(0, 2)
    ),
    "zero or undefined variance"
  )
  expect_true(is.na(standardized[[1]]))
  expect_equal(standardized[[2]], 1)
})

test_that("node impact summaries keep all-missing z statistics as NA", {
  impact <- data.frame(
    knocked_targets = c("KO", "KO"),
    affected_protein = c("A", "A"),
    metric = c("AD", "CC"),
    distance_to_knockout = c(1, 1),
    is_direct_neighbor = c(TRUE, TRUE),
    P_adjust = c(NA_real_, NA_real_),
    z_delta = c(NA_real_, NA_real_),
    stringsAsFactors = FALSE
  )

  expect_no_warning(summary <- .node_impact_protein_summary(impact))
  expect_true(is.na(summary$max_abs_z))
  expect_true(is.na(summary$mean_abs_z))
  expect_false(summary$affected)
})
