test_that("prepare_disease_weights requires scored disease target sources by default", {
  expect_error(
    prepare_disease_weights(DOSE = c("IL6", "TNF")),
    "no target-level score"
  )
})

test_that("prepare_disease_weights combines scored sources", {
  otp <- data.frame(
    gene_symbol = c("IL6", "VEGFA"),
    score = c(0.9, 0.5),
    stringsAsFactors = FALSE
  )
  gcds <- data.frame(
    symbol = c("IL6", "TNF"),
    score = c(67.2, 43.1),
    stringsAsFactors = FALSE
  )

  res <- prepare_disease_weights(OpenTargets = otp, GeneCards = gcds)

  expect_true(all(c("IL6", "TNF", "VEGFA") %in% res$symbol))
  expect_true(all(res$disease_weight >= 0 & res$disease_weight <= 1))
  expect_gt(res$disease_weight[res$symbol == "IL6"], res$disease_weight[res$symbol == "TNF"])
  expect_match(res$sources[res$symbol == "IL6"], "GeneCards")
  expect_match(res$sources[res$symbol == "IL6"], "OpenTargets")
})

test_that("DEG-like input is unscored unless explicitly enabled", {
  deg <- data.frame(
    names = c("IL6", "TNF"),
    log2FoldChange = c(4, 1),
    padj = c(1e-6, 0.05),
    stringsAsFactors = FALSE
  )

  expect_error(prepare_disease_weights(DEG = deg), "no score column")

  unscored <- prepare_disease_weights(DEG = deg, require_score = FALSE)
  scored <- prepare_disease_weights(DEG = deg, deg_score = TRUE)

  expect_equal(unscored$disease_weight[unscored$symbol == "IL6"],
               unscored$disease_weight[unscored$symbol == "TNF"])
  expect_gt(scored$disease_weight[scored$symbol == "IL6"],
            scored$disease_weight[scored$symbol == "TNF"])
})

test_that("prepare_herb_target_weights counts compound support", {
  herb_df <- data.frame(
    herb = c("H1", "H1", "H1", "H2"),
    molecule = c("M1", "M2", "M2", "M3"),
    target = c("IL6", "IL6", "TNF", "TNF"),
    stringsAsFactors = FALSE
  )

  res <- prepare_herb_target_weights(herb_df, method = "compound_count")

  expect_equal(res$herb_target_weight[res$target == "IL6"], 2)
  expect_equal(res$herb_target_weight[res$target == "TNF"], 2)
  expect_equal(res$n_herbs[res$target == "TNF"], 2)
})

test_that("rank_tcm_targets_by_ppi ranks candidates on a weighted graph", {
  edges <- data.frame(
    from = c("IL6", "TNF", "AKT1", "VEGFA", "JUN"),
    to = c("AKT1", "AKT1", "EGFR", "EGFR", "EGFR"),
    score = c(0.9, 0.85, 0.8, 0.7, 0.6),
    stringsAsFactors = FALSE
  )
  edges$distance <- 1 / edges$score
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)

  herb <- data.frame(
    target = c("AKT1", "EGFR", "JUN"),
    herb_target_weight = c(2, 1, 1),
    stringsAsFactors = FALSE
  )
  disease <- prepare_disease_weights(
    OpenTargets = data.frame(gene_symbol = c("IL6", "TNF"), score = c(0.9, 0.6))
  )

  res <- rank_tcm_targets_by_ppi(
    herb_targets = herb,
    disease_targets = disease,
    ppi = g,
    n_perm = 100,
    degree_bins = 3,
    seed = 1
  )

  expect_s3_class(res, "tcm_ppi_rank")
  expect_equal(nrow(res$result), 3)
  expect_true(all(c(
    "Score_final", "z_proximity", "z_proximity_no_self",
    "p_empirical", "p_adjust", "p_empirical_no_self", "p_adjust_no_self",
    "nearest_disease_target", "nearest_distance", "random_pool_size"
  ) %in% names(res$result)))
  expect_equal(res$result$Score_final, res$result$z_proximity_no_self)
  expect_true(all(c("IL6", "TNF") %in% res$disease_module$symbol))
})

test_that("rank_tcm_targets_by_ppi can load PPI from an RDS path", {
  edges <- data.frame(
    from = c("IL6", "AKT1", "EGFR"),
    to = c("AKT1", "EGFR", "JUN"),
    score = c(0.9, 0.8, 0.7),
    stringsAsFactors = FALSE
  )
  edges$distance <- 1 / edges$score
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  ppi_path <- tempfile(fileext = ".rds")
  saveRDS(g, ppi_path)

  res <- rank_tcm_targets_by_ppi(
    herb_targets = c("AKT1", "EGFR"),
    disease_targets = data.frame(symbol = "IL6", score = 0.9),
    ppi_path = ppi_path,
    n_perm = 100,
    seed = 1
  )

  expect_s3_class(res, "tcm_ppi_rank")
  expect_equal(nrow(res$result), 2)
})

test_that("evaluate_tcm_target_sets_by_ppi compares target sets", {
  edges <- data.frame(
    from = c("IL6", "TNF", "AKT1", "VEGFA", "JUN", "MAPK1"),
    to = c("AKT1", "AKT1", "EGFR", "EGFR", "EGFR", "JUN"),
    score = c(0.9, 0.85, 0.8, 0.7, 0.6, 0.75),
    stringsAsFactors = FALSE
  )
  edges$distance <- 1 / edges$score
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)

  sets <- list(
    near = c("AKT1", "EGFR", "JUN"),
    far = c("MAPK1", "JUN")
  )
  disease <- prepare_disease_weights(
    OpenTargets = data.frame(gene_symbol = c("IL6", "TNF"), score = c(0.9, 0.6))
  )

  res <- evaluate_tcm_target_sets_by_ppi(
    target_sets = sets,
    disease_targets = disease,
    ppi = g,
    n_perm = 100,
    degree_bins = 3,
    seed = 1
  )

  expect_s3_class(res, "tcm_ppi_set_eval")
  expect_equal(nrow(res$result), 2)
  expect_true(all(c(
    "overlap_n", "directlink_edges", "jaccard",
    "weighted_proximity", "z_proximity", "p_empirical",
    "Score_set", "Rank_set"
  ) %in% names(res$result)))
  expect_true(all(res$result$directlink_edges >= 0))
  expect_equal(res$result$Score_set, res$result$z_proximity)
})

test_that("evaluate_tcm_target_sets_by_ppi splits grouped data frames", {
  edges <- data.frame(
    from = c("IL6", "AKT1", "EGFR"),
    to = c("AKT1", "EGFR", "JUN"),
    score = c(0.9, 0.8, 0.7),
    stringsAsFactors = FALSE
  )
  edges$distance <- 1 / edges$score
  g <- igraph::graph_from_data_frame(edges, directed = FALSE)

  grouped <- data.frame(
    formula = c("F1", "F1", "F2"),
    target = c("AKT1", "EGFR", "JUN"),
    stringsAsFactors = FALSE
  )

  res <- evaluate_tcm_target_sets_by_ppi(
    target_sets = grouped,
    disease_targets = data.frame(symbol = "IL6", score = 0.9),
    ppi = g,
    set_col = "formula",
    n_perm = 100,
    seed = 1
  )

  expect_equal(sort(res$result$set), c("F1", "F2"))
  expect_true(all(res$result$n_targets_in_ppi >= 1))
})

test_that("select_tcm_targets filters ranking results", {
  ranking <- data.frame(
    target = c("A", "B", "C"),
    Score_final = c(0.9, 0.7, 0.2),
    Rank_final = c(1, 2, 3),
    z_proximity = c(3, 2, -1),
    p_empirical = c(0.01, 0.2, 0.8),
    p_adjust = c(0.03, 0.3, 0.8),
    direct_overlap = c(TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )

  selected <- select_tcm_targets(ranking, top_n = 2, min_z = 1)
  expect_equal(selected$target, c("A", "B"))

  selected_targets <- select_tcm_targets(
    ranking,
    top_n = NULL,
    max_p_adjust = 0.05,
    return = "targets"
  )
  expect_equal(selected_targets, "A")
})
