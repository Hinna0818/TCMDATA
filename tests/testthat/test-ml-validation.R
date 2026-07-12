.make_ml_fixture <- function(n_positive = 10L, n_negative = 30L, p = 6L) {
  set.seed(100)
  n <- n_positive + n_negative
  group <- factor(
    c(rep("Disease", n_positive), rep("Control", n_negative)),
    levels = c("Disease", "Control")
  )
  expr <- matrix(stats::rnorm(p * n), nrow = p)
  rownames(expr) <- paste0("G", seq_len(p))
  colnames(expr) <- paste0("S", seq_len(n))
  expr[1, group == "Disease"] <- expr[1, group == "Disease"] + 2
  list(expr = expr, group = group)
}

test_that("prepare_ml_data records balanced class distributions", {
  dat <- .make_ml_fixture(n_positive = 10, n_negative = 10)

  expect_no_warning(
    ml_data <- suppressMessages(prepare_ml_data(
      dat$expr, dat$group, imbalance_action = "warn"
    ))
  )

  expect_false(ml_data$class_balance$is_imbalanced)
  expect_equal(ml_data$class_balance$counts, c(Disease = 10L, Control = 10L))
  expect_equal(ml_data$class_balance$minority_fraction, 0.5)
  expect_true(all(c(
    "train_x", "train_y", "test_x", "test_y", "gene_names", "levels", "full_cv"
  ) %in% names(ml_data)))
})

test_that("prepare_ml_data warns for and records imbalanced training data", {
  dat <- .make_ml_fixture(n_positive = 4, n_negative = 36)

  expect_warning(
    ml_data <- suppressMessages(prepare_ml_data(
      dat$expr,
      dat$group,
      imbalance_threshold = 0.25,
      imbalance_action = "warn"
    )),
    "Class imbalance detected"
  )

  expect_true(ml_data$class_balance$is_imbalanced)
  expect_equal(ml_data$class_balance$minority_fraction, 0.1)
  expect_equal(ml_data$class_balance$imbalance_ratio, 4 / 36)
})

test_that("downsampling is reproducible and leaves held-out data untouched", {
  dat <- .make_ml_fixture(n_positive = 6, n_negative = 18)
  train_idx <- c(1:4, 7:18)
  ml_data <- suppressMessages(prepare_ml_data(
    dat$expr,
    dat$group,
    split = TRUE,
    train_idx = train_idx,
    imbalance_action = "none"
  ))
  test_counts <- table(ml_data$test_y)

  down1 <- .downsample_training_data(ml_data$train_x, ml_data$train_y, seed = 9)
  down2 <- .downsample_training_data(ml_data$train_x, ml_data$train_y, seed = 9)

  expect_equal(down1$index, down2$index)
  expect_equal(as.integer(table(down1$y)), c(4L, 4L))
  expect_equal(table(ml_data$test_y), test_counts)
})

test_that("balanced case weights use training frequencies and average one", {
  y <- factor(c(rep("Disease", 2), rep("Control", 8)),
              levels = c("Disease", "Control"))
  weights <- .balanced_case_weights(y)

  expect_equal(mean(weights), 1)
  expect_gt(unique(weights[y == "Disease"]), unique(weights[y == "Control"]))
  expect_equal(length(weights), length(y))
})

test_that("model balance interfaces preserve none as the default", {
  functions <- list(ml_enet, ml_lasso, ml_ridge, ml_rf, ml_svm_rfe,
                    ml_xgboost, run_ml_screening)
  defaults <- lapply(functions, function(fun) {
    eval(formals(fun)$balance_method, envir = environment(fun))
  })
  expect_true(all(vapply(defaults, identical, logical(1),
                         c("none", "auto", "weights", "down"))))
})

test_that("SVM-RFE rejects an unsupported explicit weight strategy", {
  dat <- .make_ml_fixture(n_positive = 5, n_negative = 15)
  ml_data <- suppressMessages(prepare_ml_data(
    dat$expr, dat$group, imbalance_action = "none"
  ))

  expect_error(
    ml_svm_rfe(ml_data, balance_method = "weights"),
    "does not support.*weights"
  )
})

test_that("LASSO auto balancing records a training-only weight strategy", {
  skip_if_not_installed("glmnet")
  dat <- .make_ml_fixture(n_positive = 8, n_negative = 32)
  ml_data <- suppressMessages(prepare_ml_data(
    dat$expr, dat$group, imbalance_action = "none"
  ))

  fit <- suppressWarnings(suppressMessages(ml_lasso(
    ml_data,
    cv_folds = 4,
    seed = 21,
    balance_method = "auto"
  )))

  expect_equal(fit$balance_info$requested, "auto")
  expect_equal(fit$balance_info$applied, "weights")
  expect_equal(fit$balance_info$train_counts_before, c(Disease = 8L, Control = 32L))
  expect_equal(mean(fit$balance_info$case_weights), 1)
  expect_equal(nrow(fit$ml_data$train_x), 40L)
})

test_that("LASSO downsampling uses class-stratified CV folds", {
  skip_if_not_installed("glmnet")
  dat <- .make_ml_fixture(n_positive = 8, n_negative = 32)
  ml_data <- suppressMessages(prepare_ml_data(
    dat$expr, dat$group, imbalance_action = "none"
  ))
  down <- .downsample_training_data(
    ml_data$train_x, ml_data$train_y, seed = 7
  )

  fit <- suppressWarnings(suppressMessages(ml_lasso(
    ml_data,
    cv_folds = 4,
    seed = 7,
    balance_method = "down"
  )))
  fold_has_both_classes <- vapply(
    split(down$y, fit$cv_fit$foldid),
    function(y) length(unique(y)) == 2L,
    logical(1)
  )

  expect_true(all(fold_has_both_classes))
  expect_equal(nrow(fit$ml_data$train_x), 40L)
})

test_that("XGBoost weighted CV uses viable stratified folds", {
  skip_if_not_installed("xgboost")
  dat <- .make_ml_fixture(n_positive = 8, n_negative = 32)
  ml_data <- suppressMessages(prepare_ml_data(
    dat$expr, dat$group, imbalance_action = "none"
  ))

  fit <- suppressWarnings(suppressMessages(ml_xgboost(
    ml_data,
    top_n = 3,
    nrounds = 8,
    max_depth = 2,
    cv_folds = 4,
    early_stopping_rounds = 3,
    seed = 7,
    balance_method = "auto"
  )))

  expect_equal(fit$balance_info$applied, "weights")
  expect_equal(fit$balance_info$scale_pos_weight, 4)
})

test_that("unified ROC metrics return AUC confidence intervals and counts", {
  truth <- factor(
    c(rep("Disease", 6), rep("Control", 6)),
    levels = c("Disease", "Control")
  )
  probability <- c(seq(0.7, 0.95, length.out = 6),
                   seq(0.05, 0.3, length.out = 6))

  metrics <- .compute_roc_metrics(
    truth, probability, positive_class = "Disease", ci_method = "delong"
  )

  expect_equal(metrics$auc, 1)
  expect_equal(metrics$n_positive, 6L)
  expect_equal(metrics$n_negative, 6L)
  expect_true(metrics$ci_lower <= metrics$auc)
  expect_true(metrics$ci_upper >= metrics$auc)
})

test_that("get_gene_auc keeps status rows for unusable genes", {
  group <- factor(
    c(rep("Disease", 6), rep("Control", 6)),
    levels = c("Disease", "Control")
  )
  expr <- rbind(
    signal = c(seq(7, 12), seq(1, 6)),
    constant = rep(1, 12),
    missing = rep(NA_real_, 12)
  )

  auc <- get_gene_auc(
    c("signal", "constant", "missing"),
    expr_mat = expr,
    group = group,
    ci_method = "delong",
    conf_level = 0.95
  )

  expect_equal(auc$status[auc$gene == "signal"], "ok")
  expect_equal(auc$status[auc$gene == "constant"], "constant_or_missing")
  expect_equal(auc$status[auc$gene == "missing"], "constant_or_missing")
  expect_true(is.na(auc$auc[auc$gene == "constant"]))
  expect_equal(auc$n_positive[auc$gene == "signal"], 6L)
  expect_equal(auc$n_negative[auc$gene == "signal"], 6L)
  expect_equal(auc$ci_method[auc$gene == "signal"], "delong")
  expect_equal(auc$evaluation_set[auc$gene == "signal"], "provided")
})

.nested_test_model_args <- function() {
  list(
    lasso = list(lambda_rule = "min"),
    rf = list(
      n_trees = 25,
      max_runs = 11,
      refit_on_selected = FALSE,
      boruta_fallback_n = 3
    ),
    svm_rfe = list(
      sizes = c(2, 3),
      top_n = 2,
      min_size = 2
    ),
    xgboost = list(
      top_n = 3,
      nrounds = 10,
      max_depth = 2,
      early_stopping_rounds = 3
    )
  )
}

test_that("nested CV returns strict OOF predictions for all supported models", {
  skip_if_not_installed("caret")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("Boruta")
  skip_if_not_installed("kernlab")
  skip_if_not_installed("xgboost")
  dat <- .make_ml_fixture(n_positive = 18, n_negative = 18)
  ml_data <- suppressMessages(prepare_ml_data(
    dat$expr, dat$group, imbalance_action = "none"
  ))

  nested <- suppressWarnings(suppressMessages(evaluate_ml_nested_cv(
    ml_data,
    methods = c("lasso", "rf", "svm_rfe", "xgboost"),
    outer_folds = 3,
    inner_folds = 3,
    inner_repeats = 1,
    seed = 31,
    model_args = .nested_test_model_args()
  )))

  expect_s3_class(nested, "tcm_nested_cv")
  expect_equal(sort(nested$performance$method),
               sort(c("lasso", "rf", "svm_rfe", "xgboost")))
  expect_equal(nrow(nested$predictions), 36L * 4L)
  per_sample <- table(nested$predictions$method, nested$predictions$sample_id)
  expect_true(all(per_sample == 1L))
  expect_true(all(c(
    "auc", "ci_lower", "ci_upper", "n_positive", "n_negative"
  ) %in% names(nested$performance)))
  expect_true(all(nested$feature_stability$selection_frequency >= 0))
  expect_true(all(nested$feature_stability$selection_frequency <= 1))
  expect_null(nested$fold_models)

  for (fold in nested$fold_indices) {
    expect_length(intersect(fold$train_sample_ids, fold$validation_sample_ids), 0L)
  }
})

test_that("repeated nested CV aggregates patient-level predictions for AUC", {
  skip_if_not_installed("glmnet")
  dat <- .make_ml_fixture(n_positive = 18, n_negative = 18)
  ml_data <- suppressMessages(prepare_ml_data(
    dat$expr, dat$group, imbalance_action = "none"
  ))

  nested <- suppressWarnings(suppressMessages(evaluate_ml_nested_cv(
    ml_data,
    methods = "lasso",
    outer_folds = 3,
    outer_repeats = 2,
    inner_folds = 3,
    seed = 32,
    model_args = list(lasso = list(lambda_rule = "min"))
  )))

  expect_equal(nrow(nested$predictions), 72L)
  expect_true(all(table(nested$predictions$sample_id) == 2L))
  expect_equal(nested$performance$n_samples, 36L)
  expect_equal(nrow(nested$aggregated_predictions), 36L)
})

test_that("nested CV leaves an existing hold-out set isolated", {
  skip_if_not_installed("glmnet")
  dat <- .make_ml_fixture(n_positive = 18, n_negative = 18)
  train_idx <- c(1:15, 19:33)
  ml_data <- suppressMessages(prepare_ml_data(
    dat$expr,
    dat$group,
    split = TRUE,
    train_idx = train_idx,
    imbalance_action = "none"
  ))
  test_x_before <- ml_data$test_x
  test_y_before <- ml_data$test_y

  nested <- suppressWarnings(suppressMessages(evaluate_ml_nested_cv(
    ml_data,
    methods = "lasso",
    outer_folds = 3,
    inner_folds = 3,
    seed = 33,
    model_args = list(lasso = list(lambda_rule = "min"))
  )))

  expect_equal(nrow(nested$predictions), nrow(ml_data$train_x))
  expect_equal(ml_data$test_x, test_x_before)
  expect_equal(ml_data$test_y, test_y_before)
  expect_equal(nested$params$isolated_test_samples, nrow(test_x_before))
})

test_that("nested CV rejects folds that cannot support class-wise inference", {
  dat <- .make_ml_fixture(n_positive = 4, n_negative = 16)
  ml_data <- suppressMessages(prepare_ml_data(
    dat$expr, dat$group, imbalance_action = "none"
  ))

  expect_error(
    evaluate_ml_nested_cv(
      ml_data,
      methods = "lasso",
      outer_folds = 3,
      inner_folds = 3
    ),
    "lower outer_folds"
  )
})

test_that("nested CV requires named model-specific argument lists", {
  dat <- .make_ml_fixture(n_positive = 18, n_negative = 18)
  ml_data <- suppressMessages(prepare_ml_data(
    dat$expr, dat$group, imbalance_action = "none"
  ))

  expect_error(
    evaluate_ml_nested_cv(
      ml_data,
      methods = "lasso",
      outer_folds = 3,
      inner_folds = 3,
      model_args = list(list(lambda_rule = "min"))
    ),
    "named by model"
  )
})
