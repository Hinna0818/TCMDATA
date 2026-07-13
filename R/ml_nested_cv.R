#' Evaluate ML feature selection with nested cross-validation
#'
#' Runs an outer stratified cross-validation loop around the existing TCMDATA
#' model functions. Feature filtering, class balancing, tuning, and feature
#' selection are performed using only the outer training fold. Existing Mode B
#' or Mode C test data in `ml_data` remain isolated.
#'
#' @param ml_data A `tcm_ml_data` object from [prepare_ml_data()].
#' @param methods Models to evaluate: `"lasso"`, `"rf"`, `"svm_rfe"`,
#'   and/or `"xgboost"`.
#' @param outer_folds Number of stratified outer folds.
#' @param outer_repeats Number of repeated outer fold partitions.
#' @param inner_folds Number of inner folds passed to model functions.
#' @param inner_repeats Number of inner repeats used by SVM-RFE.
#' @param balance_method Class-balancing strategy passed to each model.
#' @param metric Optimization metric passed to SVM-RFE. Nested performance is
#'   always summarized with ROC AUC.
#' @param seed Random seed.
#' @param keep_fold_models Whether to retain fitted models from every outer fold.
#' @param model_args Optional named list of additional arguments for each model,
#'   for example `list(rf = list(n_trees = 300))`. Fold-control arguments
#'   are managed by this function and cannot be overridden here.
#'
#' @return A `tcm_nested_cv` object containing patient-level performance,
#'   outer-fold predictions, feature stability, fold summaries, optional fold
#'   models, parameters, and auditable fold indices.
#'
#' @examples
#' \dontrun{
#' nested <- evaluate_ml_nested_cv(
#'   ml_data,
#'   methods = c("lasso", "rf"),
#'   outer_folds = 5,
#'   inner_folds = 5
#' )
#' nested$performance
#' }
#'
#' @export
evaluate_ml_nested_cv <- function(
    ml_data,
    methods = c("lasso", "rf", "svm_rfe", "xgboost"),
    outer_folds = 5L,
    outer_repeats = 1L,
    inner_folds = 5L,
    inner_repeats = 1L,
    balance_method = c("none", "auto", "weights", "down"),
    metric = "ROC",
    seed = 2025L,
    keep_fold_models = FALSE,
    model_args = list()) {
  if (!inherits(ml_data, "tcm_ml_data")) {
    stop("ml_data must be a tcm_ml_data object.", call. = FALSE)
  }
  .check_ml_deps(c("caret", "pROC"))
  methods <- match.arg(
    methods,
    c("lasso", "rf", "svm_rfe", "xgboost"),
    several.ok = TRUE
  )
  methods <- unique(methods)
  balance_method <- match.arg(balance_method)
  outer_folds <- .nested_positive_integer(outer_folds, "outer_folds", minimum = 2L)
  outer_repeats <- .nested_positive_integer(outer_repeats, "outer_repeats")
  inner_folds <- .nested_positive_integer(inner_folds, "inner_folds", minimum = 2L)
  inner_repeats <- .nested_positive_integer(inner_repeats, "inner_repeats")
  if (!is.logical(keep_fold_models) || length(keep_fold_models) != 1L ||
      is.na(keep_fold_models)) {
    stop("keep_fold_models must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.list(model_args)) {
    stop("model_args must be a named list.", call. = FALSE)
  }
  if (length(model_args) > 0L &&
      (is.null(names(model_args)) || any(!nzchar(names(model_args))))) {
    stop("Non-empty model_args must be named by model.", call. = FALSE)
  }
  unknown_model_args <- setdiff(names(model_args), methods)
  if (length(unknown_model_args) > 0L) {
    stop(
      "model_args contains entries not requested in methods: ",
      paste(unknown_model_args, collapse = ", "),
      call. = FALSE
    )
  }

  x <- as.data.frame(ml_data$train_x)
  y <- factor(ml_data$train_y, levels = ml_data$levels)
  if (nlevels(y) != 2L || any(table(y) == 0L)) {
    stop("Nested CV requires both outcome classes in training data.",
         call. = FALSE)
  }
  if (min(table(y)) < 2L * outer_folds) {
    stop(
      "Each outer validation fold needs at least two samples per class; lower outer_folds.",
      call. = FALSE
    )
  }
  sample_ids <- rownames(x)
  if (is.null(sample_ids) || any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
    sample_ids <- paste0("sample_", seq_len(nrow(x)))
  }

  prediction_rows <- list()
  fold_rows <- list()
  selection_rows <- list()
  fold_indices <- list()
  fold_models <- if (isTRUE(keep_fold_models)) list() else NULL
  prediction_index <- 0L
  fold_row_index <- 0L
  selection_index <- 0L
  fold_index <- 0L

  for (outer_repeat in seq_len(outer_repeats)) {
    set.seed(seed + outer_repeat - 1L)
    folds <- caret::createFolds(y, k = outer_folds, list = TRUE,
                                returnTrain = FALSE)

    for (outer_fold in seq_along(folds)) {
      validation_index <- sort(as.integer(folds[[outer_fold]]))
      train_index <- setdiff(seq_len(nrow(x)), validation_index)
      train_y <- factor(y[train_index], levels = ml_data$levels)
      validation_y <- factor(y[validation_index], levels = ml_data$levels)
      if (any(table(train_y) == 0L)) {
        stop(
          "An outer training fold contains only one class; lower outer_folds.",
          call. = FALSE
        )
      }
      if (min(table(train_y)) < inner_folds) {
        stop(
          "An outer training fold has fewer class samples than inner_folds; lower inner_folds.",
          call. = FALSE
        )
      }

      fold_data <- suppressMessages(prepare_ml_data(
        expr_mat = t(as.matrix(x[train_index, , drop = FALSE])),
        group = train_y,
        positive_class = ml_data$levels[[1]],
        test_expr = t(as.matrix(x[validation_index, , drop = FALSE])),
        test_group = validation_y,
        seed = seed + outer_repeat * 1000L + outer_fold,
        imbalance_threshold = if (!is.null(ml_data$class_balance$threshold)) {
          ml_data$class_balance$threshold
        } else {
          0.25
        },
        imbalance_action = "none"
      ))

      fold_index <- fold_index + 1L
      fold_indices[[fold_index]] <- list(
        outer_repeat = outer_repeat,
        outer_fold = outer_fold,
        train_sample_ids = sample_ids[train_index],
        validation_sample_ids = sample_ids[validation_index]
      )

      for (method in methods) {
        fold_seed <- seed + outer_repeat * 10000L + outer_fold * 100L +
          match(method, methods)
        fit <- .fit_nested_method(
          method = method,
          fold_data = fold_data,
          inner_folds = inner_folds,
          inner_repeats = inner_repeats,
          balance_method = balance_method,
          metric = metric,
          seed = fold_seed,
          model_args = model_args[[method]]
        )
        probability <- .nested_test_probability(fit, ml_data$levels[[1]])
        predicted <- as.character(fit$test_performance$predictions)
        if (length(probability) != length(validation_index) ||
            length(predicted) != length(validation_index)) {
          stop(
            "Model ", method,
            " did not return one prediction per outer validation sample.",
            call. = FALSE
          )
        }

        prediction_index <- prediction_index + 1L
        prediction_rows[[prediction_index]] <- data.frame(
          sample_id = sample_ids[validation_index],
          truth = as.character(validation_y),
          probability = as.numeric(probability),
          predicted = predicted,
          method = method,
          outer_fold = as.integer(outer_fold),
          outer_repeat = as.integer(outer_repeat),
          stringsAsFactors = FALSE
        )

        selected <- unique(as.character(fit$selected_features))
        selected <- selected[!is.na(selected) & nzchar(selected)]
        if (length(selected) > 0L) {
          selection_index <- selection_index + 1L
          selection_rows[[selection_index]] <- data.frame(
            method = method,
            gene = selected,
            outer_repeat = outer_repeat,
            outer_fold = outer_fold,
            stringsAsFactors = FALSE
          )
        }

        fold_row_index <- fold_row_index + 1L
        fold_rows[[fold_row_index]] <- data.frame(
          method = method,
          outer_repeat = outer_repeat,
          outer_fold = outer_fold,
          n_train = length(train_index),
          n_validation = length(validation_index),
          n_positive_train = sum(train_y == ml_data$levels[[1]]),
          n_negative_train = sum(train_y == ml_data$levels[[2]]),
          balance_applied = fit$balance_info$applied,
          n_selected_features = length(selected),
          inner_metric = .nested_inner_metric(fit),
          stringsAsFactors = FALSE
        )

        if (isTRUE(keep_fold_models)) {
          model_name <- paste(method, outer_repeat, outer_fold, sep = "_")
          fold_models[[model_name]] <- fit
        }
      }
    }
  }

  predictions <- do.call(rbind, prediction_rows)
  rownames(predictions) <- NULL
  aggregated_predictions <- stats::aggregate(
    probability ~ method + sample_id + truth,
    data = predictions,
    FUN = mean
  )
  aggregated_predictions$predicted <- ifelse(
    aggregated_predictions$probability >= 0.5,
    ml_data$levels[[1]],
    ml_data$levels[[2]]
  )
  aggregated_predictions <- aggregated_predictions[
    order(aggregated_predictions$method, aggregated_predictions$sample_id),
    ,
    drop = FALSE
  ]
  rownames(aggregated_predictions) <- NULL

  performance <- do.call(rbind, lapply(methods, function(method) {
    method_predictions <- aggregated_predictions[
      aggregated_predictions$method == method,
      ,
      drop = FALSE
    ]
    roc_metrics <- .compute_roc_metrics(
      truth = factor(method_predictions$truth, levels = ml_data$levels),
      probability = method_predictions$probability,
      positive_class = ml_data$levels[[1]],
      ci_method = "delong",
      seed = seed
    )
    data.frame(
      method = method,
      auc = roc_metrics$auc,
      ci_lower = roc_metrics$ci_lower,
      ci_upper = roc_metrics$ci_upper,
      sensitivity = roc_metrics$sensitivity,
      specificity = roc_metrics$specificity,
      n_samples = nrow(method_predictions),
      n_positive = roc_metrics$n_positive,
      n_negative = roc_metrics$n_negative,
      evaluation_set = "nested_oof",
      stringsAsFactors = FALSE
    )
  }))
  rownames(performance) <- NULL

  selected_long <- if (length(selection_rows) > 0L) {
    do.call(rbind, selection_rows)
  } else {
    data.frame(
      method = character(), gene = character(), outer_repeat = integer(),
      outer_fold = integer(), stringsAsFactors = FALSE
    )
  }
  feature_stability <- .nested_feature_stability(
    selected_long,
    methods = methods,
    evaluated_folds = outer_folds * outer_repeats
  )
  fold_summary <- do.call(rbind, fold_rows)
  rownames(fold_summary) <- NULL

  out <- list(
    performance = performance,
    predictions = predictions,
    aggregated_predictions = aggregated_predictions,
    feature_stability = feature_stability,
    fold_summary = fold_summary,
    fold_models = fold_models,
    fold_indices = fold_indices,
    params = list(
      methods = methods,
      outer_folds = outer_folds,
      outer_repeats = outer_repeats,
      inner_folds = inner_folds,
      inner_repeats = inner_repeats,
      balance_method = balance_method,
      metric = metric,
      seed = seed,
      keep_fold_models = keep_fold_models,
      model_args = model_args,
      training_samples = nrow(x),
      isolated_test_samples = if (is.null(ml_data$test_x)) 0L else nrow(ml_data$test_x)
    )
  )
  class(out) <- c("tcm_nested_cv", "list")
  out
}

#' @export
print.tcm_nested_cv <- function(x, ...) {
  cat("Nested cross-validation\n")
  cat("Methods:", paste(x$params$methods, collapse = ", "), "\n")
  cat(
    "Outer folds/repeats:", x$params$outer_folds, "/",
    x$params$outer_repeats, "\n"
  )
  print(x$performance, row.names = FALSE)
  invisible(x)
}

.nested_positive_integer <- function(x, name, minimum = 1L) {
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be an integer of at least ", minimum, ".", call. = FALSE)
  }
  value
}

.fit_nested_method <- function(method,
                               fold_data,
                               inner_folds,
                               inner_repeats,
                               balance_method,
                               metric,
                               seed,
                               model_args = NULL) {
  model_args <- if (is.null(model_args)) list() else model_args
  if (!is.list(model_args)) {
    stop("Each model_args entry must be a list.", call. = FALSE)
  }
  protected <- c(
    "ml_data", "seed", "balance_method", "cv_folds", "cv_repeats", "metric"
  )
  conflicts <- intersect(names(model_args), protected)
  if (length(conflicts) > 0L) {
    stop(
      "Nested CV controls these model argument(s): ",
      paste(conflicts, collapse = ", "),
      call. = FALSE
    )
  }

  fixed_args <- switch(
    method,
    lasso = list(
      ml_data = fold_data,
      cv_folds = inner_folds,
      seed = seed,
      balance_method = balance_method
    ),
    rf = list(
      ml_data = fold_data,
      seed = seed,
      balance_method = balance_method
    ),
    svm_rfe = list(
      ml_data = fold_data,
      cv_folds = inner_folds,
      cv_repeats = inner_repeats,
      metric = metric,
      seed = seed,
      balance_method = balance_method
    ),
    xgboost = list(
      ml_data = fold_data,
      cv_folds = inner_folds,
      seed = seed,
      balance_method = balance_method
    )
  )
  model_function <- switch(
    method,
    lasso = ml_lasso,
    rf = ml_rf,
    svm_rfe = ml_svm_rfe,
    xgboost = ml_xgboost
  )
  do.call(model_function, c(fixed_args, model_args))
}

.nested_test_probability <- function(fit, positive_class) {
  performance <- fit$test_performance
  if (is.null(performance) || is.null(performance$probabilities)) {
    stop("Nested model did not return outer validation probabilities.",
         call. = FALSE)
  }
  probabilities <- as.data.frame(performance$probabilities)
  if (!positive_class %in% names(probabilities)) {
    stop(
      "Positive-class probability column not found: ", positive_class,
      call. = FALSE
    )
  }
  as.numeric(probabilities[[positive_class]])
}

.nested_inner_metric <- function(fit) {
  candidates <- c(
    fit$cv_performance$auc,
    fit$cv_performance$metric_value,
    fit$cv_performance$accuracy
  )
  candidates <- suppressWarnings(as.numeric(candidates))
  candidates <- candidates[is.finite(candidates)]
  if (length(candidates) == 0L) NA_real_ else candidates[[1L]]
}

.nested_feature_stability <- function(selected_long, methods, evaluated_folds) {
  rows <- lapply(methods, function(method) {
    genes <- selected_long$gene[selected_long$method == method]
    if (length(genes) == 0L) return(NULL)
    counts <- sort(table(genes), decreasing = TRUE)
    data.frame(
      method = method,
      gene = names(counts),
      selected_n = as.integer(counts),
      evaluated_folds = as.integer(evaluated_folds),
      selection_frequency = as.integer(counts) / evaluated_folds,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) {
    return(data.frame(
      method = character(), gene = character(), selected_n = integer(),
      evaluated_folds = integer(), selection_frequency = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
