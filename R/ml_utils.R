#' Prepare expression data for ML feature selection
#'
#' Supports three modes:
#' - **Mode A** (default): full-data CV, no hold-out (`split = FALSE`).
#' - **Mode B**: internal train/test split (`split = TRUE`).
#' - **Mode C**: external validation (`test_expr` + `test_group`).
#'
#' @param expr_mat Numeric matrix (genes x samples). Row names = gene symbols.
#' @param group Factor/character of length `ncol(expr_mat)`, two levels.
#' @param positive_class Which level to treat as the positive class (`levels[[1]]`).
#'   Default `NULL`: if `group` is already an ordered factor, its level order is
#'   preserved; otherwise levels are sorted alphabetically.
#' @param genes Optional character vector of candidate genes to keep. Default is NULL.
#' @param split Logical. `TRUE` = Mode B. Default is `FALSE`.
#' @param train_ratio Fraction for training when `split = TRUE`. Default is 0.7.
#' @param train_idx Optional integer indices overriding `train_ratio`.
#' @param test_expr External validation matrix (genes x samples) for Mode C.
#' @param test_group Labels for `test_expr` (required if `test_expr` given).
#' @param seed Random seed. Default is 2025.
#' @param imbalance_threshold Minority-class fraction below which the training
#'   data are flagged as imbalanced. Default is 0.25.
#' @param imbalance_action How to report detected imbalance: a warning, a
#'   message, or no condition.
#'
#' @return A `tcm_ml_data` list containing the original train/test fields plus
#'   training and optional test class-balance summaries.
#' @importFrom stats var
#' @examples
#' \dontrun{
#'   ## Mode A: full CV
#'   ml_data <- prepare_ml_data(expr_mat, group, positive_class = "Disease")
#'   ## Mode B: internal split
#'   ml_data <- prepare_ml_data(expr_mat, group, split = TRUE, train_ratio = 0.7)
#' }
#' @export
prepare_ml_data <- function(expr_mat,
                            group,
                            positive_class = NULL,
                            genes = NULL,
                            split = FALSE,
                            train_ratio = 0.7,
                            train_idx = NULL,
                            test_expr = NULL,
                            test_group = NULL,
                            seed = 2025,
                            imbalance_threshold = 0.25,
                            imbalance_action = c("warn", "message", "none")) {

  imbalance_action <- match.arg(imbalance_action)

  if (isTRUE(split)) {
     .check_ml_deps("caret")
   }

  if (is.data.frame(group) || is.matrix(group)) {
    if (ncol(group) == 1) {
      group <- group[[1]]
    } else if ("group" %in% colnames(group)) {
      group <- group$group
    } else {
      stop("`group` must be a vector or factor, not a data.frame or matrix.")
    }
  }

  expr_mat <- as.matrix(expr_mat)
  group <- factor(group)
  stopifnot(is.numeric(expr_mat),
            ncol(expr_mat) == length(group),
            nlevels(group) == 2L)

  ## Re-order levels so that positive_class is levels[1]
  if (!is.null(positive_class)) {
    if (!positive_class %in% levels(group))
      stop(sprintf("positive_class '%s' not found in group levels: %s",
                   positive_class, paste(levels(group), collapse = ", ")))
    other <- setdiff(levels(group), positive_class)
    group <- factor(group, levels = c(positive_class, other))
  }
  levels(group) <- make.names(levels(group))

  if (!is.null(genes)) {
    genes <- intersect(genes, rownames(expr_mat))
    if (length(genes) == 0L)
      stop("None of the supplied genes found in rownames(expr_mat).")
    expr_mat <- expr_mat[genes, , drop = FALSE]
  }

  rv <- apply(expr_mat, 1, stats::var, na.rm = TRUE)
  expr_mat <- expr_mat[rv > 0, , drop = FALSE]
  dat <- as.data.frame(t(expr_mat))
  colnames(dat) <- make.names(colnames(dat), unique = TRUE)

  ## Mode C: external validation
  if (!is.null(test_expr)) {
    if (is.null(test_group))
      stop("test_group is required when test_expr is provided.")
    test_expr <- as.matrix(test_expr)
    test_group <- factor(test_group, levels = levels(group))
    common <- intersect(rownames(expr_mat), rownames(test_expr))
    if (length(common) == 0L)
      stop("No shared genes between train and test data.")
    dat <- as.data.frame(t(expr_mat[common, , drop = FALSE]))
    test_dat <- as.data.frame(t(test_expr[common, , drop = FALSE]))
    colnames(dat) <- make.names(colnames(dat), unique = TRUE)
    colnames(test_dat) <- make.names(colnames(test_dat), unique = TRUE)

    out <- list(train_x = dat, train_y = group,
                test_x = test_dat, test_y = test_group,
                gene_names = colnames(dat), levels = levels(group),
                full_cv = FALSE)
    mode_msg <- sprintf("Mode C: Train %d | Test %d | %d genes",
                        nrow(dat), nrow(test_dat), ncol(dat))

  ## Mode B: internal split
  } else if (isTRUE(split)) {
    if (!is.null(train_idx)) {
      idx <- train_idx
    } else {
      set.seed(seed)
      idx <- caret::createDataPartition(group, p = train_ratio, list = FALSE)[, 1]
    }
    out <- list(train_x = dat[idx, , drop = FALSE],
                train_y = group[idx],
                test_x  = dat[-idx, , drop = FALSE],
                test_y  = group[-idx],
                gene_names = colnames(dat), levels = levels(group),
                full_cv = FALSE)
    mode_msg <- sprintf("Mode B: Train %d | Test %d | %d genes",
                        nrow(out$train_x), nrow(out$test_x), ncol(dat))

  ## Mode A: full-data CV (default)
  } else {
    out <- list(train_x = dat, train_y = group,
                test_x = NULL, test_y = NULL,
                gene_names = colnames(dat), levels = levels(group),
                full_cv = TRUE)
    mode_msg <- sprintf("Mode A: %d samples | %d genes",
                        nrow(dat), ncol(dat))
  }

  out$class_balance <- .class_balance_summary(
    out$train_y,
    threshold = imbalance_threshold
  )
  out$test_class_balance <- if (!is.null(out$test_y)) {
    .class_balance_summary(out$test_y, threshold = imbalance_threshold)
  } else {
    NULL
  }
  if (isTRUE(out$class_balance$is_imbalanced) && imbalance_action != "none") {
    condition_message <- sprintf(
      "Class imbalance detected in training data: %s (minority fraction %.3f < %.3f).",
      paste(names(out$class_balance$counts), out$class_balance$counts,
            sep = "=", collapse = ", "),
      out$class_balance$minority_fraction,
      imbalance_threshold
    )
    if (imbalance_action == "warn") {
      warning(condition_message, call. = FALSE)
    } else {
      message(condition_message)
    }
  }

  class(out) <- "tcm_ml_data"
  message(mode_msg)
  return(out)
}


## Internal: check required packages for ML modules
#' @keywords internal
#' @noRd
.check_ml_deps <- function(pkgs) {
  missing <- vapply(pkgs,
                    function(p) !requireNamespace(p, quietly = TRUE),
                    logical(1))
  if (any(missing))
    stop(sprintf("Package(s) required: %s\n  install.packages(c(%s))",
                 paste(pkgs[missing], collapse = ", "),
                 paste(shQuote(pkgs[missing]), collapse = ", ")),
         call. = FALSE)
}


## Internal: caret trainControl for AUC-based repeated CV
#' @param method CV method passed to [caret::trainControl()]. Default \code{"repeatedcv"}.
#' @param number Number of CV folds. Default \code{5}.
#' @param repeats Number of CV repeats. Default \code{5}.
#' @param seed Ignored (kept for signature compatibility; callers should call
#'   [set.seed()] immediately before the \code{caret::train()} / \code{caret::rfe()} call).
#' @param class_probs Logical. Whether to compute class probabilities. Default \code{TRUE}.
#' @param summary_func Summary function for [caret::trainControl()]. \code{NULL} (default)
#'   uses \code{caret::twoClassSummary} when \code{class_probs = TRUE}, otherwise
#'   \code{caret::defaultSummary}.
#' @keywords internal
#' @noRd
.make_train_ctrl <- function(method = "repeatedcv",
                             number = 5,
                             repeats = 5,
                             seed = 2025,
                             class_probs = TRUE,
                             summary_func = NULL,
                             sampling = NULL) {
  if (is.null(summary_func)) {
    summary_func <- if (class_probs) caret::twoClassSummary else caret::defaultSummary
  }
  
  args <- list(
    method = method,
    number = number,
    classProbs = class_probs,
    summaryFunction = summary_func,
    savePredictions = "final",
    allowParallel = TRUE
  )
  
  if (method == "repeatedcv") {
    args$repeats <- repeats
  }
  if (!is.null(sampling)) args$sampling <- sampling
  
  return(do.call(caret::trainControl, args))
}


## Internal: evaluate model on hold-out test set
#' @keywords internal
#' @noRd
.eval_test <- function(model, test_x, test_y, positive_class) {
  pred_class <- stats::predict(model, newdata = test_x)
  pred_prob <- stats::predict(model, newdata = test_x, type = "prob")
  cm <- caret::confusionMatrix(pred_class, test_y, positive = positive_class)

  roc_metrics <- NULL
  if (requireNamespace("pROC", quietly = TRUE)) {
    roc_metrics <- .compute_roc_metrics(
      truth = test_y,
      probability = pred_prob[[positive_class]],
      positive_class = positive_class
    )
  }

  return(list(
    predictions = pred_class,
    probabilities = pred_prob,
    confusion = cm,
    accuracy = as.numeric(cm$overall["Accuracy"]),
    auc = if (is.null(roc_metrics)) NA_real_ else roc_metrics$auc,
    ci_lower = if (is.null(roc_metrics)) NA_real_ else roc_metrics$ci_lower,
    ci_upper = if (is.null(roc_metrics)) NA_real_ else roc_metrics$ci_upper,
    ci_method = if (is.null(roc_metrics)) NA_character_ else roc_metrics$ci_method,
    conf_level = if (is.null(roc_metrics)) NA_real_ else roc_metrics$conf_level,
    roc = if (is.null(roc_metrics)) NULL else roc_metrics$roc_object,
    sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    specificity = as.numeric(cm$byClass["Specificity"])
  ))
}

.class_balance_summary <- function(y, threshold = 0.25) {
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold <= 0 || threshold > 0.5) {
    stop("imbalance_threshold must be in (0, 0.5].", call. = FALSE)
  }
  y <- factor(y)
  counts <- table(y)
  counts <- stats::setNames(as.integer(counts), names(counts))
  total <- sum(counts)
  proportions <- if (total > 0L) counts / total else rep(NA_real_, length(counts))
  minority_index <- which.min(counts)
  majority_index <- which.max(counts)
  minority_fraction <- if (total > 0L) counts[[minority_index]] / total else NA_real_
  majority_count <- counts[[majority_index]]
  imbalance_ratio <- if (majority_count > 0L) {
    counts[[minority_index]] / majority_count
  } else {
    NA_real_
  }
  list(
    counts = counts,
    proportions = stats::setNames(as.numeric(proportions), names(counts)),
    minority_class = names(counts)[[minority_index]],
    majority_class = names(counts)[[majority_index]],
    minority_fraction = unname(minority_fraction),
    imbalance_ratio = unname(imbalance_ratio),
    threshold = threshold,
    is_imbalanced = is.finite(minority_fraction) && minority_fraction < threshold
  )
}

.balanced_case_weights <- function(y) {
  y <- factor(y)
  counts <- table(y)
  if (any(counts == 0L)) {
    stop("Case weights require at least one sample in each class.", call. = FALSE)
  }
  raw <- length(y) / (length(counts) * counts)
  weights <- unname(raw[as.character(y)])
  weights / mean(weights)
}

.downsample_training_data <- function(x, y, seed = 2025L) {
  y <- factor(y)
  indices <- split(seq_along(y), y)
  target_n <- min(lengths(indices))
  if (target_n < 1L) {
    stop("Downsampling requires at least one sample in each class.", call. = FALSE)
  }
  set.seed(seed)
  selected <- unlist(lapply(indices, function(idx) sample(idx, target_n)),
                     use.names = FALSE)
  selected <- sort(selected)
  list(
    x = x[selected, , drop = FALSE],
    y = factor(y[selected], levels = levels(y)),
    index = selected,
    counts_before = .class_balance_summary(y)$counts,
    counts_after = .class_balance_summary(y[selected])$counts
  )
}

.resolve_balance_method <- function(balance_method,
                                    y,
                                    model,
                                    imbalance_threshold = 0.25) {
  requested <- match.arg(
    balance_method,
    c("none", "auto", "weights", "down")
  )
  model <- match.arg(model, c("glmnet", "rf", "svm_rfe", "xgboost"))
  if (requested == "weights" && model == "svm_rfe") {
    stop(
      "ml_svm_rfe() does not support balance_method = 'weights'; use 'down' or 'auto'.",
      call. = FALSE
    )
  }
  if (requested != "auto") return(requested)
  summary <- .class_balance_summary(y, threshold = imbalance_threshold)
  if (!summary$is_imbalanced) return("none")
  if (model == "svm_rfe") "down" else "weights"
}

.prepare_model_balance <- function(ml_data,
                                   balance_method,
                                   model,
                                   seed = 2025L) {
  threshold <- if (!is.null(ml_data$class_balance$threshold)) {
    ml_data$class_balance$threshold
  } else {
    0.25
  }
  requested <- match.arg(
    balance_method,
    c("none", "auto", "weights", "down")
  )
  applied <- .resolve_balance_method(
    requested,
    ml_data$train_y,
    model = model,
    imbalance_threshold = threshold
  )
  x <- ml_data$train_x
  y <- ml_data$train_y
  counts_before <- .class_balance_summary(y, threshold)$counts
  case_weights <- NULL
  class_weights <- NULL
  scope <- "training_data"

  if (applied == "down" && model != "svm_rfe") {
    down <- .downsample_training_data(x, y, seed = seed)
    x <- down$x
    y <- down$y
  } else if (applied == "weights") {
    case_weights <- .balanced_case_weights(y)
    class_weights <- tapply(case_weights, y, mean)
  } else if (applied == "down" && model == "svm_rfe") {
    scope <- "resample_training_folds"
  }

  counts_after <- if (applied == "down" && model == "svm_rfe") {
    stats::setNames(rep(min(counts_before), length(counts_before)), names(counts_before))
  } else {
    .class_balance_summary(y, threshold)$counts
  }
  positive <- ml_data$levels[[1]]
  negative <- ml_data$levels[[2]]
  scale_pos_weight <- if (applied == "weights") {
    unname(counts_before[[negative]] / counts_before[[positive]])
  } else {
    NULL
  }

  list(
    x = x,
    y = factor(y, levels = ml_data$levels),
    case_weights = case_weights,
    class_weights = class_weights,
    scale_pos_weight = scale_pos_weight,
    info = list(
      requested = requested,
      applied = applied,
      scope = scope,
      train_counts_before = counts_before,
      train_counts_after = counts_after,
      class_weights = class_weights,
      case_weights = case_weights,
      seed = seed
    )
  )
}

.stratified_foldid <- function(y, nfolds, seed = 2025L) {
  y <- factor(y)
  nfolds <- as.integer(nfolds)
  if (nfolds < 2L) stop("At least 2 folds are required.", call. = FALSE)
  set.seed(seed)
  foldid <- integer(length(y))
  for (level in levels(y)) {
    idx <- which(y == level)
    foldid[idx] <- sample(rep(seq_len(nfolds), length.out = length(idx)))
  }
  foldid
}

.compute_roc_metrics <- function(truth,
                                 probability,
                                 positive_class,
                                 ci_method = c("delong", "bootstrap"),
                                 conf_level = 0.95,
                                 boot_n = 2000L,
                                 seed = 2025L) {
  .check_ml_deps("pROC")
  ci_method <- match.arg(ci_method)
  if (!is.numeric(conf_level) || length(conf_level) != 1L ||
      !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1) {
    stop("conf_level must be between 0 and 1.", call. = FALSE)
  }
  truth <- factor(truth)
  if (nlevels(truth) != 2L || !positive_class %in% levels(truth)) {
    stop("truth must contain two classes including positive_class.", call. = FALSE)
  }
  probability <- suppressWarnings(as.numeric(probability))
  keep <- !is.na(truth) & is.finite(probability)
  truth <- factor(truth[keep], levels = levels(truth))
  probability <- probability[keep]
  negative_class <- setdiff(levels(truth), positive_class)
  counts <- table(truth)
  if (any(counts < 2L)) {
    stop("At least two complete samples per class are required for ROC inference.",
         call. = FALSE)
  }

  roc_object <- pROC::roc(
    response = truth,
    predictor = probability,
    levels = c(negative_class, positive_class),
    direction = "auto",
    quiet = TRUE
  )
  if (ci_method == "bootstrap") set.seed(seed)
  ci_args <- list(
    roc = roc_object,
    conf.level = conf_level,
    method = if (ci_method == "bootstrap") "bootstrap" else "delong"
  )
  if (ci_method == "bootstrap") {
    ci_args$boot.n <- as.integer(boot_n)
    ci_args$boot.stratified <- TRUE
  }
  ci_auc <- suppressWarnings(do.call(pROC::ci.auc, ci_args))
  best <- pROC::coords(
    roc_object,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  if (is.data.frame(best) || is.matrix(best)) best <- best[1L, , drop = FALSE]
  value <- function(name) as.numeric(best[[name]][[1L]])

  list(
    auc = as.numeric(pROC::auc(roc_object)),
    ci_lower = as.numeric(ci_auc[[1L]]),
    ci_upper = as.numeric(ci_auc[[3L]]),
    ci_method = ci_method,
    conf_level = conf_level,
    threshold = value("threshold"),
    sensitivity = value("sensitivity"),
    specificity = value("specificity"),
    direction = roc_object$direction,
    n_positive = unname(as.integer(counts[[positive_class]])),
    n_negative = unname(as.integer(counts[[negative_class]])),
    roc_object = roc_object
  )
}


## Internal: S3 constructor
#' @keywords internal
#' @noRd
.new_tcm_ml <- function(method, model, importance, selected_features,
                        cv_performance, test_performance = NULL,
                        ml_data) {
  return(structure(
    list(method = method,
         model = model,
         importance = importance,
         selected_features = selected_features,
         cv_performance = cv_performance,
         test_performance = test_performance,
         ml_data = ml_data),
    class = "tcm_ml"
  ))
}


#' @export
print.tcm_ml <- function(x, ...) {
  tag <- if (!is.null(x$test_performance)) "Test" else "CV"
  perf <- if (tag == "Test") x$test_performance else x$cv_performance
  .pv <- function(v) if (is.null(v) || is.na(v)) "NA" else sprintf("%.4f", v)

  cat(sprintf("=== tcm_ml: %s ===\n", toupper(x$method)))
  cat(sprintf("  Features: %d | %s AUC: %s | Sens: %s | Spec: %s\n",
              length(x$selected_features),
              tag, .pv(perf$auc), .pv(perf$sensitivity), .pv(perf$specificity)))
  invisible(x)
}

#' @export
print.tcm_ml_list <- function(x, ...) {
  cat(sprintf("=== ML Screening (%d methods) ===\n", length(x)))
  .pv <- function(v) if (is.null(v) || is.na(v)) "NA" else sprintf("%.4f", v)
  for (m in x) {
    perf <- if (!is.null(m$test_performance)) m$test_performance else m$cv_performance
    cat(sprintf("  [%s] %d features | AUC = %s\n",
                toupper(m$method), length(m$selected_features), .pv(perf$auc)))
  }
  invisible(x)
}


#' Summarise ML screening results
#'
#' Returns a tidy data.frame with one row per method, showing the number of
#' selected features and performance metrics (CV or Test AUC, Sensitivity,
#' Specificity).
#'
#' @param object A `tcm_ml_list` produced by [run_ml_screening()].
#' @param ... Ignored.
#'
#' @return A data.frame with columns: `method`, `n_features`, `auc_type`,
#'   `auc`, `auc_sd`, `sensitivity`, `specificity`.
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res_list <- run_ml_screening(ml_data)
#'   summary(res_list)
#' }
#' @export
summary.tcm_ml_list <- function(object, ...) {
  .safe <- function(v) if (is.null(v) || length(v) == 0 || is.na(v)) NA_real_ else v

  rows <- lapply(object, function(m) {
    has_test <- !is.null(m$test_performance)
    perf <- if (has_test) m$test_performance else m$cv_performance
    cv_sd <- if (!has_test && !is.null(m$cv_performance$auc_sd))
               m$cv_performance$auc_sd else NA_real_

    data.frame(
      method = m$method,
      n_features = length(m$selected_features),
      auc_type = if (has_test) "Test" else "CV",
      auc = .safe(perf$auc),
      auc_sd = cv_sd,
      sensitivity = .safe(perf$sensitivity),
      specificity = .safe(perf$specificity),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


#' Re-select top features from a fitted ML model
#' After running a model (especially Ridge or XGBoost where all features
#' are initially retained), inspect \code{result$importance} to decide
#' how many to keep, then call this function to trim the selection.
#'
#' @param ml_obj A \code{tcm_ml} object from any \code{ml_*} function.
#' @param top_n Integer; the number of top features to keep (ranked by importance).
#' @return A modified \code{tcm_ml} object with updated
#'   \code{$selected_features} and \code{$genes}.
#' @details The underlying model and importance table are unchanged.
#'   Only the gene selection is trimmed. Performance metrics still
#'   reflect the original run; re-run the model with \code{top_n}
#'   if you need updated metrics.
#'
#' @examples
#' \dontrun{
#'   xgb <- ml_xgboost(ml_data)        # keep all features as default
#'   head(xgb$importance, 20)           # inspect ranking
#'   xgb <- select_features(xgb, 15)   # keep top 15
#' }
#' @export
select_features <- function(ml_obj, top_n) {
  stopifnot(inherits(ml_obj, "tcm_ml"))
  stopifnot(is.numeric(top_n), length(top_n) == 1L, top_n >= 1)

  imp <- ml_obj$importance
  if (is.null(imp) || nrow(imp) == 0L)
    stop("No importance data found in this model object.")

  top_n <- min(as.integer(top_n), nrow(imp))
  new_genes <- imp$gene[seq_len(top_n)]

  ml_obj$selected_features <- new_genes
  ml_obj$genes <- new_genes

  message(sprintf("[%s] Re-selected top %d / %d features",
                  toupper(ml_obj$method), top_n, nrow(imp)))
  return(ml_obj)
}
