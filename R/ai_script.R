# ai_script.R
# Deterministic replay scripts generated from recorded aisdk tool calls.

#' Print a TCM analysis script
#'
#' @param x A \code{tcm_analysis_script} object.
#' @param ... Additional arguments (ignored).
#'
#' @return \code{x}, invisibly.
#' @export
print.tcm_analysis_script <- function(x, ...) {
  cat(as.character(x), "\n")
  invisible(x)
}

.build_tcm_analysis_script <- function(result,
                                       task = NULL,
                                       turn = NULL,
                                       model = NULL) {
  tool_calls <- .tcm_result_tool_calls(result)
  tool_results <- result$all_tool_results %||% list()
  if (length(tool_calls) == 0L) {
    return(NULL)
  }

  generated_at <- Sys.time()
  task_lines <- strsplit(as.character(task %||% ""), "\n", fixed = TRUE)[[1L]]
  task_comment <- if (length(task_lines) > 0L && any(nzchar(task_lines))) {
    paste0("# Task: ", task_lines)
  } else {
    "# Task: not recorded"
  }

  lines <- c(
    "# TCMDATA agent analysis replay script",
    sprintf("# Generated: %s", format(generated_at, "%Y-%m-%d %H:%M:%S %z")),
    if (!is.null(turn)) sprintf("# Chat turn: %s", turn) else NULL,
    if (!is.null(model) && nzchar(model)) sprintf("# Model: %s", model) else NULL,
    task_comment,
    "# Generated deterministically from the executed tool-call trace.",
    "",
    "library(TCMDATA)",
    "",
    ".tcm_tools <- create_tcm_tools()",
    ".tcm_get_tool <- function(name) {",
    "  matches <- Filter(function(x) identical(x$name, name), .tcm_tools)",
    "  if (length(matches) != 1L) stop(sprintf(\"Tool '%s' not found.\", name))",
    "  matches[[1L]]",
    "}",
    ""
  )

  artifact_variables <- character(0)
  step_variables <- character(length(tool_calls))

  for (i in seq_along(tool_calls)) {
    tool_call <- tool_calls[[i]]
    tool_name <- .tcm_tool_call_name(tool_call)
    arguments <- .tcm_tool_call_arguments(tool_call)
    step_variable <- sprintf("step_%02d", i)
    step_variables[[i]] <- step_variable

    lines <- c(
      lines,
      sprintf("# Step %d: %s", i, tool_name),
      sprintf(
        "%s <- .tcm_get_tool(%s)$run(",
        step_variable,
        .deparse_tcm_script_value(tool_name)
      ),
      .render_tcm_script_arguments(arguments, artifact_variables, indent = 2L),
      ")"
    )

    tool_result <- .match_tcm_tool_result(tool_call, tool_results, i)
    artifact_id <- .tcm_tool_result_artifact_id(tool_result)
    if (!is.null(artifact_id) && nzchar(artifact_id)) {
      artifact_variables[[artifact_id]] <- paste0(step_variable, "$artifact_id")
      lines <- c(lines, sprintf("# Output artifact: %s", artifact_id))
    }
    if (is.list(tool_result) && isTRUE(tool_result$is_error)) {
      lines <- c(lines, "# Original execution status: failed")
    }
    lines <- c(lines, "")
  }

  lines <- c(
    lines,
    sprintf(
      "# Step results are available as: %s",
      paste(step_variables, collapse = ", ")
    )
  )
  script <- paste(lines, collapse = "\n")

  structure(
    script,
    class = c("tcm_analysis_script", "character"),
    task = task,
    turn = turn,
    model = model,
    generated_at = generated_at,
    tool_calls = tool_calls,
    tool_results = tool_results,
    step_variables = step_variables
  )
}

.export_tcm_analysis_script <- function(script, envir = globalenv()) {
  if (is.null(script)) {
    return(NULL)
  }
  script_name <- .next_tcm_script_name(envir)
  attr(script, "global_name") <- script_name
  assign(script_name, script, envir = envir)
  script_name
}

.next_tcm_script_name <- function(envir = globalenv()) {
  existing <- ls(envir = envir, pattern = "^tcm_script_[0-9]+$")
  existing_numbers <- suppressWarnings(as.integer(sub(
    "^tcm_script_", "", existing
  )))
  existing_numbers <- existing_numbers[!is.na(existing_numbers)]
  next_number <- if (length(existing_numbers) == 0L) {
    1L
  } else {
    max(existing_numbers) + 1L
  }
  sprintf("tcm_script_%03d", next_number)
}

.tcm_tool_call_name <- function(tool_call) {
  name <- tool_call$name %||% tool_call$tool %||% NULL
  if (is.null(name) || length(name) != 1L || !nzchar(name)) {
    stop("A recorded tool call is missing its tool name.", call. = FALSE)
  }
  as.character(name)
}

.tcm_tool_call_arguments <- function(tool_call) {
  arguments <- tool_call$arguments %||% tool_call$args %||% list()
  if (is.character(arguments) && length(arguments) == 1L) {
    arguments <- tryCatch(
      jsonlite::fromJSON(arguments, simplifyVector = FALSE),
      error = function(e) arguments
    )
  }
  if (!is.list(arguments)) {
    stop("Recorded tool arguments must be a list or JSON object.", call. = FALSE)
  }
  arguments
}

.match_tcm_tool_result <- function(tool_call, tool_results, index) {
  if (length(tool_results) == 0L) {
    return(NULL)
  }
  call_id <- tool_call$id %||% tool_call$call_id %||% NULL
  if (!is.null(call_id)) {
    result_ids <- vapply(tool_results, function(result) {
      as.character(result$id %||% result$call_id %||% "")
    }, character(1))
    matched <- which(result_ids == as.character(call_id))
    if (length(matched) > 0L) {
      return(tool_results[[matched[[1L]]]])
    }
  }
  if (index <= length(tool_results)) tool_results[[index]] else NULL
}

.tcm_tool_result_artifact_id <- function(tool_result) {
  if (!is.list(tool_result)) {
    return(NULL)
  }
  raw_result <- tool_result$raw_result %||% tool_result$result %||% NULL
  if (is.character(raw_result) && length(raw_result) == 1L) {
    raw_result <- tryCatch(
      jsonlite::fromJSON(raw_result, simplifyVector = FALSE),
      error = function(e) NULL
    )
  }
  artifact_id <- if (is.list(raw_result)) raw_result$artifact_id else NULL
  if (is.null(artifact_id) || length(artifact_id) != 1L) {
    return(NULL)
  }
  as.character(artifact_id)
}

.render_tcm_script_arguments <- function(arguments,
                                         artifact_variables,
                                         indent = 0L) {
  prefix <- strrep(" ", indent)
  if (length(arguments) == 0L) {
    return(paste0(prefix, "list()"))
  }

  argument_names <- names(arguments)
  if (is.null(argument_names) || any(!nzchar(argument_names))) {
    rendered <- .deparse_tcm_script_value(arguments)
    return(.indent_tcm_script_text(rendered, indent))
  }

  entries <- vapply(seq_along(arguments), function(i) {
    value <- .render_tcm_script_value(
      arguments[[i]],
      artifact_variables = artifact_variables
    )
    value <- .indent_tcm_script_text(value, indent + 4L)
    value <- sub(paste0("^", strrep(" ", indent + 4L)), "", value)
    sprintf(
      "%s%s = %s",
      strrep(" ", indent + 2L),
      .render_tcm_script_name(argument_names[[i]]),
      value
    )
  }, character(1))

  c(
    paste0(prefix, "list("),
    paste0(entries, ifelse(seq_along(entries) < length(entries), ",", "")),
    paste0(prefix, ")")
  )
}

.render_tcm_script_value <- function(value, artifact_variables) {
  if (is.character(value) && length(value) == 1L &&
      value %in% names(artifact_variables)) {
    return(unname(artifact_variables[[value]]))
  }
  .deparse_tcm_script_value(value)
}

.deparse_tcm_script_value <- function(value) {
  paste(deparse(value, width.cutoff = 500L), collapse = "\n")
}

.render_tcm_script_name <- function(name) {
  if (identical(make.names(name), name) &&
      !name %in% c("if", "else", "repeat", "while", "function", "for",
                   "in", "next", "break", "TRUE", "FALSE", "NULL",
                   "Inf", "NaN", "NA")) {
    return(name)
  }
  paste0("`", gsub("`", "\\`", name, fixed = TRUE), "`")
}

.indent_tcm_script_text <- function(text, indent) {
  prefix <- strrep(" ", indent)
  paste0(prefix, gsub("\n", paste0("\n", prefix), text, fixed = TRUE))
}
