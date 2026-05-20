#!/usr/bin/env Rscript
# Generate a compact Markdown report from artifacts available in the current
# TCMDATA R session. This script is best sourced from R after an agent run.
#
# Usage:
#   source("inst/skills/tcm-network-pharmacology/scripts/generate_report.R")
#   generate_tcm_artifact_report("network_pharmacology_report.md")

generate_tcm_artifact_report <- function(output_file = "network_pharmacology_report.md",
                                         max_lines = 8L) {
  if (!requireNamespace("TCMDATA", quietly = TRUE)) {
    stop("TCMDATA package is required.", call. = FALSE)
  }

  artifacts <- TCMDATA::list_tcm_artifacts()
  if (is.null(artifacts) || nrow(artifacts) == 0L) {
    stop("No artifacts found in the current TCMDATA session.", call. = FALSE)
  }

  lines <- c(
    "# Network Pharmacology Analysis Report",
    "",
    paste("Generated:", as.character(Sys.time())),
    "",
    "## Artifact Summary",
    ""
  )

  for (i in seq_len(nrow(artifacts))) {
    artifact_id <- artifacts$artifact_id[i]
    artifact_type <- artifacts$artifact_type[i]
    lines <- c(
      lines,
      paste0("### ", artifact_id, " (", artifact_type, ")"),
      "",
      "```text",
      TCMDATA::summarize_tcm_artifact(artifact_id, max_lines = max_lines),
      "```",
      ""
    )
  }

  writeLines(lines, output_file)
  message(sprintf("Report written to: %s", normalizePath(output_file, mustWork = FALSE)))
  invisible(output_file)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output_file <- if (length(args) >= 1L) args[[1L]] else "network_pharmacology_report.md"
  generate_tcm_artifact_report(output_file)
}
