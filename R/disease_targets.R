# Open Targets Platform GraphQL API utilities
# API endpoint: https://api.platform.opentargets.org/api/v4/graphql

.OT_GRAPHQL_URL <- "https://api.platform.opentargets.org/api/v4/graphql"

#' @keywords internal
#' @noRd
.check_httr <- function() {
  if (!requireNamespace("httr", quietly = TRUE))
    stop("Package 'httr' is required. Install with: install.packages('httr')", call. = FALSE)
}

#' @keywords internal
#' @noRd
.ot_post <- function(query, variables = list()) {
  .check_httr()
  resp <- httr::POST(
    .OT_GRAPHQL_URL,
    httr::content_type_json(),
    body = jsonlite::toJSON(list(query = query, variables = variables), auto_unbox = TRUE)
  )
  httr::stop_for_status(resp)
  httr::content(resp, as = "parsed", simplifyVector = TRUE)
}


#' Search Disease EFO IDs from Open Targets Platform
#'
#' Queries the Open Targets Platform GraphQL API to search for diseases
#' matching a given name and returns their EFO/MONDO/Orphanet IDs.
#'
#' @param disease_name A character string. Disease name to search for
#'   (English, e.g. \code{"breast cancer"}).
#' @param size An integer. Maximum number of results to return (default \code{10}).
#'
#' @return A \code{data.frame} with columns \code{id} (EFO/MONDO/Orphanet ID)
#'   and \code{name} (disease name). Returns an empty \code{data.frame} if no
#'   match is found.
#'
#' @examples
#' \dontrun{
#' search_disease_efo("breast cancer")
#' search_disease_efo("diabetes", size = 5)
#' }
#'
#' @export
#'
search_disease_efo <- function(disease_name, size = 10) {
  query <- '
  query SearchDisease($query: String!, $size: Int!) {
    search(queryString: $query, entityNames: ["disease"], page: {index: 0, size: $size}) {
      hits {
        id
        name
        entity
      }
    }
  }'

  parsed <- .ot_post(query, list(query = disease_name, size = as.integer(size)))
  hits <- parsed$data$search$hits

  if (is.null(hits) || length(hits) == 0)
    return(data.frame(id = character(), name = character(), stringsAsFactors = FALSE))

  result <- hits[hits$entity == "disease", c("id", "name"), drop = FALSE]
  rownames(result) <- NULL
  result
}


#' Retrieve Disease-Associated Targets from Open Targets Platform
#'
#' Queries the Open Targets Platform GraphQL API for targets associated with
#' a disease, identified by its EFO/MONDO/Orphanet ID.
#'
#' @param efo_id A character string. Disease ontology ID (e.g.
#'   \code{"EFO_0000305"} for breast carcinoma, \code{"EFO_0000400"} for
#'   diabetes mellitus).
#' @param size An integer. Maximum number of targets to return (default \code{200}).
#' @param score_threshold A numeric value between 0 and 1. Only targets with an
#'   overall association score at or above this threshold are returned
#'   (default \code{0}).
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{ensembl_id}{Ensembl gene ID.}
#'     \item{gene_symbol}{HGNC-approved gene symbol.}
#'     \item{gene_name}{Full gene name.}
#'     \item{biotype}{Gene biotype (e.g. \code{"protein_coding"}).}
#'     \item{score}{Overall association score (0-1) from Open Targets.}
#'   }
#'   Rows are sorted by \code{score} in descending order.
#'
#' @examples
#' \dontrun{
#' # Breast carcinoma targets with score >= 0.5
#' get_disease_targets("EFO_0000305", size = 100, score_threshold = 0.5)
#'
#' # All diabetes mellitus targets
#' get_disease_targets("EFO_0000400", size = 500)
#' }
#'
#' @export
#'
get_disease_targets <- function(efo_id, size = 200, score_threshold = 0) {
  query <- '
  query DiseaseTargets($efoId: String!, $size: Int!) {
    disease(efoId: $efoId) {
      id
      name
      associatedTargets(page: {index: 0, size: $size}) {
        count
        rows {
          target {
            id
            approvedSymbol
            approvedName
            biotype
          }
          score
        }
      }
    }
  }'

  parsed <- .ot_post(query, list(efoId = efo_id, size = as.integer(size)))
  disease_data <- parsed$data$disease

  if (is.null(disease_data))
    stop("No disease found for EFO ID: ", efo_id, call. = FALSE)

  rows <- disease_data$associatedTargets$rows
  if (is.null(rows) || nrow(rows) == 0)
    return(data.frame(
      ensembl_id = character(), gene_symbol = character(),
      gene_name  = character(), biotype = character(),
      score      = numeric(),
      stringsAsFactors = FALSE
    ))

  result <- data.frame(
    ensembl_id  = rows$target$id,
    gene_symbol = rows$target$approvedSymbol,
    gene_name   = rows$target$approvedName,
    biotype     = rows$target$biotype,
    score       = rows$score,
    stringsAsFactors = FALSE
  )

  result <- result[result$score >= score_threshold, , drop = FALSE]
  result <- result[order(-result$score), , drop = FALSE]
  rownames(result) <- NULL
  result
}


#' Query Disease Targets by Disease Name from Open Targets Platform
#'
#' A convenience wrapper that combines disease name search and target retrieval.
#' Searches for a disease by name, selects a match, and returns its associated
#' targets with overall association scores.
#'
#' @param disease_name A character string. Disease name to query
#'   (English, e.g. \code{"breast cancer"}, \code{"diabetes"}).
#' @param size An integer. Maximum number of targets to return (default \code{200}).
#' @param score_threshold A numeric value between 0 and 1. Minimum association
#'   score for included targets (default \code{0}).
#' @param efo_index An integer. When the search returns multiple disease matches,
#'   specifies which result to use (default \code{1}, i.e. the top hit). Use
#'   \code{search_disease_efo()} to inspect all candidates first.
#'
#' @return A \code{data.frame} with columns \code{ensembl_id}, \code{gene_symbol},
#'   \code{gene_name}, \code{biotype}, and \code{score}, sorted by \code{score}
#'   descending. Returns \code{NULL} if the disease name yields no results.
#'
#' @seealso \code{\link{search_disease_efo}}, \code{\link{get_disease_targets}}
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' targets <- query_disease_targets("breast cancer")
#' head(targets)
#'
#' # Stricter filter, more results
#' targets <- query_disease_targets("diabetes", size = 500, score_threshold = 0.1)
#'
#' # Inspect candidates, then select the correct one
#' search_disease_efo("diabetes")
#' targets <- query_disease_targets("diabetes", efo_index = 2)
#'
#' # Extract gene symbol vector for downstream analysis
#' gene_set <- targets$gene_symbol
#' }
#'
#' @export
#'
query_disease_targets <- function(disease_name, size = 200, score_threshold = 0, efo_index = 1) {
  diseases <- search_disease_efo(disease_name)

  if (nrow(diseases) == 0) {
    message("No matching disease found for: ", disease_name)
    return(NULL)
  }

  if (efo_index > nrow(diseases))
    stop("efo_index (", efo_index, ") exceeds number of matches found (",
         nrow(diseases), ").", call. = FALSE)

  message(sprintf("Using disease: %s (%s)  [match %d of %d]",
    diseases$name[efo_index], diseases$id[efo_index],
    efo_index, nrow(diseases)))

  get_disease_targets(diseases$id[efo_index], size = size, score_threshold = score_threshold)
}


# Disease-target search (DisGeNET via DOSE)
# Package-level cache
.disease_env <- new.env(parent = emptyenv())

#' @keywords internal
#' @noRd
.check_dose <- function() {
  if (!requireNamespace("DOSE", quietly = TRUE))
    stop("Package 'DOSE' is required. Install: BiocManager::install('DOSE')", call. = FALSE)
}

#' @keywords internal
#' @noRd
.get_disease_data <- function() {
  if (!is.null(.disease_env$df)) return(.disease_env$df)
  .check_dose()

  env <- new.env(parent = emptyenv())
  utils::data("DGN_PATHID2EXTID", package = "DOSE", envir = env)
  utils::data("DGN_PATHID2NAME", package = "DOSE", envir = env)

  id2gene <- env$DGN_PATHID2EXTID
  id2name <- env$DGN_PATHID2NAME

  df <- data.frame(
    disease_id = rep(names(id2gene), lengths(id2gene)),
    gene_id = as.character(unlist(id2gene, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  df$disease_name <- id2name[df$disease_id]

  .disease_env$df <- df
  .disease_env$id2name <- id2name
  df
}

#' @keywords internal
#' @noRd
.add_gene_symbol <- function(df) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("Package 'org.Hs.eg.db' required. Install: BiocManager::install('org.Hs.eg.db')", call. = FALSE)
  mapping <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                          keys = unique(df$gene_id),
                          columns = "SYMBOL", keytype = "ENTREZID")
  )
  sym_map <- stats::setNames(mapping$SYMBOL, mapping$ENTREZID)
  df$symbol <- sym_map[df$gene_id]
  df
}

#' @keywords internal
#' @noRd
.symbol_to_entrez <- function(symbols) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("Package 'org.Hs.eg.db' required. Install: BiocManager::install('org.Hs.eg.db')", call. = FALSE)
  suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                          keys = symbols, columns = "ENTREZID", keytype = "SYMBOL")
  )
}

#' Search disease targets (disease -> genes)
#'
#' Query DisGeNET (via DOSE) to find genes associated with a disease.
#' Supports UMLS CUI IDs, exact name match, or fuzzy name search.
#'
#' @param disease Character. Disease name or UMLS CUI (e.g. "sepsis" or
#'   "C0243026"). Supports a vector for multiple diseases.
#' @param readable Logical. Convert Entrez IDs to gene symbols (default TRUE).
#'
#' @return A \code{data.frame} with columns \code{disease_id}, \code{disease_name},
#'   \code{gene_id}, and optionally \code{symbol}. Returns NULL if no match.
#'
#' @examples
#' \dontrun{
#'   search_disease("sepsis")
#'   search_disease(c("sepsis", "asthma"))
#'   search_disease("C0243026")
#' }
#' @export
search_disease <- function(disease, readable = TRUE) {
  df <- .get_disease_data()
  id2name <- .disease_env$id2name

  # Match: exact ID -> exact name (case-insensitive) -> grep
  all_matched_ids <- character(0)
  for (q in disease) {
    if (q %in% names(id2name)) {
      all_matched_ids <- c(all_matched_ids, q)
    } else {
      exact <- names(id2name)[tolower(id2name) == tolower(q)]
      if (length(exact) > 0) {
        all_matched_ids <- c(all_matched_ids, exact)
      } else {
        fuzzy <- names(id2name)[grepl(q, id2name, ignore.case = TRUE)]
        if (length(fuzzy) > 0) {
          all_matched_ids <- c(all_matched_ids, fuzzy)
        } else {
          message("No disease found for: ", q)
        }
      }
    }
  }

  if (length(all_matched_ids) == 0) return(NULL)

  res <- df[df$disease_id %in% all_matched_ids, , drop = FALSE]
  if (readable) res <- .add_gene_symbol(res)
  res <- res[order(res$disease_name, res$gene_id), , drop = FALSE]
  rownames(res) <- NULL
  res
}


#' Search gene-associated diseases (gene -> diseases)
#'
#' Reverse lookup: given gene symbols or Entrez IDs, find associated diseases
#' from DisGeNET (via DOSE).
#'
#' @param gene Character. Gene symbols (e.g. "TNF") or Entrez IDs.
#'   Supports a vector for multiple genes.
#' @param readable Logical. Attach gene symbol column (default TRUE).
#'
#' @return A \code{data.frame} with columns \code{disease_id}, \code{disease_name},
#'   \code{gene_id}, and optionally \code{symbol}. Returns NULL if no match.
#'
#' @examples
#' \dontrun{
#'   search_gene_disease("TNF")
#'   search_gene_disease(c("IL6", "TNF", "PPARG"))
#'   search_gene_disease("7124")
#' }
#' @export
search_gene_disease <- function(gene, readable = TRUE) {
  df <- .get_disease_data()

  # Try as Entrez ID first
  hits <- df[df$gene_id %in% gene, , drop = FALSE]

  # If no hits, try Symbol -> Entrez conversion
  if (nrow(hits) == 0) {
    gene_map <- .symbol_to_entrez(gene)
    entrez_ids <- gene_map$ENTREZID[!is.na(gene_map$ENTREZID)]
    if (length(entrez_ids) > 0) {
      hits <- df[df$gene_id %in% entrez_ids, , drop = FALSE]
    }
  }

  if (nrow(hits) == 0) {
    message("No diseases found for: ", paste(gene, collapse = ", "))
    return(NULL)
  }

  if (readable) hits <- .add_gene_symbol(hits)
  hits <- hits[order(hits$gene_id, hits$disease_name), , drop = FALSE]
  rownames(hits) <- NULL
  hits
}