#' @keywords internal
.pug_base <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/" 


#' Generic PUG-REST request with retry and backoff
#'
#' A lightweight helper to send GET or POST requests to the PubChem PUG-REST API with automatic retry (exponential backoff) and JSON parsing.
#' @param path Character scalar. API endpoint path under
#'   \code{https://pubchem.ncbi.nlm.nih.gov/rest/pug/},
#'   e.g. \code{"compound/cid/2244/property/MolecularWeight/JSON"}.
#' @param query Named list of query parameters to be appended to the URL.
#'   Used in both GET and POST requests.
#' @param body Named list of form data to send in the POST body (ignored for GET).
#' @param method HTTP method to use, one of \code{"GET"} or \code{"POST"}.
#' @param pause Numeric. Base pause (seconds) between retries.
#' @param max_times Integer. Maximum number of retry attempts before giving up.
#' 
#' @importFrom httr GET POST status_code content user_agent
#' @importFrom jsonlite fromJSON
#' @keywords internal

pug_request <- function(path, query = list(), body = NULL, method = c("GET","POST"),
                        pause = 0.25, max_times = 5) {
  method <- match.arg(method)
  url <- paste0(.pug_base, sub("^/+","", path))
  do_req <- function() {
    if (method == "GET") {
      httr::GET(url, query = query, httr::user_agent("pubchem-helper/0.1"))
    } else {
      httr::POST(url, query = query, body = body, encode = "form",
                 httr::user_agent("pubchem-helper/0.1"))
    }
  }
  delay <- 0.5
  for (i in seq_len(max_times)) {
    resp <- do_req()
    if (httr::status_code(resp) == 200) {
      txt <- httr::content(resp, as = "text", encoding = "UTF-8")
      return(jsonlite::fromJSON(txt, simplifyVector = TRUE))
    }
    Sys.sleep(min(delay, 4)); delay <- delay * 2
  }
  stop(sprintf("PUG-REST request failed after %d tries: %s", max_times, url))
}



#' Resolve arbitrary identifiers to PubChem CIDs
#'
#' Resolve user-provided identifiers (CID/SMILES/InChI/InChIKey/Name)
#' to PubChem Compound IDs (CIDs). For structural inputs (SMILES/InChI),
#' it prefers POSTing the value in the request body to avoid URL length/encoding issues,
#' and falls back to a GET path-style endpoint if needed.
#'
#' @param x Character vector. Values can be CIDs, SMILES, InChI, InChIKey, or names.
#' @param from One of \code{c("cid","smiles","inchi","inchikey","name")}.
#'
#' @return Character vector of CIDs (NA when not found).
#'
#' @importFrom purrr map_chr
#' @export
resolve_cid <- function(x, from = c("cid","smiles","inchi","inchikey","name")) {
  from <- match.arg(from)
  x <- as.character(x)
  
  if (from == "cid") {
    return(ifelse(grepl("^[0-9]+$", x), x, NA_character_))
  }
  
  param_key <- switch(from,
                      smiles   = "smiles",
                      inchi    = "inchi",
                      inchikey = "inchikey",
                      name     = "name"
  )
  
  .extract_first_cid <- function(res) {
    ids <- try(res$IdentifierList$CID, silent = TRUE)
    if (inherits(ids, "try-error") || is.null(ids) || !length(ids)) return(NA_character_)
    as.character(ids[1])
  }
  
  purrr::map_chr(x, function(q) {
    if (is.na(q) || q == "") return(NA_character_)
    
    cid_first <- NA_character_
    
    if (from %in% c("smiles","inchi")) {
      res <- try(
        pug_request(
          path   = sprintf("compound/identity/%s/cids/JSON", param_key),
          method = "POST",
          body   = stats::setNames(list(q), param_key)
        ), silent = TRUE
      )
      if (!inherits(res, "try-error")) {
        cid_first <- .extract_first_cid(res)
      }
    }
    
    if (is.na(cid_first)) {
      res2 <- try(
        pug_request(
          path   = sprintf("compound/%s/%s/cids/JSON",
                           tolower(param_key),
                           utils::URLencode(q, reserved = TRUE)),
          method = "GET"
        ), silent = TRUE
      )
      if (!inherits(res2, "try-error")) {
        cid_first <- .extract_first_cid(res2)
      }
    }
    
    cid_first
  })
}


#' Similarity search (Tanimoto) on PubChem
#'
#' @param query identifier (CID/SMILES/InChIKey/Name)
#' @param from one of c("cid","smiles","inchikey","name")
#' @param threshold integer 0–100 (percent), e.g. 90
#' @param topn max records to return. Default 10.
#' @param ... Additional arguments passed to internal helper functions.
#' 
#' @return tibble: query, hit_cid, score
#' @export
similarity_search <- function(query, from = c("smiles","cid","inchikey","name"),
                              threshold = 90, topn = 10, ...) {
  from <- match.arg(from)
  
  if (from != "smiles") {
    if (from == "cid") {
      smi <- getprops(query, properties = "CanonicalSMILES")$CanonicalSMILES[1]
    } else if (from == "inchikey") {
      cid <- resolve_cid(query, from = "inchikey")
      smi <- getprops(cid, properties = "CanonicalSMILES")$CanonicalSMILES[1]
    } else { # name
      cid <- resolve_cid(query, from = "name")
      smi <- getprops(cid, properties = "CanonicalSMILES")$CanonicalSMILES[1]
    }
    query <- smi
  }
  if (is.na(query) || !nzchar(query)) {
    return(tibble::tibble(query = NA_character_, hit_cid = character(), score = numeric()))
  }
  
  res <- pug_request(
    path  = "compound/fastsimilarity_2d/smiles/cids/JSON",
    query = list(Threshold = threshold, MaxRecords = topn, ...),
    body  = list(smiles = query),
    method = "POST"
  )
  
  cids <- res$IdentifierList$CID
  tibble::tibble(query = query,
                 hit_cid = as.character(if (is.null(cids)) character() else cids),
                 score   = NA_real_) 
}


#' Subset PPI network by edge score and top-n node degree
#'
#' @param ppi_obj An igraph object with edge attribute 'score'
#' @param n Number of top-degree nodes to keep. Default is NULL.
#' @param score_cutoff Minimum edge score to keep. Default is 0.7.
#'
#' @return A subgraph of the original PPI network
#' @importFrom igraph degree E subgraph_from_edges vcount induced_subgraph
#' @export

ppi_subset <- function(ppi_obj, n = NULL, score_cutoff = 0.7) {

  stopifnot(inherits(ppi_obj, "igraph"))
  score <- E(ppi_obj)$score

  # edge score filter
  if (is.null(score)) {
    stop("Edges must have a 'score' attribute.")
  }

  ppi_filtered <- igraph::subgraph_from_edges(ppi_obj, eids = E(ppi_obj)[score >= score_cutoff], delete.vertices = TRUE)

  if (vcount(ppi_filtered) == 0) {
    warning("No nodes left after edge score filtering.")
    return(ppi_filtered)
  }

  # degree filter
  if (!is.null(n)){
    deg <- igraph::degree(ppi_filtered, mode = "all")
    deg <- deg[!is.na(deg)]
    top_nodes <- names(sort(deg, decreasing = TRUE))[1:min(n, length(deg))]
    ppi_filtered <- igraph::induced_subgraph(ppi_filtered, vids = top_nodes)
  }

  return(ppi_filtered)
}


#' Generate Pie Chart Data for Network Visualization from Enrichment Results
#'
#' This function extracts gene-term relationships from a clusterProfiler enrichment result
#' and formats them into a matrix suitable for scatterpie-based node pie chart visualization.
#'
#' @param enrich_obj An object of class `enrichResult`, typically from `enrichGO()`, `enrichKEGG()`, or `enricher()`.
#' @param ppi_genes A character vector of gene symbols or IDs present in the PPI network (i.e., `V(ppi)$name`).
#' @param top_n Integer. The number of top enrichment terms (e.g., GO terms) to include. Default is 5.
#' @param use_weight Logical. If `TRUE`, uses enrichment significance as weights (e.g., -log10(p.adjust)); otherwise binary (0/1). Default is `FALSE`.
#' @param weight_scale Character. One of `"logp"` (default) or `"invp"`. Defines the weighting method if `use_weight = TRUE`:
#'  - `"logp"`: use `-log10(p.adjust)`
#'  - `"invp"`: use `1 / p.adjust`
#'
#' @return A data frame with genes in rows and selected enrichment terms in columns. Values represent either binary membership or weighted scores.
#' @importFrom dplyr arrange slice_head filter
#' @export

getPieData <- function(
    enrich_obj,
    ppi_genes,
    top_n = 5,
    use_weight = FALSE,
    weight_scale = c("logp", "invp")) {

  stopifnot(inherits(enrich_obj, "enrichResult"))
  weight_scale <- match.arg(weight_scale)

  # 1. extract top n pathways
  enrich_df <- enrich_obj@result %>%
    dplyr::arrange(.data$p.adjust) %>%
    dplyr::slice_head(n = top_n)

  # 2. extract the name of term and gene id
  enrich_terms <- enrich_df$Description
  gene_lists <- strsplit(enrich_df$geneID, "/")
  names(gene_lists) <- enrich_terms

  all_genes <- unique(unlist(gene_lists))

  # 3. set gene-term df
  score_df <- data.frame(name = all_genes, stringsAsFactors = FALSE)

  for (i in seq_along(enrich_terms)) {
    term <- enrich_terms[i]
    genes_in_term <- gene_lists[[i]]

    if (use_weight) {
      if (weight_scale[1] == "logp") {
        score_df[[term]] <- ifelse(score_df$name %in% genes_in_term,
                                   -log10(enrich_df$p.adjust[i] + 1e-10), 0)
      } else if (weight_scale[1] == "invp") {
        score_df[[term]] <- ifelse(score_df$name %in% genes_in_term,
                                   1 / (enrich_df$p.adjust[i] + 1e-10), 0)
      } else {
        stop("Unknown weight_scale. Choose 'logp' or 'invp'")
      }
    } else {
      score_df[[term]] <- as.integer(score_df$name %in% genes_in_term)
    }
  }

  # 4. filter the nodes in ppi
  pie_data <- score_df %>%
    dplyr::filter(.data$name %in% ppi_genes)

  if (nrow(pie_data) == 0) {
    warning("No overlapping genes found between enrichment results and PPI network.")
  }

  return(pie_data)
}

