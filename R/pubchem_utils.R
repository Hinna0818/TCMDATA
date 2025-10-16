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


#' Get CID for compounds from PubChem
#'
#' Retrieve PubChem Compound IDs (CIDs) for a vector of compound names.
#' This function wraps `webchem::get_cid()` with built-in rate limiting
#' (≤5 requests per second as recommended by PubChem).
#'
#' @param compound A character vector of compound names.
#' @param from Source type for lookup, e.g. "name", "smiles", "inchi".
#' @param match Match mode, one of "first", "best", or "all".
#' @param pause Pause time (in seconds) between requests. Default is 0.25s (≈4 req/s).
#' @param quiet Logical, whether to suppress messages. Default TRUE.
#' @param ... Additional arguments passed to internal helper functions.
#'
#' @return A tibble with two columns: `compound` and `cid`.
#' @importFrom webchem get_cid
#' @importFrom purrr map_chr slowly
#' @importFrom tibble tibble
#' @export
#'
getcid <- function(compound,
                   from = c("name", "smiles", "inchi"),
                   match = c("first", "best", "all"),
                   pause = 0.25,
                   quiet = TRUE,
                   ...) {

  from  <- match.arg(from)
  match <- match.arg(match)

  one_query <- function(x) {
    out <- suppressWarnings(
      webchem::get_cid(x, from = from, match = match, verbose = !quiet, ...)$cid)
    if (length(out)) as.character(out[[1]]) else NA_character_
  }

  ## time limitation
  safe_query <- purrr::slowly(one_query, rate = purrr::rate_delay(pause = pause))

  cid_vec <- purrr::map_chr(compound, safe_query)

  return(tibble::tibble(cid = cid_vec, compound = compound))
}


#' get properties for compounds from PubChem
#' @param cid Integer/numeric/character vector of PubChem CIDs.
#' @param properties Character vector of PubChem property keys to request.
#' @param ... Additional arguments passed to internal helper functions.

#' @importFrom purrr slowly rate_delay insistently rate_backoff possibly map
#' @importFrom webchem pc_prop
#' @importFrom dplyr mutate
#' @importFrom tidyr unnest_wider
#' @importFrom stats na.omit setNames
#' @importFrom rlang .data
getprops <- function(cid,
                    properties = c(
                      "MolecularFormula",
                      "MolecularWeight",
                      "IUPACName",
                      "CanonicalSMILES",
                      "InChIKey",
                      "XLogP"),
                    ...){
  cid <- stats::na.omit(cid)
  slow_prop <- purrr::slowly(
    function(cid) webchem::pc_prop(cid, properties = properties, ...),
    rate = purrr::rate_delay(pause = 0.25)
  )
  prop_retry <- purrr::insistently(
    slow_prop,
    rate = purrr::rate_backoff(pause_base = 0.5, pause_cap = 4, max_times = 5),
    quiet = TRUE
  )
  prop_safe <- purrr::possibly(prop_retry, otherwise = stats::setNames(as.list(rep(NA, length(properties))), properties))

  data.frame(cid = cid) |>
    dplyr::mutate(prop = purrr::map(.data$cid, prop_safe)) |>
    tidyr::unnest_wider("prop")
}

