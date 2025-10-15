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
  cid <- na.omit(cid)
  slow_prop <- purrr::slowly(
    function(cid) webchem::pc_prop(cid, properties = properties, ...),
    rate = purrr::rate_delay(pause = 0.25)
  )
  prop_retry <- purrr::insistently(
    slow_prop,
    rate = purrr::rate_backoff(pause_base = 0.5, pause_cap = 4, max_times = 5),
    quiet = TRUE
  )
  prop_safe <- purrr::possibly(prop_retry, otherwise = setNames(as.list(rep(NA, length(properties))), properties))
  
  data.frame(cid = cid) |>
    dplyr::mutate(prop = purrr::map(.data$cid, prop_safe)) |>
    tidyr::unnest_wider(prop)
}

