# ai_core.R
# Internal helpers: aisdk dependency check, model resolution, result wrappers.
# aisdk is an optional (Suggests) dependency - missing aisdk never breaks
# non-AI TCMDATA functions.

#' Check whether aisdk is available
#' @return TRUE invisibly if aisdk is installed; otherwise stops with a
#'   user-friendly message.
#' @keywords internal
#' @noRd
.check_aisdk <- function() {
  if (!requireNamespace("aisdk", quietly = TRUE)) {
    stop(
      "The 'aisdk' package is required for AI interpretation functions.\n",
      "Install it with: remotes::install_github('YuLab-SMU/aisdk')",
      call. = FALSE
    )
  }
  if (utils::packageVersion("aisdk") < "1.5.0") {
    stop(
      "TCMDATA requires aisdk >= 1.5.0 for AI functions.\n",
      "Update it with: remotes::install_github('YuLab-SMU/aisdk')",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Resolve a model argument for AI functions
#'
#' If \code{model} is NULL, delegates to \code{aisdk::get_model()} which
#' reads the package-wide default (set via \code{aisdk::set_model()}).
#'
#' @param model NULL, a string ID like \code{"openai:gpt-4o"}, or a
#'   LanguageModelV1 object.
#' @return A resolved model suitable for \code{aisdk::generate_object()}.
#' @keywords internal
#' @noRd
.resolve_model <- function(model = NULL) {
  .check_aisdk()
  if (is.null(model)) {
    model <- aisdk::get_model()
  }
  if (is.null(model)) {
    stop("No model configured. Run tcm_setup() first.", call. = FALSE)
  }
  return(model)
}

#' Build a tcm_ai_analysis S3 object
#' @keywords internal
#' @noRd
.new_tcm_ai_analysis <- function(input, context, output, metadata) {
  structure(
    list(
      input    = input,
      context  = context,
      output   = output,
      metadata = metadata
    ),
    class = c("tcm_ai_analysis", "tcm_ai_result")
  )
}

#' Build a tcm_ai_draft S3 object
#' @keywords internal
#' @noRd
.new_tcm_ai_draft <- function(input, context, draft, metadata) {
  structure(
    list(
      input    = input,
      context  = context,
      draft    = draft,
      metadata = metadata
    ),
    class = c("tcm_ai_draft", "tcm_ai_result")
  )
}

#' Build a tcm_ai_custom S3 object
#' @keywords internal
#' @noRd
.new_tcm_ai_custom <- function(input, context, output, metadata, schema) {
  structure(
    list(
      input    = input,
      context  = context,
      output   = output,
      metadata = metadata,
      schema   = schema
    ),
    class = c("tcm_ai_custom", "tcm_ai_result")
  )
}

#' Extract the object from a generate_object() result for a custom schema,
#' with fallback parsing. Returns list(output = ..., output_mode = ...) where
#' output_mode is "structured" or "fallback_text".
#' @keywords internal
#' @noRd
.extract_custom_object <- function(result) {
  if (!is.null(result$object) && !identical(result$valid, FALSE)) {
    return(list(output = result$object, output_mode = "structured"))
  }

  raw <- trimws(result$raw_text %||% "")
  if (!nzchar(raw)) {
    return(list(output = list(raw_text = ""), output_mode = "fallback_text"))
  }

  # Strip markdown ```json ... ``` fences if present
  clean <- gsub("^```(?:json)?[[:space:]]*|[[:space:]]*```$", "",
                raw, perl = TRUE)

  parsed <- tryCatch(
    jsonlite::fromJSON(clean, simplifyVector = FALSE),
    error = function(e) NULL
  )

  if (!is.null(parsed) && !identical(result$valid, FALSE)) {
    return(list(output = parsed, output_mode = "structured"))
  }
  list(output = list(raw_text = raw), output_mode = "fallback_text")
}

#' Build standard metadata list
#' @keywords internal
#' @noRd
.build_metadata <- function(model, language, audience, input_class,
                            prompt_version = "1.0",
                            output_mode    = "structured",
                            generation_result = NULL) {
  model_id <- if (is.character(model)) {
    model
  } else if (inherits(model, "LanguageModelV1")) {
    prov <- if (is.null(model$provider)) "unknown" else model$provider
    mid  <- if (is.null(model$model_id)) "unknown" else model$model_id
    paste0(prov, ":", mid)
  } else {
    "unknown"
  }

  metadata <- list(
    model          = model_id,
    language       = language,
    audience       = audience,
    input_class    = input_class,
    generated_at   = Sys.time(),
    prompt_version = prompt_version,
    output_mode    = output_mode
  )

  if (!is.null(generation_result)) {
    metadata$structured_output_mode <- attr(
      generation_result,
      "tcm_structured_output_mode",
      exact = TRUE
    ) %||% "unknown"
    metadata$structured_output_valid <- generation_result$valid %||%
      !is.null(generation_result$object)
    metadata$structured_output_attempts <- generation_result$attempts %||% 1L
    metadata$finish_reason <- generation_result$finish_reason %||% NULL
    metadata$usage <- generation_result$usage %||% NULL
  }

  metadata
}

# Hardcoded fallback list - used when aisdk is not yet installed
# (e.g. during tcm_config() which does not require aisdk).
# "custom" is intentionally excluded: aisdk provides create_custom_provider(),
# not create_custom(), so advertising it causes misleading failures in tcm_setup().
.tcm_providers_fallback <- c(
  "openai", "anthropic", "gemini", "deepseek", "deepseek_anthropic",
  "volcengine", "stepfun", "openrouter", "xai", "nvidia", "bailian",
  "aihubmix", "aihubmix_anthropic", "aihubmix_gemini"
)

.tcm_core_providers <- c("openai", "anthropic", "gemini")

# Resolve provider factories from aisdk 1.5's split package layout. Core
# providers remain in aisdk; additional providers live in aisdk.providers.
.resolve_tcm_provider_factory <- function(provider) {
  package_name <- if (provider %in% .tcm_core_providers) {
    "aisdk"
  } else {
    "aisdk.providers"
  }

  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(
      sprintf(
        paste0(
          "Provider '%s' requires the optional package '%s'.\n",
          "Install it with: remotes::install_github('YuLab-SMU/%s')"
        ),
        provider, package_name, package_name
      ),
      call. = FALSE
    )
  }

  factory_name <- paste0("create_", provider)
  exports <- getNamespaceExports(package_name)
  if (!factory_name %in% exports) {
    stop(
      sprintf(
        "Provider '%s' is not available in %s (missing %s()).",
        provider, package_name, factory_name
      ),
      call. = FALSE
    )
  }

  getExportedValue(package_name, factory_name)
}

.normalize_tcm_api_format <- function(api_format = NULL) {
  value <- tolower(trimws(api_format %||% "auto"))
  aliases <- c(
    auto = "auto",
    chat = "chat_completions",
    chat_completions = "chat_completions",
    responses = "responses",
    anthropic = "anthropic_messages",
    anthropic_messages = "anthropic_messages"
  )
  if (!value %in% names(aliases)) {
    stop(
      "api_format must be one of: auto, chat_completions, responses, ",
      "or anthropic_messages.",
      call. = FALSE
    )
  }
  unname(aliases[[value]])
}

.resolve_tcm_relay_format <- function(provider, api_format) {
  if (!identical(api_format, "auto")) {
    return(api_format)
  }
  if (provider == "anthropic" || grepl("_anthropic$", provider)) {
    "anthropic_messages"
  } else {
    "chat_completions"
  }
}

.validate_tcm_native_api_format <- function(provider, api_format) {
  allowed <- if (provider == "openai") {
    c("auto", "chat_completions", "responses")
  } else if (provider %in% c("anthropic", "deepseek_anthropic",
                             "aihubmix_anthropic")) {
    c("auto", "anthropic_messages")
  } else if (provider %in% c("gemini", "aihubmix_gemini")) {
    "auto"
  } else {
    c("auto", "chat_completions")
  }

  if (!api_format %in% allowed) {
    stop(
      sprintf(
        "api_format='%s' is not supported by the native '%s' provider.",
        api_format, provider
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.as_optional_logical <- function(x, name) {
  if (is.null(x) || (is.character(x) && length(x) == 1L && !nzchar(x))) {
    return(NULL)
  }
  if (is.logical(x) && length(x) == 1L && !is.na(x)) {
    return(x)
  }
  if (is.character(x) && length(x) == 1L) {
    value <- tolower(trimws(x))
    if (value %in% c("true", "1", "yes")) return(TRUE)
    if (value %in% c("false", "0", "no")) return(FALSE)
  }
  stop(name, " must be TRUE, FALSE, or NULL.", call. = FALSE)
}

# Return the validated provider whitelist.
# Dynamic scanning of aisdk create_*() functions is intentionally avoided:
# it would also match create_agent(), create_skill_registry(), etc.
.available_providers <- function() {
  .tcm_providers_fallback
}

# Providers whose create_* function does not accept base_url.
.tcm_no_base_url <- c(
  "deepseek_anthropic", "aihubmix_anthropic", "aihubmix_gemini"
)

#' Write AI provider credentials to .env
#'
#' Saves provider, model, endpoint, and optional API protocol settings to a
#' \code{.env} file. Existing \code{TCM_*} settings managed by this function
#' are overwritten; all other lines are preserved.
#'
#' Supported \code{provider} values: \code{"openai"}, \code{"anthropic"},
#' \code{"gemini"}, \code{"deepseek"}, \code{"deepseek_anthropic"},
#' \code{"volcengine"}, \code{"stepfun"}, \code{"openrouter"}, \code{"xai"},
#' \code{"nvidia"}, \code{"bailian"}, \code{"aihubmix"},
#' \code{"aihubmix_anthropic"}, \code{"aihubmix_gemini"}.
#'
#' @param provider Character. Provider name (see Details).
#' @param api_key Character. Your API key.
#' @param model Character. Model name, e.g. \code{"gpt-4o-mini"},
#'   \code{"claude-3-5-haiku-20241022"}, \code{"gemini-2.0-flash"}.
#' @param base_url Character or NULL. Override the default API endpoint.
#'   Required for proxies or self-hosted endpoints.
#' @param path Character. Path to the \code{.env} file. Default \code{".env"}.
#' @param api_format Character or NULL. API protocol used by a custom
#'   \code{base_url}: \code{"chat_completions"}, \code{"responses"}, or
#'   \code{"anthropic_messages"}. NULL uses provider-specific defaults.
#' @param supports_native_tools Logical or NULL. Whether a relay endpoint
#'   accepts native tool calls. The relay default is TRUE.
#' @param disable_stream_options Logical or NULL. Whether an OpenAI-compatible
#'   relay should omit \code{stream_options}. The relay default is TRUE.
#' @param responses_state_mode Character or NULL. Conversation state mode for
#'   Responses-compatible endpoints: \code{"stateless"}, \code{"auto"}, or
#'   \code{"server"}.
#'
#' @return The path to the \code{.env} file, invisibly.
#' @examples
#' \dontrun{
#'   tcm_config("openai",    "sk-xxx",     "gpt-4o-mini")
#'   tcm_config("anthropic", "sk-ant-xxx", "claude-3-5-haiku-20241022")
#'   tcm_config("gemini",    "AIza-xxx",   "gemini-2.0-flash")
#'   tcm_config("deepseek",  "sk-xxx",     "deepseek-chat")
#'   tcm_config("openai",    "sk-xxx",     "gpt-5-minimal",
#'              base_url = "https://www.packyapi.com/v1")
#' }
#' @export
tcm_config <- function(provider,
                       api_key,
                       model,
                       base_url = NULL,
                       path = ".env",
                       api_format = NULL,
                       supports_native_tools = NULL,
                       disable_stream_options = NULL,
                       responses_state_mode = NULL) {
  if (!provider %in% .available_providers()) {
    stop(sprintf("Unknown provider '%s'. Supported: %s",
                 provider, paste(.available_providers(), collapse = ", ")),
         call. = FALSE)
  }

  env_path <- if (dirname(path) == ".") file.path(getwd(), path) else path
  existing <- if (file.exists(env_path)) readLines(env_path) else character(0)
  managed <- c(
    "TCM_PROVIDER", "TCM_API_KEY", "TCM_MODEL", "TCM_BASE_URL",
    "TCM_API_FORMAT", "TCM_SUPPORTS_NATIVE_TOOLS",
    "TCM_DISABLE_STREAM_OPTIONS", "TCM_RESPONSES_STATE_MODE"
  )
  kept <- existing[
    !grepl(paste0("^(", paste(managed, collapse = "|"), ")="), existing)
  ]

  new_lines <- c(
    kept,
    sprintf("TCM_PROVIDER=%s", provider),
    sprintf("TCM_API_KEY=%s",  api_key),
    sprintf("TCM_MODEL=%s",    model)
  )
  if (!is.null(base_url) && nzchar(base_url)) {
    new_lines <- c(new_lines, sprintf("TCM_BASE_URL=%s", base_url))
  }
  if (!is.null(api_format) && nzchar(api_format)) {
    new_lines <- c(new_lines, sprintf(
      "TCM_API_FORMAT=%s", .normalize_tcm_api_format(api_format)
    ))
  }
  supports_native_tools <- .as_optional_logical(
    supports_native_tools, "supports_native_tools"
  )
  disable_stream_options <- .as_optional_logical(
    disable_stream_options, "disable_stream_options"
  )
  if (!is.null(supports_native_tools)) {
    new_lines <- c(new_lines, sprintf(
      "TCM_SUPPORTS_NATIVE_TOOLS=%s", tolower(supports_native_tools)
    ))
  }
  if (!is.null(disable_stream_options)) {
    new_lines <- c(new_lines, sprintf(
      "TCM_DISABLE_STREAM_OPTIONS=%s", tolower(disable_stream_options)
    ))
  }
  if (!is.null(responses_state_mode) && nzchar(responses_state_mode)) {
    responses_state_mode <- match.arg(
      responses_state_mode, c("stateless", "auto", "server")
    )
    new_lines <- c(new_lines, sprintf(
      "TCM_RESPONSES_STATE_MODE=%s", responses_state_mode
    ))
  }

  writeLines(new_lines, env_path)
  message(sprintf("tcm_config: saved [provider=%s, model=%s] -> %s",
                  provider, model, env_path))
  invisible(env_path)
}

#' Initialise the AI model from .env or explicit arguments
#'
#' Loads \code{.env} (if present), resolves \code{TCM_*} variables, creates the
#' requested provider, and registers the model via \code{aisdk::set_model()}.
#' OpenAI, Anthropic, and Gemini are provided by \pkg{aisdk}; DeepSeek and
#' other additional providers are resolved from \pkg{aisdk.providers}.
#'
#' @param provider Character or NULL. Overrides \code{TCM_PROVIDER}.
#' @param api_key Character or NULL. Overrides \code{TCM_API_KEY}.
#' @param model Character or NULL. Overrides \code{TCM_MODEL}.
#' @param base_url Character or NULL. Overrides \code{TCM_BASE_URL}.
#' @param .env Logical. Load \code{.env} before reading env vars (default TRUE).
#' @param save Logical. If TRUE, also calls \code{\link{tcm_config}()} to
#'   persist the resolved credentials to \code{.env}. Useful for first-time
#'   setup when you want a single call to both initialise and save.
#'   Default FALSE.
#' @param test Logical. If TRUE, sends a minimal test request after setup to
#'   verify the API key and endpoint are reachable. Warnings (not errors) are
#'   issued on failure so the model is still registered. Default FALSE.
#' @param force_json_schema Logical. Default \code{TRUE}. When \code{TRUE},
#'   the schema is also passed as \code{response_format} when structured output
#'   must use JSON mode. With native tool support, TCMDATA uses the aisdk 1.5
#'   forced-tool structured-output mode and this option has no effect. Set
#'   \code{FALSE} only when a JSON-only relay rejects \code{response_format}.
#' @param skip_internet_check Logical. Default \code{TRUE}. Sets
#'   \code{options(aisdk.skip_internet_check = TRUE)} before live requests so
#'   \code{curl::has_internet()} false negatives in proxy/VPN environments do
#'   not block otherwise reachable API endpoints.
#' @param api_format Character or NULL. API protocol for a custom
#'   \code{base_url}. Supported values are \code{"chat_completions"},
#'   \code{"responses"}, and \code{"anthropic_messages"}; common aliases
#'   \code{"chat"} and \code{"anthropic"} are accepted. NULL selects the
#'   provider-specific default.
#' @param supports_native_tools Logical or NULL. Whether a relay endpoint
#'   accepts native tool calls. Defaults to TRUE for relay endpoints; set FALSE
#'   to use aisdk's text-embedded tool-call fallback.
#' @param disable_stream_options Logical or NULL. Whether an OpenAI-compatible
#'   relay should omit \code{stream_options}. Defaults to TRUE for relays.
#' @param responses_state_mode Character or NULL. One of \code{"stateless"},
#'   \code{"auto"}, or \code{"server"}. Relay endpoints default to
#'   \code{"stateless"}.
#'
#' @return The model object, invisibly.
#' @examples
#' \dontrun{
#'   # Standard two-step workflow
#'   tcm_config("openai", "sk-xxx", "gpt-4o-mini")
#'   tcm_setup()
#'
#'   # One-step: configure + initialise in a single call
#'   tcm_setup("openai", "sk-xxx", "gpt-4o-mini", save = TRUE)
#'
#'   # Verify connectivity after setup
#'   tcm_setup(test = TRUE)
#'
#'   # Override at runtime without touching .env
#'   tcm_setup("deepseek", api_key = "sk-xxx", model = "deepseek-chat")
#'
#'   # OpenAI-compatible relay (also works with a relayed DeepSeek model)
#'   tcm_setup("deepseek", "sk-xxx", "deepseek-chat",
#'             base_url = "https://relay.example.com/v1")
#'
#'   # Relay implementing the OpenAI Responses API
#'   tcm_setup("openai", "sk-xxx", "gpt-5-mini",
#'             base_url = "https://relay.example.com/v1",
#'             api_format = "responses")
#' }
#' @export
tcm_setup <- function(provider         = NULL,
                      api_key          = NULL,
                      model            = NULL,
                      base_url         = NULL,
                      .env             = TRUE,
                      save             = FALSE,
                      test             = FALSE,
                      force_json_schema = TRUE,
                      skip_internet_check = TRUE,
                      api_format = NULL,
                      supports_native_tools = NULL,
                      disable_stream_options = NULL,
                      responses_state_mode = NULL) {
  .check_aisdk()

  if (isTRUE(skip_internet_check)) {
    options(aisdk.skip_internet_check = TRUE)
  }

  if (.env && requireNamespace("dotenv", quietly = TRUE)) {
    env_file <- file.path(getwd(), ".env")
    if (file.exists(env_file)) dotenv::load_dot_env(env_file)
  }

  provider <- .env_or(provider, Sys.getenv("TCM_PROVIDER", "openai"))
  api_key  <- .env_or(api_key,  Sys.getenv("TCM_API_KEY",  ""))
  model    <- .env_or(model,    Sys.getenv("TCM_MODEL",    "gpt-4o-mini"))
  base_url <- .env_or(base_url, Sys.getenv("TCM_BASE_URL", ""))
  api_format <- .normalize_tcm_api_format(
    .env_or(api_format, Sys.getenv("TCM_API_FORMAT", "auto"))
  )
  supports_native_tools <- .as_optional_logical(
    supports_native_tools %||% Sys.getenv("TCM_SUPPORTS_NATIVE_TOOLS", ""),
    "supports_native_tools"
  )
  disable_stream_options <- .as_optional_logical(
    disable_stream_options %||% Sys.getenv("TCM_DISABLE_STREAM_OPTIONS", ""),
    "disable_stream_options"
  )
  responses_state_mode <- .env_or(
    responses_state_mode,
    Sys.getenv("TCM_RESPONSES_STATE_MODE", "")
  )
  if (nzchar(responses_state_mode)) {
    responses_state_mode <- match.arg(
      responses_state_mode, c("stateless", "auto", "server")
    )
  } else {
    responses_state_mode <- NULL
  }

  if (!provider %in% .available_providers()) {
    stop(
      sprintf(
        "Unknown provider '%s'. Supported: %s",
        provider, paste(.available_providers(), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!nzchar(api_key)) {
    env_path <- file.path(getwd(), ".env")
    stop(
      "No API key found.\n",
      "  Option A: pass it directly -\n",
      "    tcm_setup(provider=\"openai\", api_key=\"sk-...\",",
      " model=\"gpt-4o-mini\")\n",
      "  Option B: save to .env first -\n",
      "    tcm_config(\"openai\", \"sk-...\", \"gpt-4o-mini\")",
      " then tcm_setup()\n",
      "  .env searched at: ", env_path, "\n",
      "  (tcm_config() writes to .env in the current working directory)",
      call. = FALSE
    )
  }

  if (save) {
    tcm_config(
      provider = provider,
      api_key  = api_key,
      model    = model,
      base_url = if (nzchar(base_url)) base_url else NULL,
      api_format = api_format,
      supports_native_tools = supports_native_tools,
      disable_stream_options = disable_stream_options,
      responses_state_mode = responses_state_mode
    )
  }

  use_relay <- nzchar(base_url) &&
    !provider %in% c("gemini", "aihubmix_gemini")
  use_openai_responses <- FALSE

  if (nzchar(base_url) && provider == "aihubmix_gemini") {
    stop(
      "provider='aihubmix_gemini' does not support a custom base_url.",
      call. = FALSE
    )
  }

  if (use_relay) {
    relay_format <- .resolve_tcm_relay_format(provider, api_format)
    provider_obj <- aisdk::create_custom_provider(
      provider_name = provider,
      base_url = base_url,
      api_key = api_key,
      api_format = relay_format,
      disable_stream_options = disable_stream_options %||% TRUE,
      supports_native_tools = supports_native_tools %||% TRUE,
      responses_state_mode = responses_state_mode %||% "stateless"
    )
  } else {
    .validate_tcm_native_api_format(provider, api_format)
    create_fn <- .resolve_tcm_provider_factory(provider)
    args <- list(api_key = api_key)
    if (nzchar(base_url) && !provider %in% .tcm_no_base_url) {
      args$base_url <- base_url
    }
    if (provider == "openai") {
      args$api_format <- switch(
        api_format,
        auto = "auto",
        chat_completions = "chat",
        responses = "responses"
      )
      args$responses_state_mode <- responses_state_mode %||% "auto"
      if (!is.null(disable_stream_options)) {
        args$disable_stream_options <- disable_stream_options
      }
      use_openai_responses <- identical(api_format, "responses")
    }
    provider_obj <- do.call(create_fn, args)
  }

  model_obj <- if (use_openai_responses) {
    provider_obj$responses_model(model)
  } else {
    provider_obj$language_model(model)
  }
  aisdk::set_model(model_obj)

  if (test) {
    tryCatch(
      aisdk::generate_text(
        model = model_obj,
        prompt = "ping",
        temperature = NULL
      ),
      error = function(e) warning(
        "Connection test failed: ", conditionMessage(e), call. = FALSE
      )
    )
  }

  options(tcm.force_json_schema = isTRUE(force_json_schema))
  options(
    tcm.supports_native_tools = if (use_relay) {
      supports_native_tools %||% TRUE
    } else {
      TRUE
    }
  )

  message(sprintf("tcm_setup: %s / %s ready.", provider, model))
  invisible(model_obj)
}

# null-coalescing operator (NULL only - consistent with rlang::`%||%`)
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Wrapper around aisdk::generate_object() with provider compatibility
#'
#' Uses aisdk 1.5 forced-tool structured output when native tool calls are
#' available. JSON mode remains available for endpoints configured with
#' \code{supports_native_tools = FALSE}; the legacy \code{response_format}
#' compatibility option is applied only on that path.
#' @keywords internal
#' @noRd
.call_generate_object <- function(model, prompt, schema, system,
                                  temperature = 0.3) {
  mode <- if (isTRUE(getOption("tcm.supports_native_tools", TRUE))) {
    "tool"
  } else {
    "json"
  }
  args <- list(
    model       = model,
    prompt      = prompt,
    schema      = schema,
    schema_name = "tcm_result",
    system      = system,
    temperature = temperature,
    mode        = mode,
    max_retries = 1L
  )
  if (identical(mode, "json") &&
      isTRUE(getOption("tcm.force_json_schema", FALSE))) {
    args$response_format <- schema
  }
  result <- do.call(aisdk::generate_object, args)
  attr(result, "tcm_structured_output_mode") <- mode
  result
}

# env-var helper: NULL *or* empty string falls back to default.
# Used exclusively inside tcm_setup() for Sys.getenv() resolution.
.env_or <- function(x, default) if (is.null(x) || !nzchar(x)) default else x

#' Extract structured output from a generate_object result
#'
#' Tries three strategies in order:
#' 1. Native structured object from aisdk (\code{result$object}).
#' 2. JSON parsed from \code{result$raw_text} (models that return JSON text
#'    but do not support the structured-output API).
#' 3. Plain-text fallback: wraps \code{raw_text} into \code{$summary} so the
#'    result is never NULL, even for models that ignore schema instructions.
#'
#' @param result A \code{GenerateObjectResult} from
#'   \code{aisdk::generate_object()}.
#' @param type One of \code{"analysis"} or \code{"draft"}, controls the shape
#'   of the plain-text fallback object.
#' @return A list with two elements: \code{$output} (a named list matching the
#'   schema, never NULL) and \code{$output_mode} ("structured" when the model
#'   returned valid JSON, "fallback_text" otherwise).
#' @keywords internal
#' @noRd
.extract_object <- function(result, type = "analysis") {
  # Case 1: model returned a proper structured object
  if (!is.null(result$object) && !identical(result$valid, FALSE)) {
    return(list(output = result$object, output_mode = "structured"))
  }

  raw <- trimws(result$raw_text %||% "")

  # Case 2: empty response
  if (!nzchar(raw)) {
    return(list(output = .fallback_object("", type),
                output_mode = "fallback_text"))
  }

  # Case 3: try stripping markdown code fences then parsing JSON
  clean <- gsub("^```(?:json)?[[:space:]]*|[[:space:]]*```$", "",
                raw, perl = TRUE)
  parsed <- tryCatch(
    jsonlite::fromJSON(clean, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (!is.null(parsed) && !identical(result$valid, FALSE)) {
    return(list(output = parsed, output_mode = "structured"))
  }

  # Case 4: plain-text fallback - model ignored the schema
  list(output = .fallback_object(raw, type), output_mode = "fallback_text")
}

#' Build a minimal fallback object from plain text
#' @keywords internal
#' @noRd
.fallback_object <- function(text, type) {
  if (type == "draft") {
    list(paragraph = text, figure_legend_hint = "")
  } else {
    list(
      summary                   = text,
      key_findings              = list(),
      biological_interpretation = "",
      tcm_relevance             = "",
      caveats                   = list()
    )
  }
}
