#' Internal publication plotting helpers
#'
#' @importFrom ggplot2 element_blank element_line element_text theme
#' @importFrom ggplot2 theme_classic theme_void
#' @importFrom grid unit
#' @keywords internal
#' @noRd
.theme_tcm_text <- function(base_size = 7, base_family = "Arial") {
  base_family <- .resolve_tcm_font_family(base_family)
  theme(
    text = element_text(
      family = base_family,
      size = base_size,
      face = "plain",
      colour = "grey10"
    ),
    axis.text = element_text(
      size = base_size * 0.9,
      face = "plain",
      colour = "grey10"
    ),
    axis.title = element_text(
      size = base_size,
      face = "plain",
      colour = "grey10"
    ),
    legend.title = element_text(size = base_size * 0.9, face = "plain"),
    legend.text = element_text(size = base_size * 0.85, face = "plain"),
    strip.text = element_text(size = base_size * 0.95, face = "plain"),
    plot.title = element_text(size = base_size * 1.05, face = "plain"),
    plot.subtitle = element_text(
      size = base_size * 0.9,
      face = "plain",
      colour = "grey35"
    )
  )
}

#' @keywords internal
#' @noRd
.theme_tcm_pub <- function(base_size = 7,
                           base_family = "Arial",
                           grid = c("none", "x", "y", "both")) {
  grid <- match.arg(grid)
  base_family <- .resolve_tcm_font_family(base_family)

  major_grid <- element_line(linewidth = 0.2, colour = "grey88")
  blank_grid <- element_blank()

  theme_classic(base_size = base_size, base_family = base_family) +
    .theme_tcm_text(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = unit(1.5, "pt"),
      legend.key.height = unit(3.5, "mm"),
      legend.key.width = unit(3.5, "mm"),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      strip.background = element_blank(),
      panel.grid.major.x = if (grid %in% c("x", "both")) major_grid else blank_grid,
      panel.grid.major.y = if (grid %in% c("y", "both")) major_grid else blank_grid,
      panel.grid.minor = element_blank()
    )
}

#' @keywords internal
#' @noRd
.theme_tcm_void <- function(base_size = 7, base_family = "Arial") {
  base_family <- .resolve_tcm_font_family(base_family)
  theme_void(base_size = base_size, base_family = base_family) +
    .theme_tcm_text(base_size = base_size, base_family = base_family)
}

#' @keywords internal
#' @noRd
.resolve_tcm_font_family <- function(base_family = "Arial") {
  if (!identical(tolower(base_family), "arial")) {
    return(base_family)
  }

  device <- names(grDevices::dev.cur())
  unsupported <- c(
    "pdf", "postscript", "xfig", "pictex"
  )
  if (length(device) > 0L && device %in% unsupported) {
    return("")
  }
  if (identical(device, "null device")) {
    default_device <- getOption("device")
    if (identical(default_device, grDevices::pdf) ||
        identical(default_device, grDevices::postscript) ||
        (is.character(default_device) &&
         default_device %in% unsupported)) {
      return("")
    }
  }
  base_family
}

#' @keywords internal
#' @noRd
.par_tcm_pub <- function(base_family = "Arial") {
  graphics::par(
    family = .resolve_tcm_font_family(base_family),
    cex = 0.75,
    cex.axis = 0.75,
    cex.lab = 0.85,
    font.lab = 1,
    font.main = 1,
    lwd = 0.8,
    bty = "l"
  )
}

#' @keywords internal
#' @noRd
.pal_tcm_methods <- function(n) {
  pal <- c(
    "#4E79A7", "#59A14F", "#E15759", "#B07AA1",
    "#F28E2B", "#76B7B2", "#9C755F", "#BAB0AC",
    "#EDC948", "#8CD17D", "#A0CBE8", "#FFBE7D"
  )
  rep_len(pal, n)
}

#' @keywords internal
#' @noRd
.pal_tcm_direction <- function() {
  c(Positive = "#B85C5C", Negative = "#4E79A7")
}

#' @keywords internal
#' @noRd
.pal_tcm_decision <- function() {
  c(Confirmed = "#4E9F7A", Tentative = "#D8A24A", Rejected = "#B85C5C")
}

#' @keywords internal
#' @noRd
.pal_tcm_groups <- function() {
  c("#B85C5C", "#4E79A7", "#59A14F", "#B07AA1")
}
