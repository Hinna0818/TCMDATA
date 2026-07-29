test_that("the shared publication theme has standardized typography", {
  th <- .theme_tcm_pub()

  expect_identical(eval(formals(.theme_tcm_pub)$base_family), "Arial")
  expect_equal(th$text$family, .resolve_tcm_font_family("Arial"))
  expect_equal(th$text$size, 7)
  expect_equal(th$axis.title$face, "plain")
  expect_equal(th$legend.title$face, "plain")
  expect_equal(th$plot.title$face, "plain")
})

test_that("major plotting interfaces use the shared font defaults", {
  expect_identical(eval(formals(ggppi_network)$base_family), "Arial")
  expect_identical(eval(formals(ggtcm_network)$base_family), "Arial")
  expect_identical(eval(formals(ggdock)$base_family), "Arial")
  expect_identical(eval(formals(radar_plot)$base_family), "Arial")
  expect_identical(eval(formals(go_barplot)$base_family), "Arial")
  expect_identical(eval(formals(gglollipop)$font.family), "Arial")
})

test_that("ggvenn text layers use standardized typography", {
  skip_if_not_installed("ggvenn")

  venn_df <- data.frame(
    id = letters[1:3],
    A = c(TRUE, TRUE, FALSE),
    B = c(TRUE, FALSE, TRUE)
  )
  p <- ggvenn_plot(venn_df)
  text_layers <- Filter(
    function(layer) inherits(layer$geom, "GeomText"),
    p$layers
  )

  expect_true(length(text_layers) > 0L)
  expect_true(all(vapply(
    text_layers,
    function(layer) identical(
      layer$aes_params$family,
      .resolve_tcm_font_family("Arial")
    ),
    logical(1)
  )))
  expect_true(all(vapply(
    text_layers,
    function(layer) identical(layer$aes_params$fontface, "plain"),
    logical(1)
  )))
})

test_that("upset panels inherit standardized typography", {
  skip_if_not_installed("aplotExtra")

  p <- upsetplot(list(A = c("a", "b"), B = c("b", "c")))

  expect_true(all(vapply(
    p$plotlist,
    function(panel) identical(
      panel$theme$text$family,
      .resolve_tcm_font_family("Arial")
    ),
    logical(1)
  )))
  expect_true(all(vapply(
    p$plotlist,
    function(panel) identical(panel$theme$axis.title$face, "plain"),
    logical(1)
  )))
})
