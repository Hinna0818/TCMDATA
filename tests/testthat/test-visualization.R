# Tests for visualization functions in TCMDATA package

test_that("gglollipop works with basic input", {
  # Create a simple mock enrichment result
  mock_enrich <- data.frame(
    ID = paste0("GO:", 1:10),
    Description = paste("Pathway", 1:10),
    GeneRatio = paste0(1:10, "/100"),
    BgRatio = paste0(5:14, "/200"),
    pvalue = runif(10, 0, 0.05),
    p.adjust = runif(10, 0, 0.05),
    qvalue = runif(10, 0, 0.05),
    geneID = paste0("Gene", 1:10),
    Count = 1:10,
    RichFactor = (1:10)/100
  )
  
  # Convert to enrichResult-like object (simplified)
  class(mock_enrich) <- c("enrichResult", "data.frame")
  attr(mock_enrich, "ontology") <- "BP"
  attr(mock_enrich, "organism") <- "hsa"
  attr(mock_enrich, "keytype") <- "ENTREZID"
  attr(mock_enrich, "pvalueCutoff") <- 0.05
  attr(mock_enrich, "qvalueCutoff") <- 0.05
  
  # Test that function runs without error
  expect_silent(p <- gglollipop(mock_enrich, top_n = 5))
  expect_s3_class(p, "ggplot")
})

test_that("ggvenn_plot works with basic input", {
  # Create a simple list of gene sets
  gene_sets <- list(
    Set1 = paste0("Gene", 1:20),
    Set2 = paste0("Gene", 10:30),
    Set3 = paste0("Gene", 15:35)
  )
  
  # Test that function runs without error
  expect_silent(p <- ggvenn_plot(gene_sets))
  expect_s3_class(p, "ggplot")
})

test_that("ggdock works with basic input", {
  # Create mock docking data
  mock_docking <- data.frame(
    herb = rep(paste0("Herb", 1:5), each = 3),
    molecule = rep(paste0("Mol", 1:3), times = 5),
    affinity = rnorm(15, -7, 1)
  )
  
  # Test that function runs without error
  expect_silent(p <- ggdock(mock_docking))
  expect_s3_class(p, "ggplot")
})

test_that("TCM_sankey works with basic input", {
  # Create mock data for sankey plot
  mock_data <- data.frame(
    herb = rep(paste0("Herb", 1:3), each = 4),
    molecule = rep(paste0("Mol", 1:4), times = 3),
    target = rep(paste0("Target", 1:6), times = 2)
  )
  
  # Test that function runs without error
  expect_silent(p <- TCM_sankey(mock_data))
  expect_s3_class(p, "ggplot")
})

test_that("PlotNodeHeatmap requires ComplexHeatmap", {
  # Skip if ComplexHeatmap not installed
  testthat::skip_if_not_installed("ComplexHeatmap")
  
  # Create mock node metrics data
  mock_metrics <- data.frame(
    name = paste0("Node", 1:10),
    centrality = runif(10),
    betweenness = runif(10),
    degree = sample(1:10, 10)
  )
  
  # Test that function runs without error when ComplexHeatmap is available
  expect_silent(p <- PlotNodeHeatmap(mock_metrics))
  expect_s3_class(p, "Heatmap")
})

test_that("gocircle_plot works with basic input", {
  # Create mock data
  mock_data <- data.frame(
    ID = paste0("GO:", 1:10),
    Description = paste("Pathway", 1:10),
    GeneRatio = paste0(1:10, "/100"),
    Count = 1:10,
    p.adjust = runif(10, 0, 0.05)
  )
  
  # Test that function runs without error
  expect_silent(p <- gocircle_plot(mock_data))
  expect_s3_class(p, "ggplot")
})

test_that("radar_plot works with basic input", {
  # Create mock data
  mock_data <- data.frame(
    Category = rep(LETTERS[1:5], each = 3),
    Group = rep(paste0("Group", 1:3), times = 5),
    Value = runif(15, 0, 1)
  )
  
  # Test that function runs without error
  expect_silent(p <- radar_plot(mock_data, category = "Category", group = "Group", value = "Value"))
  expect_s3_class(p, "ggplot")
})

test_that("visualization functions handle invalid input gracefully", {
  # Test gglollipop with non-enrichment input
  expect_error(gglollipop(data.frame()), "enrichResult")
  
  # Test ggvenn_plot with empty list
  expect_error(ggvenn_plot(list()), "length")
  
  # Test ggdock with missing columns
  expect_error(ggdock(data.frame()), "herb")
})