# TCMDATA Network Pharmacology Workflow

Use this file to plan and execute TCMDATA tasks. It covers every package documentation section, but the scope rule remains: execute only what the user requested.

## 0. Setup and Documentation Coverage

Install and load TCMDATA before analysis:

```r
library(TCMDATA)
```

Use `package-sections.md` as the section checklist. Use `function-reference.md` for exact function parameters and examples.

## 1. Load TCM Data

### Search by herb name

Use `search_herb()` for herb-compound-target records.

```r
huangqi <- search_herb("黄芪", type = "Herb_cn_name")
targets <- unique(huangqi$target)
```

Name-type rule:

- Chinese characters: `type = "Herb_cn_name"`
- Continuous Pinyin, such as `Huangqi`: `type = "Herb_pinyin_name"`
- English/common lowercase names: `type = "Herb_en_name"`

Avoid spaces or hyphens in pinyin queries unless the database explicitly stores them.

### Search by target gene

Use `search_target()` when the user starts from genes.

```r
tnf_records <- search_target(c("TNF", "IL6", "AKT1"))
```

## 2. Pre-Analysis

### PubMed literature mining

Use PubMed before or after the computational analysis when the user asks for evidence, references, publication trends, or result interpretation.

```r
pubmed <- get_pubmed_data(tcm_name = "Huangqi", disease = "diabetic nephropathy")
plot(pubmed)
tbl <- get_pubmed_table(pubmed, n = 10)
```

### Herb enrichment analysis

Use `herb_enricher()` when the user provides genes and asks which herbs are enriched.

```r
genes <- c("TNF", "IL6", "AKT1", "VEGFA", "TP53")
herb_or <- herb_enricher(genes)
```

## 3. Molecule Annotation

### Resolve compound identifiers

Use `resolve_cid()` for direct PubChem identifier resolution and `getcid()` for the webchem-backed alternative.

```r
cid_tbl <- resolve_cid(c("quercetin", "kaempferol"), from = "name")
cid_alt <- getcid(c("quercetin", "kaempferol"))
```

### Retrieve chemical properties

Use `getprops()` after resolving CIDs.

```r
props <- getprops(cid_tbl$cid)
```

### Similarity searching

Use `compound_similarity()` for 2D structural similarity.

```r
similar <- compound_similarity(query = 5280343, threshold = 0.8)
```

### Download structures and convert formats

Use PubChem for ligands and RCSB PDB for receptors. Convert structures before docking when needed.

```r
sdf <- download_ligand_structure(cid = 5280343, format = "SDF")
pdb <- download_receptor_structure(pdb_id = "1A2B")
ligand_pdb <- convert_structure(sdf, output_format = "pdb")
```

## 4. Target Retrieval and DEG Analysis

### Programmatic disease target retrieval

Use `search_disease()` for disease to gene queries and `search_gene_disease()` for gene to disease queries.

```r
dn <- search_disease("diabetic nephropathy")
tnf_diseases <- search_gene_disease("TNF")
```

If the disease name is ambiguous, try a more specific term first. If fewer than 50 genes return, try a broader disease term or exact CUI if available.

### External disease target databases

When the user brings GeneCards, Open Targets, CTD, SwissTargetPrediction, or SEA results, import them as additional gene sets and keep their source labels. Do not mix evidence sources without naming them.

```r
gene_sets <- list(
  TCMDATA_herb = unique(huangqi$target),
  DisGeNET_disease = unique(dn$symbol),
  GeneCards = genecards_genes
)
```

### DEG visualization

When the user provides DEG output, filter and visualize by the documented cutoffs unless they specify different thresholds.

```r
deg_sig <- subset(deg, abs(log2FoldChange) > 1 & padj < 0.05)
```

### Intersection analysis

Use `getvenndata()` and `getvennresult()` for Venn-style intersections, `ggvenn_plot()` for Venn plots, and `upsetplot()` for larger set combinations.

```r
venn_df <- getvenndata(
  Herb = unique(huangqi$target),
  Disease = unique(dn$symbol)
)
common_targets <- getvennresult(venn_df, category = c("Herb", "Disease"))
p_intersection <- ggvenn_plot(venn_df)
```

Decision points:

- Fewer than 5 common targets: warn and suggest broader terms or relaxed filters.
- More than 500 common targets: suggest narrowing the disease, herb, or evidence source.
- For herb-disease mechanism claims, use common targets downstream.

## 5. TCM Network Construction

### Prepare herb-compound-target data

Use a clear edge table with source, target, and node type.

```r
graph_data <- prepare_herb_graph(huangqi)
```

### Sankey diagram

Use `tcm_sankey()` for herb-compound-target flows.

```r
p_sankey <- tcm_sankey(huangqi)
```

### Network graph

When the user asks for graph visualization, use `ggtangle` if available and keep node type, degree, or centrality visible.

Recommended mapping:

- Herb: muted green
- Compound: muted amber
- Target: muted blue
- Key/hub target: restrained red accent

## 6. Enrichment Analysis

### Over-representation analysis

For TCMDATA herb enrichment, use `herb_enricher()`. For GO/KEGG ORA, use `clusterProfiler` with explicit gene ID conversion.

```r
library(clusterProfiler)
library(org.Hs.eg.db)

gene_map <- bitr(common_targets, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
kegg <- enrichKEGG(gene = gene_map$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
go_bp <- enrichGO(gene = gene_map$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
```

### Advanced visualization

Use:

- `gglollipop()` for ranked enrichment terms.
- `ggdot_sankey()` for term-gene dot-sankey views.
- `go_barplot()` for GO bar plots.
- `gocircle_plot()` after preparing fold-change data with `getGores()`.

```r
p_kegg <- gglollipop(kegg, showCategory = 10)
p_go <- go_barplot(go_bp, showCategory = 10)
```

### GSEA

Use GSEA only when the user provides a ranked gene list, such as named logFC scores.

```r
gene_list <- sort(setNames(deg$log2FoldChange, deg$SYMBOL), decreasing = TRUE)
```

## 7. PPI Network Analysis

### Retrieve and filter PPI data

Use STRING through `get_ppi()` and filter with `ppi_subset()` when the user asks for confidence filtering or top interactions.

```r
ppi <- get_ppi(common_targets, taxID = 9606)
ppi_filtered <- ppi_subset(ppi, score_threshold = 400)
```

### Topological metrics and integrated ranking

Use `compute_nodeinfo()` and `rank_ppi_nodes()`.

```r
ppi_metrics <- compute_nodeinfo(ppi_filtered, weight_attr = "score")
hub_rank <- rank_ppi_nodes(ppi_metrics, top_n = 10)
```

### Radar and heatmap

Use `get_node_profile()` + `radar_plot()` for individual hub genes and `plot_node_heatmap()` for topological metrics.

```r
profile <- get_node_profile(hub_rank, node_name = hub_rank$name[1])
p_radar <- radar_plot(profile)
p_heatmap <- plot_node_heatmap(hub_rank)
```

### Community detection

Use community methods based on the question:

- `run_louvain()` for modularity-based communities.
- `run_mcode()` / `runMCODE()` for dense molecular complexes.
- `run_MCL()` for flow-based modules.
- `run_fastgreedy()` for fast modularity clustering.
- `add_cluster_score()` to rank clusters.

### Robustness analysis

Use `ppi_knock()` to simulate removal of hub genes or target sets.

```r
robust <- ppi_knock(ppi_filtered, targets = hub_rank$name[1])
```

### PPI visualization

Use topology mapping, community coloring, hub-centric layouts, or layout comparisons according to the user request. Keep labels sparse; label only top hubs unless the network is small.

## 8. Machine Learning-Based Key Target Screening

### Data preparation

Use expression matrix rows as genes and columns as samples. Group labels must match sample order.

```r
ml_data <- prepare_ml_data(expr_mat = expr, group = group, gene_filter = common_targets)
```

### Mode A: full cross-validation

Run all supported methods when the user asks for comprehensive screening.

```r
ml_all <- run_ml_screening(ml_data)
consensus <- get_ml_consensus(ml_all, min_methods = 2)
```

### Individual methods

Use individual methods when requested:

```r
fit_lasso <- ml_lasso(ml_data)
fit_enet <- ml_enet(ml_data)
fit_rf <- ml_rf(ml_data)
fit_svm <- ml_svm_rfe(ml_data)
fit_xgb <- ml_xgboost(ml_data)
```

### ML visualization

Use:

- `plot_enet_cv()`, `plot_enet_path()`, `plot_enet_coefs()`
- `plot_rf_boruta()`, `plot_rf_importance()`
- `plot_svm_rfe_curve()`
- `plot_xgb_importance()`
- `plot_ml_venn()`, `plot_ml_roc()`

### Modes B and C

Use internal train/test split when the user has one dataset and asks for validation. Use external validation when the user provides separate training and validation datasets.

### Post-hoc feature trimming and manual assembly

Use `select_features()`, `create_tcm_ml_list()`, and `get_ml_gene_sets()` to assemble or trim models.

### Single-gene diagnostic analysis

Use ROC and group expression plots for candidate biomarkers.

```r
auc_tbl <- get_gene_auc(genes = consensus$gene, expr_mat = expr, group = group)
p_roc <- plot_gene_roc(genes = consensus$gene[1:4], expr_mat = expr, group = group)
p_box <- plot_gene_boxplot(genes = consensus$gene[1:4], expr_mat = expr, group = group)
```

## 9. AI Module

### Setup

Use `tcm_setup()` for model configuration. Set `skip_internet_check = TRUE` in proxy or TUN environments when endpoint checks pass but `curl::has_internet()` misfires.

```r
tcm_setup(
  provider = "openai",
  api_key = Sys.getenv("OPENAI_API_KEY"),
  model = "gpt-5.4-mini",
  base_url = "https://www.packyapi.com/v1",
  save = TRUE,
  test = TRUE
)
```

### Interpretation layer

Use `tcm_interpret()`, `interpret_enrichment()`, `interpret_ppi()`, `interpret_table()`, and `draft_result_paragraph()` for grounded text.

```r
interpret_enrichment(kegg, language = "zh")
interpret_ppi(ppi_metrics, language = "zh")
draft_result_paragraph(kegg, result_type = "KEGG enrichment", language = "zh")
```

### Agent layer

Use:

- `tcm_agent()` for one-shot natural-language tasks.
- `tcm_chat()` for interactive sessions.
- `create_tcm_task_agent()` for programmatic sessions.
- `create_tcm_tools()` to expose a bounded set of tools.

### Artifact management

Use `save_tcm_artifact()`, `load_tcm_artifact()`, `list_tcm_artifacts()`, `summarize_tcm_artifact()`, `artifact_exists()`, and `clear_tcm_artifacts()`.

```r
art <- save_tcm_artifact(kegg, artifact_type = "enrichment_result")
kegg2 <- load_tcm_artifact(art$artifact_id)
```

### Skill management

Use `tcm_init_skills()`, `tcm_use_skills()`, `tcm_skill_dir()`, `tcm_reset_skills()`, and `tcm_aisdk_skill()` when creating or switching custom skills.

## 10. Other Resources

### Additional databases

When the user asks for transcription factor regulation, microbiota-host interactions, or external validation, keep these sources separate from TCMDATA-derived evidence and label them clearly.

### TF-gene interactions

Use external TF resources such as DoRothEA only when requested. Treat TF-target edges as regulatory evidence, not direct herb-target evidence.

### Gut microbiota-host interactions

Use GutMGene-style bacterium-metabolite-target relationships for gut-TCM integration requests. Visualize as Bacteria -> Metabolite -> Target networks.

### Docking heatmap

Use `ggdock()` for docking score matrices.

```r
p_dock <- ggdock(docking_matrix)
```

## 11. Full Pipeline Decision Tree

```text
User request
  -> Herb/target lookup only: Section 1
  -> Literature/evidence only: Section 2
  -> Molecules/structures/docking: Sections 3 and 10
  -> Disease targets or intersections: Section 4
  -> Network figure: Section 5
  -> GO/KEGG/Herb enrichment: Section 6
  -> PPI/hub genes/modules/robustness: Section 7
  -> Expression ML/biomarkers: Section 8
  -> AI interpretation/chat/report: Section 9
  -> Complete network pharmacology: Sections 1, 4, 7, 6, 5, 2, 9
```

For a complete herb-disease mechanism study, always compute common targets before PPI and enrichment.
