# TCMDATA Plotting Style

Use this file when producing report-ready R figures from TCMDATA results. The goal is a clean scientific figure with a clear biological message, not decorative complexity.

## Figure Contract

Before plotting, define:

1. Claim: one sentence describing what the figure should show.
2. Evidence: the object being plotted, such as KEGG terms, PPI hub ranking, ML features, or docking scores.
3. Audience: exploratory notebook, report, manuscript, or slide.
4. Output: object only, current R environment assignment, or saved PNG/PDF/SVG.

## Global R Theme

Use this as the default style unless the user supplies a theme.

```r
library(ggplot2)

tcm_pub_theme <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      legend.title = element_text(size = base_size * 0.9),
      legend.text = element_text(size = base_size * 0.85),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(colour = "grey30"),
      panel.grid = element_blank()
    )
}
```

## Palette

Use restrained, low-saturation colors with one accent for the main signal.

```r
tcm_cols <- c(
  herb = "#5D8A66",
  compound = "#C49A4A",
  target = "#4E79A7",
  disease = "#B85C5C",
  neutral = "#6F6F6F",
  light = "#D9DEE7"
)
```

Rules:

- Do not use rainbow palettes for enrichment or hub rankings.
- Use red only for highlighted key targets, risk, or stronger signal.
- Use color and size redundantly for enrichment plots: color for adjusted p-value, size for gene count or ratio.
- Keep labels sparse. Label top terms, hub genes, or user-requested genes only.

## Enrichment Plots

### Dot or scatter plot

Use a dot plot when the user asks for scatter/dot visualization of KEGG or GO results.

```r
df <- as.data.frame(kegg)
df <- df[order(df$p.adjust), ][seq_len(min(10, nrow(df))), ]
df$Description <- factor(df$Description, levels = rev(df$Description))

p_kegg <- ggplot(df, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, colour = p.adjust), alpha = 0.9) +
  scale_colour_gradient(low = "#B85C5C", high = "#4E79A7", name = "Adjusted P") +
  scale_size_continuous(name = "Gene count", range = c(2, 6)) +
  labs(x = "Gene ratio", y = NULL, title = "KEGG enrichment") +
  tcm_pub_theme()
```

If `GeneRatio` is stored as `"x/y"`, convert it before plotting:

```r
ratio_to_num <- function(x) vapply(strsplit(as.character(x), "/"), function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
df$GeneRatio <- ratio_to_num(df$GeneRatio)
```

### Lollipop plot

Use `gglollipop()` for compact ranked pathway views.

```r
p_go <- gglollipop(go_bp, showCategory = 10) + tcm_pub_theme()
```

### Bar plot

Use bars for a single metric such as `-log10(p.adjust)`.

```r
df$neg_log10_p <- -log10(df$p.adjust)
p_bar <- ggplot(df, aes(x = neg_log10_p, y = Description)) +
  geom_col(fill = tcm_cols["target"], width = 0.68) +
  labs(x = "-log10 adjusted P", y = NULL) +
  tcm_pub_theme()
```

## PPI Plots

### Hub heatmap

Use a heatmap for comparing centrality metrics across top hub genes.

```r
p_heat <- plot_node_heatmap(hub_rank, top_n = 10)
```

### Radar plot

Use radar plots only for one to three genes. For more genes, prefer heatmaps.

```r
profile <- get_node_profile(hub_rank, node_name = "TNF")
p_radar <- radar_plot(profile)
```

### Network plot

For network graphs, map:

- node size to degree or hub score
- node color to type or community
- label to top hubs only
- edge alpha to confidence score

Avoid plotting hundreds of labels. Filter or facet by community when a network is dense.

## Machine Learning Plots

Use the model-specific plotting functions first:

```r
p_lasso_cv <- plot_enet_cv(fit_lasso)
p_lasso_coef <- plot_enet_coefs(fit_lasso, top_n = 20)
p_rf <- plot_rf_importance(fit_rf, top_n = 20)
p_svm <- plot_svm_rfe_curve(fit_svm)
p_xgb <- plot_xgb_importance(fit_xgb, top_n = 20)
p_venn <- plot_ml_venn(ml_all)
p_roc <- plot_ml_roc(ml_all)
```

For diagnostic gene plots:

```r
p_gene_roc <- plot_gene_roc(genes = c("TNF", "IL6"), expr_mat = expr, group = group)
p_gene_box <- plot_gene_boxplot(genes = c("TNF", "IL6"), expr_mat = expr, group = group)
```

## Docking Heatmaps

Use `ggdock()` for ligand-target score matrices. Cluster rows/columns only when the figure remains interpretable.

```r
p_dock <- ggdock(docking_matrix)
```

## Saving

When the user asks to save a plot to the current R environment, assign it exactly:

```r
assign("p_kegg", p_kegg, envir = .GlobalEnv)
```

When the user asks for files, save both PNG and PDF unless they specify one format.

```r
ggsave("p_kegg.png", p_kegg, width = 7, height = 4.8, dpi = 300)
ggsave("p_kegg.pdf", p_kegg, width = 7, height = 4.8)
```

For manuscript figures, prefer PDF/SVG with editable text and 300-600 dpi raster output when required.

## Quality Checks

- Axis text must be readable without rotation unless long pathway names require wrapping.
- Pathway descriptions should be ordered by adjusted p-value or gene ratio.
- Legends must explain color and size encodings.
- Empty enrichment results should not be plotted; report the empty result and suggest ID conversion or relaxed thresholds.
- Figures should be reproducible from named objects, not only displayed interactively.
