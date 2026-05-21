---
name: tcm-network-pharmacology
description: >
  TCMDATA network pharmacology analysis skill for Traditional Chinese Medicine workflows in R.
  Use when the user asks to query herbs, compounds, targets, disease genes, PubMed evidence,
  molecule annotation, compound similarity, herb/disease target intersection, GO/KEGG/Herb
  enrichment, GSEA, PPI network construction, topology calculation and ranking, community detection, PPI robustness
  analysis, machine-learning key target screening, diagnostic ROC/boxplots, AI interpretation,
  report drafting, or publication-style visualization for TCM network pharmacology. Also use
  for Chinese requests such as 中药网络药理学分析, 靶点查询, 疾病靶点, 富集分析, PPI分析, 机器学习筛选,
  分子注释, 文献挖掘, and AI模块.
---

# TCMDATA Network Pharmacology

Use this skill as an execution guide for TCMDATA-based network pharmacology tasks. Execute exactly the analysis scope requested by the user. Do not expand a target query into a full pipeline unless the user explicitly asks for a complete, systematic, or network pharmacology analysis.

## First Move

1. Identify the requested scope: lookup, molecule annotation, intersection, enrichment, PPI, machine learning, visualization, AI interpretation, or full pipeline.
2. Load only the reference file needed for that scope.
3. Prefer TCMDATA functions and built-in agent tools over ad hoc code.
4. Preserve object names requested by the user, especially plot names such as `p_kegg`.
5. Report exact counts, thresholds, and artifact IDs from tool outputs. Never invent genes, terms, p-values, pathways, or evidence.

## Reference Files

| File | Open when |
|------|-----------|
| `references/workflow.md` | Need a full or partial analysis workflow covering all package documentation sections |
| `references/function-reference.md` | Need exact function usage, parameter interpretation, or examples for any exported TCMDATA function |
| `references/package-sections.md` | Need to ensure every documentation section is represented in a plan or answer |
| `references/plotting-style.md` | Need publication-style R plots for enrichment, PPI, ML, network, docking, or report figures |
| `references/quality-checkpoints.md` | Need thresholds, warnings, or sanity checks for target counts, intersections, PPI, enrichment, hubs, or literature |
| `references/data-sources.md` | Need source grounding, evidence hierarchy, anti-hallucination rules, or cross-validation guidance |

## Scope Rules

- For simple lookup requests, run only the relevant query function and summarize the result.
- For enrichment requests, use the gene set explicitly requested. If the user says herb targets, use herb targets; if the user says herb-disease intersection, compute the intersection first.
- For full network pharmacology requests, run target collection, intersection, PPI, enrichment, visualization, evidence validation, and report synthesis in that order.
- For expression, DEG, WGCNA, ML, or single-cell validation, run those sections only when the user provides suitable data or explicitly asks to incorporate them.
- For AI chat/tool workflows, use artifact IDs returned by tools; if executing R code in the same turn, load artifacts through `load_tcm_artifact()` or use exported artifact variables when available.

## Core Pipeline

Use this order for a complete herb-disease network pharmacology analysis:

1. Herb targets: `search_herb()`
2. Disease targets: `search_disease()`
3. Intersection: `getvenndata()` / `getvennresult()` or the agent intersection tool
4. PPI: `get_ppi()` -> `compute_nodeinfo()` -> `rank_ppi_nodes()`
5. Community and robustness when requested: `run_louvain()`, `run_mcode()`, `run_MCL()`, `run_fastgreedy()`, `ppi_knock()`
6. Enrichment: `herb_enricher()` and standard GO/KEGG workflows from `workflow.md`
7. Machine learning when expression data is present: `prepare_ml_data()` -> `run_ml_screening()` -> `get_ml_consensus()`
8. Visualization: select plots from `plotting-style.md` and `function-reference.md`
9. Evidence and interpretation: `get_pubmed_data()`, `get_pubmed_table()`, `tcm_interpret()`, `draft_result_paragraph()`

Key rule: use intersection genes, not all herb or disease genes, for downstream PPI/enrichment when the user asks for herb-disease mechanism analysis.

## Reporting Rules

- State source and exact counts: herb records, unique targets, disease genes, intersection size, PPI nodes/edges, significant terms, selected hub genes.
- Use adjusted p-values for enrichment claims.
- Distinguish computational predictions from experimental, literature, or clinical evidence.
- Mention empty results honestly and suggest parameter relaxation only after reporting the failure.
- Keep Chinese answers in Chinese unless the user requests English.

## Minimal Report Template

```markdown
## [Herb] - [Disease] Network Pharmacology

### Target Retrieval
- Herb targets: [N]
- Disease genes: [N]
- Common targets: [N]

### PPI and Hub Genes
- STRING threshold: [score]
- Network: [nodes] nodes, [edges] edges
- Top hub genes: [genes]

### Functional Enrichment
- GO-BP top terms: [terms with adjusted p-values]
- KEGG top pathways: [pathways with adjusted p-values]

### Evidence and Interpretation
- PubMed evidence: [N] records
- Mechanistic summary: [grounded interpretation]
```
