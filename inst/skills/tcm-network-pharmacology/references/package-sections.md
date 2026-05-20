# TCMDATA Documentation Section Map

Use this file to ensure every package documentation section is represented when planning or answering network pharmacology tasks. Only execute sections explicitly requested by the user.

## 01-Load-TCM-Data.Rmd - Load TCM Data {#load-data}

- Search by herb name
- Using Pinyin format
- Search by target gene
- Session information

## 02-Pre-analysis.Rmd - Pre-analysis {#pre-analysis}

- PubMed literature mining
- Visualization
- Exporting results
- Herb enrichment analysis {#herb-ora}
- Example: Diabetic Nephropathy
- References
- Session information

## 03-Molecule-detection.Rmd - Molecule annotation {#molecule-detection}

- Resolving compound identifiers
- Using `resolve_cid()`
- Using `getcid()`
- Retrieving chemical properties
- Similarity searching
- Download molecular structures and format conversion
- Downloading ligand structures from PubChem
- Downloading receptor structures from the RCSB PDB
- Converting structure file formats
- References
- Session information

## 04-Target-retrieval-analysis.Rmd - Target retrieval and DEG analysis {#target-retrieval}

- Programmatic disease target retrieval with TCMDATA {#disease-target-tcmdata}
- `search_disease()`: disease → genes
- `search_gene_disease()`: gene → diseases
- Disease target databases {#disease-db}
- GeneCards
- Open Targets Platform
- CTD (Comparative Toxicogenomics Database)
- Other databases
- Compound target prediction {#compound-target}
- SwissTargetPrediction {#compound-target-swisstarget}
- Similarity Ensemble Approach (SEA)
- DEG visualization{#volcano}
- Load and inspect data
- Volcano plot with `ivolcano`
- Intersection analysis {#intersection}
- Prepare gene sets
- Venn diagram
- UpSet plot
- Extract intersection results
- References
- Session information

## 05-TCM-network-construction.Rmd - TCM network construction {#network-construction}

- Preparing the data
- Sankey diagram visualization
- Basic usage
- Customizing colors and appearance
- Network graph with `ggtangle`
- Building the graph
- Basic network visualization
- Adding node labels
- Alternative layouts
- Enriching the network with external data
- References
- Session information

## 06-Enrichment-analysis.Rmd - Enrichment analysis {#enrichment}

- Over Representation Analysis (ORA)
- KEGG enrichment analysis
- GO enrichment analysis
- Advanced visualization
- Lollipop plot
- Dot-sankey plot
- GO bar plot
- GO circle plot
- Gene Set Enrichment Analysis (GSEA)
- Prepare ranked gene list
- KEGG GSEA analysis
- GO GSEA analysis
- Session information

## 07-PPI-analysis.Rmd - PPI network analysis {#ppi-analysis}

- Retrieving PPI data from STRING
- Network filtering
- Topological metrics {#topology}
- Integrated ranking
- Radar plot of topological metrics
- Heatmap of topological metrics
- Community detection {#clustering}
- Louvain modularity optimization
- MCODE (Molecular Complex Detection)
- MCL (Markov Clustering)
- Fast greedy modularity optimization
- PPI network robustness analysis {#robustness}
- Algorithm overview
- Example: knocking out IL6
- Visualization: before vs. after knockout
- PPI network visualization {#ppi-vis}
- Basic network with topological mapping
- Community-colored network
- Hub-centric star layout
- Layout algorithm comparison
- References
- Session information

## 08-Machine-Learning-analysis.Rmd - Machine learning–based key target screening {#ml-analysis}

- Introduction
- Example dataset {#ml-data}
- Data preparation {#ml-prepare}
- Mode A: full cross-validation {#ml-mode-a}
- LASSO regression {#ml-lasso}
- CV error curve
- Coefficient path
- Coefficient bar chart
- Random Forest and Boruta {#ml-rf}
- Boruta importance plot
- Variable importance (Gini)
- SVM-RFE {#ml-svm-rfe}
- RFE accuracy profile
- XGBoost {#ml-xgboost}
- Feature importance (Gain)
- Consensus analysis {#ml-consensus}
- Running all methods together
- Venn diagram
- Extracting consensus genes
- ROC curves
- Mode B: internal train/test split {#ml-mode-b}
- Mode C: external validation {#ml-mode-c}
- Post-hoc feature trimming {#ml-select}
- Assembling models manually {#ml-manual-list}
- Single-gene diagnostic analysis {#ml-gene-diag}
- AUC summary table
- Single-gene ROC curves
- Two-group expression boxplot
- Summary {#ml-summary}
- References
- Session Information {#ml-session-info}

## 09-AI-module.Rmd - AI Module {#ai-module}

- Prerequisites
- Interpretation layer
- Free-text interpretation
- Interpreting enrichment objects
- Interpreting PPI objects
- Drafting result paragraphs
- Custom structured output
- Agent layer
- One-shot task: `tcm_agent()`
- Interactive session: `tcm_chat()`
- Programmatic agent: `create_tcm_task_agent()`
- Artifact management
- Skill layer
- Built-in skills
- How skills work
- Creating custom skills
- Managing skills
- Session information

## 10-Other-resources.Rmd - Other resources {#other-resources}

- Additional databases
- TF–gene interactions (DoRothEA)
- Application example
- Gut microbiota–host interactions (GutMGene)
- Application example: Gut–TCM integration
- Network visualization with ggtangle
- Additional visualization
- Heatmap for molecular docking results
- To be continued
- References
- Session information

## index.Rmd - Preface {-}

- Installation {-}

