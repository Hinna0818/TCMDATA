# TCMDATA Function Reference

Use this file when the user asks for exact function parameters or examples. Each exported function below includes usage, parameter interpretation, and a runnable or fill-in usage pattern. Prefer examples from package documentation when present; otherwise adapt the minimal pattern to the current objects.

## PPI network and topology

### `add_cluster_score()`

Score and Rank Network Clusters

This function evaluates the clusters/modules identified in the graph.
It calculates a score for each cluster based on density and size (MCODE style),
and returns a ranked data frame.

Usage:
```r
add_cluster_score(g, cluster_attr = "louvain_cluster", min_size = 3)

Addclusterscore(...)
```

参数解析:
- `g`: An igraph object. The graph must have a vertex attribute containing cluster labels.
- `cluster_attr`: Character. The name of the vertex attribute that stores cluster labels. Default is louvain_cluster.
- `min_size`: Integer. Clusters smaller than this size will be ignored. Default is 3.
- `...`: Additional arguments passed to add_cluster_score when using the deprecated alias Addclusterscore.

使用示例:
```r
data(demo_ppi)
ppi <- run_louvain(demo_ppi, resolution = 1)
louvain_score <- add_cluster_score(ppi, cluster_attr = "louvain_cluster", min_size = 3)
head(louvain_score)
```

### `Addclusterscore()`

Score and Rank Network Clusters

This function evaluates the clusters/modules identified in the graph.
It calculates a score for each cluster based on density and size (MCODE style),
and returns a ranked data frame.

Usage:
```r
add_cluster_score(g, cluster_attr = "louvain_cluster", min_size = 3)

Addclusterscore(...)
```

参数解析:
- `g`: An igraph object. The graph must have a vertex attribute containing cluster labels.
- `cluster_attr`: Character. The name of the vertex attribute that stores cluster labels. Default is louvain_cluster.
- `min_size`: Integer. Clusters smaller than this size will be ignored. Default is 3.
- `...`: Additional arguments passed to add_cluster_score when using the deprecated alias Addclusterscore.

使用示例:
```r
data(demo_ppi)
ppi <- run_louvain(demo_ppi, resolution = 1)
louvain_score <- add_cluster_score(ppi, cluster_attr = "louvain_cluster", min_size = 3)
head(louvain_score)
```

### `compute_nodeinfo()`

Compute some useful metrics for PPI nodes

Compute some useful metrics for PPI nodes

Usage:
```r
compute_nodeinfo(g, weight_attr = "score", normalize = FALSE, seed = 42)
```

参数解析:
- `g`: An igraph object containing PPI information.
- `weight_attr`: Character; The attribute of weight in this PPI object.
- `normalize`: Logical; whether to normalize betweenness and closeness centrality to the range [0, 1]. When TRUE, each value is divided by the theoretical maximum.
- `seed`: Numeric; Random seed for function compute_EPC. Default is 42.

使用示例:
```r
data(demo_ppi)
library(igraph)
ppi <- compute_nodeinfo(demo_ppi, weight_attr = "score")
str(ppi)
```

### `get_mcode_res()`

Extract MCODE analysis results as a data frame

Extract MCODE analysis results as a data frame

Usage:
```r
get_mcode_res(g, only_clusters = FALSE)

getMCODE_res(...)
```

参数解析:
- `g`: An igraph object processed by runMCODE.
- `only_clusters`: Logical. If TRUE, returns only nodes that belong to a cluster. Default FALSE.
- `...`: Additional arguments passed to get_mcode_res when using the deprecated alias getMCODE_res.

使用示例:
```r
data(demo_ppi)
ppi <- runMCODE(demo_ppi, max_depth = 100)
MCODE_res <- get_mcode_res(ppi)
print(head(MCODE_res))
```

### `get_node_profile()`

Extract normalized centrality profile for a node

This function extracts selected metrics for a given node and normalizes
them to the range [0, 1] for use in plots such as radar charts.

Usage:
```r
get_node_profile(
tab,
node_name,
metrics = c("degree", "betweenness", "closeness", "MCC", "MNC", "DMNC", "coreness",
"EPC")
)
```

参数解析:
- `tab`: A data frame containing at least a name column and the specified metric columns.
- `node_name`: Character; node name to extract.
- `metrics`: Character vector; names of metric columns to use. Defaults to common centrality measures.

使用示例:
```r
# Fill required arguments according to the parameter list above
get_node_profile(
tab,
node_name,
metrics = c("degree", "betweenness", "closeness", "MCC", "MNC", "DMNC", "coreness",
"EPC")
)
```

### `get_ppi()`

Retrieve a STRING PPI network

Thin wrapper around clusterProfiler::getPPI() so the function is
available directly from TCMDATA.

Usage:
```r
get_ppi(x, taxID = "auto", ...)
```

参数解析:
- `x`: An enrichResult object or a character vector of proteins.
- `taxID`: NCBI taxonomy identifier. Use "auto" to let
clusterProfiler::getPPI() infer it from an enrichResult
when possible.
- `...`: Additional arguments passed to clusterProfiler::getPPI(),
such as ID, required_score, network_type,
add_nodes, show_query_node_labels, and output.

使用示例:
```r
genes <- c("TP53", "BRCA1", "AKT1", "EGFR")
g <- get_ppi(genes, taxID = 9606)
```

### `interpret_ppi()`

Interpret PPI network analysis results with AI

Interpret PPI network analysis results with AI

Usage:
```r
interpret_ppi(
x,
top_n = 10,
audience = "researcher",
language = c("zh", "en"),
role = NULL,
system = NULL,
prompt = NULL,
model = NULL
)
```

参数解析:
- `x`: An igraph object with vertex attributes from
compute_nodeinfo(), a node-metric data.frame, or the list
output of rank_ppi_nodes().
- `top_n`: Integer. Number of top nodes. Default 10.
- `audience`: Character. Default "researcher".
Also accepts free-text descriptions.
- `language`: Character. Default "zh".
- `role`: Character or NULL. The AI's identity description. Replaces the
default opening line; other constraints are preserved. Ignored when
system is provided.
- `system`: Character or NULL. Full system prompt override.
- `prompt`: Character or NULL. Custom user prompt instruction.
- `model`: Model identifier or NULL.

使用示例:
```r
ppi <- compute_nodeinfo(demo_ppi)
interpret_ppi(ppi)
interpret_ppi(ppi, role = "You are a systems biology expert focusing on hub gene identification.")
```

### `plot_node_heatmap()`

Plot PPI nodes metrics Heatmap

Create a Z-score normalized heatmap of PPI node topological metrics using
ComplexHeatmap.  Provides convenient parameters for the most common
styling adjustments while still accepting any additional argument via
....

Usage:
```r
plot_node_heatmap(
data,
id_col = "name",
select_cols = NULL,
colors = c("#2166AC", "white", "#B2182B"),
cluster_rows = TRUE,
cluster_cols = FALSE,
row_fontsize = 12,
row_fontface = "italic",
col_fontsize = 10,
col_fontface = "bold",
col_rotation = 45,
show_row_names = TRUE,
show_column_names = TRUE,
legend_title = "Z-score",
border_color = "white",
border_width = 1,
...
)

PlotNodeHeatmap(...)
```

参数解析:
- `data`: A data frame containing the node names and metric values.
- `id_col`: Character. The column name serving as row identifiers.
Default is "name".
- `select_cols`: Character vector. Columns to include in the heatmap.
If NULL (default), all numeric columns except id_col are
used.
- `colors`: Character vector of length 3. Colours mapped to Z-score
values -2, 0, and 2. Default is
c("#2166AC", "white", "#B2182B").
- `cluster_rows`: Logical. Perform hierarchical clustering on rows?
Default TRUE.
- `cluster_cols`: Logical. Perform hierarchical clustering on columns?
Default FALSE.
- `row_fontsize`: Numeric. Font size for row (gene) names.
Default 12.
- `row_fontface`: Character. Font face for row names. One of
"plain", "italic", "bold", or
"bold.italic". Default "italic".
- `col_fontsize`: Numeric. Font size for column (metric) names.
Default 10.
- `col_fontface`: Character. Font face for column names. Default
"bold".
- `col_rotation`: Numeric. Rotation angle (degrees) for column labels.
Default 45.
- `show_row_names`: Logical. Show row names? Default TRUE.
- `show_column_names`: Logical. Show column names? Default TRUE.
- `legend_title`: Character. Title for the colour legend. Default
"Z-score".
- `border_color`: Character. Cell border colour. Default "white".
Set to NA to remove borders.
- `border_width`: Numeric. Cell border line width. Default 1.
- `...`: Additional arguments passed to
Heatmap.

使用示例:
```r
data(demo_ppi)
ppi <- compute_nodeinfo(demo_ppi)
rk_res <- rank_ppi_nodes(ppi)[[2]]
selected_cols <- colnames(rk_res)[c(2, 3, 4, 6, 9, 12, 14, 15, 16)]
p1 <- plot_node_heatmap(rk_res, select_cols = selected_cols)
print(p1)
```

### `ppi_knock()`

Evaluate PPI network robustness via drug attack model

Simulates a targeted node-knockout attack on a weighted PPI network and
evaluates the statistical significance of network disruption using a
permutation test.  Four topological metrics are tracked (ASPL, AD, DC, CC),
exactly following the methodology of Xi et al. (2022).

Usage:
```r
ppi_knock(
g,
targets,
n_perm = 100L,
weight_attr = "score",
rewire_niter = 10L,
seed = 42L
)
```

参数解析:
- `g`: An undirected igraph object representing the PPI network.
- `targets`: Character vector of node names to knock out (drug targets).
- `n_perm`: Integer. Number of permutation iterations for the null
distribution.  Default is 100, matching the original paper.
- `weight_attr`: Character. Edge attribute name for confidence weights.
Default is "score".
- `rewire_niter`: Integer. Multiplier for edge-swap attempts per rewiring
step (ecount(g) * rewire_niter swaps).  Default is 10.
- `seed`: Integer. Random seed for reproducibility. Default is 42.

使用示例:
```r
data(demo_ppi)
targets <- "IL6"
res <- ppi_knock(demo_ppi, targets, n_perm = 100)
print(res$Summary)
cat("Total Score:", res$Total_Score, "\n")
cat("Total P-value:", res$Total_Pvalue, "\n")
```

### `ppi_subset()`

Subset PPI network by edge score and top-n node degree

This function filters a PPI network in two steps:

Removes edges with a score below score_cutoff.
(Optional) Keeps only the top n nodes with the highest degree.

Usage:
```r
ppi_subset(
ppi_obj,
n = NULL,
score_cutoff = 0.7,
edge_attr = "score",
rm_isolates = TRUE
)
```

参数解析:
- `ppi_obj`: An igraph object of PPI network.
- `n`: Integer. Number of top-degree nodes to keep. If NULL (default), no degree filtering is performed.
- `score_cutoff`: Numeric. Minimum edge score to keep. Default is 0.7.
- `edge_attr`: Character. The name of the edge attribute representing confidence score. Default is "score" (keep the same as STRING).
- `rm_isolates`: Logical. Whether to remove isolated nodes (degree = 0) after edge filtering. Default is TRUE.

使用示例:
```r
data(demo_ppi)
library(igraph)
ppi.subset <- ppi_subset(demo_ppi, score_cutoff = 0.7)
str(ppi.subset)
```

### `radar_plot()`

Plot a radar chart for node centrality profile

This function draws a radar-style polygon plot for a single node,
showing its values on multiple centrality metrics.

Usage:
```r
radar_plot(
profile_df,
category_col = "metric",
value_col = "value",
fill_color = "#B6B0D4",
line_color = "#6D65A6",
title = NULL
)
```

参数解析:
- `profile_df`: A data frame containing at least one column of metric
names and one column of numeric values.
- `category_col`: Character string; column name in profile_df
for metric labels. Defaults to "metric".
- `value_col`: Character string; column name in profile_df
for metric values (usually scaled to [0, 1]). Defaults to "value".
- `fill_color`: Polygon fill color. "#A3BEDD", "#A7CBA9", and "#D59390" are recommended.
- `line_color`: Polygon border color.
- `title`: Character; plot title. Default is NULL.

使用示例:
```r
# Fill required arguments according to the parameter list above
radar_plot(
profile_df,
category_col = "metric",
value_col = "value",
fill_color = "#B6B0D4",
line_color = "#6D65A6",
title = NULL
)
```

### `rank_ppi_nodes()`

Rank PPI nodes by integrated network centrality

Rank PPI nodes by integrated network centrality

Usage:
```r
rank_ppi_nodes(
g,
metrics = c("degree", "betweenness", "closeness", "eccentricity", "radiality",
"Stress", "MCC", "MNC", "DMNC", "BN", "EPC"),
weights = NULL,
use_weight = TRUE,
na_rm = TRUE
)
```

参数解析:
- `g`: An igraph object that has already been processed by compute_nodeinfo() or have added node metrics manually.
- `metrics`: Character vector; which vertex attributes to use for scoring. Defaults are the same as Cytohubba.
- `weights`: Numeric vector of the same length as metrics; relative weights for each metric. If NULL (default), all metrics are equally weighted.
- `use_weight`: Logical; whether use weighted metrics(beweenness and closeness) instead. Default is TRUE.
- `na_rm`: Logical; if TRUE, NAs are ignored in normalization (set to 0.5).

使用示例:
```r
data(demo_ppi)
library(igraph)
ppi <- compute_nodeinfo(demo_ppi)
rank_res <- rank_ppi_nodes(ppi)

# update igraph object and extract rank info
ppi <- rank_res[[1]]
rank_df <- rank_res[[2]]
print(rank_df)
```

### `run_fastgreedy()`

Compute Fast Greedy Clustering for PPI Network

This function performs fast greedy modularity optimization on an igraph object.
It is suitable for detecting communities in undirected PPI networks.

Usage:
```r
run_fastgreedy(g, weights = NULL)
```

参数解析:
- `g`: An igraph object.
- `weights`: Numeric vector, character, NULL, or NA. Edge weights to use for
clustering. If NULL (default), the function attempts to use the weight
edge attribute first, then the score edge attribute. If a character
value is supplied, it is treated as the edge attribute name. Set to NA to
perform unweighted clustering.

使用示例:
```r
data(demo_ppi)
library(igraph)
ppi <- run_fastgreedy(demo_ppi)
head(V(ppi)$fastgreedy_cluster)
```

### `run_louvain()`

Compute Louvain Clustering for PPI Network

This function performs Louvain clustering (multi-level modularity optimization)
on an igraph object. It is suitable for detecting larger, macro-scale functional
modules in PPI networks.

Usage:
```r
run_louvain(g, resolution = 1, weights = NULL)
```

参数解析:
- `g`: An igraph object.
- `resolution`: Numeric. The resolution parameter controls the granularity of clusters. Default is 1.0 (standard modularity).
- `weights`: Numeric vector or NULL. Edge weights to use for clustering.
If NULL (default), the function attempts to use the 'weight' or 'score' edge attribute.
Set to NA to perform unweighted clustering.

使用示例:
```r
data(demo_ppi)
library(igraph)
ppi <- run_louvain(demo_ppi, resolution = 1)
head(V(ppi)$louvain_cluster)
```

### `run_MCL()`

Perform Markov Clustering (MCL) on a Graph

This function implements the Markov Clustering (MCL) algorithm for detecting
communities (clusters) in a graph. MCL simulates random walks within the graph
by alternating between two operations: expansion and inflation. It is particularly
efficient for biological networks.

Usage:
```r
run_MCL(g, inflation = 2.5, max_iter = 100, pruning = 1e-05, allow1 = FALSE)
```

参数解析:
- `g`: An igraph object. The graph to be clustered. It can be directed or undirected.
- `inflation`: Numeric. Controls cluster granularity. Higher values (e.g., > 2) yield smaller, tighter clusters; lower values yield larger clusters. Default is 2.5.
- `max_iter`: Integer. The maximum number of iterations to perform if convergence is not reached. Default is 100.
- `pruning`: Numeric. A threshold for pruning small values in the matrix to zero.
This preserves the sparsity of the matrix and significantly speeds up computation
while saving memory. Default is 1e-5.
- `allow1`: Logical. If TRUE, clusters with only 1 node are kept as unique clusters. If FALSE,
cluster of size 1 are interpreted as background noise and grouped in one cluster. Default is FALSE.

使用示例:
```r
library(igraph)
g <- make_graph("Zachary")
g <- run_MCL(g, inflation = 2.5)
print(head(V(g)$MCL_cluster))

# Visualize
plot(g,
vertex.color = V(g)$MCL_cluster,
vertex.size = 15,
vertex.label = V(g)$name)
```

### `run_mcode()`

Molecular Complex Detection (MCODE) method for identifying PPI clusters.

Molecular Complex Detection (MCODE) method for identifying PPI clusters.

Usage:
```r
run_mcode(
g,
vwp = 0.2,
degree_cutoff = 2,
k_core_threshold = 2,
haircut = TRUE,
fluff = FALSE,
fdt = 0.1,
loops = FALSE,
max_depth = 100
)

runMCODE(...)
```

参数解析:
- `g`: An igraph object (PPI network).
- `vwp`: Numeric. Node Score Cutoff (Vertex Weight Percentage). Default 0.2.
- `degree_cutoff`: Numeric. Nodes with degree < cutoff will get a score of 0. Default 2.
- `k_core_threshold`: Numeric. Filters out clusters that do not contain a k-core of at least this level. Default 2.
- `haircut`: Logical. Whether to prune singly-connected nodes. Default TRUE.
- `fluff`: Logical. Whether to expand cluster by adding dense neighbors. Default FALSE.
- `fdt`: Numeric. Fluff Node Density Cutoff. Used if fluff=TRUE. Default 0.1.
- `loops`: Logical. Whether to include self-loops in scoring. Default FALSE.
- `max_depth`: Numeric. Maximum recursion depth for cluster finding (to prevent stack overflow on huge networks). Default 100.
- `...`: Additional arguments passed to run_mcode when using the deprecated alias runMCODE.

使用示例:
```r
data(demo_ppi)
ppi <- run_mcode(demo_ppi, max_depth = 100)
str(ppi)
```

## Artifact management

### `artifact_exists()`

Check if artifact exists

Check if artifact exists

Usage:
```r
artifact_exists(artifact_id)
```

参数解析:
- `artifact_id`: Character. The artifact identifier.

使用示例:
```r
# Fill required arguments according to the parameter list above
artifact_exists(artifact_id)
```

### `clear_tcm_artifacts()`

Clear all artifacts from registry

Removes all stored artifacts. Use with caution.

Usage:
```r
clear_tcm_artifacts()
```

参数解析:
- No user-facing parameters; call exactly as shown in Usage.

使用示例:
```r
# Fill required arguments according to the parameter list above
clear_tcm_artifacts()
```

### `list_tcm_artifacts()`

List all TCM artifacts in the registry

Returns a data.frame summarizing all stored artifacts.

Usage:
```r
list_tcm_artifacts()
```

参数解析:
- No user-facing parameters; call exactly as shown in Usage.

使用示例:
```r
list_tcm_artifacts()
```

### `load_tcm_artifact()`

Load a TCM artifact from the registry

Retrieves the actual R object by its artifact_id.

Usage:
```r
load_tcm_artifact(artifact_id)
```

参数解析:
- `artifact_id`: Character. The artifact identifier.

使用示例:
```r
enrich_obj <- load_tcm_artifact("enrich_001")
```

### `save_tcm_artifact()`

Save a TCM artifact to the registry

Stores a complex R object (igraph, enrichResult, data.frame, etc.) and
returns a lightweight handle that can be passed to the LLM.

Usage:
```r
save_tcm_artifact(
object,
artifact_type,
artifact_id = NULL,
summary = NULL,
provenance = NULL
)
```

参数解析:
- `object`: The R object to store.
- `artifact_type`: Character. Type label, e.g. "enrichment_result",
"herb_graph", "ppi_graph", "ranking_table".
- `artifact_id`: Character or NULL. If NULL, auto-generated.
- `summary`: Character or NULL. Human-readable summary for LLM context.
- `provenance`: List or NULL. Metadata about how this artifact was created.

使用示例:
```r
# Save an enrichment result
enrich_res <- herb_enricher(genes = my_genes)
artifact <- save_tcm_artifact(
object = enrich_res,
artifact_type = "enrichment_result",
summary = "Herb enrichment on 126 genes; 14 terms passed cutoff."
)
# artifact$artifact_id can be passed to next tool call
```

### `summarize_tcm_artifact()`

Summarize a TCM artifact for LLM context

Generates a concise text summary suitable for passing to the LLM.

Usage:
```r
summarize_tcm_artifact(artifact_id, max_lines = 5L)
```

参数解析:
- `artifact_id`: Character. The artifact identifier.
- `max_lines`: Integer. Maximum lines for preview tables.

使用示例:
```r
summarize_tcm_artifact("enrich_001")
```

## Molecule annotation and docking support

### `compound_similarity()`

Similarity search (Tanimoto) on PubChem

Similarity search (Tanimoto) on PubChem

Usage:
```r
compound_similarity(
query,
from = c("smiles", "cid", "inchikey", "name"),
threshold = 90,
topn = 10,
fetch_factor = 3,
compute_score = TRUE,
...
)
```

参数解析:
- `query`: identifier (CID/SMILES/InChIKey/Name)
- `from`: one of c("cid","smiles","inchikey","name")
- `threshold`: integer 0–100 (percent), e.g. 90
- `topn`: max records to return. Default 10.
- `fetch_factor`: The multiplier for overfetching candidate compounds (default = 5).
- `compute_score`: Logical. If TRUE, use rcdk to compute Tanimoto score locally. Set to FALSE if rJava causes crashes.
- `...`: Additional arguments passed to internal helper functions.

使用示例:
```r
sim_hits <- compound_similarity(
query = "996",
from = "cid",
compute_score = TRUE)
```

### `convert_structure()`

Convert molecular structure files between formats

Convert molecular structure format via Open Babel

Usage:
```r
convert_structure(
input_file,
output_file = NULL,
input_type = NULL,
output_type = "pdb",
add_hydrogens = FALSE,
gen3d = FALSE,
overwrite = FALSE,
quiet = FALSE
)
```

参数解析:
- `input_file`: Path to input structure file.
- `output_file`: Output path, or NULL to auto-derive.
- `input_type`: Input format (e.g. "sdf"), or NULL
to detect from extension.
- `output_type`: Output format. Default "pdb".
- `add_hydrogens`: Add hydrogens? Default FALSE.
- `gen3d`: Generate 3D coordinates? Default FALSE.
- `overwrite`: Overwrite existing file? Default FALSE.
- `quiet`: Suppress messages? Default FALSE.

使用示例:
```r
sdf <- download_ligand_structure(2244, destdir = tempdir())
convert_structure(sdf)                          # SDF -> PDB
convert_structure(sdf, output_type = "mol2")    # SDF -> MOL2
```

### `download_ligand_structure()`

Download ligand SDF from PubChem

Download ligand SDF from PubChem

Usage:
```r
download_ligand_structure(
cid,
type = c("auto", "3d", "2d"),
destdir = ".",
filename = NULL,
overwrite = FALSE,
quiet = FALSE
)
```

参数解析:
- `cid`: PubChem Compound ID.
- `type`: "auto" (default, 3D then 2D), "3d", or "2d".
- `destdir`: Save directory. Default ".".
- `filename`: Output name, or NULL for auto-naming.
- `overwrite`: Overwrite existing file? Default FALSE.
- `quiet`: Suppress messages? Default FALSE.

使用示例:
```r
download_ligand_structure(2244)
download_ligand_structure(2244, type = "3d", destdir = tempdir())
```

### `download_receptor_structure()`

Download receptor structure from RCSB PDB

Download receptor structure from RCSB PDB

Usage:
```r
download_receptor_structure(
pdb_id,
format = c("pdb", "cif"),
destdir = ".",
filename = NULL,
overwrite = FALSE,
quiet = FALSE
)
```

参数解析:
- `pdb_id`: 4-character PDB ID (e.g. "4hhb"). Case-insensitive.
- `format`: "pdb" (default) or "cif".
- `destdir`: Save directory. Default ".".
- `filename`: Output name, or NULL for auto-naming.
- `overwrite`: Overwrite existing file? Default FALSE.
- `quiet`: Suppress messages? Default FALSE.

使用示例:
```r
download_receptor_structure("4hhb")
download_receptor_structure("4hhb", format = "cif", destdir = tempdir())
```

### `getcid()`

Get CID for compounds from PubChem

Retrieve PubChem Compound IDs (CIDs) for a vector of compound names.
This function wraps webchem::get_cid() with built-in rate limiting
(≤5 requests per second as recommended by PubChem).

Usage:
```r
getcid(
compound,
from = c("name", "smiles", "inchi"),
match = c("first", "best", "all"),
pause = 0.25,
quiet = TRUE,
...
)
```

参数解析:
- `compound`: A character vector of compound names.
- `from`: Source type for lookup, e.g. "name", "smiles", "inchi".
- `match`: Match mode, one of "first", "best", or "all".
- `pause`: Pause time (in seconds) between requests. Default is 0.25s (≈4 req/s).
- `quiet`: Logical, whether to suppress messages. Default TRUE.
- `...`: Additional arguments passed to internal helper functions.

使用示例:
```r
res_smiles <- getcid("CCO", from = "smiles")
print(res_smiles)
```

### `getprops()`

get properties for compounds from PubChem

get properties for compounds from PubChem

Usage:
```r
getprops(
cid,
properties = c("MolecularFormula", "MolecularWeight", "IUPACName", "CanonicalSMILES",
"InChIKey", "XLogP"),
...
)
```

参数解析:
- `cid`: Integer/numeric/character vector of PubChem CIDs.
- `properties`: Character vector of PubChem property keys to request.
- `...`: Additional arguments passed to internal helper functions.

使用示例:
```r
herbs <- c("灵芝")
lz <- search_herb(herb = herbs, type = "Herb_cn_name")
lz_mol <- sample(unique(lz$molecule), 5, replace = FALSE)
lz_mol_cid <- resolve_cid(lz_mol, from = "name")
props <- getprops(lz_mol_cid)
print(props)
```

### `ggdock()`

Plot molecule–target docking affinity heatmaps

Plot molecule–target docking affinity heatmaps

Usage:
```r
ggdock(
dock_data,
order = NULL,
type = "dot",
point_size = 8,
base_size = 12,
angle = 50,
hjust = 1,
vjust = 1,
palette = "ggthemes::Orange-Blue-White Diverging",
label = FALSE,
label_digits = 2,
label_size = 3,
label_color = "black",
label_family = "sans",
label_fontface = "plain",
...
)
```

参数解析:
- `dock_data`: A numeric matrix or data frame.
- `order`: Reordering method for factor levels (optional).
- `type`: Type of visualization, either "dot" or "tile".
- `point_size`: Numeric. Point size used in dot plots. Default is 6.
- `base_size`: Numeric. Base font size for the theme. Default is 12.
- `angle`: Numeric. Rotation angle for x-axis text labels. Default is 50.
- `hjust`: Numeric. Horizontal justification for x-axis labels.
- `vjust`: Numeric. Vertical justification for x-axis labels.
- `palette`: Character. Name of the continuous color palette to use. Default is "ggthemes::Orange-Blue-White Diverging".
- `label`: Logical. If TRUE, print affinity numbers on the plot. Default FALSE.
- `label_digits`: Integer. Number of digits to show for the affinity text. Default 2.
- `label_size`: Numeric. Text size for affinity labels. Default 3.5.
- `label_color`: Character. Text color for affinity labels. Default "black".
- `label_family`: Character. Font family for labels (e.g., "sans", "Times"). Default "sans".
- `label_fontface`: Character. Font face for labels: "plain", "bold", "italic", "bold.italic". Default "plain".
- `...`: Additional parameters passed to geom_point() or geom_tile().

使用示例:
```r
# Fill required arguments according to the parameter list above
ggdock(
dock_data,
order = NULL,
type = "dot",
point_size = 8,
base_size = 12,
angle = 50,
hjust = 1,
vjust = 1,
palette = "ggthemes::Orange-Blue-White Diverging",
label = FALSE,
label_digits = 2,
label_size = 3,
label_color = "black",
label_family = "sans",
label_fontface = "plain",
...
)
```

### `resolve_cid()`

Resolve arbitrary identifiers to PubChem CIDs

Resolve user-provided identifiers (CID/SMILES/InChI/InChIKey/Name)
to PubChem Compound IDs (CIDs). For structural inputs (SMILES/InChI),
it prefers POSTing the value in the request body to avoid URL length/encoding issues,
and falls back to a GET path-style endpoint if needed.

Usage:
```r
resolve_cid(x, from = c("cid", "smiles", "inchi", "inchikey", "name"))
```

参数解析:
- `x`: Character vector. Values can be CIDs, SMILES, InChI, InChIKey, or names.
- `from`: One of c("cid","smiles","inchi","inchikey","name").

使用示例:
```r
herbs <- c("灵芝")
lz <- search_herb(herb = herbs, type = "Herb_cn_name")
lz_mol <- sample(unique(lz$molecule), 5, replace = FALSE)
lz_mol_cid <- resolve_cid(lz_mol, from = "name")
print(lz_mol_cid)
```

### `sdf_to_pdb()`

Convert SDF file to PDB format

Backward-compatible wrapper around convert_structure for
the common SDF-to-PDB conversion workflow.

Usage:
```r
sdf_to_pdb(
sdf_file,
pdb_file = NULL,
add_hydrogens = TRUE,
gen3d = FALSE,
overwrite = FALSE,
quiet = FALSE
)
```

参数解析:
- `sdf_file`: Path to input SDF file.
- `pdb_file`: Output PDB path, or NULL to auto-derive.
- `add_hydrogens`: Add hydrogens? Default TRUE.
- `gen3d`: Generate 3D coordinates? Default FALSE.
- `overwrite`: Overwrite existing file? Default FALSE.
- `quiet`: Suppress messages? Default FALSE.

使用示例:
```r
sdf <- download_ligand_structure(2244, destdir = tempdir())
pdb <- sdf_to_pdb(sdf)
```

## AI and agent layer

### `create_tcm_agent()`

Create a generic aisdk agent

A thin wrapper around aisdk::create_agent().
Unlike create_tcm_task_agent, this creates a generic agent
without TCM-specific tools or system prompt.
Pass the result to run_tcm_agent to call it.

Usage:
```r
create_tcm_agent(name, description, system_prompt)
```

参数解析:
- `name`: Character. A short identifier for the agent.
- `description`: Character. One-line description of the agent's role.
- `system_prompt`: Character. The system prompt defining agent behaviour.

使用示例:
```r
bio_agent <- create_tcm_agent(
name = "BioInterpreter",
description = "Bioinformatics result interpreter",
system_prompt = "You are a bioinformatics expert. Explain results
concisely in 3 sentences for a research report."
)
```

### `create_tcm_ml_list()`

Assemble individually fitted ML models into a tcm_ml_list

Assemble individually fitted ML models into a tcm_ml_list

Usage:
```r
create_tcm_ml_list(...)
```

参数解析:
- `...`: Named or unnamed tcm_ml objects. Names become the method labels used in downstream plots. At least one object must be supplied.

使用示例:
```r
ml_data <- prepare_ml_data(expr_mat, group)
res_lasso <- ml_lasso(ml_data)
res_rf <- ml_rf(ml_data)
res_svm <- ml_svm_rfe(ml_data)

## assemble into a tcm_ml_list
ml_list <- create_tcm_ml_list(
lasso   = res_lasso,
rf      = res_rf,
svm_rfe = res_svm)
```

### `create_tcm_task_agent()`

Create a TCM task agent

Creates an aisdk agent configured for TCM network pharmacology analysis.
The agent can call TCMDATA tools natively.

Usage:
```r
create_tcm_task_agent(tools = NULL, system_prompt = NULL, skills = NULL)
```

参数解析:
- `tools`: List of Tool objects. If NULL, uses default tools from
create_tcm_tools.
- `system_prompt`: Character. Custom system prompt. If NULL, uses default.
- `skills`: Character vector of skill directories. If NULL, uses
the active TCMDATA skill directory plus aisdk's
skill-creator skill (if available). Use character(0) to
disable skills entirely.

使用示例:
```r
agent <- create_tcm_task_agent()
result <- run_tcm_task(agent, "Search targets of Astragalus")
```

### `create_tcm_tools()`

Create all TCM analysis tools

Returns a list of aisdk Tool objects that wrap TCMDATA analysis functions.
These tools can be passed to an aisdk agent for native tool calling.

Usage:
```r
create_tcm_tools(task_type = NULL, tool_names = NULL)
```

参数解析:
- `task_type`: Character or NULL. If provided, only returns tools relevant
to this task type (for example, "enrichment", "ppi_analysis",
or "ml_screening"). If NULL, returns all tools.
- `tool_names`: Character vector or NULL. If provided, overrides task_type
and returns only tools whose names are in this vector. Useful when the
router merges tools from multiple matched task types.

使用示例:
```r
tools <- create_tcm_tools()
tools <- create_tcm_tools(task_type = "enrichment")
```

### `create_tcm_workflow()`

Create a TCM analysis workflow

Defines a multi-step analysis workflow that can be executed automatically.

Usage:
```r
create_tcm_workflow(name, steps)
```

参数解析:
- `name`: Character. Workflow name.
- `steps`: List of step definitions, each with 'tool' and 'params'.

使用示例:
```r
workflow <- create_tcm_workflow(
name = "herb_to_hub",
steps = list(
list(tool = "search_herb_records",
params = list(herb = "{{herb}}", type = "Herb_pinyin_name")),
list(tool = "run_herb_enrichment",
params = list(genes = "{{from_previous}}")),
list(tool = "interpret_artifact",
params = list(artifact_id = "{{from_previous}}"))
)
)
run_tcm_workflow(workflow, herb = "Astragalus")
```

### `make_tcm_function()`

Create a reusable AI function from an instruction

A function factory: supply an instruction (system prompt) once, get back a
plain R function. The returned function accepts a query string and behaves
like any ordinary R function -- it works directly in sapply(),
purrr::map(), mutate(), and other vectorised workflows.

Usage:
```r
make_tcm_function(instruction, name = "CustomFn")
```

参数解析:
- `instruction`: Character. The system prompt that defines the function's
behaviour (its "role" and task description).
- `name`: Character. A short identifier used for the underlying agent.
Useful for debugging. Default "CustomFn".

使用示例:
```r
explain_pathway <- make_tcm_function(
"You are an MSigDB/KEGG pathway expert. Explain the biological
function of the pathway in 2 sentences for a research report."
)
explain_pathway("HALLMARK_INFLAMMATORY_RESPONSE")

pathways <- c("HALLMARK_E2F_TARGETS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
explanations <- sapply(pathways, explain_pathway)
```

### `route_tcm_task()`

Route a task to the appropriate tool set

Analyzes the user's task description and determines which tools should be
made available to the agent. Uses rule-based routing first, then optionally
falls back to LLM routing for ambiguous cases.

Usage:
```r
route_tcm_task(task, use_llm = FALSE, model = NULL)
```

参数解析:
- `task`: Character. The user's task description.
- `use_llm`: Logical. Whether to use LLM for ambiguous routing.
- `model`: Model object for LLM routing when use_llm = TRUE.

使用示例:
```r
route_tcm_task("Run GO enrichment on the target genes")
route_tcm_task("Retrieve a PPI network and rank hub genes")
```

### `run_tcm_agent()`

Run an agent as a plain function

Wraps agent$run() so the agent behaves like an ordinary R function.
Prints the response via cat() and returns it invisibly.

Usage:
```r
run_tcm_agent(agent, query, model = NULL, verbose = TRUE)
```

参数解析:
- `agent`: An agent object from create_tcm_agent.
- `query`: Character. The query or task to send to the agent.
- `model`: A model identifier, LanguageModelV1 object, or NULL (uses the
package-wide default set via aisdk::set_model()).
- `verbose`: Logical. If TRUE (default), print the response to the
console. Set FALSE for programmatic / batch use.

使用示例:
```r
bio_agent <- create_tcm_agent(
name = "BioInterpreter",
description = "Bioinformatics result interpreter",
system_prompt = "You are a bioinformatics expert."
)
run_tcm_agent(bio_agent, "Explain HALLMARK_E2F_TARGETS")
```

### `run_tcm_task()`

Run a task with a TCM agent

Executes a natural language task using the provided agent.

Usage:
```r
run_tcm_task(agent, task, model = NULL, verbose = TRUE)
```

参数解析:
- `agent`: An Agent object from create_tcm_task_agent.
- `task`: Character. The task description in natural language.
- `model`: Model object or NULL (uses package default).
- `verbose`: Logical. Print progress and results (default TRUE).

使用示例:
```r
agent <- create_tcm_task_agent()
result <- run_tcm_task(agent, "Run herb enrichment on Astragalus targets")
```

### `run_tcm_workflow()`

Run a TCM analysis workflow

Executes a predefined workflow step by step.

Usage:
```r
run_tcm_workflow(workflow, ..., verbose = TRUE)
```

参数解析:
- `workflow`: A workflow object from create_tcm_workflow.
- `...`: Named parameters to substitute into step definitions.
- `verbose`: Logical. Print progress (default TRUE).

使用示例:
```r
# Fill required arguments according to the parameter list above
run_tcm_workflow(workflow, ..., verbose = TRUE)
```

### `tcm_agent()`

One-step TCM task execution

The main entry point for natural language TCM analysis. Automatically
routes the task, creates appropriate agent, and executes the task.

Usage:
```r
tcm_agent(task, model = NULL, verbose = TRUE, use_router = TRUE)
```

参数解析:
- `task`: Character. Natural language task description.
- `model`: Model object or NULL (uses package default).
- `verbose`: Logical. Print progress and results (default TRUE).
- `use_router`: Logical. Use task router to select tools (default TRUE).

使用示例:
```r
# Basic usage
tcm_agent("Search targets of Astragalus")

# Herb enrichment workflow
tcm_agent("Run herb enrichment on Astragalus targets and find top hub genes")

# PPI analysis
tcm_agent("Analyse PPI network for diabetes-related genes and find hub genes")

# Interpretation
tcm_agent("Interpret the previous enrichment result")
```

### `tcm_aisdk_skill()`

Get path to a bundled aisdk skill

Returns the absolute path to a skill bundled with the aisdk
package, such as "skill-creator". This is useful when combining
TCMDATA's own skills with upstream aisdk skills without copying them into
the TCMDATA package.

Usage:
```r
tcm_aisdk_skill(name = "skill-creator")
```

参数解析:
- `name`: Character. Skill name in aisdk/inst/skills.
Default "skill-creator".

使用示例:
```r
skill_creator <- tcm_aisdk_skill()
agent <- create_tcm_task_agent(
skills = c(tcm_skill_dir(), skill_creator)
)
```

### `tcm_chat()`

Interactive TCM analysis chat

Starts a rich interactive chat session for multi-turn TCM analysis
conversations in the terminal. Each exchange is formatted with clear
visual structure including routing information, tool call logs, artifact
updates, and the agent response.

Usage:
```r
tcm_chat(model = NULL, verbose = TRUE, stream = TRUE, skills = NULL)
```

参数解析:
- `model`: Model object or NULL.
- `verbose`: Logical. Print agent responses (default TRUE).
- `stream`: Logical. Enable streaming output (default TRUE). Toggle
at runtime with the /stream command.
- `skills`: Character vector of skill directories to load for the chat
session. Default NULL disables aisdk skills for chat stability and
uses TCMDATA's built-in tool-oriented system prompt. Pass explicit skill
paths, for example c(tcm_skill_dir(), tcm_aisdk_skill()), to enable
skill loading.

使用示例:
```r
tcm_chat()
# With streaming disabled
tcm_chat(stream = FALSE)
# Add aisdk's skill-creator alongside TCMDATA skills
tcm_chat(skills = c(tcm_skill_dir(), tcm_aisdk_skill()))
```

### `tcm_config()`

Write AI provider credentials to .env

Saves TCM_PROVIDER, TCM_API_KEY, TCM_MODEL, and
optionally TCM_BASE_URL to a .env file. Existing values for
these four keys are overwritten; all other lines are preserved.

Usage:
```r
tcm_config(provider, api_key, model, base_url = NULL, path = ".env")
```

参数解析:
- `provider`: Character. Provider name (see Details).
- `api_key`: Character. Your API key.
- `model`: Character. Model name, e.g. "gpt-4o-mini",
"claude-3-5-haiku-20241022", "gemini-2.0-flash".
- `base_url`: Character or NULL. Override the default API endpoint.
Required for proxies or self-hosted endpoints.
- `path`: Character. Path to the .env file. Default ".env".

使用示例:
```r
tcm_config("openai",    "sk-xxx",     "gpt-4o-mini")
tcm_config("anthropic", "sk-ant-xxx", "claude-3-5-haiku-20241022")
tcm_config("gemini",    "AIza-xxx",   "gemini-2.0-flash")
tcm_config("deepseek",  "sk-xxx",     "deepseek-chat")
tcm_config("openai",    "sk-xxx",     "gpt-5-minimal",
base_url = "https://www.packyapi.com/v1")
```

### `tcm_field_array()`

Build a custom output schema for tcm_interpret_schema()

A family of thin wrappers around aisdk::z_object() and the
associated field constructors. These functions let you define structured
output contracts without writing aisdk:: anywhere in your code.

Usage:
```r
tcm_schema(...)

tcm_field_string(description = "")

tcm_field_number(description = "")

tcm_field_boolean(description = "")

tcm_field_array(description = "")

tcm_field_enum(choices, description = "")
```

参数解析:
- `...`: Named field definitions created with tcm_field_*().
- `description`: Character. Guidance for the model about this field.
Default "" (no description).
- `choices`: Character vector of allowed values (for tcm_field_enum
only).

使用示例:
```r
my_schema <- tcm_schema(
summary     = tcm_field_string("2-3 sentence overview"),
mechanism   = tcm_field_string("Key molecular mechanism"),
key_targets = tcm_field_array("Top gene or protein targets"),
score       = tcm_field_number("Confidence score 0-1"),
is_tcm      = tcm_field_boolean("Whether TCM herbs are involved"),
confidence  = tcm_field_enum(c("high", "medium", "low"))
)

res <- tcm_interpret_schema(enrich_res, schema = my_schema,
language = "zh")
res$output$summary
res$output$key_targets   # character vector
res$output$confidence    # one of "high" / "medium" / "low"
```

### `tcm_field_boolean()`

Build a custom output schema for tcm_interpret_schema()

A family of thin wrappers around aisdk::z_object() and the
associated field constructors. These functions let you define structured
output contracts without writing aisdk:: anywhere in your code.

Usage:
```r
tcm_schema(...)

tcm_field_string(description = "")

tcm_field_number(description = "")

tcm_field_boolean(description = "")

tcm_field_array(description = "")

tcm_field_enum(choices, description = "")
```

参数解析:
- `...`: Named field definitions created with tcm_field_*().
- `description`: Character. Guidance for the model about this field.
Default "" (no description).
- `choices`: Character vector of allowed values (for tcm_field_enum
only).

使用示例:
```r
my_schema <- tcm_schema(
summary     = tcm_field_string("2-3 sentence overview"),
mechanism   = tcm_field_string("Key molecular mechanism"),
key_targets = tcm_field_array("Top gene or protein targets"),
score       = tcm_field_number("Confidence score 0-1"),
is_tcm      = tcm_field_boolean("Whether TCM herbs are involved"),
confidence  = tcm_field_enum(c("high", "medium", "low"))
)

res <- tcm_interpret_schema(enrich_res, schema = my_schema,
language = "zh")
res$output$summary
res$output$key_targets   # character vector
res$output$confidence    # one of "high" / "medium" / "low"
```

### `tcm_field_enum()`

Build a custom output schema for tcm_interpret_schema()

A family of thin wrappers around aisdk::z_object() and the
associated field constructors. These functions let you define structured
output contracts without writing aisdk:: anywhere in your code.

Usage:
```r
tcm_schema(...)

tcm_field_string(description = "")

tcm_field_number(description = "")

tcm_field_boolean(description = "")

tcm_field_array(description = "")

tcm_field_enum(choices, description = "")
```

参数解析:
- `...`: Named field definitions created with tcm_field_*().
- `description`: Character. Guidance for the model about this field.
Default "" (no description).
- `choices`: Character vector of allowed values (for tcm_field_enum
only).

使用示例:
```r
my_schema <- tcm_schema(
summary     = tcm_field_string("2-3 sentence overview"),
mechanism   = tcm_field_string("Key molecular mechanism"),
key_targets = tcm_field_array("Top gene or protein targets"),
score       = tcm_field_number("Confidence score 0-1"),
is_tcm      = tcm_field_boolean("Whether TCM herbs are involved"),
confidence  = tcm_field_enum(c("high", "medium", "low"))
)

res <- tcm_interpret_schema(enrich_res, schema = my_schema,
language = "zh")
res$output$summary
res$output$key_targets   # character vector
res$output$confidence    # one of "high" / "medium" / "low"
```

### `tcm_field_number()`

Build a custom output schema for tcm_interpret_schema()

A family of thin wrappers around aisdk::z_object() and the
associated field constructors. These functions let you define structured
output contracts without writing aisdk:: anywhere in your code.

Usage:
```r
tcm_schema(...)

tcm_field_string(description = "")

tcm_field_number(description = "")

tcm_field_boolean(description = "")

tcm_field_array(description = "")

tcm_field_enum(choices, description = "")
```

参数解析:
- `...`: Named field definitions created with tcm_field_*().
- `description`: Character. Guidance for the model about this field.
Default "" (no description).
- `choices`: Character vector of allowed values (for tcm_field_enum
only).

使用示例:
```r
my_schema <- tcm_schema(
summary     = tcm_field_string("2-3 sentence overview"),
mechanism   = tcm_field_string("Key molecular mechanism"),
key_targets = tcm_field_array("Top gene or protein targets"),
score       = tcm_field_number("Confidence score 0-1"),
is_tcm      = tcm_field_boolean("Whether TCM herbs are involved"),
confidence  = tcm_field_enum(c("high", "medium", "low"))
)

res <- tcm_interpret_schema(enrich_res, schema = my_schema,
language = "zh")
res$output$summary
res$output$key_targets   # character vector
res$output$confidence    # one of "high" / "medium" / "low"
```

### `tcm_field_string()`

Build a custom output schema for tcm_interpret_schema()

A family of thin wrappers around aisdk::z_object() and the
associated field constructors. These functions let you define structured
output contracts without writing aisdk:: anywhere in your code.

Usage:
```r
tcm_schema(...)

tcm_field_string(description = "")

tcm_field_number(description = "")

tcm_field_boolean(description = "")

tcm_field_array(description = "")

tcm_field_enum(choices, description = "")
```

参数解析:
- `...`: Named field definitions created with tcm_field_*().
- `description`: Character. Guidance for the model about this field.
Default "" (no description).
- `choices`: Character vector of allowed values (for tcm_field_enum
only).

使用示例:
```r
my_schema <- tcm_schema(
summary     = tcm_field_string("2-3 sentence overview"),
mechanism   = tcm_field_string("Key molecular mechanism"),
key_targets = tcm_field_array("Top gene or protein targets"),
score       = tcm_field_number("Confidence score 0-1"),
is_tcm      = tcm_field_boolean("Whether TCM herbs are involved"),
confidence  = tcm_field_enum(c("high", "medium", "low"))
)

res <- tcm_interpret_schema(enrich_res, schema = my_schema,
language = "zh")
res$output$summary
res$output$key_targets   # character vector
res$output$confidence    # one of "high" / "medium" / "low"
```

### `tcm_init_skills()`

Initialize a local skills directory

Copies all bundled skills from the TCMDATA package to a local directory.
Once initialized, the agent will use the local directory instead of the
package defaults, allowing full customization.

Usage:
```r
tcm_init_skills(path = "tcm_skills", overwrite = FALSE)
```

参数解析:
- `path`: Character. Directory to create. Default "tcm_skills".
- `overwrite`: Logical. Overwrite existing directory? Default FALSE.

使用示例:
```r
# Copy all skills to ./tcm_skills/ for customization
tcm_init_skills()

# Then edit the preferences:
# file.edit("tcm_skills/analysis-preferences/SKILL.md")

# Or place new custom skills under tcm_skills/
# aisdk's skill-creator remains available by default
```

### `tcm_interpret()`

Interpret analysis results or free text with AI

A unified entry point that dispatches to the appropriate
interpret_*() function based on the class of x.
When x is a character string (or vector), it routes to
an agent-based text interpreter instead.

Usage:
```r
tcm_interpret(
x,
type = NULL,
top_n = NULL,
max_genes = 5L,
audience = "researcher",
language = c("zh", "en"),
role = NULL,
system = NULL,
prompt = NULL,
model = NULL,
verbose = TRUE
)
```

参数解析:
- `x`: Analysis result object, or a character string / vector
for free-text interpretation.
- `type`: Character or NULL. Overrides auto-detection.
One of "enrichment", "ppi", "table".
Takes priority over class-based dispatch, so a data.frame
enrichment table can be forced with type = "enrichment".
- `top_n`: Integer. Passed to the underlying interpret
function (ignored for text).
- `max_genes`: Integer. Max genes shown per term when
type = "enrichment". Ignored for text, PPI, and generic tables.
- `audience`: Character. Preset: "researcher",
"wetlab", "paper". Or any free-text description
like "clinical doctor, no bioinformatics background".
- `language`: Character. "zh" or "en".
Default "zh".
- `role`: Character or NULL. The AI's identity description. Passed to
the underlying interpret_*() function. Ignored for free-text
input and when system is provided.
- `system`: Character or NULL. Full system prompt override.
- `prompt`: Character or NULL. Custom user prompt instruction.
- `model`: Model identifier or NULL.
- `verbose`: Logical. If TRUE (default), print the response to
the console. Set FALSE for programmatic / batch use to suppress
cat() output. Only affects the free-text path; structured results
are always returned silently and printed via the S3 print method.

使用示例:
```r
tcm_interpret(enrich_res)
tcm_interpret(ppi_graph,
role = "You are a systems biologist focusing on drug target identification.")
# Force enrichment path for a plain data.frame enrichment table
tcm_interpret(enrich_df, type = "enrichment")
tcm_interpret(my_df, type = "table")
# Suppress printing for batch use
results <- tcm_interpret("IL6 logFC=3.2, TNF logFC=2.8", verbose = FALSE)
```

### `tcm_interpret_schema()`

Interpret analysis results with a user-defined output schema

A low-level entry point that keeps all existing context-compression and
prompt logic intact, but replaces the fixed output contract with a schema
supplied by the caller. Useful when you need a different set of fields than
the default tcm_interpret produces.

Usage:
```r
tcm_interpret_schema(
x,
schema,
type = NULL,
top_n = NULL,
max_genes = 5L,
audience = "researcher",
language = c("zh", "en"),
role = NULL,
system = NULL,
prompt = NULL,
model = NULL
)
```

参数解析:
- `x`: Analysis result object. Supported: enrichResult,
enrichment data.frame, igraph, any data.frame.
- `schema`: An aisdk::z_object() defining the output structure.
- `type`: Character or NULL. Override auto-detection.
One of "enrichment", "ppi", "table".
- `top_n`: Integer. Rows / nodes to include in the compressed context.
Default 10 for enrichment/ppi, 20 for table.
- `max_genes`: Integer. Max genes per term (enrichment only). Default 5.
- `audience`: Character. Passed to the default system prompt builder.
Ignored when system is provided.
- `language`: Character. "zh" or "en". Default "zh".
- `role`: Character or NULL. Replaces the AI identity line in the default
system prompt. Ignored when system is provided.
- `system`: Character or NULL. Full system prompt override.
- `prompt`: Character or NULL. Instruction prepended to the context block.
Defaults to "Interpret the following <type> data:".
- `model`: Model identifier or NULL.

使用示例:
```r
# Define a custom schema
my_schema <- aisdk::z_object(
summary     = aisdk::z_string("Brief overview"),
mechanism   = aisdk::z_string("Molecular mechanism"),
key_targets = aisdk::z_array(aisdk::z_string(), "Top targets"),
confidence  = aisdk::z_enum(c("high", "medium", "low"))
)

res <- tcm_interpret_schema(
enrich_res,
schema   = my_schema,
language = "zh",
prompt   = "Focus on anti-inflammatory pathways:"
)

# Access fields by name
res$output$summary
res$output$key_targets   # character vector
res$output$confidence

print(res)               # shows all fields generically
```

### `tcm_reset_skills()`

Reset skills to package defaults

Clears the user skills directory setting so the agent uses the bundled
package skills again.

Usage:
```r
tcm_reset_skills()
```

参数解析:
- No user-facing parameters; call exactly as shown in Usage.

使用示例:
```r
tcm_reset_skills()
```

### `tcm_sankey()`

TCM Sankey Plot for herb-molecule-target

TCM Sankey Plot for herb-molecule-target

Usage:
```r
tcm_sankey(
data,
axis_order = c("herb", "molecule", "target"),
herb_cols = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
"#e377c2", "#7f7f7f"),
mol_cols = c("#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
"#c5b0d5", "#E41A1C", "#377EB8"),
target_cols = c("#c49c94", "#f7b6d2", "#dbdb8d", "#c7e9c0", "#f4cae4", "#e6f598",
"#ffeda0", "#BB5234", "#BB7813", "#FF6158"),
plot_font = "sans",
font_face = "plain",
target_fontface = "italic",
font_size = 3.6,
width = 0.05,
alpha = 0.3,
knot.pos = 0.3
)

TCM_sankey(...)
```

参数解析:
- `data`: A data frame containing at least three columns: herb, molecule, and target.
- `axis_order`: Character vector specifying the order of axes in the Sankey diagram. Default is c("herb", "molecule", "target").
- `herb_cols`: Character vector defining the base color palette for the herb layer.
- `mol_cols`: Character vector defining the base color palette for the molecule.
- `target_cols`: Character vector defining the base color palette for the target.
- `plot_font`: Character string specifying the font family used for text. Default is "sans".
- `font_face`: Character string specifying the font face. Default is "plain".
- `target_fontface`: Character string specifying the font face for target labels (rightmost axis). Default is "italic".
- `font_size`: Numeric value controlling the size of node labels. Default is 3.5.
- `width`: Numeric value controlling the width of both nodes and flows. Default is 0.05.
- `alpha`: Numeric value controlling the transparency of the flows. Default is 0.3.
- `knot.pos`: Numeric value (between 0 and 1) determining the curvature position of flow lines. Default is 0.3.
- `...`: Additional arguments passed to tcm_sankey when using the deprecated alias TCM_sankey.

使用示例:
```r
# Fill required arguments according to the parameter list above
tcm_sankey(
data,
axis_order = c("herb", "molecule", "target"),
herb_cols = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
"#e377c2", "#7f7f7f"),
mol_cols = c("#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
"#c5b0d5", "#E41A1C", "#377EB8"),
target_cols = c("#c49c94", "#f7b6d2", "#dbdb8d", "#c7e9c0", "#f4cae4", "#e6f598",
"#ffeda0", "#BB5234", "#BB7813", "#FF6158"),
plot_font = "sans",
font_face = "plain",
target_fontface = "italic",
font_size = 3.6,
width = 0.05,
alpha = 0.3,
knot.pos = 0.3
)

TCM_sankey(...)
```

### `tcm_schema()`

Build a custom output schema for tcm_interpret_schema()

A family of thin wrappers around aisdk::z_object() and the
associated field constructors. These functions let you define structured
output contracts without writing aisdk:: anywhere in your code.

Usage:
```r
tcm_schema(...)

tcm_field_string(description = "")

tcm_field_number(description = "")

tcm_field_boolean(description = "")

tcm_field_array(description = "")

tcm_field_enum(choices, description = "")
```

参数解析:
- `...`: Named field definitions created with tcm_field_*().
- `description`: Character. Guidance for the model about this field.
Default "" (no description).
- `choices`: Character vector of allowed values (for tcm_field_enum
only).

使用示例:
```r
my_schema <- tcm_schema(
summary     = tcm_field_string("2-3 sentence overview"),
mechanism   = tcm_field_string("Key molecular mechanism"),
key_targets = tcm_field_array("Top gene or protein targets"),
score       = tcm_field_number("Confidence score 0-1"),
is_tcm      = tcm_field_boolean("Whether TCM herbs are involved"),
confidence  = tcm_field_enum(c("high", "medium", "low"))
)

res <- tcm_interpret_schema(enrich_res, schema = my_schema,
language = "zh")
res$output$summary
res$output$key_targets   # character vector
res$output$confidence    # one of "high" / "medium" / "low"
```

### `tcm_setup()`

Initialise the AI model from .env or explicit arguments

Loads .env (if present), resolves TCM_* variables, calls
the matching aisdk::create_*() function, and registers the model via
aisdk::set_model(). All subsequent AI functions then work without
further setup.

Usage:
```r
tcm_setup(
provider = NULL,
api_key = NULL,
model = NULL,
base_url = NULL,
.env = TRUE,
save = FALSE,
test = FALSE,
force_json_schema = TRUE,
skip_internet_check = TRUE
)
```

参数解析:
- `provider`: Character or NULL. Overrides TCM_PROVIDER.
- `api_key`: Character or NULL. Overrides TCM_API_KEY.
- `model`: Character or NULL. Overrides TCM_MODEL.
- `base_url`: Character or NULL. Overrides TCM_BASE_URL.
- `.env`: Logical. Load .env before reading env vars (default TRUE).
- `save`: Logical. If TRUE, also calls tcm_config() to
persist the resolved credentials to .env. Useful for first-time
setup when you want a single call to both initialise and save.
Default FALSE.
- `test`: Logical. If TRUE, sends a minimal test request after setup to
verify the API key and endpoint are reachable. Warnings (not errors) are
issued on failure so the model is still registered. Default FALSE.
- `force_json_schema`: Logical. Default TRUE. When TRUE,
every generate_object() call internally re-passes the output schema
as response_format, which (a) enables native JSON schema output on
models that support the OpenAI structured-output API and (b) forces
OpenAI-compatible proxies to preserve the system message that
contains the schema instruction (some 3rd-party relay APIs silently strip
the system message when response_format is absent).
Set FALSE only when targeting a provider whose API rejects an
unknown response_format field entirely.
- `skip_internet_check`: Logical. Default TRUE. Sets
options(aisdk.skip_internet_check = TRUE) before live requests so
curl::has_internet() false negatives in proxy/VPN environments do
not block otherwise reachable API endpoints.

使用示例:
```r
# Standard two-step workflow
tcm_config("openai", "sk-xxx", "gpt-4o-mini")
tcm_setup()

# One-step: configure + initialise in a single call
tcm_setup("openai", "sk-xxx", "gpt-4o-mini", save = TRUE)

# Verify connectivity after setup
tcm_setup(test = TRUE)

# Override at runtime without touching .env
tcm_setup("deepseek", api_key = "sk-xxx", model = "deepseek-chat")
```

### `tcm_skill_dir()`

Get the active skills directory

Returns the current skills directory path. If a user-local directory has been
set via tcm_init_skills or tcm_use_skills, that
path is returned. Otherwise, returns the package default.

Usage:
```r
tcm_skill_dir()
```

参数解析:
- No user-facing parameters; call exactly as shown in Usage.

使用示例:
```r
tcm_skill_dir()
```

### `tcm_use_skills()`

Set the active skills directory

Points the agent to an existing skills directory. All subsequent calls to
tcm_agent, tcm_chat, and
create_tcm_task_agent will use skills from this directory.

Usage:
```r
tcm_use_skills(path)
```

参数解析:
- `path`: Character. Path to an existing skills directory containing
skill subdirectories with SKILL.md files.

使用示例:
```r
tcm_use_skills("my_project/skills")
```

## Other utilities

### `draft_result_paragraph()`

Draft a publication-ready result paragraph

Accepts either a tcm_ai_analysis object (from any
interpret_* function) or a raw analysis object. When a
raw object is supplied, the context is compressed on the fly.

Usage:
```r
draft_result_paragraph(
x,
type = c("enrichment", "ppi", "table"),
language = c("zh", "en"),
prompt = NULL,
model = NULL
)
```

参数解析:
- `x`: A tcm_ai_analysis object, or a raw analysis
object (enrichResult, igraph, data.frame).
- `type`: Character. Required when x is a raw object.
One of "enrichment", "ppi", "table".
- `language`: Character. "zh" or "en".
Default "zh".
- `prompt`: Character or NULL. Custom drafting instruction. Replaces the
default "Draft a result paragraph from this {type} data:"
opening; the compressed data context is always appended regardless.
Use this to steer tone, focus, or word count, e.g.
"Summarise the following enrichment in no more than 3 sentences,
focusing on immune-related pathways:".
- `model`: Model identifier or NULL.

使用示例:
```r
ai_res <- interpret_enrichment(enrich_res)
draft  <- draft_result_paragraph(ai_res)
cat(draft$draft$paragraph)

# Custom instruction
draft2 <- draft_result_paragraph(ai_res,
prompt = "用不超过3句话总结，聚焦免疫相关通路，语气适合高水平期刊：")
```

### `getGores()`

Prepare input data for enrich_circle_plot

Prepare input data for enrich_circle_plot

Usage:
```r
getGores(x, up_genes = NULL, down_genes = NULL, top = 10)
```

参数解析:
- `x`: enrichResult object or data.frame from clusterProfiler
with columns: ID, ONTOLOGY, BgRatio, p.adjust, geneID, Count, GeneRatio and RichFactor.
- `up_genes`: Character vector of up-regulated genes (same ID space as geneID).
- `down_genes`: Character vector of down-regulated genes (same ID space as geneID).
- `top`: Numeric. Number of pathways to keep per category (by ascending p.adjust). Default 10.

使用示例:
```r
# Fill required arguments according to the parameter list above
getGores(x, up_genes = NULL, down_genes = NULL, top = 10)
```

### `getMCODE_res()`

Extract MCODE analysis results as a data frame

Extract MCODE analysis results as a data frame

Usage:
```r
get_mcode_res(g, only_clusters = FALSE)

getMCODE_res(...)
```

参数解析:
- `g`: An igraph object processed by runMCODE.
- `only_clusters`: Logical. If TRUE, returns only nodes that belong to a cluster. Default FALSE.
- `...`: Additional arguments passed to get_mcode_res when using the deprecated alias getMCODE_res.

使用示例:
```r
data(demo_ppi)
ppi <- runMCODE(demo_ppi, max_depth = 100)
MCODE_res <- get_mcode_res(ppi)
print(head(MCODE_res))
```

### `interpret_table()`

Interpret a data.frame (table) with AI

Accepts any data.frame and serialises it row-by-row into a compact
key=value format before sending to the model.

Usage:
```r
interpret_table(
x,
top_n = 20,
audience = "researcher",
language = c("zh", "en"),
role = NULL,
system = NULL,
prompt = NULL,
model = NULL
)
```

参数解析:
- `x`: A data.frame.
- `top_n`: Integer. Number of top rows to include. Default 20.
- `audience`: Character. Default "researcher".
Also accepts free-text descriptions.
- `language`: Character. Default "zh".
- `role`: Character or NULL. The AI's identity description. Replaces the
default opening line; other constraints are preserved. Ignored when
system is provided.
- `system`: Character or NULL. Full system prompt override.
- `prompt`: Character or NULL. Custom user prompt instruction.
- `model`: Model identifier or NULL.

使用示例:
```r
interpret_table(my_df)
interpret_table(my_df, role = "You are a clinical data analyst.")
```

### `PlotNodeHeatmap()`

Plot PPI nodes metrics Heatmap

Create a Z-score normalized heatmap of PPI node topological metrics using
ComplexHeatmap.  Provides convenient parameters for the most common
styling adjustments while still accepting any additional argument via
....

Usage:
```r
plot_node_heatmap(
data,
id_col = "name",
select_cols = NULL,
colors = c("#2166AC", "white", "#B2182B"),
cluster_rows = TRUE,
cluster_cols = FALSE,
row_fontsize = 12,
row_fontface = "italic",
col_fontsize = 10,
col_fontface = "bold",
col_rotation = 45,
show_row_names = TRUE,
show_column_names = TRUE,
legend_title = "Z-score",
border_color = "white",
border_width = 1,
...
)

PlotNodeHeatmap(...)
```

参数解析:
- `data`: A data frame containing the node names and metric values.
- `id_col`: Character. The column name serving as row identifiers.
Default is "name".
- `select_cols`: Character vector. Columns to include in the heatmap.
If NULL (default), all numeric columns except id_col are
used.
- `colors`: Character vector of length 3. Colours mapped to Z-score
values -2, 0, and 2. Default is
c("#2166AC", "white", "#B2182B").
- `cluster_rows`: Logical. Perform hierarchical clustering on rows?
Default TRUE.
- `cluster_cols`: Logical. Perform hierarchical clustering on columns?
Default FALSE.
- `row_fontsize`: Numeric. Font size for row (gene) names.
Default 12.
- `row_fontface`: Character. Font face for row names. One of
"plain", "italic", "bold", or
"bold.italic". Default "italic".
- `col_fontsize`: Numeric. Font size for column (metric) names.
Default 10.
- `col_fontface`: Character. Font face for column names. Default
"bold".
- `col_rotation`: Numeric. Rotation angle (degrees) for column labels.
Default 45.
- `show_row_names`: Logical. Show row names? Default TRUE.
- `show_column_names`: Logical. Show column names? Default TRUE.
- `legend_title`: Character. Title for the colour legend. Default
"Z-score".
- `border_color`: Character. Cell border colour. Default "white".
Set to NA to remove borders.
- `border_width`: Numeric. Cell border line width. Default 1.
- `...`: Additional arguments passed to
Heatmap.

使用示例:
```r
data(demo_ppi)
ppi <- compute_nodeinfo(demo_ppi)
rk_res <- rank_ppi_nodes(ppi)[[2]]
selected_cols <- colnames(rk_res)[c(2, 3, 4, 6, 9, 12, 14, 15, 16)]
p1 <- plot_node_heatmap(rk_res, select_cols = selected_cols)
print(p1)
```

### `prepare_herb_graph()`

Prepare an igraph object from a Herb–Molecule–Target data frame

Prepare an igraph object from a Herb–Molecule–Target data frame

Usage:
```r
prepare_herb_graph(
df,
n = NULL,
herb_col = "herb",
molecule_col = "molecule",
target_col = "target",
compute_metrics = TRUE
)
```

参数解析:
- `df`: A data.frame with three columns representing herb–molecule–target relationships.
- `n`: Integer. Optional. The number of rows to use from df. Default is 'nrow(df)'
- `herb_col`: Character. The column name corresponding to herbs. Default is "herb".
- `molecule_col`: Character. The column name corresponding to molecules/compounds. Default is "molecule".
- `target_col`: Character. The column name corresponding to targets. Default is "target".
- `compute_metrics`: Logical. Whether to compute and attach node-level metrics. Default is TRUE.

使用示例:
```r
library(igraph)
library(ggplot2)
herbs <- c("灵芝")
lz <- search_herb(herb = herbs, type = "Herb_cn_name")[1:40, ]
res <- prepare_herb_graph(lz)

# visualization
p <- ggtangle::ggplot(res, layout = "kk") +
ggtangle::geom_edge(alpha = 0.4) +
ggrepel::geom_text_repel(aes(label = name), size = 3.5, max.overlaps = 20) +
geom_point(aes(color = type, size = centrality)) +
scale_color_manual(values = c(
"Herb" = "#E41A1C", "Molecule" = "#377EB8", "Target" = "#4DAF4A")) +
theme_void()
print(p)
```

### `runMCODE()`

Molecular Complex Detection (MCODE) method for identifying PPI clusters.

Molecular Complex Detection (MCODE) method for identifying PPI clusters.

Usage:
```r
run_mcode(
g,
vwp = 0.2,
degree_cutoff = 2,
k_core_threshold = 2,
haircut = TRUE,
fluff = FALSE,
fdt = 0.1,
loops = FALSE,
max_depth = 100
)

runMCODE(...)
```

参数解析:
- `g`: An igraph object (PPI network).
- `vwp`: Numeric. Node Score Cutoff (Vertex Weight Percentage). Default 0.2.
- `degree_cutoff`: Numeric. Nodes with degree < cutoff will get a score of 0. Default 2.
- `k_core_threshold`: Numeric. Filters out clusters that do not contain a k-core of at least this level. Default 2.
- `haircut`: Logical. Whether to prune singly-connected nodes. Default TRUE.
- `fluff`: Logical. Whether to expand cluster by adding dense neighbors. Default FALSE.
- `fdt`: Numeric. Fluff Node Density Cutoff. Used if fluff=TRUE. Default 0.1.
- `loops`: Logical. Whether to include self-loops in scoring. Default FALSE.
- `max_depth`: Numeric. Maximum recursion depth for cluster finding (to prevent stack overflow on huge networks). Default 100.
- `...`: Additional arguments passed to run_mcode when using the deprecated alias runMCODE.

使用示例:
```r
data(demo_ppi)
ppi <- run_mcode(demo_ppi, max_depth = 100)
str(ppi)
```

### `select_features()`

Re-select top features from a fitted ML model
After running a model (especially Ridge or XGBoost where all features
are initially retained), inspect result$importance to decide
how many to keep, then call this function to trim the selection.

Re-select top features from a fitted ML model
After running a model (especially Ridge or XGBoost where all features
are initially retained), inspect result$importance to decide
how many to keep, then call this function to trim the selection.

Usage:
```r
select_features(ml_obj, top_n)
```

参数解析:
- `ml_obj`: A tcm_ml object from any ml_* function.
- `top_n`: Integer; the number of top features to keep (ranked by importance).

使用示例:
```r
xgb <- ml_xgboost(ml_data)        # keep all features as default
head(xgb$importance, 20)           # inspect ranking
xgb <- select_features(xgb, 15)   # keep top 15
```

## Machine learning and diagnostic plots

### `get_enet_coefs()`

Extract ranked coefficients from LASSO / Elastic Net / Ridge

Extract ranked coefficients from LASSO / Elastic Net / Ridge

Usage:
```r
get_enet_coefs(
ml_obj,
top_n = NULL,
min_coef = 0,
sort_by = c("abs", "coef"),
include_zero = FALSE
)
```

参数解析:
- `ml_obj`: A tcm_ml object from ml_lasso(), ml_enet(), or ml_ridge().
- `top_n`: Return only top-n genes (optional).
- `min_coef`: Minimum absolute coefficient. Default 0.
- `sort_by`: "abs" or "coef". Default "abs".
- `include_zero`: Logical. Whether to include zero-coefficient genes. Default FALSE.

使用示例:
```r
ml_data <- prepare_ml_data(expr_mat, group)
res <- ml_lasso(ml_data)
get_enet_coefs(res, top_n = 20)
```

### `get_gene_auc()`

Compute single-gene AUC for diagnostic evaluation

For each gene the function treats its expression as a univariate predictor,
computes a ROC curve via pROC, and returns a summary table with AUC,
95 \% DeLong confidence interval, optimal cutoff (Youden's J), sensitivity,
specificity and a one-sided test against AUC = 0.5.

Usage:
```r
get_gene_auc(
genes,
ml_data = NULL,
expr_mat = NULL,
group = NULL,
use = c("auto", "train", "test"),
ci = TRUE,
p_adjust_method = "BH"
)
```

参数解析:
- `genes`: Character vector of gene symbols to evaluate.
- `ml_data`: A tcm_ml_data object from prepare_ml_data
(preferred input).  When provided, expr_mat / group are
ignored.
- `expr_mat`: Numeric matrix (genes \times samples) with
rownames = gene symbols.  Used only when ml_data is
NULL.
- `group`: Factor or character vector of sample labels (length =
ncol(expr_mat), two levels).
- `use`: Which data split to evaluate when using ml_data:
"auto" (default; test set if available, else train),
"train", or "test".
- `ci`: Logical; compute 95 \% DeLong confidence interval?
Default TRUE.
- `p_adjust_method`: Method for multiple-testing correction passed to
p.adjust.  Default "BH" (Benjamini--Hochberg
FDR).  Other common choices: "bonferroni", "holm",
"none".  When only a single gene is supplied no correction is
needed and p_adj equals p_value regardless of method.

使用示例:
```r
consensus <- get_ml_consensus(ml_list, min_methods = 2)
auc_tbl <- get_gene_auc(consensus, ml_data = ml_data)
auc_tbl
```

### `get_ml_consensus()`

Consensus features across ML methods
Genes selected by >= min_methods models.

Consensus features across ML methods
Genes selected by >= min_methods models.

Usage:
```r
get_ml_consensus(ml_list, min_methods = 2)
```

参数解析:
- `ml_list`: A tcm_ml_list or plain list of tcm_ml objects.
- `min_methods`: Minimum agreement count. Default 2.

使用示例:
```r
ml_data <- prepare_ml_data(expr_mat, group)
res_list <- run_ml_screening(ml_data)
get_ml_consensus(res_list, min_methods = 2)
```

### `get_ml_gene_sets()`

Extract gene sets from ML results

Flexibly accepts the output of run_ml_screening(), individually fitted
tcm_ml objects, or plain character vectors, and returns a named list
of character vectors ready for downstream plotting.

Usage:
```r
get_ml_gene_sets(..., set_names = NULL)
```

参数解析:
- `...`: One of the following:

A single tcm_ml_listfrom run_ml_screening().
Multiple tcm_ml objectse.g.
get_ml_gene_sets(res_lasso, res_rf, res_xgb).
A single named listof character vectors or tcm_ml objects.
Multiple character vectorsgene sets directly.
- `set_names`: Optional character vector of names for the sets.
Ignored when names can be inferred from the input.

使用示例:
```r
## From run_ml_screening()
gene_sets <- get_ml_gene_sets(ml_list)

## From individual models
gene_sets <- get_ml_gene_sets(res_lasso, res_rf, res_xgb)

## Venn diagram (2-4 sets)
venn_df <- do.call(getvenndata,
c(gene_sets, list(set_names = names(gene_sets))))
ggvenn_plot(venn_df)

## UpSet plot (usually more than 4 sets)
upsetplot(gene_sets)
```

### `ml_enet()`

Elastic Net feature selection via cv.glmnet
Penalised logistic regression with configurable alpha.

alpha = 1 : LASSO
alpha = 0 : Ridge
0 < alpha < 1 : Elastic Net

Caution: For Ridge (alpha = 0), all coefficients are non-zero.
When top_n = NULL (default), all features are retained so that
users can inspect $importance first, then call
select_features() to trim to the desired number.

Elastic Net feature selection via cv.glmnet
Penalised logistic regression with configurable alpha.

alpha = 1 : LASSO
alpha = 0 : Ridge
0 < alpha < 1 : Elastic Net

Caution: For Ridge (alpha = 0), all coefficients are non-zero.
When top_n = NULL (default), all features are retained so that
users can inspect $importance first, then call
select_features() to trim to the desired number.

Usage:
```r
ml_enet(
ml_data,
alpha = 0.5,
lambda_rule = c("1se", "min"),
top_n = NULL,
type_measure = "auc",
cv_folds = 10,
seed = 2025,
...
)
```

参数解析:
- `ml_data`: A tcm_ml_data object from prepare_ml_data().
- `alpha`: Mixing parameter, 0 to 1. Default 0.5.
- `lambda_rule`: "1se" (parsimonious) or "min" (best AUC).
- `top_n`: For Ridge: max number of top genes to keep (by
abs(coefficient)). Default NULL (keep all).
Ignored for LASSO / Elastic Net.
- `type_measure`: Loss function for CV in glmnet::cv.glmnet().
Default "auc". Other options: "deviance",
"class", "mse", "mae".
- `cv_folds`: Number of CV folds. Default 10.
- `seed`: Random seed. Default 2025.
- `...`: Passed to glmnet::cv.glmnet().

使用示例:
```r
ml_data <- prepare_ml_data(expr_mat, group)
res <- ml_enet(ml_data)
```

### `ml_lasso()`

LASSO feature selection (a wrapper function from ml_enet())
Calls ml_enet() with alpha = 1.

LASSO feature selection (a wrapper function from ml_enet())
Calls ml_enet() with alpha = 1.

Usage:
```r
ml_lasso(
ml_data,
lambda_rule = c("1se", "min"),
cv_folds = 10,
seed = 2025,
...
)
```

参数解析:
- `ml_data`: A tcm_ml_data object from prepare_ml_data().
- `lambda_rule`: "1se" (parsimonious) or "min" (best AUC).
- `cv_folds`: Number of CV folds. Default 10.
- `seed`: Random seed. Default 2025.
- `...`: Passed to glmnet::cv.glmnet().

使用示例:
```r
ml_data <- prepare_ml_data(expr_mat, group)
res <- ml_lasso(ml_data)
```

### `ml_rf()`

Random Forest + Boruta feature selection for key-gene screening

Supervised feature-selection pipeline designed for network-pharmacology
workflows: starting from PPI-derived candidate genes and RNA-seq expression
profiles, this function identifies a stable, interpretable set of key genes
rather than optimising prediction accuracy.

Pipeline overview

Train a Random Forest on all candidate genes to obtain
importance scores (MeanDecreaseGini / MeanDecreaseAccuracy).
Run Boruta all-relevant feature selection; tentative features are
resolved via TentativeRoughFix.
(Optional, default) Re-fit a Random Forest on the selected
genes only so that OOB and test-set metrics honestly reflect the
chosen subset (refit_on_selected = TRUE).

Usage:
```r
ml_rf(
ml_data,
n_trees = 500,
top_n = NULL,
boruta_fallback_n = 20,
max_runs = 100,
refit_on_selected = TRUE,
seed = 2025,
...
)
```

参数解析:
- `ml_data`: A tcm_ml_data object produced by prepare_ml_data. Must contain exactly two classes.
- `n_trees`: Number of trees for both the initial RF and the refit RF. Default 500.
- `top_n`: Hard cap: keep at most top_n Boruta-confirmed genes (ranked by MeanDecreaseGini). NULL = no cap (default).
- `boruta_fallback_n`: When Boruta confirms zero features, fall back to the top-boruta_fallback_n genes ranked by MeanDecreaseGini.
Default 20.
- `max_runs`: Maximum Boruta iterations passed to Boruta. Default 100.
- `refit_on_selected`: Logical. Re-fit a RF on selected genes only so that OOB / test metrics reflect the actual feature subset. Default TRUE (recommended).
- `seed`: Random seed. Default 2025.
- `...`: Additional arguments forwarded to randomForest.

使用示例:
```r
ml_data <- prepare_ml_data(expr_mat, group)
res <- ml_rf(ml_data)
```

### `ml_ridge()`

Ridge regression feature selection (a wrapper function from ml_enet())
Calls ml_enet() with alpha = 0.
Caution: Ridge does not shrink coefficients to zero; when top_n = NULL
(default), all features are kept. Use select_features() after inspecting importance to trim.

Ridge regression feature selection (a wrapper function from ml_enet())
Calls ml_enet() with alpha = 0.
Caution: Ridge does not shrink coefficients to zero; when top_n = NULL
(default), all features are kept. Use select_features() after inspecting importance to trim.

Usage:
```r
ml_ridge(
ml_data,
lambda_rule = c("1se", "min"),
top_n = NULL,
cv_folds = 10,
seed = 2025,
...
)
```

参数解析:
- `ml_data`: A tcm_ml_data object from prepare_ml_data().
- `lambda_rule`: "1se" (parsimonious) or "min" (best AUC).
- `top_n`: For Ridge: max number of top genes to keep (by
abs(coefficient)). Default NULL (keep all).
Ignored for LASSO / Elastic Net.
- `cv_folds`: Number of CV folds. Default 10.
- `seed`: Random seed. Default 2025.
- `...`: Passed to glmnet::cv.glmnet().

使用示例:
```r
ml_data <- prepare_ml_data(expr_mat, group)
res <- ml_ridge(ml_data)
```

### `ml_svm_rfe()`

SVM-RFE feature selection
Recursive Feature Elimination with SVM via caret::rfe().
Designed for binary classification of gene expression data (e.g. RNA-seq
profiles of PPI candidate genes vs. disease/control labels).

SVM-RFE feature selection
Recursive Feature Elimination with SVM via caret::rfe().
Designed for binary classification of gene expression data (e.g. RNA-seq
profiles of PPI candidate genes vs. disease/control labels).

Usage:
```r
ml_svm_rfe(
ml_data,
sizes = NULL,
top_n = NULL,
min_size = 5,
kernel = c("svmLinear", "svmRadial"),
metric = "ROC",
cv_folds = 5,
cv_repeats = 5,
seed = 2025,
...
)
```

参数解析:
- `ml_data`: A tcm_ml_data object. Must contain exactly two classes.
- `sizes`: Integer vector of feature-subset sizes to evaluate in RFE.
Default NULL: sizes are derived automatically from the data by
log-spacing 10 values from 2 to min(p - 1, max(p %/% 2, 50)),
where p is the number of input features. When provided manually,
values outside [2, p] are silently dropped.
- `top_n`: Override the auto-selected optimal size: keep the top
top_n genes by aggregated importance. Default NULL
(use the size that maximises metric in the RFE profile).
- `min_size`: Minimum number of features to retain. The profile search
is restricted to Variables >= min_size. Default 5.
Useful to avoid biologically trivial 2-3 gene subsets.
- `kernel`: SVM kernel to use. "svmLinear" (default, recommended
for gene expression data — weights |w_j|^2 are directly used for
recursive elimination per Guyon 2002) or "svmRadial".
- `metric`: Optimization metric for RFE. Default "ROC".
Other options include "Accuracy", "Kappa".
When "ROC" / "Sens" / "Spec" is used,
twoClassSummary is applied automatically; otherwise
defaultSummary is used.
- `cv_folds`: Outer CV folds. Default 5.
- `cv_repeats`: Outer CV repeats. Default 5.
- `seed`: Random seed. Default 2025.
- `...`: Passed to caret::rfe().

使用示例:
```r
ml_data <- prepare_ml_data(expr_mat, group)
res <- ml_svm_rfe(ml_data)
```

### `ml_xgboost()`

XGBoost feature selection
Trains an XGBoost classifier with nrounds-round CV, then ranks features by information gain.
Since XGBoost does not inherently perform variable selection, top_n controls how many top features to keep.
When top_n = NULL (default), all features with Gain > 0 are retained; use select_features() afterwards to trim.

XGBoost feature selection
Trains an XGBoost classifier with nrounds-round CV, then ranks features by information gain.
Since XGBoost does not inherently perform variable selection, top_n controls how many top features to keep.
When top_n = NULL (default), all features with Gain > 0 are retained; use select_features() afterwards to trim.

Usage:
```r
ml_xgboost(
ml_data,
top_n = NULL,
eval_metric = "auc",
nrounds = 200,
max_depth = 6,
eta = 0.1,
cv_folds = 5,
early_stopping_rounds = 20,
seed = 2025,
...
)
```

参数解析:
- `ml_data`: A tcm_ml_data object.
- `top_n`: Number of top features to select by Gain. Default NULL (keep all).
- `eval_metric`: Evaluation metric for XGBoost CV. Default "auc".
Other options include "error", "logloss", "aucpr".
See xgboost documentation for full list.
- `nrounds`: Maximum boosting rounds. Default 200.
- `max_depth`: Maximum tree depth. Default 6.
- `eta`: Learning rate. Default 0.1.
- `cv_folds`: CV folds for early stopping. Default 5.
- `early_stopping_rounds`: Stop if no improvement for this many rounds. Default 20.
- `seed`: Random seed. Default 2025.
- `...`: Passed to xgboost::xgb.cv() and xgboost::xgb.train().

使用示例:
```r
ml_data <- prepare_ml_data(expr_mat, group)
res <- ml_xgboost(ml_data)
```

### `plot_enet_coefs()`

Bar chart of LASSO / Elastic Net / Ridge coefficients

Horizontal bar chart coloured by sign (positive / negative).

Usage:
```r
plot_enet_coefs(result, top_n = 20)
```

参数解析:
- `result`: A tcm_ml object from ml_enet(), ml_lasso() or
ml_ridge().
- `top_n`: Number of top features to show. Default 20.

使用示例:
```r
# Fill required arguments according to the parameter list above
plot_enet_coefs(result, top_n = 20)
```

### `plot_enet_cv()`

Plot Elastic Net / LASSO / Ridge CV curve

Wrapper around the default plot.cv.glmnet() method from the
glmnet package. Displays the CV metric (AUC / deviance) as a
function of \log(\lambda), with ± 1 SE error bars and vertical
dashed lines at lambda.min and lambda.1se.

Usage:
```r
plot_enet_cv(result, ...)
```

参数解析:
- `result`: A tcm_ml object from ml_enet(), ml_lasso() or ml_ridge().
- `...`: Additional arguments passed to plot.cv.glmnet().

使用示例:
```r
# Fill required arguments according to the parameter list above
plot_enet_cv(result, ...)
```

### `plot_enet_path()`

Plot Elastic Net / LASSO / Ridge coefficient path

Wrapper around the default plot.glmnet() method.
Each curve traces one feature's coefficient across the \lambda grid.
When top_n is set, only the features with the largest absolute
coefficients at the selected \lambda are shown.

Usage:
```r
plot_enet_path(result, top_n = 20, xvar = "lambda", ...)
```

参数解析:
- `result`: A tcm_ml object from ml_enet(), ml_lasso() or ml_ridge().
- `top_n`: Integer or NULL. If non-NULL and fewer features
than the total, only the top_n features with largest absolute
coefficient at the selected \lambda are drawn. Default 20.
- `xvar`: Variable on x-axis: "lambda" (default), "norm", or "dev".
- `...`: Additional arguments passed to plot.glmnet().

使用示例:
```r
# Fill required arguments according to the parameter list above
plot_enet_path(result, top_n = 20, xvar = "lambda", ...)
```

### `plot_gene_boxplot()`

Two-group expression boxplots for diagnostic genes

Draws a faceted boxplot comparing expression levels between two groups for
each gene, with optional jittered points, violin overlay and statistical
test annotation.  Designed for publication-quality figures with zero extra
dependencies beyond ggplot2 and base R stats.

Usage:
```r
plot_gene_boxplot(
genes,
ml_data = NULL,
expr_mat = NULL,
group = NULL,
use = c("auto", "train", "test"),
test_method = c("wilcox", "t.test"),
palette = c("#E74C3C", "#3498DB"),
show_points = TRUE,
point_alpha = 0.35,
violin = FALSE,
ncol = NULL,
scales = "free_y",
p_label = c("p.format", "p.signif", "p.adj"),
p_adjust_method = "BH",
base_size = 12
)
```

参数解析:
- `genes`: Character vector of gene symbols to evaluate.
- `ml_data`: A tcm_ml_data object from prepare_ml_data
(preferred input).  When provided, expr_mat / group are
ignored.
- `expr_mat`: Numeric matrix (genes \times samples) with
rownames = gene symbols.  Used only when ml_data is
NULL.
- `group`: Factor or character vector of sample labels (length =
ncol(expr_mat), two levels).
- `use`: Which data split to evaluate when using ml_data:
"auto" (default; test set if available, else train),
"train", or "test".
- `test_method`: Statistical test: "wilcox" (Wilcoxon rank-sum,
default) or "t.test".
- `palette`: Character(2); fill colours for the two groups.
Default c("#E74C3C", "#3498DB") (red / blue).
- `show_points`: Logical; overlay jittered data points? Default TRUE.
- `point_alpha`: Alpha transparency for jittered points. Default 0.35.
- `violin`: Logical; add a violin layer behind the boxes?
Default FALSE.
- `ncol`: Integer; number of columns in facet_wrap.
Default NULL (auto-computed).
- `scales`: Passed to facet_wrap.
Default "free_y".
- `p_label`: How to annotate p-values: "p.format" (numeric,
default), "p.signif" (stars: ns / * / ** / *** / ****), or
"p.adj" (adjusted, numeric).
- `p_adjust_method`: Method for multiple-testing correction passed to
p.adjust.  Default "BH".  Applies both to the
p_adj column in attr(, "test_table") and to the annotation
when p_label = "p.adj".  When only a single gene is tested,
p_adj equals p_value regardless of method.
- `base_size`: Base font size for theme_bw.
Default 12.

使用示例:
```r
consensus <- get_ml_consensus(ml_list, min_methods = 2)
plot_gene_boxplot(consensus, ml_data = ml_data)
plot_gene_boxplot(consensus, expr_mat = expr_mat, group = group,
violin = TRUE, p_label = "p.signif")
```

### `plot_gene_roc()`

ROC curves for individual diagnostic genes

Overlays one ROC curve per gene on a single panel (combine = TRUE,
default), or draws one ROC per facet (combine = FALSE).
This is the natural follow-up after get_ml_consensus: the
consensus genes are fed back to evaluate their individual diagnostic power.

Usage:
```r
plot_gene_roc(
genes,
ml_data = NULL,
expr_mat = NULL,
group = NULL,
use = c("auto", "train", "test"),
combine = TRUE,
palette = NULL,
show_ci = FALSE,
ncol = NULL,
label_size = 9,
p_adjust_method = "BH"
)
```

参数解析:
- `genes`: Character vector of gene symbols to evaluate.
- `ml_data`: A tcm_ml_data object from prepare_ml_data
(preferred input).  When provided, expr_mat / group are
ignored.
- `expr_mat`: Numeric matrix (genes \times samples) with
rownames = gene symbols.  Used only when ml_data is
NULL.
- `group`: Factor or character vector of sample labels (length =
ncol(expr_mat), two levels).
- `use`: Which data split to evaluate when using ml_data:
"auto" (default; test set if available, else train),
"train", or "test".
- `combine`: Logical; if TRUE (default) all curves are overlaid on one panel.  If FALSE each gene gets its own facet.
- `palette`: Character vector of colours (one per gene), or NULL for the built-in 12-colour publication palette.
- `show_ci`: Logical; Whether to add a translucent 95 percent CI band around each curve. Default FALSE.  Requires pROC::ci.se().
- `ncol`: Integer; number of columns when combine = FALSE. Default NULL (auto).
- `label_size`: Legend text size. Default 9.
- `p_adjust_method`: Method for multiple-testing correction passed to
p.adjust.  Default "BH" (Benjamini--Hochberg
FDR).  Other common choices: "bonferroni", "holm",
"none".  When only a single gene is supplied no correction is
needed and p_adj equals p_value regardless of method.

使用示例:
```r
consensus <- get_ml_consensus(ml_list, min_methods = 2)
plot_gene_roc(consensus, ml_data = ml_data)
plot_gene_roc(consensus, ml_data = ml_data, combine = FALSE)
```

### `plot_ml_roc()`

ROC curves for ML methods

Supports all modes.

Mode B / C: ROC is computed on the held-out test set.
Mode A: uses out-of-fold (LASSO / Enet / XGBoost / SVM-RFE) or OOB (RF)
predictions.
SVM-RFE falls back to resubstitution when OOF is unavailable.

Usage:
```r
plot_ml_roc(ml_list)
```

参数解析:
- `ml_list`: A tcm_ml_list from run_ml_screening().

使用示例:
```r
# Fill required arguments according to the parameter list above
plot_ml_roc(ml_list)
```

### `plot_ml_venn()`

Venn diagram of selected genes

Convenience wrapper: calls get_ml_gene_sets() to extract gene sets,
then draws a Venn diagram via getvenndata() + ggvenn_plot().
Accepts the same flexible inputs as get_ml_gene_sets().
For > 4 sets consider using upsetplot(get_ml_gene_sets(...)).

Usage:
```r
plot_ml_venn(..., set_names = NULL)
```

参数解析:
- `...`: One of the following:

A single tcm_ml_listfrom run_ml_screening().
Multiple tcm_ml objectse.g.
get_ml_gene_sets(res_lasso, res_rf, res_xgb).
A single named listof character vectors or tcm_ml objects.
Multiple character vectorsgene sets directly.
- `set_names`: Optional character vector of names for the sets.
Ignored when names can be inferred from the input.

使用示例:
```r
# Fill required arguments according to the parameter list above
plot_ml_venn(..., set_names = NULL)
```

### `plot_rf_boruta()`

Plot Boruta feature selection result

Plot Boruta feature selection result

Usage:
```r
plot_rf_boruta(result, top_n = 20)
```

参数解析:
- `result`: A tcm_ml object from ml_rf().
- `top_n`: Integer; show at most this many features.

使用示例:
```r
# Fill required arguments according to the parameter list above
plot_rf_boruta(result, top_n = 20)
```

### `plot_rf_importance()`

Plot RF variable importance

Plot RF variable importance

Usage:
```r
plot_rf_importance(result, top_n = 20)
```

参数解析:
- `result`: A tcm_ml object from ml_rf().
- `top_n`: Integer; number of top features to show.

使用示例:
```r
# Fill required arguments according to the parameter list above
plot_rf_importance(result, top_n = 20)
```

### `plot_svm_rfe_curve()`

Plot SVM-RFE accuracy profile

Plot SVM-RFE accuracy profile

Usage:
```r
plot_svm_rfe_curve(result)
```

参数解析:
- `result`: A tcm_ml object from ml_svm_rfe().

使用示例:
```r
# Fill required arguments according to the parameter list above
plot_svm_rfe_curve(result)
```

### `plot_xgb_importance()`

Plot XGBoost feature importance (Gain)

Plot XGBoost feature importance (Gain)

Usage:
```r
plot_xgb_importance(result, top_n = 20)
```

参数解析:
- `result`: A tcm_ml object from ml_xgboost().
- `top_n`: Integer; number of top features to show.

使用示例:
```r
# Fill required arguments according to the parameter list above
plot_xgb_importance(result, top_n = 20)
```

### `prepare_ml_data()`

Prepare expression data for ML feature selection

Supports three modes:

Mode A (default): full-data CV, no hold-out (split = FALSE).
Mode B: internal train/test split (split = TRUE).
Mode C: external validation (test_expr + test_group).

Usage:
```r
prepare_ml_data(
expr_mat,
group,
positive_class = NULL,
genes = NULL,
split = FALSE,
train_ratio = 0.7,
train_idx = NULL,
test_expr = NULL,
test_group = NULL,
seed = 2025
)
```

参数解析:
- `expr_mat`: Numeric matrix (genes x samples). Row names = gene symbols.
- `group`: Factor/character of length ncol(expr_mat), two levels.
- `positive_class`: Which level to treat as the positive class (levels[[1]]).
Default NULL: if group is already an ordered factor, its level order is
preserved; otherwise levels are sorted alphabetically.
- `genes`: Optional character vector of candidate genes to keep. Default is NULL.
- `split`: Logical. TRUE = Mode B. Default is FALSE.
- `train_ratio`: Fraction for training when split = TRUE. Default is 0.7.
- `train_idx`: Optional integer indices overriding train_ratio.
- `test_expr`: External validation matrix (genes x samples) for Mode C.
- `test_group`: Labels for test_expr (required if test_expr given).
- `seed`: Random seed. Default is 2025.

使用示例:
```r
## Mode A: full CV
ml_data <- prepare_ml_data(expr_mat, group, positive_class = "Disease")
## Mode B: internal split
ml_data <- prepare_ml_data(expr_mat, group, split = TRUE, train_ratio = 0.7)
```

### `run_ml_screening()`

Run multiple ML methods for feature selection
Convenience wrapper that runs LASSO / Elastic Net / Ridge / RF / SVM-RFE / XGBoost.
By default top_n = NULL: Ridge and XGBoost retain all features, SVM-RFE auto-selects the optimal size.
After inspecting results, call select_features() on individual models to trim to the desired count.

Run multiple ML methods for feature selection
Convenience wrapper that runs LASSO / Elastic Net / Ridge / RF / SVM-RFE / XGBoost.
By default top_n = NULL: Ridge and XGBoost retain all features, SVM-RFE auto-selects the optimal size.
After inspecting results, call select_features() on individual models to trim to the desired count.

Usage:
```r
run_ml_screening(
ml_data,
methods = c("lasso", "rf", "svm_rfe", "xgboost"),
seed = 2025,
cv_folds = 5,
cv_repeats = 5,
alpha = 0.5,
top_n = NULL,
...
)
```

参数解析:
- `ml_data`: A tcm_ml_data object.
- `methods`: Subset of c("lasso", "enet", "ridge", "rf", "svm_rfe", "xgboost").
- `seed`: Random seed. Default 2025.
- `cv_folds`: CV folds. Default 5.
- `cv_repeats`: SVM-RFE outer CV repeats. Default 5.
- `alpha`: Elastic Net alpha. Default 0.5 (ignored for lasso/ridge).
- `top_n`: Top-n features for Ridge / XGBoost / SVM-RFE. Default NULL (keep all or auto-select).
- `...`: Forwarded to individual model functions.

使用示例:
```r
ml_data <- prepare_ml_data(expr_mat, group)
res_list <- run_ml_screening(ml_data, methods = c("lasso", "rf"))
summary(res_list)
```

## Search and evidence retrieval

### `get_pubmed_data()`

PubMed Literature Scraper for TCMDATA

Fetches and summarizes PubMed literature based on specific Traditional Chinese Medicine (TCM)
and disease keywords within a defined time frame.

Usage:
```r
get_pubmed_data(
tcm_name,
disease_name,
year_range = NULL,
email = NULL,
retmax = 1000,
...
)
```

参数解析:
- `tcm_name`: Character. Keywords for the TCM (e.g., "Ginseng").
- `disease_name`: Character. Keywords for the disease (e.g., "Fatigue").
- `year_range`: Integer vector of length 2. The publication year range (e.g., c(2015, 2026)).
- `email`: Character. User email address required by NCBI E-utils.
- `retmax`: Integer. Maximum number of PMIDs to retrieve. Default is 1000.
- `...`: Additional arguments passed to other methods.

使用示例:
```r
result <- get_pubmed_data(tcm_name = "Artemisinin", 
disease_name = "Malaria", 
year_range = c(2020, 2025),
retmax = 100,
email = "your email")

print(str(result))
```

### `get_pubmed_table()`

Extract and Sort PubMed Literature Table

Accesses the literature data within a 'tcm_pubmed' object and returns it as
a data frame sorted by publication year in descending order.

Usage:
```r
get_pubmed_table(x, n = NULL, file = NULL)
```

参数解析:
- `x`: An object of class 'tcm_pubmed' generated by get_pubmed_data.
- `n`: Integer. Number of top records to return. If NULL (default), returns all records.
- `file`: Character. Optional file path (e.g., "results.csv"). If provided, the table will be saved to this location.

使用示例:
```r
## get result from pubmed
result <- get_pubmed_data(tcm_name = "Artemisinin", 
disease_name = "Malaria", 
year_range = c(2020, 2025),
retmax = 100,
email = "your email")

# Get the sorted table
lit_table <- get_pubmed_table(result)

# Get top 10 records and save to CSV
lit_table <- get_pubmed_table(result, n = 10, file = "top10_lit.csv")
```

### `search_disease()`

Search disease targets (disease -> genes)

Query DisGeNET (via DOSE) to find genes associated with a disease.
Supports UMLS CUI IDs, exact name match, or fuzzy name search.

Usage:
```r
search_disease(disease, readable = TRUE)
```

参数解析:
- `disease`: Character. Disease name or UMLS CUI (e.g. "sepsis" or
"C0243026"). Supports a vector for multiple diseases.
- `readable`: Logical. Convert Entrez IDs to gene symbols (default TRUE).

使用示例:
```r
search_disease("sepsis")
search_disease(c("sepsis", "asthma"))
search_disease("C0243026")
```

### `search_gene_disease()`

Search gene-associated diseases (gene -> diseases)

Reverse lookup: given gene symbols or Entrez IDs, find associated diseases
from DisGeNET (via DOSE).

Usage:
```r
search_gene_disease(gene, readable = TRUE)
```

参数解析:
- `gene`: Character. Gene symbols (e.g. "TNF") or Entrez IDs.
Supports a vector for multiple genes.
- `readable`: Logical. Attach gene symbol column (default TRUE).

使用示例:
```r
search_gene_disease("TNF")
search_gene_disease(c("IL6", "TNF", "PPARG"))
search_gene_disease("7124")
```

### `search_geo_datasets()`

Search GEO for datasets related to a disease or keyword

Queries the NCBI GEO database (via E-utilities) to find relevant
expression datasets (GDS/GSE) for a given search term. Useful for
discovering RNA-seq, microarray, or single-cell datasets before analysis.

Usage:
```r
search_geo_datasets(
query,
organism = "Homo sapiens",
dataset_type = NULL,
email = NULL,
retmax = 20L
)
```

参数解析:
- `query`: Character. Search query (e.g. disease name, tissue + disease).
- `organism`: Character. Organism filter. Default "Homo sapiens".
- `dataset_type`: Character. Optional type filter: "expression profiling by array",
"expression profiling by high throughput sequencing", or "single-cell".
Pass NULL to search all types.
- `email`: Character. Email address required by NCBI E-utilities.
- `retmax`: Integer. Maximum number of results. Default 20.

使用示例:
```r
res <- search_geo_datasets("sepsis", email = "user@example.com")
head(res)
```

### `search_herb()`

Retrieve Herbs and Retrieve Associated Molecules and Targets

This function retrieves herb, molecule, and target information from the internal dataset, tcm_data,
based on specified herb names and their corresponding name types (Chinese, Pinyin, or English).

Usage:
```r
search_herb(herb, type)
```

参数解析:
- `herb`: A character vector containing the names of herbs to be queried.
- `type`: A string indicating the type of herb names provided. Must be one of "Herb_cn_name",
"Herb_pinyin_name", or "Herb_en_name".

使用示例:
```r
# Example usage with Chinese herb names
herbs <- c("灵芝")
lz <- search_herb(herb = herbs, type = "Herb_cn_name")
print(lz)

# Example usage with pinyin herb names
herbs <- c("Ginseng", "Licorice")
comp <- search_herb(herb = herbs, type = "Herb_pinyin_name")
print(comp)
```

### `search_target()`

Retrieve Herbs, Molecules, and Targets Based on Gene List

This function retrieves herb, molecule, and target information from the internal dataset, tcm_data,
based on a provided list of target genes.

Usage:
```r
search_target(gene_list)
```

参数解析:
- `gene_list`: A character vector containing gene symbols which can be determined freely by users.

使用示例:
```r
# Example usage with a list of gene symbols
genes <- c("TP53", "EGFR", "BRCA1")
herbs_targets <- search_target(genes)
print(herbs_targets)
```

## General visualization

### `getvenndata()`

Construct a Logical Matrix for Venn Diagram Visualization

This function takes 2 to 4 named character vectors as input, removes duplicates,
and returns a data frame in logical format that is compatible with ggvenn for Venn diagram plotting.

Usage:
```r
getvenndata(..., set_names = NULL)
```

参数解析:
- `...`: 2 to 4 character vectors representing different gene sets or element sets.
- `set_names`: Optional character vector to assign custom names to the sets.
If not provided, names will be automatically assigned as "Set1", "Set2", etc.

使用示例:
```r
# Usage with 3 sets
targets_db1 <- c("TP53", "EGFR", "KRAS", "MYC", "AKT1")
targets_db2 <- c("TP53", "KRAS", "PTEN", "BRCA1")
targets_db3 <- c("EGFR", "MYC", "BRAF", "PTEN", "TP53")

df_venn3 <- getvenndata(targets_db1, targets_db2, targets_db3,
set_names = c("Database_1", "Database_2", "Database_3"))

head(df_venn3)
```

### `getvennresult()`

Get all Venn set intersections from getvenndata-style logical matrix

Get all Venn set intersections from getvenndata-style logical matrix

Usage:
```r
getvennresult(venn_df, col_names = NULL, drop_empty = TRUE)
```

参数解析:
- `venn_df`: A data.frame from getvenndata(), with logical columns
- `col_names`: Character vector of columns to use (default: all except 1st)
- `drop_empty`: Logical; whether to drop combinations with 0 genes

使用示例:
```r
# Usage with 3 sets
targets_db1 <- c("TP53", "EGFR", "KRAS", "MYC", "AKT1")
targets_db2 <- c("TP53", "KRAS", "PTEN", "BRCA1")
targets_db3 <- c("EGFR", "MYC", "BRAF", "PTEN", "TP53")

df_venn3 <- getvenndata(targets_db1, targets_db2, targets_db3,
set_names = c("Database_1", "Database_2", "Database_3"))

venn_res <- getvennresult(df_venn3)
print(venn_res)
```

### `ggvenn_plot()`

Plot a Venn Diagram for 2–4 Sets Using ggvenn

This is a wrapper around ggvenn() to draw customizable Venn diagrams
for 2–4 sets with automatic color matching and error handling.

Usage:
```r
ggvenn_plot(
venn_df,
col_names = NULL,
set.color = c("#E41A1C", "#1E90FF", "#FF8C00", "#4DAF4A", "#75cbdc"),
set.name.color = "black",
use.color.as.text = TRUE,
name.size = 5,
text.size = 4,
stroke.color = "black",
stroke.size = 0.6,
show.percentage = FALSE,
show.elements = FALSE,
digits = 1,
expand_ratio = 0.2
)
```

参数解析:
- `venn_df`: A data.frame where the first column is element ID and the remaining 2–4 columns are logical (TRUE/FALSE) indicators for set membership.
- `col_names`: Character vector of column names to be used as sets. If NULL, all columns except the first are used.
- `set.color`: Fill colors for Venn sets.
- `set.name.color`: Color for set name labels (ignored if use.color.as.text = TRUE).
- `use.color.as.text`: Logical; if TRUE, use fill color as set name text color.
- `name.size`: Font size for set names.
- `text.size`: Font size for intersection counts.
- `stroke.color`: Circle border color. Default is 'black'.
- `stroke.size`: Circle border thickness. Default is 0.6.
- `show.percentage`: Logical; whether to show percentages instead of raw counts.
- `show.elements`: Logical; whether to display individual elements in each region.
- `digits`: Number of decimal places for percentage display.
- `expand_ratio`: Plot margin expansion ratio for aesthetics.

使用示例:
```r
# Fill required arguments according to the parameter list above
ggvenn_plot(
venn_df,
col_names = NULL,
set.color = c("#E41A1C", "#1E90FF", "#FF8C00", "#4DAF4A", "#75cbdc"),
set.name.color = "black",
use.color.as.text = TRUE,
name.size = 5,
text.size = 4,
stroke.color = "black",
stroke.size = 0.6,
show.percentage = FALSE,
show.elements = FALSE,
digits = 1,
expand_ratio = 0.2
)
```

### `TCM_sankey()`

TCM Sankey Plot for herb-molecule-target

TCM Sankey Plot for herb-molecule-target

Usage:
```r
tcm_sankey(
data,
axis_order = c("herb", "molecule", "target"),
herb_cols = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
"#e377c2", "#7f7f7f"),
mol_cols = c("#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
"#c5b0d5", "#E41A1C", "#377EB8"),
target_cols = c("#c49c94", "#f7b6d2", "#dbdb8d", "#c7e9c0", "#f4cae4", "#e6f598",
"#ffeda0", "#BB5234", "#BB7813", "#FF6158"),
plot_font = "sans",
font_face = "plain",
target_fontface = "italic",
font_size = 3.6,
width = 0.05,
alpha = 0.3,
knot.pos = 0.3
)

TCM_sankey(...)
```

参数解析:
- `data`: A data frame containing at least three columns: herb, molecule, and target.
- `axis_order`: Character vector specifying the order of axes in the Sankey diagram. Default is c("herb", "molecule", "target").
- `herb_cols`: Character vector defining the base color palette for the herb layer.
- `mol_cols`: Character vector defining the base color palette for the molecule.
- `target_cols`: Character vector defining the base color palette for the target.
- `plot_font`: Character string specifying the font family used for text. Default is "sans".
- `font_face`: Character string specifying the font face. Default is "plain".
- `target_fontface`: Character string specifying the font face for target labels (rightmost axis). Default is "italic".
- `font_size`: Numeric value controlling the size of node labels. Default is 3.5.
- `width`: Numeric value controlling the width of both nodes and flows. Default is 0.05.
- `alpha`: Numeric value controlling the transparency of the flows. Default is 0.3.
- `knot.pos`: Numeric value (between 0 and 1) determining the curvature position of flow lines. Default is 0.3.
- `...`: Additional arguments passed to tcm_sankey when using the deprecated alias TCM_sankey.

使用示例:
```r
# Fill required arguments according to the parameter list above
tcm_sankey(
data,
axis_order = c("herb", "molecule", "target"),
herb_cols = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
"#e377c2", "#7f7f7f"),
mol_cols = c("#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
"#c5b0d5", "#E41A1C", "#377EB8"),
target_cols = c("#c49c94", "#f7b6d2", "#dbdb8d", "#c7e9c0", "#f4cae4", "#e6f598",
"#ffeda0", "#BB5234", "#BB7813", "#FF6158"),
plot_font = "sans",
font_face = "plain",
target_fontface = "italic",
font_size = 3.6,
width = 0.05,
alpha = 0.3,
knot.pos = 0.3
)

TCM_sankey(...)
```

### `upsetplot()`

UpSet plot for gene sets

Convenience wrapper around aplotExtra::upset_plot() so the plot is
directly available from TCMDATA. This is particularly useful for the output
of get_ml_gene_sets().

Usage:
```r
upsetplot(
list,
nintersects = NULL,
order.intersect.by = c("size", "name"),
order.set.by = c("size", "name"),
color.intersect.by = "none",
color.set.by = "none",
remove_empty_intersects = TRUE
)
```

参数解析:
- `list`: A named list of gene sets.
- `nintersects`: Number of intersections to show. NULL shows all.
- `order.intersect.by`: One of "size" or "name".
- `order.set.by`: One of "size" or "name".
- `color.intersect.by`: Color scheme for intersection bars.
- `color.set.by`: Color scheme for set bars.
- `remove_empty_intersects`: Whether to remove empty intersections.

使用示例:
```r
gene_sets <- list(
LASSO = c("TP53", "BRCA1", "EGFR"),
RF = c("TP53", "AKT1"),
XGBOOST = c("EGFR", "AKT1", "MTOR")
)
upsetplot(gene_sets)
```

## Enrichment and pathway visualization

### `ggdot_sankey()`

Dot–Sankey plot for enrichment results

Dot–Sankey plot for enrichment results

Usage:
```r
ggdot_sankey(
enrich_obj,
n = 10,
axis_order = c("Pathway", "Gene"),
id_colors = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
"#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"),
desc_colors = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
"#91D1C2", "#C8DC96", "#7E6148", "#B09C85", "#6A5ACD", "#A0522D"),
id_y_pos = 1.1,
desc_y_pos = 1,
pathway_wrap = 50,
sankey_text_size = 4,
bubble_size_range = c(3, 8),
dot_palette = "RdBu",
dot_x_var = c("GeneRatio", "RichFactor", "FoldEnrichment"),
bubble_p_label = "p.adjust",
sankey_width = 2,
dot_width = 1,
font_family = "sans",
font_face = "plain",
gene_fontface = "italic",
sankey_lab = "Gene-Pathway",
seed = 2025,
...
)
```

参数解析:
- `enrich_obj`: enrichResult object or data.frame containing enrichment results.
- `n`: Integer. The number of pathways to visualize (top n ranked by significance).
- `axis_order`: Character vector of length 2, specifying the order of axes in
the Sankey diagram. Default is c("Pathway", "Gene").
- `id_colors`: Character vector of base colors for gene nodes.
- `desc_colors`: Character vector of base colors for pathway nodes.
- `id_y_pos`: Numeric constant controlling the y-position of “Gene” strata in the Sankey diagram. Default is 1.1.
- `desc_y_pos`: Numeric constant controlling the y-position of “Pathway” strata in the Sankey diagram. Default is 1.0.
- `pathway_wrap`: Integer. The maximum line width for pathway label wrapping. Default is 50.
- `sankey_text_size`: Numeric. Font size for text labels in the Sankey diagram. Default is 4.
- `bubble_size_range`: Numeric vector of length 2. Range of point sizes in the dot plot. Default is c(3, 8).
- `dot_palette`: Character. Name of the RColorBrewer palette used for color gradients in the dot plot. Default is "RdBu".
- `dot_x_var`: Character. Variable used for the x-axis in the dot plot.
- `bubble_p_label`: Character. The column name for significance values used to color dotsd. Default is p.adjust.
- `sankey_width`: Numeric. Relative width of the Sankey panel in the combined plot. Default is 2.
- `dot_width`: Numeric. Relative width of the dot plot panel in the combined plot. Default is 1.
- `font_family`: Character. Font family for all text elements. Default is "Arial".
- `font_face`: Character. Font face for all text. Default is "plain".
- `gene_fontface`: Character. Font face for gene labels (rightmost axis). Default is "italic".
- `sankey_lab`: Character. Label for the x-axis of the Sankey diagram. Default is "Gene-Pathway".
- `seed`: Integer. Random seed for reproducibility of layout. Default is 2025.
- `...`: Additional arguments passed to internal helper functions.

使用示例:
```r
# Fill required arguments according to the parameter list above
ggdot_sankey(
enrich_obj,
n = 10,
axis_order = c("Pathway", "Gene"),
id_colors = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
"#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"),
desc_colors = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
"#91D1C2", "#C8DC96", "#7E6148", "#B09C85", "#6A5ACD", "#A0522D"),
id_y_pos = 1.1,
desc_y_pos = 1,
pathway_wrap = 50,
sankey_text_size = 4,
bubble_size_range = c(3, 8),
dot_palette = "RdBu",
dot_x_var = c("GeneRatio", "RichFactor", "FoldEnrichment"),
bubble_p_label = "p.adjust",
sankey_width = 2,
dot_width = 1,
font_family = "sans",
font_face = "plain",
gene_fontface = "italic",
sankey_lab = "Gene-Pathway",
seed = 2025,
...
)
```

### `gglollipop()`

Lollipop Plot for Enrichment Results

Lollipop Plot for Enrichment Results

Usage:
```r
gglollipop(
enrich_obj,
x = "RichFactor",
top_n = 10,
orderBy = NULL,
text.col = "black",
text.size = 8,
text.width = 35,
palette = "RdBu",
line.col = "grey60",
line.type = "solid",
line.size = 0.9,
plot_title = NULL,
show_count = TRUE,
...
)
```

参数解析:
- `enrich_obj`: An enrichment result object from clusterProfiler.
- `x`: Character. The variable used for the x-axis. Default is "RichFactor".
- `top_n`: Integer. Number of top enriched terms to display. Default is 10.
- `orderBy`: Character. Variable used to order the y-axis terms. Default is "x".
- `text.col`: Character. Colors for text. Default is black.
- `text.size`: Numeric. Font size for axis text and title. Default is 8.
- `text.width`: Numeric. Font width for axis text and title. Default is 35.
- `palette`: Character. Color palette name from RColorBrewer to use for dot color. Default is "RdBu".
- `line.col`: Character. Color of the segment lines. Default is "grey60".
- `line.type`: Character. Line type for segments. Default is "solid".
- `line.size`: Numeric. Line width for segments. Default is 0.9.
- `plot_title`: Character. Optional plot title. Default is NULL.
- `show_count`: Logical. Whether to display the count value as a text label next to each dot. Default is TRUE
- `...`: Additional arguments passed to internal helper functions.

使用示例:
```r
# Fill required arguments according to the parameter list above
gglollipop(
enrich_obj,
x = "RichFactor",
top_n = 10,
orderBy = NULL,
text.col = "black",
text.size = 8,
text.width = 35,
palette = "RdBu",
line.col = "grey60",
line.type = "solid",
line.size = 0.9,
plot_title = NULL,
show_count = TRUE,
...
)
```

### `go_barplot()`

Bar Plot for GO Enrichment Results

Create a bar plot for GO enrichment analysis results, with bars colored
by ontology category (BP, CC, MF). Based on barplot method from enrichplot.

Usage:
```r
go_barplot(
enrich_obj,
x = "Count",
order = TRUE,
top_n = 10,
colors = NULL,
label.size = 3.2,
show_count = TRUE,
label.bold = FALSE,
x.angle = 60,
legend.position = "top",
plot_title = NULL,
...
)
```

参数解析:
- `enrich_obj`: An enrichResult object from clusterProfiler::enrichGO().
- `x`: Character. The variable to be plotted on the x-axis. Default is "Count".
- `order`: Logical. Whether to sort x. Default is TRUE.
- `top_n`: Integer. Number of top terms to display per ontology category. Default is 10.
- `colors`: Named character vector for BP, CC, MF colors. Default is c(BP = "#E64B35", CC = "#4DBBD5", MF = "#00A087").
- `label.size`: Numeric. Font size for count labels on bars. Default is 3.2.
- `show_count`: Logical. Whether to show count labels on bars. Default is TRUE.
- `label.bold`: Logical. Whether count labels should be bold. Default is FALSE.
- `x.angle`: Numeric. Rotation angle of x-axis text (pathway descriptions). Default is 60.
- `legend.position`: Character. Position of the legend. Default is "top".
- `plot_title`: Character or NULL. Title of the plot. Default is NULL.
- `...`: Additional parameters if neccessary.

使用示例:
```r
## demo data
herbs <- c("灵芝")
lz <- search_herb(herb = herbs, type = "Herb_cn_name")
set.seed(2025)
g <- sample(lz$target, 200)

## GO enrichment analysis
library(clusterProfiler)
x <- enrichGO(g, ont="all", OrgDb='org.Hs.eg.db', keyType="SYMBOL")
p1 <- go_barplot(x)
print(p1)
```

### `gocircle_plot()`

GO circle plot

Visualize GO enrichment results as a circular plot highlighting p.adjust,
up/down-regulated genes, and RichFactor.

Usage:
```r
gocircle_plot(
x,
top,
up_genes = NULL,
down_genes = NULL,
cat_col = c(BP = "#E64B35", CC = "#4DBBD5", MF = "#00A087"),
up_col = "#AE2A8A",
down_col = "#7B87BF",
padjust_col = c("#FF906F", "#861D30"),
bg.col = "gray95",
max_width = 2.5,
fontsize = 8,
fontface = "bold",
fontfamily = "sans",
...
)
```

参数解析:
- `x`: enrichResult object or data.frame from clusterProfiler.
- `top`: Number of top terms to show per category.
- `up_genes`: Vector of up-regulated genes (optional).
- `down_genes`: Vector of down-regulated genes (optional).
- `cat_col`: Named colors for GO categories (BP/CC/MF).
- `up_col`: Colors for up-regulated gene bars.
- `down_col`: Colors for down-regulated gene bars.
- `padjust_col`: Gradient colors for -log10(p.adjust).
- `bg.col`: Background color of outer ring.
- `max_width`: Maximum width of up/down bars.
- `fontsize`: Numeric. Font size for labels and legends.
- `fontface`: Character. Font style, e.g. "plain", "bold", "italic".
- `fontfamily`: Character. Font family, e.g. "sans", "serif", "mono".
- `...`: Additional parameters for flexibility.

使用示例:
```r
# Fill required arguments according to the parameter list above
gocircle_plot(
x,
top,
up_genes = NULL,
down_genes = NULL,
cat_col = c(BP = "#E64B35", CC = "#4DBBD5", MF = "#00A087"),
up_col = "#AE2A8A",
down_col = "#7B87BF",
padjust_col = c("#FF906F", "#861D30"),
bg.col = "gray95",
max_width = 2.5,
fontsize = 8,
fontface = "bold",
fontfamily = "sans",
...
)
```

### `herb_enricher()`

A wrapper function of clusterProfiler::enricher() to do herb-target enrichment analysis.

A wrapper function of clusterProfiler::enricher() to do herb-target enrichment analysis.

Usage:
```r
herb_enricher(
genes,
bg_universe = NULL,
type = c("Herb_pinyin_name", "Herb_cn_name", "Herb_en_name"),
pvalueCutoff = 0.05,
qvalueCutoff = 0.2,
pAdjustMethod = "BH",
minGSSize = 10,
maxGSSize = 5000
)
```

参数解析:
- `genes`: Character vector. Query genes such as DEGs or disease targets.
- `bg_universe`: Character vector or NULL. Background gene set. If NULL, defaults to all unique targets in TCMDATA.
- `type`: Character. Specifies which column of tcm_data to use as herb names. Default is "Herb_pinyin_name".
- `pvalueCutoff`: Numeric. P-value threshold for significance (default 0.05).
- `qvalueCutoff`: Numeric. Q-value (FDR) threshold for significance (default 0.2).
- `pAdjustMethod`: Character. Multiple testing correction method (default "BH").
- `minGSSize`: Integer. Minimum herb gene-set size to include (default 10).
- `maxGSSize`: Integer. Maximum herb gene-set size to include (default 5000).

使用示例:
```r
library(enrichplot)
data("dn_gcds")
enrich_res <- herb_enricher(genes = dn_gcds, type = "Herb_pinyin_name")
head(enrich_res)

p <- dotplot(enrich_res)
print(p)
```

### `interpret_enrichment()`

Interpret enrichment analysis results with AI
Uses a large language model to generate a structured interpretation of enrichment results (GO, KEGG, MSigDB, herb enrichment, etc.).

Interpret enrichment analysis results with AI
Uses a large language model to generate a structured interpretation of enrichment results (GO, KEGG, MSigDB, herb enrichment, etc.).

Usage:
```r
interpret_enrichment(
x,
top_n = 10,
max_genes = 5L,
audience = "researcher",
language = c("zh", "en"),
role = NULL,
system = NULL,
prompt = NULL,
model = NULL
)
```

参数解析:
- `x`: An enrichResult object or a data.frame with enrichment
results.
- `top_n`: Integer. Number of top terms to include. Default 10.
- `max_genes`: Integer. Max genes shown per term in the compressed
context. Default 5. Increase (e.g. 20) when you want the model
to reason about the full gene list of specific terms.
- `audience`: Character. Preset: "researcher",
"wetlab", "paper". Or any free-text description
like "clinical doctor unfamiliar with bioinformatics".
- `language`: Character. "zh" or "en".
Default "zh".
- `role`: Character or NULL. The AI's identity / expertise description.
Replaces the default "You are a bioinformatics expert..." opening
line while keeping all other system-prompt constraints intact. Ignored
when system is also provided (full override takes precedence).
Example: "You are a TCM pharmacologist specialising in
herb-target interactions.".
- `system`: Character or NULL. Full system prompt override. Replaces
the entire generated system prompt when provided (including role,
audience, and language instructions).
- `prompt`: Character or NULL. Custom user prompt instruction.
Replaces the default instruction; the compressed data context
is always appended.
- `model`: A model identifier, LanguageModelV1 object, or NULL.

使用示例:
```r
enrich_res <- herb_enricher(genes = my_genes)
interpret_enrichment(enrich_res)

# Custom role — keeps audience/language/constraints
interpret_enrichment(enrich_res,
role = "You are a TCM pharmacologist specialising in herb-target networks.")

# Show more genes per term
interpret_enrichment(enrich_res, max_genes = 20)
```
