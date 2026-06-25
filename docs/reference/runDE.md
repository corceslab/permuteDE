# Run differential expression analysis to compare two groups

This function identifies differentially expressed features between two
groups using indicated differential expression analysis methods.

## Usage

``` r
runDE(
  object,
  metadata = NULL,
  replicate_labels = NULL,
  group_labels,
  split_labels = NULL,
  reference_group = NULL,
  design = NULL,
  use_cells = NULL,
  pseudobulk = "generate",
  de_method = "edgeR",
  de_test = "LRT",
  de_params = list(),
  return_raw_de = FALSE,
  normalize_prefilter = FALSE,
  p_adjust_method = "fdr",
  min_cells_per_split = 100,
  min_cells_per_replicate = 10,
  min_replicates_per_split = 6,
  min_replicates_per_group = 3,
  min_cells_per_feature = 10,
  min_prop_cells_per_feature = 0.1,
  force_balance = FALSE,
  use_assay = NULL,
  use_layer = NULL,
  random_seed = 1,
  n_cores = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  An object of class `Seurat`, `SingleCellExperiment`, or `matrix`. Data
  supplied as class `matrix` may be either a feature x cell matrix or a
  pre-computed pseudobulk feature x replicate matrix. Note that raw
  counts are expected, and the normalization method applied during
  differential expression analysis differs across the methods and tests.

- metadata:

  An optional dataframe containing relevant metadata columns
  corresponding to the data provided to parameter `object`. Default =
  `NULL` looks for metadata in `object` or other provided inputs.

- replicate_labels:

  A string indicating the name of the metadata column containing the
  biological replicate labels or a vector containing the biological
  replicate labels in order. For pseudobulk DE analysis, the biological
  replicate labels are used to construct/define the pseudobulks. Input
  is not required for cell-level DE analysis.

- group_labels:

  A string indicating the name of the metadata column containing the two
  comparison group labels or a vector containing the comparison labels
  in order.

- split_labels:

  A string indicating the name of a metadata column by which to split
  the cells prior to pseudobulking and performing differential
  expression (e.g., cell types). Alternately, a vector containing the
  split labels for each cell in order. Results will be returned for each
  unique value indicated by `split_labels`. Default = `NULL` will run
  pseudobulk differential expression on all cells together.

- reference_group:

  A string specifying the reference group. Defaults to `NULL`, in which
  case the first value alphabetically is used as the reference.

- design:

  An optional string specifying a model formula for more complex
  designs. Last term in formula must correspond to group labels. Default
  = `NULL` will run a pairwise group comparison (~ group) based on the
  input provided to parameter `group_labels`.

- use_cells:

  A vector of cell names to subset the object to prior to subsequent
  pseudobulk and differential expression steps. Default = `NULL` will
  use all cells.

- pseudobulk:

  A string indicating pseudobulk handling. Permitted values are:
  "generate" (pseudobulk matrices will be generated), "supplied"
  (pseudobulk matrix was supplied by the user to parameter `object`), or
  "none" (pseudobulking will not be used, cell-level differential
  expression analysis will be run). Defaults to "generate".

- de_method:

  Which tool to use for differential expression analysis. Permitted
  values are "edgeR", "DESeq2", "limma", "presto", and "BPCells".
  Defaults to "edgeR".

- de_test:

  Which test to use for differential expression analysis. Available
  values are dependent on the `de_method`: "edgeR" ("LRT", "QLF",
  "exact"), "DESeq2" ("LRT", "Wald"), "limma" ("trend", "voom",
  "wilcox_cpm", "wilcox_log_cpm"), "presto" ("wilcox_cpm",
  "wilcox_log_cpm"), and "BPCells" ("wilcox_cpm", "wilcox_log_cpm").
  Defaults to "LRT".

- de_params:

  A list of lists containing additional parameters to be passed to
  specific DE functions. The name of each element must be the specific
  DE function to which those parameters are passed. Defaults to an empty
  list. The special top-level option `return_all_coefficients = TRUE`
  can be used with `return_raw_de = TRUE` for coefficient-based model
  tests to include raw results for all model coefficients.

- return_raw_de:

  A Boolean value indicating whether to also return the raw output from
  the selected DE method/test. Defaults to `FALSE`.

- normalize_prefilter:

  A Boolean value indicating whether normalization should be applied
  before (`TRUE`) or after (`FALSE`) filtering out features with low
  counts. Defaults to `FALSE`.

- p_adjust_method:

  A string indicating which multiple comparison adjustment to use. For
  permitted values, see
  [`stats::p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html).
  Defaults to "fdr" (Benjamini & Hochberg, 1995). For advanced users,
  this parameter can also be set to "fdrtool" to use the `fdrtool`
  package (applied to raw p-values by default, or for DESeq2 Wald, apply
  to z-scores by setting parameter `de_params` to
  `list(fdrtools = list(statistic = "zscore"))`).

- min_cells_per_split:

  A numeric value indicating the minimum number of cells within one
  split. Pseudobulk and differential expression steps will not be
  performed for splits with fewer cells. Defaults to 100.

- min_cells_per_replicate:

  A numeric value indicating the minimum number of cells within one
  replicate for one split. Pseudobulk steps will not be performed for
  replicates with fewer cells for that split. Defaults to 10.

- min_replicates_per_split:

  A numeric value indicating the minimum number of distinct replicates
  represented within one split. Pseudobulk expression and differential
  expression will not be performed for splits with fewer replicates.
  Defaults to 6.

- min_replicates_per_group:

  A numeric value indicating the minimum number of distinct replicates
  represented within each of the two comparison groups. Pseudobulk and
  differential expression steps will not be performed for splits with
  fewer replicates. Defaults to 3.

- min_cells_per_feature:

  A numeric value indicating the minimum number of cells (within a
  split) with expression of a feature. Pseudobulk and differential
  expression will not be calculated for features expressed in fewer
  cells. Defaults to 10.

- min_prop_cells_per_feature:

  A numeric value indicating the minimum proportion of cells (within a
  split) with expression of a feature. Pseudobulk and differential
  expression will not be calculated for features expressed in fewer
  cells. Defaults to 0.1.

- force_balance:

  A boolean indicating whether to force the two comparison groups to
  have the same sample size. Defaults to `FALSE`. If `TRUE`, the larger
  group will be randomly downsampled to the size of the smaller group.

- use_assay:

  A string indicating the assay to use in the provided object. Default =
  `NULL` will choose the current active assay for `Seurat` objects and
  the `counts` assay for `SingleCellExperiment` objects.

- use_layer:

  For `Seurat` objects, a string or vector indicating the layer
  (previously known as slot) to use in the provided object. Default =
  `NULL` will use the `counts` layer.

- random_seed:

  A numerical value indicating the random seed to be used. Defaults
  to 1. Only relevant in this function when parameter
  `force_balance = TRUE`.

- n_cores:

  A numeric value indicating the number of cores to use for
  parallelization. Default = `NULL` will use the number of available
  cores minus 2.

- verbose:

  A Boolean value indicating whether to use verbose output during the
  execution of this function. Defaults to `TRUE`. Can be set to `FALSE`
  for a cleaner output.

## Value

Returns a list containing the following elements:

- DE_results:

  Dataframe containing DE results for each feature, by split

- PB_values:

  If using pseudobulk data, a list of feature x replicate matri(ces)
  containing pseudobulk values for each feature, one matrix per split

- cell_values:

  Alternately, if using cell-level data, a list of feature x cell
  matri(ces) containing counts for each feature, one matrix per split

- metadata:

  List recording characteristics of the data and runtime

- parameters:

  List recording parameter values used

## Details

By default, pseudobulk matri(ces) are generated or supplied by the user,
then used to run pseudobulk differential expression. The following
existing tools are supported: `edgeR`, `DESeq2`, `limma`, and the
Wilcoxon rank-sum test. Alternately, users may skip pseudobulking and
run cell-level differential expression (not recommended in most cases).
