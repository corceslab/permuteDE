# Calculate pseudobulk values

This function will generate pseudobulk expression matri(ces) from a
provided Seurat object, SingleCellExperiment object, or feature x cell
matrix. Multiple pseudobulk matrices can be generated across different
divisions of a dataset, such as cell types.

## Usage

``` r
getPseudobulk(
  object,
  metadata = NULL,
  replicate_labels = NULL,
  split_labels = NULL,
  use_cells = NULL,
  min_cells_per_split = 100,
  min_cells_per_replicate = 10,
  min_replicates_per_split = 6,
  min_cells_per_feature = 10,
  min_prop_cells_per_feature = 0.1,
  filter = TRUE,
  pseudobulk = "generate",
  use_assay = NULL,
  use_layer = NULL,
  n_cores = NULL,
  verbose = TRUE
)
```

## Arguments

- object:

  An object of class `Seurat`, `SingleCellExperiment`, or `matrix`. Data
  supplied as class `matrix` should be a feature x cell matrix.

- metadata:

  An optional dataframe containing relevant metadata columns
  corresponding to the data provided to parameter `object`. Default =
  `NULL` looks for metadata in `object` or other provided inputs.

- replicate_labels:

  A string indicating the name of the metadata column containing the
  biological replicate labels or a character vector containing the
  biological replicate labels in order. The biological replicate labels
  are used to construct/define the pseudobulks. Default = `NULL` allowed
  only when `pseudobulk` is "none".

- split_labels:

  A string indicating the name of a column by which to split the cells
  prior to pseudobulking and performing differential expression (e.g.,
  cell types). Alternately, a character vector containing the split
  labels for each cell in order. Results will be returned for each
  unique value indicated by `split_labels`. Default = `NULL` will
  generate a single pseudobulk matrix that includes all cells.

- use_cells:

  A vector of cell names to subset the object to prior to subsequent
  pseudobulk steps. Default = `NULL` will use all cells.

- min_cells_per_split:

  A numeric value indicating the minimum number of cells within one
  split. Pseudobulk steps will not be performed for splits with fewer
  cells. Defaults to 100.

- min_cells_per_replicate:

  A numeric value indicating the minimum number of cells within one
  replicate for one split. Pseudobulk steps will not be performed for
  replicates with fewer cells for that split. Defaults to 10.

- min_replicates_per_split:

  A numeric value indicating the minimum number of distinct replicates
  represented within one split. Pseudobulk steps will not be performed
  for splits with fewer replicates. Defaults to 6.

- min_cells_per_feature:

  A numeric value indicating the minimum number of cells (within a
  split) with expression of a feature. Pseudobulk expression will not be
  calculated for features expressed in fewer cells. Defaults to 10.

- min_prop_cells_per_feature:

  A numeric value indicating the minimum proportion of cells (within a
  split) with expression of a feature. Pseudobulk expression will not be
  calculated for features expressed in fewer cells. Defaults to 0.1.

- filter:

  A Boolean value indicating whether to remove features from the
  matrices in the output (`TRUE`) or simply list features that do not
  meet the criteria as part of the metadata (`FALSE`). Defaults to
  `TRUE`.

- pseudobulk:

  A string indicating whether to actually pseudobulk ("generate"), or
  simply apply the filtering thresholds ("none"). Defaults to
  "generate".

- use_assay:

  A string indicating the assay to use in the provided object. Default =
  `NULL` will choose the current active assay for `Seurat` objects and
  the `counts` assay for `SingleCellExperiment` objects.

- use_layer:

  For `Seurat` objects, a string or vector indicating the layer
  (previously known as slot) to use in the provided object. Default =
  `NULL` will use the `counts` layer.

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

- PB_values:

  A list of feature x replicate matri(ces) containing pseudobulk values
  for each feature, one matrix per split

- metadata:

  Dataframe record of quality control metrics for each split

- parameters:

  List recording parameter values used
