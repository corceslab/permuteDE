# Generate combinations

For a set of replicates that belong to two groups, this function will
generate a set of randomized group combinations.

## Usage

``` r
getCombinations(
  object = NULL,
  metadata = NULL,
  replicate_labels = NULL,
  group_labels = NULL,
  use_cells = NULL,
  n_replicates = NULL,
  n_group1 = NULL,
  n_combinations = 1000,
  confound_check = NULL,
  progress_message = "",
  random_seed = 1,
  verbose = TRUE
)
```

## Arguments

- object:

  An optional 'Seurat' or 'SingleCellExperiment' object. If `NULL`,
  `countCombinations` will expect either (1) vector input to parameters
  `replicate_label` and `group_label` or (2) numeric input to parameters
  `n_replicates` and `n_group1`.

- metadata:

  An optional dataframe containing relevant metadata columns
  corresponding to the data provided to parameter `object`. Default =
  `NULL` looks for metadata in `object` or other provided inputs.

- replicate_labels:

  A string indicating the name of the metadata column containing the
  biological replicate labels or a character vector containing the
  biological replicate labels in order.

- group_labels:

  A string indicating the name of the column containing the two
  comparison group labels or a character vector containing the
  comparison labels in order.

- use_cells:

  A vector of cell names to subset prior to calculating possible group
  comibinations. Default = `NULL` will use all cells.

- n_replicates:

  A numeric value indicating the total number of replicates.
  Alternately, a vector can be provided to generate combinations when
  shuffling separately within multiple sets. Defaults to `NULL`.

- n_group1:

  A numeric value indicating the number of replicates in one group
  (doesn't matter which). Alternately, a vector can be provided to
  generate combinations when shuffling separately within multiple sets.
  Defaults to `NULL`.

- n_combinations:

  A numeric value indicating the number of combinations to generate.
  Defaults to 1000.

- confound_check:

  An optional dataframe of covariates for which to exclude permutations
  that are perfectly confounded. Defaults to `NULL`.

- progress_message:

  A character string indicating additional progress messaging (internal
  use). Defaults to "".

- random_seed:

  A numeric value indicating the random seed to be used. Defaults to 1.

- verbose:

  A Boolean value indicating whether to use verbose output during the
  execution of this function. Defaults to `TRUE`. Can be set to `FALSE`
  for a cleaner output.

## Value

Returns a matrix where each column contains a combination of index
values indicating which replicates to assign to the first group.

## Details

Users may provide input in three ways:

\(1\) A vector of replicate labels and a vector of group labels (in
order)

\(2\) The total number of replicates (n) and the number of replicates in
one group (k) (doesn't matter which group).

\(3\) Column names indicating replicate and group metadata columns in a
provided 'Seurat' or 'SingleCellExperiment' object or a provided
metadata dataframe.
