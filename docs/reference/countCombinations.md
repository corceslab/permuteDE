# Calculate the number of possible group combinations

Given two groups of replicates, this function calculates the number of
unique ways to divide the replicates into two groups while preserving
group sizes. The formula is: C = n!/(k!(n-k)!). (see base function
`choose`)

## Usage

``` r
countCombinations(
  object = NULL,
  metadata = NULL,
  replicate_labels = NULL,
  group_labels = NULL,
  use_cells = NULL,
  n_replicates = NULL,
  n_group1 = NULL
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

  A numeric value indicating the total number of replicates. Defaults to
  `NULL`.

- n_group1:

  A numeric value indicating the number of replicates in one group
  (doesn't matter which). Defaults to `NULL`.

## Value

Returns a numeric value indicating the number of all possible group
combinations.

## Details

Users may provide input in three ways:

\(1\) A vector of replicate labels and a vector of group labels (in
order)

\(2\) The total number of replicates (n) and the number of replicates in
one group (k) (doesn't matter which group).

\(3\) Column names indicating replicate and group metadata columns in a
provided 'Seurat' or 'SingleCellExperiment' object or a provided
metadata dataframe.
