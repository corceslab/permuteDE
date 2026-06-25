# Perform permutation test

This function performs a permutation test by shuffling the group labels
for the given set of replicates and performing pseudobulk differential
expression between those permuted groups. The metric used by the
permutation test is the number of significantly differentially expressed
features for a given significance level `alpha` and log fold change
threshold `lfc_threshold`.

## Usage

``` r
permuteDE(
  input,
  alpha = 0.05,
  lfc_threshold = 0,
  n_iterations = 1000,
  use_splits = NULL,
  permute_by = NULL,
  permute_within = NULL,
  min_DE = 1,
  return_all = FALSE,
  random_seed = 1,
  n_cores = NULL,
  verbose = TRUE
)
```

## Arguments

- input:

  Output from function `runDE`: a list containing differential
  expression results, pseudobulk values, and parameters used.

- alpha:

  A numeric value indicating the significance level used for permutation
  test comparisons of the number of differentially expressed features.
  Defaults to 0.05.

- lfc_threshold:

  A numeric value indicating the minimum absolute value log fold change
  for a feature to be counted as a "hit". Default = 0 disregards log
  fold change when counting hits.

- n_iterations:

  A numeric value indicating the number of iterations run for the
  permutation test. Defaults to 1000. Computational time increases
  approximately linearly with the number of iterations.

- use_splits:

  A vector containing the names of splits to use. Defaults to `NULL`,
  which will try all splits.

- permute_by:

  An optional character string indicating the name of the column in
  `input$metadata$group_key` indicating partitions of the data within
  which group labels should be shuffled together as a unit. Alternately,
  a vector ordered according to the values in
  `input$metadata$group_key`. Only for specific use cases such as
  cell-level tests. Default = `NULL` will shuffle group labels across
  all replicates.

- permute_within:

  An optional character string indicating the name of the column in
  `input$metadata$group_key` indicating partitions of the data within
  which group labels should be shuffled separately. Alternately, a
  vector ordered according to the values in `input$metadata$group_key`.
  Only for specific use cases such as cell-level tests. Default = `NULL`
  will shuffle group labels across all replicates.

- min_DE:

  A numeric value indicating the minimum number of differentially
  expressed features between the true group labels for a split, below
  which permutations will not be run. Defaults to 1. Set to 0 to run
  permutation test for all splits, regardless of the true number of
  DEGs.

- return_all:

  A Boolean value indicating whether to store and return all DE results
  (log fold changes and p-values per feature per split) for every single
  permutation. Defaults = `FALSE` will return only high-level
  permutation test results. Note that setting this to `TRUE` will
  substantially increase the size of the returned output.

- random_seed:

  A numeric value indicating the random seed to be used. Defaults to 1.

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

- permutation_test_results:

  Dataframe containing the permutation test results by split

- permutation_DE_summary:

  Dataframe containing the permutation DE summary metrics by split

- permutation_DE_results:

  If parameter 'return_all' is TRUE, dataframe DE results for each
  feature, by split, for each iteration

- metadata:

  List recording additional characteristics of the data: the number of
  significant DE features from runDE, the group indices for each
  iteration of the permutation test, and runtime

- parameters:

  List recording parameter values used

## Details

As input, this function requires the output from function `runDE`
containing the pseudobulk (or cell-level) values for each feature and
the differential expression results for the true group labels.

Permutation is always performed at the replicate level, rather than the
cell level, and is performed without replacement, such that each
iteration is a unique permutation of the group labels.

Notably, this permutation test does not pass judgment on any individual
feature, rather, it is intended to assess how many false positive
significant differentially expressed features can be expected by chance.
In addition, it can be used to characterize the log fold change and
significance observed for such false positives to help users prioritize
reliable DE results.
