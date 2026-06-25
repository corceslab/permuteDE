# Generate volcano plot(s) from differential expression results

This function takes the output from function `runDE` and creates volcano
plots for each indicated split.

## Usage

``` r
plotVolcano(
  input,
  alpha = 0.05,
  lfc_threshold = 0,
  use_splits = NULL,
  title = NULL,
  subtitle = NULL,
  n_max_label = 10,
  label_features = NULL,
  center = TRUE
)
```

## Arguments

- input:

  Output from function `runDE` or a list containing (at minimum) a
  dataframe named "DE_results" including columns:

  feature

  :   Feature identifiers (character).

  lfc

  :   Log2 fold change values (numeric).

  padj

  :   Adjusted p-values (numeric).

  split

  :   Grouping variable used to subset results (character or factor).

- alpha:

  A numeric value indicating the significance level used for permutation
  test comparisons of the number of differentially expressed features.
  Defaults to 0.05.

- lfc_threshold:

  A numeric value indicating the minimum absolute value log fold change
  for a feature to be counted as a "hit". Default = 0 disregards log
  fold change when counting hits.

- use_splits:

  A character string or vector containing the names of splits to use.
  Defaults to `NULL`, which will try all splits.

- title:

  Character string indicating the plot title. Default = `NULL` sets a
  title automatically for each split.

- subtitle:

  Character string indicating the plot subtitle. Default = `NULL`
  automatically generates a subtitle describing the DE method used.

- n_max_label:

  A numeric value indicating how many of the top significant DE features
  to label. Defaults to 10.

- label_features:

  An optional vector containing feature names to label. Default = `NULL`
  will label top significant DE features.

- center:

  A Boolean value indicating whether to center the x-axis at 0. Defaults
  to `TRUE`.

## Value

If only one split is provided, a single `ggplot2` object, otherwise a
named list of `ggplot2` objects, where each element corresponds to a
split.
