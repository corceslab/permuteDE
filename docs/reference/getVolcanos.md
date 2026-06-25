# Generate volcano plots from differential expression results

This function takes the output from function `runDE` and creates a list
containing a volcano plot for each subset of differential expression
(DE) results. The DE results are subset by the "split" column, and each
split produces its own volcano plot.

## Usage

``` r
getVolcanos(
  input,
  alpha = 0.05,
  lfc_threshold = 0.5,
  title = NULL,
  subtitle = NULL,
  n_max_label = 10,
  center = TRUE
)
```

## Arguments

- input:

  Output from function `runDE` or a list containing (at minimum) a
  dataframe named "DE_results" including columns:

  gene

  :   Gene identifiers (character).

  lfc

  :   Log2 fold change values (numeric).

  padj

  :   Adjusted p-values (numeric).

  split

  :   (Optionally) Grouping variable used to subset results (character
      or factor).

- alpha:

  A numeric value indicating the significance level used for permutation
  test comparisons of the number of differentially expressed features.
  Defaults to 0.05.

- lfc_threshold:

  A numeric value indicating the minimum absolute value log fold change
  for a gene to be counted as a "hit". Defaults to 0.5. Set to 0 to
  disregard log fold change when counting hits.

- title:

  Character string indicating the plot title. Default = `NULL` sets a
  title automatically for each split.

- subtitle:

  Character string indicating the plot subtitle. Default = `NULL`
  automatically generates a subtitle describing the DE method used.

- n_max_label:

  A numeric value indicating how many of the top significant DE features
  to label. Defaults to 10.

- center:

  A Boolean value indicating whether to center the x-axis at 0. Defaults
  to `TRUE`.

## Value

A named list of `ggplot2` objects, where each element corresponds to a
volcano plot for one split of the DE results.
