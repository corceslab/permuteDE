# Generate histogram plot(s) of permutation test results

This function takes the output from function `permuteDE` and creates a
histogram plot showing the permutation test result for each indicated
split. Each histogram includes a vertical line at the observed number of
DE features in the unpermuted "true" comparison.

## Usage

``` r
plotHistogram(
  input,
  use_splits = NULL,
  title = NULL,
  subtitle = NULL,
  label_pvalue = TRUE
)
```

## Arguments

- input:

  Output from function `permuteDE`.

- use_splits:

  A character string or vector containing the names of splits to use.
  Defaults to `NULL`, which will try all splits.

- title:

  Character string indicating the plot title. Default = `NULL` sets a
  title automatically for each split.

- subtitle:

  Character string indicating the plot subtitle. Default = `NULL`
  automatically generates a subtitle.

- label_pvalue:

  A Boolean value indicating whether to label the permutation test
  p-value. Defaults to TRUE.

## Value

If only one split is provided, a single `ggplot2` object, otherwise a
named list of `ggplot2` objects, where each element corresponds to a
split.
