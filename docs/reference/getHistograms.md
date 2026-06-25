# Generate histogram plots of permutation test results

This function takes the output from function `permuteDE` and creates a
list containing a histogram plot showing the permutation test result for
each split. Each histogram includes a vertical line at the observed
number of DE features in the unpermuted "true" comparison.

## Usage

``` r
getHistograms(input, title = NULL, subtitle = NULL, label_pvalue = TRUE)
```

## Arguments

- input:

  Output from function `permuteDE`.

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

A named list of `ggplot2` objects, where each element corresponds to a
split.
