# Plot levels of a single feature across groups

This functions takes the output from function `runDE` and plots the
expression of a single feature across the replicates in each group for
each indicated split.

## Usage

``` r
plotFeature(
  input,
  feature_name,
  use_splits = NULL,
  normalization_method = "cpm",
  plot_type = "boxplot",
  title = NULL,
  subtitle = NULL,
  label_replicates = FALSE,
  label_statistics = TRUE
)
```

## Arguments

- input:

  Output from function `runDE`.

- feature_name:

  A character string indicating which feature to plot.

- use_splits:

  A character string or vector containing the names of splits to use.
  Defaults to `NULL`, which will try all splits.

- normalization_method:

  A character string indicating which normalization method to apply to
  the replicate x feature matrix. Permitted values are "cpm", "log_cpm",
  and "none". Defaults to "cpm".

- plot_type:

  A character string indicating the type of plot. Permitted values are
  "boxplot", "bar_se", "bar_sd", "beeswarm". Defaults to "boxplot".

- title:

  A character string indicating the plot title. Default = `NULL` sets a
  title automatically for each split.

- subtitle:

  Character string indicating the plot subtitle. Default = `NULL`
  automatically generates a subtitle.

- label_replicates:

  A Boolean value indicating whether to label the replicates. Defaults
  to `FALSE`.

- label_statistics:

  A Boolean value indicating whether to label the differential
  expression analysis LFC and p-value. Defaults to TRUE.

## Value

If only one split is provided, a single `ggplot2` object, otherwise a
named list of `ggplot2` objects, where each element corresponds to a
split.
