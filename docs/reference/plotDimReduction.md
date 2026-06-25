# Plot dimensionality reduction

This function will generate a dimensionality reduction plot colored
according to the selected metric.

## Usage

``` r
plotDimReduction(
  reduction,
  input = NULL,
  split_labels = NULL,
  feature_values = NULL,
  feature_name = NULL,
  use_cells = NULL,
  color_by = NULL,
  permutation_test_alpha = 1,
  label_splits = FALSE,
  label_statistics = FALSE,
  palette = NULL,
  color_limits = NULL,
  fix_coords = TRUE,
  ...
)
```

## Arguments

- reduction:

  Dimensionality reduction cell coordinates. If there are more than two
  columns, only the first two will be used.

- input:

  Output from function `permuteDE`.

- split_labels:

  An character vector containing the split labels for each cell in
  order. Default = `NULL` will assume all cells belong to the same
  split.

- feature_values:

  A vector containing the value for a feature for each cell in order,
  used if `color_by` = "feature".

- feature_name:

  A string indicating the name of the plotted feature if `color_by` =
  "feature". Used for legend title.

- use_cells:

  A vector of cell names to subset the dimensionality reduction cell
  coordinates to. Default = `NULL` will use all cells.

- color_by:

  A character vector indicating what metric to use to color each split.
  Permitted values are "n_sig" (number of significant DE features),
  "pvalue" (permutation test p-value), "split" (values provided to
  parameter `split_labels`), and "feature" (values provided to parameter
  `feature_values`). Default = `NULL` will not apply a color scheme.

- permutation_test_alpha:

  A numeric value indicating the significance level to apply to the
  permutation test results. Splits that do not pass this threshold will
  be grayed out. Default = 1 applies no threshold.

- label_splits:

  A Boolean value indicating whether to label the splits. Defaults to
  `FALSE`.

- label_statistics:

  A Boolean value indicating whether to label the values of the selected
  metric. Defaults to `FALSE`.

- palette:

  A character string indicating the palette name. Permitted values are
  "choir", "archr", "inferno", and "frozen". Default = `NULL` will use
  "choir" for discrete colors and "inferno" for gradient colors.
  Alternately, provide a vector of color values to use as starting
  values for the color palette.

- color_limits:

  A vector with the minimum and maximum values indicated by `color_by`
  ("n_sig" or "pvalue") for display on the color bar legend. Default =
  `NULL` sets limits automatically.

- fix_coords:

  A Boolean value indicating whether the aspect ratio of the x and y
  axis should be 1.

- ...:

  Extra parameters passed to
  [`Seurat::DimPlot()`](https://satijalab.org/seurat/reference/DimPlot.html)
  or
  [`Seurat::FeaturePlot()`](https://satijalab.org/seurat/reference/FeaturePlot.html).

## Value

A dimensionality reduction plot.
