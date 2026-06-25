# Generate color palette

Generate a color palette. For discrete palettes, hex values are standard
up to n = 100, but for larger values of n, they are generated using
[`Polychrome::createPalette()`](https://rdrr.io/pkg/Polychrome/man/createPalette.html).

## Usage

``` r
permuteDEpalette(type = "discrete", n = NULL, palette = NULL, swatch = FALSE)
```

## Arguments

- type:

  A character string indicating the palette type. Permitted values are
  "discrete" and "gradient". Defaults to discrete.

- n:

  Number of colors. Default = `NULL` will return all of the pre-set
  colors in the palette.

- palette:

  A character string indicating the palette name. Permitted values are
  "choir", "archr", "inferno", and "frozen". Default = `NULL` will use
  "choir" for discrete colors and "inferno" for gradient colors.
  Alternately, provide a vector of color values to use as starting
  values for the color palette.

- swatch:

  A Boolean value indicating whether to plot a swatch of the palette.

## Value

Returns a vector of n hex values.

## Details

The "choir" palette is adapted from function `CHOIR::CHOIRpalette` in R
package `CHOIR` (Sant et al. 2025). The "archr" palette is adapted from
palette "stallion" in R package `ArchR` (Granja & Corces et al. 2020).
