# Generate color palette

Generate a color palette. For discrete palettes, hex values are standard
up to n = 100, but for larger values of n, they are generated using
[`Polychrome::createPalette()`](https://rdrr.io/pkg/Polychrome/man/createPalette.html).

## Usage

``` r
palette_permuteDE(
  type = "discrete",
  n = NULL,
  palette_name = NULL,
  swatch = FALSE
)
```

## Arguments

- type:

  A character string indicating the palette type. Permitted values are
  "discrete" and "gradient". Defaults to discrete.

- n:

  Number of colors. Default = `NULL` will return all of the pre-set
  colors in the palette.

- palette_name:

  A character string indicating the palette name. Permitted values are
  "choir", "archr", "corces_cold", and "corces_warm". Default = `NULL`
  will use "choir" when `type` is "discrete" and "corces_cold" when
  `type` is "gradient.

- swatch:

  A Boolean value indicating whether to plot a swatch of the palette.

## Value

Returns a vector of n hex values.

## Details

The "choir" palette is adapted from function
[`CHOIR::CHOIRpalette`](https://rdrr.io/pkg/CHOIR/man/CHOIRpalette.html)
in R package `CHOIR` (Sant et al. 2025). The "archr" palette is adapted
from palette "stallion" in R package `ArchR` (Granja & Corces et al.
2020).
