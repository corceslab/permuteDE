# permuteDE ggplot theme

This function returns the ggplot2 theme used throughout the permuteDE
plotting functions.

## Usage

``` r
theme_permuteDE(
  color = "black",
  base_size = 10,
  base_line_size = 0.5,
  base_rect_size = 0.5,
  axis_title_size = 12,
  plot_title_size = 14,
  plot_margin_cm = 1,
  legend_text_size = 10,
  legend_position = "bottom",
  axis_tick_length_mm = 1,
  rotate_x_axis_text_90 = FALSE,
  rotate_y_axis_text_90 = FALSE
)
```

## Arguments

- color:

  Color of text and lines of the plot

- base_size:

  The size of the font for the axis labels

- base_line_size:

  The line width for most lines

- base_rect_size:

  The line width for rectangular boxes

- axis_title_size:

  The font size of the axis title

- plot_title_size:

  The font size of the plot title

- plot_margin_cm:

  The margin around the plot in centimeters

- legend_text_size:

  0.75\*base_size

- legend_position:

  The placement of the legend in the plot

- axis_tick_length_mm:

  Axis tick length in mm

- rotate_x_axis_text_90:

  Boolean value indicating whether to rotate the x-axis text by 90
  degrees

- rotate_y_axis_text_90:

  Boolean value indicating whether to rotate the y-axis text by 90
  degrees
