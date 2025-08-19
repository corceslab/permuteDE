#' GG Histogram
#'
#' Create a customized histogram with optional density overlay and vertical reference lines.
#'
#' @param x Numeric vector of values to plot.
#' @param xlabel Label for the x-axis.
#' @param ylabel Label for the y-axis.
#' @param addDensity Logical; if `TRUE`, adds a density curve (default: `FALSE`).
#' @param bins Number of bins in the histogram (default: 20).
#' @param baseSize Base font size used in `themeCorcesRegular()` (default: 14).
#' @param histFill Fill color of histogram bars (default: "#3B9AB2").
#' @param histColor Outline color of histogram bars (default: `NA`).
#' @param densityFill Fill color of the density area (default: "#3B9AB2").
#' @param densityColor Line color of the density curve (default: "#3B9AB2").
#' @param size Line width for histogram and density elements (default: 1.25).
#' @param histAlpha Transparency (alpha) of the histogram fill (default: 0.85).
#' @param densityAlpha Transparency (alpha) of the density fill (default: 0.15).
#' @param title Title for the plot (default: empty string).
#' @param ratioYX Desired aspect ratio of y over x for `coord_equal()` (default: 0.8).
#' @param vline Numeric value (or vector) for vertical reference lines (default: `NULL`).
#' @param vlineColor Color for vertical lines (default: "red").
#' @param vlineType Line type for vertical lines (e.g., "solid", "dashed") (default: "solid").
#' @param vlineWidth Line width for vertical lines (default: 1).
#'
#' @return A ggplot2 object.
#' @export

ggPlotHistogram <- function(x, xlabel = "values", ylabel = "Count", 
                            addDensity = FALSE, bins = 20, baseSize = 14,
                            histFill = "#3B9AB2", histColor = NA,
                            densityFill = "#3B9AB2", densityColor = "#3B9AB2", title = "",
                            size = 1.25, histAlpha = 0.85, densityAlpha = 0.15, ratioYX = 0.8,
                            vline = NULL, vlineColor = "red", vlineType = "solid", vlineWidth = 1) {
  
  stopifnot(is.numeric(x))
  
  df <- data.frame(y = x)
  
  p <- ggplot(df, aes(y)) + 
    geom_histogram(bins = bins, fill = histFill, color = histColor, size = size, alpha = histAlpha)
  
  if (addDensity) {
    p <- p +
      geom_density(aes(y = ..count..), fill = densityFill, color = NA, alpha = densityAlpha, size = size) +
      stat_density(aes(y = ..count..), color = densityColor, geom = "line", size = size)
  }
  
  if (!is.null(vline)) {
    p <- p + geom_vline(xintercept = vline, color = vlineColor, linewidth = vlineWidth, linetype = vlineType)
  }
  
  p <- p +
    themeCorcesRegular(base_size = baseSize) +
    xlab(xlabel) + ylab(ylabel)
  
  # Include vline in xlim computation
  xlim_raw <- range(c(df$y, vline), na.rm = TRUE)
  xlim <- extendrange(xlim_raw, f = 0.1)
  
  data <- ggplot_build(p)$data[[1]]
  ylim <- c(0, max(data$y, na.rm = TRUE) * 1.1)
  ratioXY <- ratioYX * diff(xlim) / diff(ylim)
  
  p <- p + coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE)
  
  return(p)
}


#' ggplot2 default theme for Corces Lab
#'
#' This function returns a ggplot2 theme that is black borded with black font.
#' 
#' @param color color of text and lines of the plot
#' @param base_size the size of the font for the axis labels
#' @param base_line_size the line width for most lines
#' @param base_rect_size the line width for rectangular boxes
#' @param axis_title_size the font size of the axis title
#' @param plot_title_size the font size of the plot title
#' @param plot_margin_cm the margin around the plot in centimeters
#' @param legend_text_size 0.75*base_size
#' @param legend_position the placement of the legend in the plot
#' @param axis_tick_length_mm axis tick length in mm
#' @param rotate_x_axis_text_90 Boolean value indicating whether to rotate the x-axis text by 90 degrees
#' @param rotate_y_axis_text_90 Boolean value indicating whether to rotate the y-axis text by 90 degrees
#' @export
themeCorcesRegular <- function(color = "black",
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
                               rotate_y_axis_text_90 = FALSE){
  theme <- theme(
    text = element_text(color = color, size = base_size),
    axis.text = element_text(color = color, size = base_size), 
    axis.title = element_text(color = color, size = axis_title_size),
    plot.title = element_text(color = color, size = plot_title_size),
    plot.margin = unit(c(plot_margin_cm,plot_margin_cm,plot_margin_cm,plot_margin_cm), "cm"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = color, size = (4/3) * base_rect_size * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
    axis.ticks.length = unit(axis_tick_length_mm, "mm"), 
    axis.ticks = element_line(color = color, size = base_line_size * (4/3) * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(color = color, size = legend_text_size),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent"),
    legend.position = legend_position,
    strip.text = element_text(size = base_size, color="black"),
    plot.background = element_rect(fill = "transparent", color = NA)
  )
  if(rotate_x_axis_text_90){
    theme <- theme %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  if(rotate_y_axis_text_90){
    theme <- theme %+replace% theme(axis.text.y = element_text(angle = 90, vjust = 1))
  }
  return(theme)
}
