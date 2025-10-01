# ---------------------------------------------------------------------------
# Plotting-related functions
# ---------------------------------------------------------------------------

#' permuteDE ggplot theme
#'
#' This function returns the ggplot2 theme used throughout the permuteDE
#' plotting functions.
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
theme_permuteDE <- function(color = "black",
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
  custom_theme <- ggplot2::theme(
    text = ggplot2::element_text(color = color, size = base_size),
    axis.text = ggplot2::element_text(color = color, size = base_size),
    axis.title = ggplot2::element_text(color = color, size = axis_title_size),
    plot.title = ggplot2::element_text(color = color, size = plot_title_size),
    plot.margin = grid::unit(c(plot_margin_cm, plot_margin_cm, plot_margin_cm, plot_margin_cm), "cm"),
    panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.border = ggplot2::element_rect(fill = NA, color = color, size = (4/3) * base_rect_size),
    axis.ticks.length = grid::unit(axis_tick_length_mm, "mm"),
    axis.ticks = ggplot2::element_line(color = color, size = (4/3) * base_line_size),
    legend.key = ggplot2::element_rect(fill = "transparent", colour = NA),
    legend.text = ggplot2::element_text(color = color, size = legend_text_size),
    legend.background = ggplot2::element_rect(fill = "transparent"),
    legend.box.background = ggplot2::element_rect(fill = "transparent"),
    legend.position = legend_position,
    strip.text = ggplot2::element_text(size = base_size, color="black"),
    plot.background = ggplot2::element_rect(fill = "transparent", color = NA)
  )

  if(rotate_x_axis_text_90){
    custom_theme <- custom_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  }
  if(rotate_y_axis_text_90){
    custom_theme <- custom_theme + ggplot2::theme(axis.text.y = ggplot2::element_text(angle = 90, vjust = 1))
  }
  return(custom_theme)
}

#' Generate volcano plots from differential expression results
#'
#' This function takes the output from function \code{runDE} and creates a list
#' containing a volcano plot for each subset of differential expression (DE)
#' results. The DE results are subset by the "split" column, and each split
#' produces its own volcano plot using the \code{EnhancedVolcano} package.
#'
#' @param input Output from function \code{runDE} or a list containing
#'  (at minimum) a dataframe named "DE_result" including columns:
#'   \describe{
#'     \item{gene}{Gene identifiers (character).}
#'     \item{lfc}{Log2 fold change values (numeric).}
#'     \item{padj}{Adjusted p-values (numeric).}
#'     \item{split}{(Optionally) Grouping variable used to subset results
#'     (character or factor).}
#'   }
#' @param alpha A numeric value indicating the significance level used for
#' permutation test comparisons of the number of differentially expressed
#' features. Defaults to 0.05.
#' @param lfc_threshold A numeric value indicating the minimum absolute value
#' log fold change for a gene to be counted as a "hit". Defaults to 0.5. Set to
#' 0 to disregard log fold change when counting hits.
#' @param title Character string indicating the plot title. Default = `NULL`
#' sets a title automatically for each split.
#' @param subtitle Character string indicating the plot subtitle. Default =
#' `NULL` automatically generates a subtitle describing the significance
#' thresholds.
#' @param significant_color Character string color name or hex code indicating
#' the color used to highlight significant differentially expressed features
#' (default = `"red2"`).
#' @param ... Additional parameters passed to
#' \code{EnhancedVolcano::EnhancedVolcano()}.
#'
#' @return A named list of \code{ggplot2} objects, where each element
#'   corresponds to a volcano plot for one split of the DE results.
#'
#' @details
#' Each volcano plot is created using \code{EnhancedVolcano::EnhancedVolcano()},
#' customized with:
#' \itemize{
#'   \item The x-axis showing log2 fold changes (`lfc`).
#'   \item The y-axis showing adjusted p-values (`padj`).
#'   \item Genes colored grey, blue, or red depending on significance status.
#'   \item A caption showing the total number of genes in that split.
#'   \item \code{theme_permuteDE()} applied, with the legend removed.
#' }
#'
#' @examples
#' \dontrun{
#' plots <- generateVolcanoPlot(runDE_output)
#' print(plots[[1]])
#' }
#'
#' @export
#'
getVolcanos <- function(input,
                        alpha = 0.05,
                        lfc_threshold = 0.5,
                        title = NULL,
                        subtitle = NULL,
                        significant_color = 'red2',
                        ...){

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(input, "input")
  .validInput(alpha, "alpha")
  .validInput(lfc_threshold, "lfc_threshold")
  .validInput(title, "title")
  .validInput(subtitle, "subtitle")
  .validInput(significant_color, "significant_color") #

  # ---------------------------------------------------------------------------
  # Generate plots
  # ---------------------------------------------------------------------------

  # If available, grab group names and reference group
  # TODO

  # Separate results from each split
  split_results <- split(runDE_output$DE_result, runDE_output$DE_result$split)

  volcano_list <- lapply(names(split_results), function(s) {
    df <- split_results[[s]]

    EnhancedVolcano::EnhancedVolcano(
      df,
      lab = df$gene,
      x = 'lfc',
      y = 'padj',
      pCutoff = alpha,
      FCcutoff = lfc_threshold,
      title = if (is.null(title)) paste('Volcano Plot of Gene Expression in', s) else title,
      subtitle = if (is.null(subtitle))
        paste("Significance: adjusted p-value <", alpha, "and |log2 fold change| >", lfc_threshold)
      else subtitle,
      caption = paste0("total = ", nrow(df), " genes"),
      col = c("grey30", "grey30", "royalblue", significant_color)
    ) +
      theme_permuteDE() +
      ggplot2::theme(legend.position = 'none')
  })

  names(volcano_list) <- names(split_results)

  return(volcano_list)
}


#' Generate histogram plots of permutation test results
#'
#' This function creates a list of `ggplot2` histogram plots showing the distribution
#' of the number of differentially expressed (DE) genes across permutations for each
#' split in the permutation DE results. Each histogram includes a vertical line at
#' the observed number of DE genes and a p-value annotation.
#'
#' @param permuteDE_result A returned list from `permuteDE`
#' @param title A string to be the title of histograms. Default to \code{NULL}.
#' @export

getHistograms <- function(permuteDE_result, title = NULL) {
  permutation_DE_results <- permuteDE_result$permutation_DE_results
  permutation_test_results <- permuteDE_result$permutation_test_results
  all_splits <- unique(permutation_DE_results$split)
  n_iterations <- permuteDE_result$parameters$n_iterations

  # generate a list of histogram(s) per split
  histogram_list <- lapply(all_splits, function(split_id) {
    df_perm <- dplyr::filter(permutation_DE_results, split == split_id)
    df_test <- dplyr::filter(permutation_test_results, split == split_id)

    x_max <- max(df_perm$n_sig, df_test$true_n_sig)
    label_x <- min(df_test$true_n_sig + 5, x_max * 0.95)


    .plotHistogram(df_perm$n_sig, xlabel = 'Number of DE genes',
                    histAlpha = 0.5, vline = df_test$true_n_sig) +
      ggplot2::geom_vline(xintercept = df_test$true_n_sig,
                 color = 'red',
                 linewidth = 1,
                 linetype = 'solid') +
      ggplot2::annotate(
        "label",
        x = label_x,
        y = Inf,
        label = ifelse(df_test$p_n_sig < 0.01,
                       "p < 0.01",
                       paste0("p = ", round(df_test$p_n_sig, 2))),
        vjust = 1.5,
        hjust = 0.44,
        fill = "white",
        color = 'red',
        fontface = "bold",
        label.size = 0.3
      ) +
      ggplot2::labs(
        title = if (is.null(title)) {
          stringr::str_wrap(
            paste0(
              'Distribution of # DE genes across ',
              n_iterations, ' Iterations in ', split_id
            ),
            width = 40
          )
        } else {
          title
        }
      )
  })
  names(histogram_list) <- all_splits
  return(histogram_list)
}

# Plot custom histogram ---------------------------
#
# x -- Numeric vector of values to plot.
# xlabel -- Label for the x-axis.
# ylabel -- Label for the y-axis.
# addDensity -- Logical; if `TRUE`, adds a density curve (default: `FALSE`).
# bins -- Number of bins in the histogram (default: 20).
# baseSize -- Base font size used in `themeCorcesRegular()` (default: 14).
# histFill -- Fill color of histogram bars (default: "#3B9AB2").
# histColor -- Outline color of histogram bars (default: `NA`).
# densityFill -- Fill color of the density area (default: "#3B9AB2").
# densityColor -- Line color of the density curve (default: "#3B9AB2").
# size -- Line width for histogram and density elements (default: 1.25).
# histAlpha -- Transparency (alpha) of the histogram fill (default: 0.85).
# densityAlpha -- Transparency (alpha) of the density fill (default: 0.15).
# title -- Title for the plot (default: empty string).
# ratioYX -- Desired aspect ratio of y over x for `coord_equal()` (default: 0.8).
# vline -- Numeric value (or vector) for vertical reference lines (default: `NULL`).
# vlineColor -- Color for vertical lines (default: "red").
# vlineType -- Line type for vertical lines (e.g., "solid", "dashed") (default: "solid").
# vlineWidth -- Line width for vertical lines (default: 1).

.plotHistogram <- function(x,
                           xlabel = "values",
                           ylabel = "Count",
                           addDensity = FALSE,
                           bins = 20,
                           baseSize = 14,
                           histFill = "#3B9AB2",
                           histColor = NA,
                           densityFill = "#3B9AB2",
                           densityColor = "#3B9AB2",
                           title = "",
                           size = 1.25,
                           histAlpha = 0.85,
                           densityAlpha = 0.15,
                           ratioYX = 0.8,
                           vline = NULL,
                           vlineColor = "red",
                           vlineType = "solid",
                           vlineWidth = 1) {

  stopifnot(is.numeric(x))

  df <- data.frame(y = x)

  p <- ggplot2::ggplot(df, ggplot2::aes(y)) +
    ggplot2::geom_histogram(bins = bins, fill = histFill, color = histColor, size = size, alpha = histAlpha)

  if (addDensity) {
    p <- p +
      ggplot2::geom_density(ggplot2::aes(y = ..count..), fill = densityFill, color = NA, alpha = densityAlpha, size = size) +
      ggplot2::stat_density(ggplot2::aes(y = ..count..), color = densityColor, geom = "line", size = size)
  }

  if (!is.null(vline)) {
    p <- p + ggplot2::geom_vline(xintercept = vline, color = vlineColor, linewidth = vlineWidth, linetype = vlineType)
  }

  p <- p +
    theme_permuteDE(base_size = baseSize) +
    ggplot2::xlab(xlabel) + ggplot2::ylab(ylabel)

  # Include vline in xlim computation
  xlim_raw <- range(c(df$y, vline), na.rm = TRUE)
  xlim <- grDevices::extendrange(xlim_raw, f = 0.1)

  data <- ggplot2::ggplot_build(p)$data[[1]]
  ylim <- c(0, max(data$y, na.rm = TRUE) * 1.1)
  ratioXY <- ratioYX * diff(xlim) / diff(ylim)

  p <- p + ggplot2::coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE)

  return(p)
}
