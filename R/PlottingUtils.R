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
#' produces its own volcano plot.
#'
#' @param input Output from function \code{runDE} or a list containing
#'  (at minimum) a dataframe named "DE_results" including columns:
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
#' @param n_max_label A numeric value indicating how many of the top
#' significant DE features to label. Defaults to 10.
#' @param center A Boolean value indicating whether to center the x-axis at 0.
#' Defaults to \code{TRUE}.
#'
#' @return A named list of \code{ggplot2} objects, where each element
#'   corresponds to a volcano plot for one split of the DE results.
#'
#' @export
#'
getVolcanos <- function(input,
                        alpha = 0.05,
                        lfc_threshold = 0.5,
                        title = NULL,
                        subtitle = NULL,
                        n_max_label = 10,
                        center = TRUE){

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .requirePackage("ggrepel", source = "cran")
  .requirePackage("ggtext", source = "cran")

  .validInput(input, "input", "getVolcanos")
  .validInput(alpha, "alpha")
  .validInput(lfc_threshold, "lfc_threshold")
  .validInput(title, "title")
  .validInput(subtitle, "subtitle")
  .validInput(n_max_label, "n_max_label")
  .validInput(center, "center")

  # ---------------------------------------------------------------------------
  # Set up
  # ---------------------------------------------------------------------------

  # If available, grab group names and reference group
  reference_group <- NULL
  non_reference_group <- "non-reference group"
  if ("metadata" %in% names(input)) {
    groups <- unique(input$metadata$group_key$group)
    if ("parameters" %in% names(input)) {
      reference_group <- input$parameters[["reference_group"]]
      non_reference_group <- groups[groups != reference_group]
    } else {
      reference_group <- sort(groups)[1]
      non_reference_group <- sort(groups)[2]
    }
  }

  # If available, grab test info
  de_method <- " "
  de_test <- " "
  p_adjust_method <- " "
  if ("parameters" %in% names(input)) {
    de_method <- paste0(" ", input$parameters[["de_method"]], " ")
    de_test <- input$parameters[["de_test"]]
    p_adjust_method <- input$parameters[["p_adjust_method"]]
  }

  # Set subtitle
  if (is.null(subtitle)) {
    subtitle <- paste0("DE using", de_method,
                       switch(de_test,
                              wilcox_cpm = "Wilcoxon Rank Sum Test with CPM normalization, ",
                              wilcox_log_cpm = "Wilcoxon Rank Sum Test with log CPM normalization, ",
                              paste0(de_test, ", ")),
                       switch(p_adjust_method,
                              holm =  paste0("Holm method, \u03b1=", alpha),
                              hochberg = paste0("Hochberg adjustment, \u03b1=", alpha),
                              hommel = paste0("Hommel procedure, \u03b1=", alpha),
                              bonferroni = paste0("Bonferroni method, \u03b1=", alpha),
                              BH = paste0("FDR=", alpha),
                              fdr = paste0("FDR=", alpha),
                              BY = paste0("Benjamini & Yekutieli, FDR=", alpha),
                              paste0("no multiple comparison correction, \u03b1=", alpha)),
                       ", |LFC|>", lfc_threshold)
    subtitle <- paste(strwrap(subtitle, 80), collapse = "\n")
  }

  # Separate results from each split
  split_results <- split(input$DE_results, input$DE_results$split)

  # ---------------------------------------------------------------------------
  # Generate plots
  # ---------------------------------------------------------------------------

  volcano_list <- lapply(names(split_results), function(s) {
    # Current split
    split_results_s <- split_results[[s]]

    # Set title
    if (is.null(title) & !is.null(reference_group)) {
      current_title <- paste0(non_reference_group, " *vs.* ", reference_group,": ", s)
      current_title <- paste(strwrap(current_title, 80), collapse = "\n")
    } else if (is.null(title)) {
      current_title <- paste0("Differential expression: ", s)
      current_title <- paste(strwrap(current_title, 80), collapse = "\n")
    } else {
      current_title <- title
    }

    # Set limits
    if (center == FALSE) {
      x_limits <- c(min(split_results_s$lfc, na.rm = TRUE)*1.1, max(split_results_s$lfc, na.rm = TRUE)*1.1)
    } else {
      x_limits <- c(max(abs(split_results_s$lfc), na.rm = TRUE)*(-1.1), max(abs(split_results_s$lfc), na.rm = TRUE)*1.1)
    }
    y_limits <- c(0, max(-log10(split_results_s$padj), na.rm = TRUE)*1.1)

    # Set color groups & label set
    split_results_s <- split_results_s |>
      dplyr::mutate(sig_group = ifelse(lfc > lfc_threshold & padj < alpha, paste0("Higher in ", non_reference_group),
                                       ifelse(lfc < lfc_threshold*(-1) & padj < alpha, paste0("Lower in ", non_reference_group),
                                              "Not significant")))
    label_features <- dplyr::arrange(dplyr::filter(split_results_s, padj < 0.05, abs(lfc) > lfc_threshold),
                                     padj, -abs(lfc))$gene[1:min(n_max_label, nrow(dplyr::filter(split_results_s, padj < 0.05)))]
    label_split_results_s <- split_results_s |>
      dplyr::filter(gene %in% label_features)

    # Plot
    ggplot2::ggplot(data = split_results_s,
                    ggplot2::aes(x = lfc,
                                 y = -log10(padj),
                                 color = sig_group)) +
      theme_permuteDE() +
      ggplot2::theme(plot.title = ggtext::element_markdown(),
                     axis.title.x = ggtext::element_markdown(),
                     axis.title.y = ggtext::element_markdown()) +
      ggplot2::geom_point(alpha = 0.5,
                          size = 2) +
      ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                          linetype = "dashed") +
      ggplot2::geom_hline(yintercept = -log10(alpha),
                          linetype = "dashed") +
      ggrepel::geom_text_repel(data = label_split_results_s,
                               ggplot2::aes(label = gene),
                               color = "black",
                               max.overlaps = Inf) +
      ggplot2::xlim(x_limits) +
      ggplot2::ylim(y_limits) +
      ggplot2::scale_color_manual(values = c("#3A4BED", "#EE3751", "#BBBBBB"),
                                  breaks = c(paste0("Lower in ", non_reference_group),
                                             paste0("Higher in ", non_reference_group))) +
      ggplot2::labs(title = current_title,
                    subtitle = subtitle,
                    color = "",
                    x = "Log<sub>2</sub> Fold Change",
                    y = "-Log<sub>10</sub> Adjusted P Value")
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
#' @param input A returned list from `permuteDE`
#' @param title A string to be the title of histograms. Default to \code{NULL}.
#' @export

getHistograms <- function(input,
                          title = NULL) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(input, "input", "getHistograms")
  .validInput(title, "title")

  # ---------------------------------------------------------------------------
  # Generate plots
  # ---------------------------------------------------------------------------

  permutation_DE_summary <- input$permutation_DE_summary
  permutation_test_results <- input$permutation_test_results
  all_splits <- unique(permutation_DE_summary$split)
  n_iterations <- input$parameters$n_iterations

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
        label = ifelse(df_test$pvalue_n_sig < 0.01,
                       "p < 0.01",
                       paste0("p = ", round(df_test$pvalue_n_sig, 2))),
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
