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
#' `NULL` automatically generates a subtitle describing the DE method used.
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

  # If available, grab group names
  if ("parameters" %in% names(input)) {
    reference_group <- input$parameters[["reference_group"]]
    non_reference_group <- input$parameters[["non_reference_group"]]
  } else {
    reference_group <- "reference group"
    non_reference_group <- "non-reference group"
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
    if (is.null(title)) {
      current_title <- paste0(non_reference_group, " *vs.* ", reference_group,": ", s)
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
#' This function takes the output from function \code{permuteDE} and creates a
#' list containing a histogram plot showing the permutation test result for
#' each split. Each histogram includes a vertical line at the observed number
#' of DE features in the unpermuted "true" comparison.
#'
#' @param input Output from function \code{permuteDE}.
#' @param title Character string indicating the plot title. Default = `NULL`
#' sets a title automatically for each split.
#' @param subtitle Character string indicating the plot subtitle. Default =
#' `NULL` automatically generates a subtitle.
#' @param label_pvalue A Boolean value indicating whether to label the
#' permutation test p-value. Defaults to TRUE.
#'
#' @return A named list of \code{ggplot2} objects, where each element
#'   corresponds to a volcano plot for one split of the DE results.
#' @export

getHistograms <- function(input,
                          title = NULL,
                          subtitle = NULL,
                          label_pvalue = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .requirePackage("ggtext", source = "cran")

  .validInput(input, "input", "getHistograms")
  .validInput(title, "title")
  .validInput(subtitle, "subtitle")
  .validInput(label_pvalue, "label_pvalue")

  # ---------------------------------------------------------------------------
  # Set up
  # ---------------------------------------------------------------------------

  # If available, grab group names
  if ("parameters" %in% names(input)) {
    reference_group <- input$parameters[["reference_group"]]
    non_reference_group <- input$parameters[["non_reference_group"]]
  } else {
    reference_group <- "reference group"
    non_reference_group <- "non-reference group"
  }

  # Separate results from each split
  split_results <- split(input$permutation_DE_summary, input$permutation_DE_summary$split)

  # ---------------------------------------------------------------------------
  # Generate plots
  # ---------------------------------------------------------------------------

  histogram_list <- lapply(names(split_results), function(s) {
    # Current split
    split_results_s <- split_results[[s]]

    # Set title
    if (is.null(title)) {
      current_title <- paste0(non_reference_group, " *vs.* ", reference_group,": ", s)
      current_title <- paste(strwrap(current_title, 80), collapse = "\n")
    } else {
      current_title <- title
    }

    # Set subtitle
    if (is.null(subtitle)) {
      current_subtitle <- paste0("Permutation test with ", (nrow(split_results_s) + 1), " iterations")
      current_subtitle <- paste(strwrap(current_subtitle, 80), collapse = "\n")
    } else {
      current_subtitle <- subtitle
    }

    # Set color groups & label
    true_n_sig_s <- dplyr::filter(input$permutation_test_results, split == s)$true_n_sig[1]
    pvalue_s <- dplyr::filter(input$permutation_test_results, split == s)$pvalue[1]
    if (pvalue_s < 0.0001) {
      pvalue_s <- paste0("*p* < 0.0001")
    } else {
      pvalue_s <- paste0("*p* = ", round(pvalue_s, 4))
    }
    split_results_s <- split_results_s |>
      dplyr::mutate(fill_group = ifelse(n_sig >= true_n_sig_s, TRUE, FALSE))
    # Set limits
    y_limits <- c(0, (max(table(split_results_s$n_sig))*1.05))
    x_max <- max(c(split_results_s$n_sig, true_n_sig_s)) + 1
    x_breaks <- floor(pretty(seq(0, x_max, 1)))
    if (true_n_sig_s/x_max < 0.5) {
      x_label_position <- true_n_sig_s + x_max*0.025
      x_label_hjust <- "left"
    } else {
      x_label_position <- true_n_sig_s - x_max*0.025
      x_label_hjust <- "right"
    }

    # Plot
    p <- ggplot2::ggplot(data = split_results_s,
                         ggplot2::aes(x = n_sig,
                                      fill = fill_group)) +
      theme_permuteDE() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::theme(plot.title = ggtext::element_markdown()) +
      ggplot2::geom_histogram(alpha = 0.6,
                              binwidth = 1,
                              breaks = seq(-0.01, x_max, 1)) +
      ggplot2::geom_vline(xintercept = true_n_sig_s, linetype = "longdash", color = "#EE3751") +
      ggplot2::scale_fill_manual(values = c("#AAAAAA", "#EE3751")) +
      ggplot2::scale_y_continuous(limits = y_limits, expand = c(0,0)) +
      ggplot2::scale_x_continuous(breaks = x_breaks) +
      ggplot2::labs(title = current_title,
                    subtitle = current_subtitle,
                    x = "Number of DE features",
                    y = "Count")
    # Add p-value label
    if (label_pvalue == TRUE) {
      p + ggtext::geom_richtext(x = x_label_position,
                                y = y_limits[2]*0.975,
                                label = pvalue_s,
                                hjust = x_label_hjust,
                                vjust = "top",
                                color = "#EE3751",
                                fill = "white",
                                label.size = 0.4)
    } else {
      p
    }
  })

  names(histogram_list) <- names(split_results)
  return(histogram_list)
}

