# ---------------------------------------------------------------------------
# Plotting-related functions
# ---------------------------------------------------------------------------

#' permuteDE ggplot theme
#'
#' This function returns the ggplot2 theme used throughout the permuteDE
#' plotting functions.
#'
#' @param color Color of text and lines of the plot
#' @param base_size The size of the font for the axis labels
#' @param base_line_size The line width for most lines
#' @param base_rect_size The line width for rectangular boxes
#' @param axis_title_size The font size of the axis title
#' @param plot_title_size The font size of the plot title
#' @param plot_margin_cm The margin around the plot in centimeters
#' @param legend_text_size 0.75*base_size
#' @param legend_position The placement of the legend in the plot
#' @param legend_title_hjust Left-right justification of legend title
#' @param legend_title_vjust Top-bottom justification of legend title
#' @param axis_tick_length_mm Axis tick length in mm
#' @param rotate_x_axis_text_90 Boolean value indicating whether to rotate the x-axis text by 90 degrees
#' @param rotate_y_axis_text_90 Boolean value indicating whether to rotate the y-axis text by 90 degrees
#'
#' @export
permuteDEtheme <- function(color = "black",
                           base_size = 10,
                           base_line_size = 0.5,
                           base_rect_size = 0.5,
                           axis_title_size = 12,
                           plot_title_size = 14,
                           plot_margin_cm = 1,
                           legend_text_size = 10,
                           legend_position = "bottom",
                           legend_title_hjust = 0,
                           legend_title_vjust = 0.5,
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
    legend.title = ggplot2::element_text(color = color, size = legend_text_size, hjust = legend_title_hjust, vjust = legend_title_vjust),
    legend.box.background = ggplot2::element_rect(fill = "transparent"),
    legend.position = legend_position,
    legend.spacing = ggplot2::unit(15, "pt"),
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

##' Generate color palette
#'
#' Generate a color palette. For discrete palettes, hex values are standard up
#' to n = 100, but for larger values of n, they are generated using
#' \code{Polychrome::createPalette()}.
#'
#' The "choir" palette is adapted from function  \code{CHOIR::CHOIRpalette}
#' in R package \code{CHOIR} (Sant et al. 2025). The "archr" palette is
#' adapted from palette "stallion" in R package \code{ArchR} (Granja & Corces
#' et al. 2020).
#'
#' @param type A character string indicating the palette type. Permitted values
#' are "discrete" and "gradient". Defaults to discrete.
#' @param n Number of colors. Default = \code{NULL} will return all of the
#' pre-set colors in the palette.
#' @param palette_name A character string indicating the palette name. Permitted
#' values are "choir", "archr", "frozen", and "inferno". Default = \code{NULL}
#' will use "choir" when \code{type} is "discrete" and "inferno" when
#' \code{type} is "gradient.
#' @param swatch A Boolean value indicating whether to plot a swatch of the
#' palette.
#'
#' @return Returns a vector of n hex values.
#'
#' @export
#'
permuteDEpalette <- function(type = "discrete",
                             n = NULL,
                             palette_name = NULL,
                             swatch = FALSE) {

  # ---------------------------------------------------------------------------
  # Check parameter input validity
  # ---------------------------------------------------------------------------

  .validInput(type, "type")
  .validInput(n, "n")
  .validInput(palette_name, "palette_name", type)
  .validInput(swatch, "swatch")

  # ---------------------------------------------------------------------------
  # Create palette
  # ---------------------------------------------------------------------------
  if (type == "discrete") {
    if (is.null(palette_name)) {
      palette_name <- "choir"
    }
    if (palette_name == "choir") {
      starting_colors <- c("#00CCE3", "#F8A100", "#E81AEF", "#F56900", "#6560FF",
                           "#00D456", "#A25AFF", "#C1D400", "#E58CCC", "#1990FF",
                           "#00DEA3", "#FF5B4B", "#F7D823", "#3BA833", "#AEA9FF",
                           "#FF9C88", "#7DA6CC", "#CB8251", "#61BBFF", "#FF4CF9",
                           "#8B7CFF", "#A99900", "#FCB467", "#FB798C", "#83DD00",
                           "#D98AF7", "#67D4D9", "#08A38E", "#8099FF", "#24B167",
                           "#DBAA00", "#FF3FAD", "#D2C0CA", "#888DB2", "#8D9D63",
                           "#E934FF", "#33C100", "#E0BFA1", "#21ADFD", "#86D582",
                           "#FFA2D9", "#FFB703", "#9DCCC6", "#C078C1", "#F88800",
                           "#9FD362", "#A876EB", "#B4908B", "#FF6C69", "#6693DE",
                           "#E2B0FD", "#AD8EA8", "#93AC00", "#E37847", "#E368FB",
                           "#ECB5D8", "#789F33", "#3AA2AA", "#F3B19C", "#CCC482",
                           "#34DBC9", "#B486E1", "#C39131", "#FF7D4A", "#7F9F87",
                           "#FF9BF3", "#0899D0", "#CCB0DB", "#7DD6B9", "#EABD7F",
                           "#FC76AE", "#D97C9E", "#75C7E4", "#5EA157", "#A5D300",
                           "#FF00FF", "#6A99C5", "#66D100", "#979177", "#E08500",
                           "#EA53C2", "#BAC0FA", "#C3C4AE", "#9785C6", "#6BAE00",
                           "#A78D60", "#EA4974", "#3AEF83", "#D168CC", "#B77E7B",
                           "#77C65D", "#E2D452", "#FFB6C9", "#A2B3CC", "#C660A2",
                           "#616DFF", "#FF9240", "#9B9CA3", "#7DE3FF", "#FF69A3")
    } else if (palette_name == "archr") {
      starting_colors <- c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B",
                           "#FEE500","#8A9FD1","#C06CAB","#D8A767","#90D5E4",
                           "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
                           "#6E4B9E","#0C727C","#7E1416","#E6C2DC","#3D3D3D")
    }
    if (is.null(n)) {
      n <- length(starting_colors)
    }
    if (n <= length(starting_colors)) {
      values <- starting_colors[1:n]
    } else {
      .requirePackage("Polychrome", source = "cran")
      values <- Polychrome::createPalette(N = n,
                                          seedcolors = starting_colors,
                                          range = c(50, 80))
      names(values) <- NULL
    }
  } else if (type == "gradient") {
    if (is.null(palette_name)) {
      palette_name <- "inferno"
    }
    if (palette_name == "inferno") {
      starting_colors <- c("#FFC715", "#FF9D33", "#FF7145", "#FF527B", "#ED35B9", "#CB27E2", "#9031FF", "#5939F7", "#2A26EA", "#1A17BA", "#0C0782", "#000A51")
    } else if (palette_name == "frozen") {
      starting_colors <- c("#61DAFF", "#48ACFF", "#337DFF", "#2E5CEF", "#5F36F2", "#511DCE", "#32217A", "#1E1551")
    }
    if (is.null(n)) {
      n <- length(starting_colors)
    }
    if (n != length(starting_colors)) {
      color_function <- grDevices::colorRampPalette(colors = starting_colors)
      values <- color_function(n)
    } else {
      values <- starting_colors
    }
  }

  # ---------------------------------------------------------------------------
  # Plot swatch
  # ---------------------------------------------------------------------------
  if (swatch == TRUE) {
    n_rows <- ceiling(n / 10)
    swatch_data <- data.frame(x = rep(1:10, n_rows)[1:n],
                              y = rep(n_rows:1, each = 10)[1:n],
                              color = values)
    p <- ggplot2::ggplot(swatch_data, ggplot2::aes(x, y)) +
      ggplot2::geom_tile(ggplot2::aes(fill = color), color = NA) +
      ggplot2::scale_fill_identity() +
      ggplot2::coord_fixed() +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.05, vjust = 1, size = 10)) +
      ggplot2::labs(title = paste0('Palette: "', palette_name, '"'))
    print(p)
  }
  return(values)
}

#' Generate volcano plot(s) from differential expression results
#'
#' This function takes the output from function \code{runDE} and creates volcano
#' plots for each indicated split.
#'
#' @param input Output from function \code{runDE} or a list containing
#'  (at minimum) a dataframe named "DE_results" including columns:
#'   \describe{
#'     \item{feature}{Feature identifiers (character).}
#'     \item{lfc}{Log2 fold change values (numeric).}
#'     \item{padj}{Adjusted p-values (numeric).}
#'     \item{split}{Grouping variable used to subset results
#'     (character or factor).}
#'   }
#' @param alpha A numeric value indicating the significance level used for
#' permutation test comparisons of the number of differentially expressed
#' features. Defaults to 0.05.
#' @param lfc_threshold A numeric value indicating the minimum absolute value
#' log fold change for a feature to be counted as a "hit". Defaults to 0.5. Set
#' to 0 to disregard log fold change when counting hits.
#' @param use_splits A character string or vector containing the names of splits
#' to use. Defaults to \code{NULL}, which will try all splits.
#' @param title Character string indicating the plot title. Default =
#' \code{NULL} sets a title automatically for each split.
#' @param subtitle Character string indicating the plot subtitle. Default =
#' \code{NULL} automatically generates a subtitle describing the DE method used.
#' @param n_max_label A numeric value indicating how many of the top
#' significant DE features to label. Defaults to 10.
#' @param center A Boolean value indicating whether to center the x-axis at 0.
#' Defaults to \code{TRUE}.
#'
#' @return If only one split is provided, a single \code{ggplot2} object,
#' otherwise a named list of \code{ggplot2} objects, where each element
#' corresponds to a split.
#'
#' @export
#'
plotVolcano <- function(input,
                        alpha = 0.05,
                        lfc_threshold = 0.5,
                        use_splits = NULL,
                        title = NULL,
                        subtitle = NULL,
                        n_max_label = 10,
                        center = TRUE){

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .requirePackage("ggrepel", source = "cran")
  .requirePackage("ggtext", source = "cran")

  .validInput(input, "input", "plotVolcano")
  .validInput(alpha, "alpha")
  .validInput(lfc_threshold, "lfc_threshold")
  .validInput(use_splits, "use_splits", list(input, "plotVolcano"))
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
                              wilcox_cpm = "Wilcoxon Rank Sum Test with CPM normalization\n",
                              wilcox_log_cpm = "Wilcoxon Rank Sum Test with log CPM normalization\n",
                              paste0(de_test, "\n")),
                       switch("fdr",
                              holm =  paste0("Holm method, \u03b1 = ", alpha),
                              hochberg = paste0("Hochberg adjustment, \u03b1 = ", alpha),
                              hommel = paste0("Hommel procedure, \u03b1 = ", alpha),
                              bonferroni = paste0("Bonferroni method, \u03b1 = ", alpha),
                              BH = paste0("FDR = ", alpha),
                              fdr = paste0("FDR = ", alpha),
                              BY = paste0("Benjamini & Yekutieli, FDR = ", alpha),
                              paste0("No multiple comparison correction, \u03b1 = ", alpha)),
                       ", |LFC| > ", lfc_threshold)
    if (nchar(subtitle) < 50) {
      subtitle <- gsub("\n", ", ", subtitle)
    }
  }

  # Separate results from each split
  split_results <- split(input$DE_results, input$DE_results$split)

  # Subset data if split(s) are provided
  if (!is.null(use_splits)) {
    split_results <- split_results[use_splits]
  }

  # ---------------------------------------------------------------------------
  # Generate plots
  # ---------------------------------------------------------------------------

  volcano_list <- lapply(names(split_results), function(s) {
    # Current split
    split_results_s <- split_results[[s]]

    # Set title
    if (is.null(title)) {
      current_title <- paste0(non_reference_group, " *vs.* ", reference_group,"<br>", s)
      if (nchar(current_title) < 50) {
        subtitle <- gsub("<br>", ": ", subtitle)
      }
    } else {
      current_title <- title
    }

    # Set limits
    if (center == FALSE) {
      x_limits <- c(min(split_results_s$lfc, na.rm = TRUE)*1.1, max(split_results_s$lfc, na.rm = TRUE)*1.1)
    } else {
      x_limits <- c(max(abs(split_results_s$lfc), na.rm = TRUE)*(-1.1), max(abs(split_results_s$lfc), na.rm = TRUE)*1.1)
    }
    y_limits <- c(0, max(c(-log10(split_results_s$padj), -log10(alpha)), na.rm = TRUE)*1.1)

    # Set color groups & label set
    split_results_s <- split_results_s |>
      dplyr::mutate(sig_group = ifelse(lfc > lfc_threshold & padj < alpha, paste0("Higher in ", non_reference_group),
                                       ifelse(lfc < lfc_threshold*(-1) & padj < alpha, paste0("Higher in ", reference_group),
                                              "Not significant")))
    label_features <- dplyr::arrange(dplyr::filter(split_results_s, padj < 0.05, abs(lfc) > lfc_threshold),
                                     padj, -abs(lfc))$feature[1:min(n_max_label, nrow(dplyr::filter(split_results_s, padj < 0.05)))]
    label_split_results_s <- split_results_s |>
      dplyr::filter(feature %in% label_features)

    # Plot
    ggplot2::ggplot(data = split_results_s,
                    ggplot2::aes(x = lfc,
                                 y = -log10(padj),
                                 color = sig_group)) +
      permuteDEtheme() +
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
                               ggplot2::aes(label = feature),
                               color = "black",
                               max.overlaps = Inf) +
      ggplot2::xlim(x_limits) +
      ggplot2::ylim(y_limits) +
      ggplot2::scale_color_manual(values = c("#3A4BED", "#EE3751", "#BBBBBB"),
                                  breaks = c(paste0("Higher in ", reference_group),
                                             paste0("Higher in ", non_reference_group))) +
      ggplot2::labs(title = current_title,
                    subtitle = subtitle,
                    color = "",
                    x = "Log<sub>2</sub> Fold Change",
                    y = "-Log<sub>10</sub> Adjusted P Value")
  })

  names(volcano_list) <- names(split_results)
  if (length(volcano_list) == 1) {
    return(volcano_list[[1]])
  } else {
    return(volcano_list)
  }
}


#' Generate histogram plot(s) of permutation test results
#'
#' This function takes the output from function \code{permuteDE} and creates a
#' histogram plot showing the permutation test result for each indicated split.
#' Each histogram includes a vertical line at the observed number of DE features
#' in the unpermuted "true" comparison.
#'
#' @param input Output from function \code{permuteDE}.
#' @param use_splits A character string or vector containing the names of splits
#' to use. Defaults to \code{NULL}, which will try all splits.
#' @param title Character string indicating the plot title. Default =
#' \code{NULL} sets a title automatically for each split.
#' @param subtitle Character string indicating the plot subtitle. Default =
#' \code{NULL} automatically generates a subtitle.
#' @param label_pvalue A Boolean value indicating whether to label the
#' permutation test p-value. Defaults to TRUE.
#'
#' @return If only one split is provided, a single \code{ggplot2} object,
#' otherwise a named list of \code{ggplot2} objects, where each element
#' corresponds to a split.
#'
#' @export

plotHistogram <- function(input,
                          use_splits = NULL,
                          title = NULL,
                          subtitle = NULL,
                          label_pvalue = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .requirePackage("ggtext", source = "cran")

  .validInput(input, "input", "plotHistogram")
  .validInput(use_splits, "use_splits", list(input, "plotHistogram"))
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

  # Subset data if split(s) are provided
  if (!is.null(use_splits)) {
    split_results <- split_results[use_splits]
  }

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
    runDE_n_sig_s <- dplyr::filter(input$permutation_test_results, split == s)$runDE_n_sig[1]
    pvalue_s <- dplyr::filter(input$permutation_test_results, split == s)$pvalue[1]
    if (pvalue_s == 1) {
      palette <- c("#EE3751")
    } else {
      palette <- c("#AAAAAA", "#EE3751")
    }
    if (pvalue_s < 0.0001) {
      pvalue_s <- paste0("*p* < 0.0001")
    } else {
      pvalue_s <- paste0("*p* = ", round(pvalue_s, 4))
    }
    split_results_s <- split_results_s |>
      dplyr::mutate(fill_group = ifelse(n_sig >= runDE_n_sig_s,
                                        "Permuted \u2265 true\ncomparison",
                                        "Permuted < true\ncomparison"))
    # Set limits
    y_limits <- c(0, (max(table(split_results_s$n_sig))*1.1))
    x_max <- max(c(split_results_s$n_sig, runDE_n_sig_s)) + 1
    x_breaks <- floor(pretty(seq(0, x_max, 1)))
    if (runDE_n_sig_s/x_max < 0.5) {
      x_label_position <- runDE_n_sig_s + x_max*0.025
      x_label_hjust <- "left"
    } else {
      x_label_position <- runDE_n_sig_s - x_max*0.025
      x_label_hjust <- "right"
    }

    # Plot
    p <- ggplot2::ggplot(data = split_results_s,
                         ggplot2::aes(x = n_sig,
                                      fill = fill_group)) +
      permuteDEtheme() +
      ggplot2::theme(plot.title = ggtext::element_markdown()) +
      ggplot2::geom_histogram(alpha = 0.6,
                              binwidth = 1,
                              breaks = seq(-0.01, x_max, 1)) +
      ggplot2::geom_vline(xintercept = runDE_n_sig_s, linetype = "longdash", color = "#EE3751") +
      ggplot2::scale_fill_manual(values = palette) +
      ggplot2::scale_y_continuous(limits = y_limits, expand = c(0,0)) +
      ggplot2::scale_x_continuous(breaks = x_breaks) +
      ggplot2::labs(title = current_title,
                    subtitle = current_subtitle,
                    x = "Number of DE features",
                    y = "Count",
                    fill = "")
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
  if (length(histogram_list) == 1) {
    return(histogram_list[[1]])
  } else {
    return(histogram_list)
  }
}


#' Plot levels of a single feature across groups
#'
#' This functions takes the output from function \code{runDE} and plots the
#' expression of a single feature across the replicates in each group for each
#' indicated split.
#'
#' @param input Output from function \code{runDE}.
#' @param feature A character string indicating which feature to plot.
#' @param use_splits A character string or vector containing the names of splits
#' to use. Defaults to \code{NULL}, which will try all splits.
#' @param normalization_method A character string indicating which normalization
#' method to apply to the replicate x feature matrix. Permitted values are
#' "cpm", "log_cpm", and "none". Defaults to "cpm".
#' @param title A character string indicating the plot title. Default =
#' \code{NULL} sets a title automatically for each split.
#' @param subtitle Character string indicating the plot subtitle. Default =
#' \code{NULL} automatically generates a subtitle.
#' @param label_replicates A Boolean value indicating whether to label the
#' replicates. Defaults to \code{FALSE}.
#' @param label_statistics A Boolean value indicating whether to label the
#' differential expression analysis LFC and p-value. Defaults to TRUE.
#'
#' @return If only one split is provided, a single \code{ggplot2} object,
#' otherwise a named list of \code{ggplot2} objects, where each element
#' corresponds to a split.
#'
#' @export
#'
plotFeature <- function(input,
                        feature,
                        use_splits = NULL,
                        normalization_method = "cpm",
                        title = NULL,
                        subtitle = NULL,
                        label_replicates = FALSE,
                        label_statistics = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .requirePackage("ggtext", source = "cran")
  .requirePackage("ggbeeswarm", source = "cran")
  if (label_replicates == TRUE) {
    .requirePackage("ggrepel", source = "cran")
  }

  .validInput(input, "input", "plotFeature")
  .validInput(feature, "feature")
  .validInput(use_splits, "use_splits", list(input, "plotFeature"))
  .validInput(normalization_method, "normalization_method")
  .validInput(title, "title")
  .validInput(subtitle, "subtitle")
  .validInput(label_replicates, "label_replicates")
  .validInput(label_statistics, "label_statistics")

  # ---------------------------------------------------------------------------
  # Set up
  # ---------------------------------------------------------------------------

  # Retrieve groups
  reference_group <- input$parameters[["reference_group"]]
  non_reference_group <- input$parameters[["non_reference_group"]]

  # Subset data if split(s) are provided
  if ("PB_values" %in% names(input)) {
    if (!is.null(use_splits)) {
      split_results <- input$PB_values[use_splits]
    } else {
      split_results <- input$PB_values
    }
  } else if ("cell_values" %in% names(input)) {
    if (!is.null(use_splits)) {
      split_results <- input$cell_values[use_splits]
    } else {
      split_results <- input$cell_values
    }
    warning("This plotting function may be less efficient for cell-level data.")
  } else {
    stop("Structure of list provided for parameter 'input' is unexpected, it should be the output returned by function 'runDE()' ",
         "or a list containing (at minimum) elements named 'PB_values' (or 'cell_values'), 'metadata', and 'parameters'. Please supply valid input!")
  }

  # ---------------------------------------------------------------------------
  # Generate plots
  # ---------------------------------------------------------------------------

  plot_list <- lapply(names(split_results), function(s) {
    # Current split
    split_results_s <- split_results[[s]]

    # Normalized data
    if (normalization_method == "cpm") {
      split_results_s <- edgeR::cpm(y = split_results_s)
      y_axis_title <- paste0("Normalized expression of *", feature, "* (CPM)")
    } else if (normalization_method == "log_cpm") {
      split_results_s <- edgeR::cpm(y = split_results_s,
                                    log = TRUE)
      y_axis_title <- paste0("Normalized expression of *", feature, "* (Log CPM)")
    } else {
      y_axis_title <- paste0("Expression of *", feature, "*")
    }

    # Set title & subtitle
    if (is.null(title)) {
      current_title <- paste0("*", feature, "*")
    } else {
      current_title <- title
    }
    if (is.null(subtitle)) {
      current_subtitle <- s
    } else {
      current_subtitle <- subtitle
    }

    # If feature is not in matrix, return NULL
    if (!(feature %in% rownames(split_results_s))) {
      warning("Feature '", feature, "' is not present in split '", s, "'.")
      return(NULL)
    } else {
      # Extract feature & proceed
      feature_s <- data.frame(feature = split_results_s[feature,],
                              replicate = input$metadata$group_key[colnames(split_results_s), "replicate"],
                              group = input$metadata$group_key[colnames(split_results_s), "group"])
      feature_s$group <- stats::relevel(factor(feature_s$group), ref = reference_group)

      # If labeling statistics, extract
      f <- feature
      feature_statistics <- dplyr::filter(input$DE_results,
                                          feature == f,
                                          split == s)
      if (nrow(feature_statistics) == 0) {
        label_statistics <- FALSE
      }
      if (label_statistics == TRUE) {
        if (!("DE_results" %in% names(input))) {
          stop("When parameter 'label_statistics' is set to TRUE, list provided to parameter 'input' must contain an elements named 'DE_results'. Please supply valid input!")
        }
        if (feature_statistics$padj[1] < 0.0001) {
          statistics_text <- paste0("LFC = ", round(feature_statistics$lfc[1], 4), ", *p* < 0.0001")
        } else {
          statistics_text <- paste0("LFC = ", round(feature_statistics$lfc[1], 4),
                                    ", *p* = ", round(feature_statistics$padj[1], 4))
        }
      }

      # Plot
      p <- ggplot2::ggplot(data = feature_s,
                           ggplot2::aes(x = group,
                                        y = feature,
                                        fill = group)) +
        permuteDEtheme() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::theme(plot.title = ggtext::element_markdown(),
                       axis.title.y = ggtext::element_markdown()) +
        ggbeeswarm::geom_beeswarm(size = 2, shape = 21) +
        ggplot2::scale_fill_manual(values = c("#AAAAAA", "#EE3751")) +
        ggplot2::labs(title = current_title,
                      subtitle = current_subtitle,
                      x = "Group",
                      y = y_axis_title)
      # Add replicate labels
      if (label_replicates == TRUE) {
        p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = replicate))
      }
      # Add LFC & p-value
      if (label_statistics == TRUE) {
        p <- p + ggtext::geom_richtext(x = 1.5,
                                       y = max(feature_s$feature)*1.05,
                                       label = statistics_text,
                                       hjust = "center",
                                       vjust = "bottom",
                                       fill = NA,
                                       label.color = NA) +
          ggplot2::geom_segment(x = 1, xend = 2,
                                y = max(feature_s$feature)*1.05, yend = max(feature_s$feature)*1.05) +
          ggplot2::ylim(c(min(feature_s$feature, na.rm = TRUE), max(feature_s$feature, na.rm = TRUE)*1.1))
      }
      return(p)
    }
  })
  names(plot_list) <- names(split_results)

  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
  }
}


#' Plot dimensionality reduction
#'
#' This function will generate a dimensionality reduction plot colored according
#' to the selected metric.
#'
#' @param reduction Dimensionality reduction cell coordinates. If there are
#' more than two columns, only the first two will be used.
#' @param input Output from function \code{permuteDE}.
#' @param split_labels An character vector containing the split labels for each
#' cell in order. Default = \code{NULL} will assume all cells belong to the same
#' split.
#' @param feature_values A vector containing the value for a feature for each
#' cell in order, used if \code{color_by} = "feature".
#' @param feature_name A string indicating the name of the plotted feature if
#' \code{color_by} = "feature". Used for legend title.
#' @param use_cells A vector of cell names to subset the dimensionality
#' reduction cell coordinates to. Default = \code{NULL} will use all cells.
#' @param color_by A character vector indicating what metric to use to color
#' each split. Permitted values are "n_sig" (number of significant DE features),
#' "pvalue" (permutation test p-value), "split" (values provided to parameter
#' \code{split_labels}), and "feature" (values provided to parameter
#' \code{feature_values}). Default = "split" will color the cells by the split
#' they belong to if input to parameter \code{split_labels} is provided.
#' @param permutation_test_alpha A numeric value indicating the significance
#' level to apply to the permutation test results. Splits that do not pass this
#' threshold will be grayed out. Default = 1 applies no threshold.
#' @param label_splits A Boolean value indicating whether to label the
#' splits. Defaults to \code{FALSE}.
#' @param label_statistics A Boolean value indicating whether to label the
#' values of the selected metric. Defaults to \code{FALSE}.
#' @param palette_name A character string indicating the palette name. Permitted
#' values are "choir", "archr", "frozen", and "inferno". Default =
#' \code{NULL} will use "choir" for discrete colors and "frozen" for
#' gradient colors.
#' @param color_limits A vector with the minimum and maximum values indicated by
#' `color_by` ("n_sig" or "pvalue") for display on the color bar legend. Default
#' = \code{NULL} sets limits automatically.
#' @param fix_coords A Boolean value indicating whether the aspect ratio of the
#' x and y axis should be 1.
#' @param ... Extra parameters passed to \code{Seurat::DimPlot()} or
#' \code{Seurat::FeaturePlot()}.
#'
#' @return A dimensionality reduction plot.
#' @export
#'
plotDimReduction <- function(reduction,
                             input = NULL,
                             split_labels = NULL,
                             feature_values = NULL,
                             feature_name = NULL,
                             use_cells = NULL,
                             color_by = NULL,
                             permutation_test_alpha = 1,
                             label_splits = FALSE,
                             label_statistics = FALSE,
                             palette_name = NULL,
                             color_limits = NULL,
                             fix_coords = TRUE,
                             ...) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .requirePackage("ggtext", source = "cran")

  .validInput(reduction, "reduction")
  .validInput(input, "input", "plotDimReduction")
  .validInput(split_labels, "split_labels", reduction)
  .validInput(use_cells, "use_cells", list(t(reduction), "none"))
  .validInput(color_by, "color_by", list(split_labels, input))
  .validInput(feature_values, "feature_values", list(reduction, use_cells, color_by))
  .validInput(feature_name, "feature_name", color_by)
  .validInput(permutation_test_alpha, "permutation_test_alpha", color_by)
  .validInput(label_splits, "label_splits")
  .validInput(label_statistics, "label_statistics")
  .validInput(color_limits, "color_limits", color_by)
  .validInput(fix_coords, "fix_coords")

  # ---------------------------------------------------------------------------
  # Set up
  # ---------------------------------------------------------------------------

  # Subset if necessary
  if(!is.null(use_cells)) {
    reduction <- reduction_[use_cells, ]
  }

  # Create temporary Seurat object
  tmp <- matrix(stats::rnorm(nrow(reduction) * 3, 10),
                ncol = nrow(reduction), nrow = 3)
  colnames(tmp) <- rownames(reduction)
  rownames(tmp) <- paste0("t",seq_len(nrow(tmp)))
  tmp_seurat <- Seurat::CreateSeuratObject(tmp, min.cells = 0, min.features = 0, assay = 'tmp')
  # Add dimensionality reduction
  tmp_seurat@reductions$dim_reduction <- suppressWarnings(Seurat::CreateDimReducObject(embeddings = reduction,
                                                                                       key = 'dimreduction',
                                                                                       assay = 'tmp'))
  # Add color groupings
  na_color <- "#BBBBBB"
  na_cutoff_low <- NA
  na_cutoff_high <- NA
  na_zero_legend <- FALSE
  add_labels <- FALSE
  if (is.null(color_by) || (is.null(split_labels) && identical(color_by, "split"))) {
    type <- "DimPlot"
    # No groups
    tmp_seurat$color_groups <- "all"
    palette <- c(na_color)
    color_legend <- ""
    if (label_splits == TRUE) {
      warning("Input value for 'label_splits' is not used when parameters 'split_labels' or 'color_by' are NULL.")
    }
    if (label_statistics == TRUE) {
      warning("Input value for 'label_statistics' is not used when parameter 'color_by' is not 'n_sig' or 'pvalue'.")
    }
  } else if (color_by == "split") {
    type <- "DimPlot"
    # Add splits to metadata
    tmp_seurat$color_groups <- split_labels
    palette <- permuteDEpalette(type = "discrete",
                                n = dplyr::n_distinct(split_labels),
                                palette_name = palette_name)
    color_legend <- ""
    if (label_splits == TRUE) {
      add_labels <- TRUE
    }
    if (label_statistics == TRUE) {
      warning("Input value for 'label_statistics' is not used when parameter 'color_by' is not 'n_sig' or 'pvalue'.")
    }
  } else if (color_by %in% c("n_sig", "pvalue")) {
    legend_barwidth <- NULL
    legend_barheight <- 1
    type <- "FeaturePlot"
    # Extract data
    key <- input$permutation_test_results[, c("split", "runDE_n_sig", "pvalue")]
    if (permutation_test_alpha < 1) {
      key$runDE_n_sig[key$pvalue >= permutation_test_alpha] <- NA
      key$pvalue[key$pvalue >= permutation_test_alpha] <- NA
      na_legend <- paste0("Permutation<br>test *p* \u2265 ", permutation_test_alpha)
    } else if (("metadata" %in% names(input) & ("runDE_values" %in% names(input$metadata)))) {
      # If not all splits were run using permuteDE, pull runDE_n_sig values
      if (!all(input$metadata$runDE_values$split %in% key$split)) {
        indices <- which(!(input$metadata$runDE_values$split %in% key$split))
        key <- rbind(key, data.frame(split = input$metadata$runDE_values$split[indices],
                                     runDE_n_sig = input$metadata$runDE_values$runDE_n_sig[indices],
                                     pvalue = NA))
      }
    }
    tmp_seurat$n_sig <- key$runDE_n_sig[match(split_labels, key$split)]
    tmp_seurat$pvalue <- key$pvalue[match(split_labels, key$split)]
    # Color palette
    palette <- permuteDEpalette(type = "gradient",
                                palette_name = palette_name)
    if (color_by == "n_sig") {
      color_legend <- "Number of<br>significant<br>DE features"
      if (all(is.na(tmp_seurat$n_sig))) {
        tmp_seurat$n_sig <- 0
      }
      tmp_seurat$color_groups <- tmp_seurat$n_sig
      # Set NA cutoff
      if (!is.null(color_limits)) {
        na_cutoff_low <- color_limits[1]
        na_cutoff_high <- color_limits[2]
        if (na_cutoff_low == 1) {
          na_zero_legend <- TRUE
        }
      } else if (min(tmp_seurat$n_sig, na.rm = TRUE) == 0) {
        na_cutoff_low <- 1
        na_zero_legend <- TRUE
      }
      # Set barwidth if only 1 value
      if ((dplyr::n_distinct(tmp_seurat$n_sig[!is.na(tmp_seurat$n_sig)]) < 2)) {
        if (is.null(color_limits) || (color_limits[1] == color_limits[2])) {
          legend_barwidth <- 1
        }
      }
    } else if (color_by == "pvalue") {
      color_legend <- "Permutation<br>test p-value"
      if (all(is.na(tmp_seurat$pvalue))) {
        tmp_seurat$n_sig <- 2
        na_cutoff_high <- 1
      }
      tmp_seurat$color_groups <- tmp_seurat$pvalue
      # Set NA cutoff
      if (!is.null(color_limits)) {
        na_cutoff_low <- max(0, color_limits[1])
        na_cutoff_high <- min(1, color_limits[2])
      }
      # Set barwidth if only 1 value
      if ((dplyr::n_distinct(tmp_seurat$pvalue[!is.na(tmp_seurat$pvalue)]) < 2)) {
        if (is.null(color_limits) || (color_limits[1] == color_limits[2])) {
          legend_barwidth <- 1
        }
      }
      # Reverse palette
      palette <- rev(palette)
    }
    # Labels
    if (label_splits == TRUE & label_statistics == TRUE) {
      add_labels <- TRUE
      tmp_seurat$labels <- paste0(split_labels, "\n", tmp_seurat$n_sig, " DE features\np = ", round(tmp_seurat$pvalue, 4))
      tmp_seurat$labels[is.na(tmp_seurat$n_sig) & is.na(tmp_seurat$pvalue)] <- split_labels[is.na(tmp_seurat$n_sig) & is.na(tmp_seurat$pvalue)]
      tmp_seurat$labels[!is.na(tmp_seurat$n_sig) & is.na(tmp_seurat$pvalue)] <- paste0(split_labels[!is.na(tmp_seurat$n_sig) & is.na(tmp_seurat$pvalue)],
                                                                                       "\n", tmp_seurat$n_sig[!is.na(tmp_seurat$n_sig) & is.na(tmp_seurat$pvalue)],
                                                                                       " DE features")

      Seurat::Idents(object = tmp_seurat) <- "labels"
    } else if (label_splits == TRUE) {
      add_labels <- TRUE
      tmp_seurat$labels <- split_labels
      Seurat::Idents(object = tmp_seurat) <- "labels"
    } else if (label_statistics == TRUE) {
      add_labels <- TRUE
      tmp_seurat$labels <- paste0(tmp_seurat$n_sig, " DE features\np = ", round(tmp_seurat$pvalue, 4))
      tmp_seurat$labels[is.na(tmp_seurat$n_sig) & is.na(tmp_seurat$pvalue)] <- NA
      tmp_seurat$labels[!is.na(tmp_seurat$n_sig) & is.na(tmp_seurat$pvalue)] <- paste0(tmp_seurat$n_sig[!is.na(tmp_seurat$n_sig) & is.na(tmp_seurat$pvalue)],
                                                                                       " DE features")
      Seurat::Idents(object = tmp_seurat) <- "labels"
    }
  } else if (color_by == "feature") {
    legend_barwidth <- NULL
    legend_barheight <- 1
    type <- "FeaturePlot"
    # Add feature values to metadata
    tmp_seurat$color_groups <- feature_values
    # Color palette
    palette <- permuteDEpalette(type = "gradient",
                                palette_name = palette_name)
    # Set NA cutoff
    if (!is.null(color_limits)) {
      na_cutoff_low <- color_limits[1]
      na_cutoff_high <- color_limits[2]
    } else if (min(tmp_seurat$color_groups, na.rm = TRUE) == 0) {
      na_cutoff_low <- 0.00001
      na_zero_legend <- TRUE
    }
    # Set barwidth if only 1 value
    if ((dplyr::n_distinct(tmp_seurat$color_groups[!is.na(tmp_seurat$color_groups)]) < 2)) {
      if (is.null(color_limits) || (color_limits[1] == color_limits[2])) {
        legend_barwidth <- 1
      }
    }
    # Color legend
    if (is.null(feature_name)) {
      color_legend <- ""
    } else {
      color_legend <- feature_name
    }
    # Labels
    if (label_splits == TRUE & !is.null(split_labels)) {
      add_labels <- TRUE
      tmp_seurat$labels <- split_labels
      Seurat::Idents(object = tmp_seurat) <- "labels"
    }
    if (label_statistics == TRUE) {
      warning("Input value for 'label_statistics' is not used when parameter 'color_by' is not 'n_sig' or 'pvalue'.")
    }
    if (permutation_test_alpha < 1) {
      permutation_test_alpha <- 1
    }
  }

  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------
  if (type == "DimPlot") {
    p <- Seurat::DimPlot(tmp_seurat,
                         reduction = "dim_reduction",
                         group.by = "color_groups",
                         label = add_labels,
                         ...) +
      ggplot2::theme_void() +
      permuteDEtheme() +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.5),
                     plot.title = ggplot2::element_blank(),
                     legend.box.margin = ggplot2::margin(5, 5, 5, 5)) +
      ggplot2::scale_color_manual(values = palette) +
      ggplot2::labs(color = color_legend) +
      ggplot2::xlab("Dim 1") +
      ggplot2::ylab("Dim 2")
  } else if (type == "FeaturePlot") {
    p <- Seurat::FeaturePlot(tmp_seurat,
                             features = "color_groups",
                             reduction = "dim_reduction",
                             label = add_labels,
                             label.size = 3.5,
                             ...) +
      ggplot2::theme_void() +
      permuteDEtheme(legend_title_vjust = 1) +
      ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.5),
                     plot.title = ggplot2::element_blank(),
                     legend.title = ggtext::element_markdown(margin = ggplot2::margin(l = 5, r = 5)),
                     legend.box.margin = ggplot2::margin(5, 12, 5, 5)) +
      ggplot2::scale_color_gradientn(colors = palette,
                                     na.value = na_color,
                                     limits = c(na_cutoff_low, na_cutoff_high),
                                     guide = ggplot2::guide_colorbar(frame.colour = "black",
                                                                     ticks.colour = "black",
                                                                     barwidth = legend_barwidth,
                                                                     barheight = legend_barheight)) +
      ggplot2::xlab("Dim 1") +
      ggplot2::ylab("Dim 2")
    # Add NA legend
    if (permutation_test_alpha < 1) {
      sample_x <- min(tmp_seurat@reductions$dim_reduction@cell.embeddings[,1])
      sample_y <- min(tmp_seurat@reductions$dim_reduction@cell.embeddings[,2])
      p <- p + ggplot2::geom_point(x = sample_x,
                                   y = sample_y,
                                   alpha = 0,
                                   ggplot2::aes(fill = 0)) +
        ggplot2::scale_fill_gradientn(colors = na_color,
                                      labels = "",
                                      guide = ggplot2::guide_colorbar(title = na_legend,
                                                                      ticks = FALSE,
                                                                      frame.colour = "black",
                                                                      barwidth = 1,
                                                                      barheight = 1,
                                                                      order = 2)) +
        ggplot2::guides(color = ggplot2::guide_colorbar(color_legend,
                                                        frame.colour = "black",
                                                        ticks.colour = "black",
                                                        barwidth = legend_barwidth,
                                                        barheight = legend_barheight,
                                                        order = 1))
    } else if (na_zero_legend == TRUE) {
      sample_x <- min(tmp_seurat@reductions$dim_reduction@cell.embeddings[,1])
      sample_y <- min(tmp_seurat@reductions$dim_reduction@cell.embeddings[,2])
      p <- p + ggplot2::geom_point(x = sample_x,
                                   y = sample_y,
                                   alpha = 0,
                                   ggplot2::aes(fill = 0)) +
        ggplot2::scale_fill_gradientn(colors = na_color,
                                      labels = "0",
                                      guide = ggplot2::guide_colorbar(title = color_legend,
                                                                      frame.colour = "black",
                                                                      ticks.colour = "black",
                                                                      barwidth = 1,
                                                                      barheight = 1,
                                                                      order = 1)) +
        ggplot2::guides(color = ggplot2::guide_colorbar("",
                                                        frame.colour = "black",
                                                        ticks.colour = "black",
                                                        barwidth = legend_barwidth,
                                                        barheight = legend_barheight,
                                                        order = 2)) +
        ggplot2::theme(legend.spacing = ggplot2::unit(0, "pt"))
    } else {
      p <- p + ggplot2::labs(color = color_legend)
    }
  }

  # Fix coords
  if (fix_coords == TRUE) {
    p <- p + ggplot2::coord_fixed()
  }

  return(p)
}
