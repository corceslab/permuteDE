
#' Generate Volcano Plots from Differential Expression Results
#'
#' This function creates a list of volcano plots for each subset of 
#' differential expression (DE) results contained in a `runDE_output` object. 
#' The DE results are split by the `split` column, and each split produces 
#' its own volcano plot using the \pkg{EnhancedVolcano} package. 
#'
#' @param runDE_output A list-like object returned from a DE analysis, 
#'   containing at least a data frame `DE_result` with columns:
#'   \describe{
#'     \item{gene}{Gene identifiers (character).}
#'     \item{lfc}{Log2 fold change values (numeric).}
#'     \item{padj}{Adjusted p-values (numeric).}
#'     \item{split}{Grouping variable used to subset results (character or factor).}
#'   }
#' @param alpha Numeric scalar. Adjusted p-value cutoff for significance 
#'   (default = 0.05).
#' @param lfc_threshold Numeric scalar. Log2 fold change threshold for 
#'   significance (default = 0.5).
#' @param title Character string or `NULL`. Plot title for all volcano plots. 
#'   If `NULL`, a title is automatically generated per split (default = `NULL`).
#' @param subtitle Character string or `NULL`. Subtitle for all volcano plots. 
#'   If `NULL`, a subtitle describing the significance thresholds is generated 
#'   automatically (default = `NULL`).
#' @param significant_color Character string. Color used to highlight 
#'   significant upregulated genes (default = `"red2"`).
#'
#' @return A named list of \pkg{ggplot2} objects, where each element 
#'   corresponds to a volcano plot for one split of the DE results. 
#'
#' @details
#' Each volcano plot is created using \code{EnhancedVolcano()}, customized with:
#' \itemize{
#'   \item The x-axis showing log2 fold changes (`lfc`).
#'   \item The y-axis showing adjusted p-values (`padj`).
#'   \item Genes colored grey, blue, or red depending on significance status.
#'   \item A caption showing the total number of genes in that split.
#'   \item \code{themeCorcesRegular()} applied, with the legend removed.
#' }
#'
#' @examples
#' \dontrun{
#' plots <- generateVolcanoPlot(runDE_output)
#' print(plots[[1]])
#' }
#'
#' @export
generateVolcanoPlot <- function(runDE_output,
                                alpha = 0.05,
                                lfc_threshold = 0.5,
                                title = NULL,
                                subtitle = NULL,
                                significant_color = 'red2'){
  
  library(EnhancedVolcano)
  
  split_results <- split(runDE_output$DE_result, runDE_output$DE_result$split)
  
  volcano_list <- lapply(names(split_results), function(split_id) {
    df <- split_results[[split_id]]
    
    EnhancedVolcano(
      df,
      lab = df$gene,
      x = 'lfc',
      y = 'padj',
      pCutoff = alpha,
      FCcutoff = lfc_threshold,
      title = if (is.null(title)) paste('Volcano Plot of Gene Expression in', split_id) else title,
      subtitle = if (is.null(subtitle)) 
        paste("Significance: adjusted p-value <", alpha, "and |log2 fold change| >", lfc_threshold) 
      else subtitle,
      caption = paste0("total = ", nrow(df), " genes"),
      col = c("grey30", "grey30", "royalblue", significant_color)
    ) + 
      themeCorcesRegular() + 
      theme(legend.position = 'none')
  })

  names(volcano_list) <- names(split_results)
  
  return(volcano_list)
  
  
  
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
    plot.margin = unit(c(plot_margin_cm, plot_margin_cm, plot_margin_cm, plot_margin_cm), "cm"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = color, size = (4/3) * base_rect_size),
    axis.ticks.length = unit(axis_tick_length_mm, "mm"), 
    axis.ticks = element_line(color = color, size = (4/3) * base_line_size),
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

