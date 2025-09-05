#' Generate Histograms of Permutation Test Results
#'
#' This function creates a list of `ggplot2` histogram plots showing the distribution 
#' of the number of differentially expressed (DE) genes across permutations for each 
#' split in the permutation DE results. Each histogram includes a vertical line at 
#' the observed number of DE genes and a p-value annotation.
#'
#' @param permuteDE_result A returned list from `permuteDE`
#' @export

generateHistogram <- function(permuteDE_result) {
  library(ggplot2)
  library(stringr)
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
    
    ggPlotHistogram(df_perm$n_sig, xlabel = 'Number of DE genes', histAlpha = 0.5, vline = df_test$true_n_sig) + 
      geom_vline(xintercept = df_test$true_n_sig, 
                 color = 'red', 
                 linewidth = 1, 
                 linetype = 'solid') +
      annotate(
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
      labs(title = str_wrap(
        paste0('Distribution of # DE genes across ', 
               n_iterations, ' Iterations in ', split_id), 
        width = 40
      ))
  })
  names(histogram_list) <- all_splits
  return(histogram_list)
}
