#' Perform permutation test
#'
#' This function performs a permutation test by shuffling the group labels for
#' the given set of replicates and performing pseudobulk differential expression
#' between those permuted groups. The metric used by the permutation test is the
#' number of significantly differentially expressed features for a given
#' significance level 'alpha'.
#'
#' As input, this function requires the output from function \code{runDE}
#' containing the pseudobulk values for each feature and the differential
#' expression results for the true group labels.
#'
#' Permutation is performed at the replicate level, rather than the cell level,
#' and is performed without replacement, such that each iteration is a unique
#' permutation of the group labels.
#'
#' Notably, this permutation test is not performed to identify the true
#' differentially expressed features, but rather to assess how many false
#' positive significant differentially expressed features can be expected by
#' chance, and to characterize the log fold change and significance observed for
#' such false positives.
#'
#' @param input Output from function 'runDE' containing differential expression
#' results, pseudobulk values, and parameters used.
#' @param alpha A numeric value indicating the significance level used for
#' permutation test comparisons of the number of differentially expressed
#' features. Defaults to 0.05.
#' @param n_iterations A numeric value indicating the number of iterations run
#' for the permutation test. Defaults to 1000.
#' @param use_splits A vector containing the names of splits to use. Defaults to
#' \code{NULL}.
#' @param random_seed A numeric value indicating the random seed to be used.
#' Defaults to 1.
#' @param n_cores A numeric value indicating the number of cores to use for
#' parallelization. Default = \code{NULL} will use the number of available cores
#' minus 2.
#' @param verbose A boolean value indicating whether to use verbose output
#' during the execution of this function. Defaults to \code{TRUE}.
#' Can be set to \code{FALSE} for a cleaner output.
#'
#' @return Returns a list containing the following elements: \describe{
#'   \item{permutation_test_results}{Dataframe containing the permutation test
#'   results by split}
#'   \item{permutation_summary}{Dataframe containing the permutation DE summary
#'   metrics by split}
#'   \item{parameters}{Dataframe record of parameter values used}
#'   }
#'
#' @export
#'
permuteDE <- function(input,
                      alpha = 0.05,
                      n_iterations = 1000,
                      use_splits = NULL,
                      random_seed = 1,
                      n_cores = NULL,
                      verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(input, "input")
  .validInput(alpha, "alpha")
  .validInput(n_iterations, "n_iterations")
  .validInput(return_all, "return_all")
  .validInput(use_splits, "use_splits", input)
  .validInput(random_seed, "random_seed")
  .validInput(n_cores, "n_cores")
  .validInput(verbose, "verbose")

  # Set defaults & fetch values
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 2
  }
  if (is.null(use_splits)) {
    use_splits <- names(input$PB_values)
  }
  n_splits <- length(use_splits)
  de_method <- input$parameters$de_method
  de_test <- input$parameters$de_test
  p_adjust_method <- input$parameters$p_adjust_method

  # ---------------------------------------------------------------------------
  # Compute permuted DE results
  # ---------------------------------------------------------------------------

  # Progress
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                       " : Differential expression results for group ",
                       unique(input$group_key$group)[1], " vs. ", unique(input$group_key$group)[2],
                       " across ", n_splits, " pseudobulk matrices..")

  # Set up
  permutation_DE_results <- data.frame(split = NULL,
                                       permutation = NULL,
                                       n_sig = NULL,
                                       min_lfc = NULL,
                                       max_lfc = NULL)
  permutation_test_results <- data.frame(split = NULL,
                                         true_n_sig = NULL,
                                         p_n_sig = NULL,
                                         n_iterations = NULL)

  # For each split
  for (s in 1:n_splits) {
    # Extract relevant values
    current_split <- use_splits[s]
    current_pb <- input$PB_values[[current_split]]
    current_replicates <- colnames(current_pb)
    n_replicates <- length(current_replicates)
    true_groups <- input$group_key[current_replicates, "group"]
    # Proceed if there are two groups present
    if (dplyr::n_distinct(true_groups) == 2) {
      # Set remaining values
      group1 <- sort(unique(true_groups))[1]
      group2 <- sort(unique(true_groups))[2]
      n_group1 <- sum(true_groups == group1)
      # Compute a set of permuted group labels
      permuted_group_labels <- getCombinations(n_replicates = n_replicates,
                                               n_group1 = n_group1,
                                               n_combinations = n_iterations,
                                               verbose = TRUE)
      current_n_iterations <- ncol(permuted_group_labels)
      # Progress
      if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                           " : Running ", current_n_iterations, " permutations for split ", current_split,
                           " (", s, "/", n_splits, ")..")
      # Run permutations
      permutation_DE_results_list <- pbmcapply::pbmclapply(seq_len(current_n_iterations),
                                                     FUN = function(i) {
                                                       # Create target dataframe with permuted group labels
                                                       targets_i <- data.frame(replicate = current_replicates,
                                                                               group = group2)
                                                       targets_i[permuted_group_labels[,i], "group"] <- group1
                                                       rownames(targets_i) <- targets_i$replicate
                                                       # Skip if true groups
                                                       if (!identical(targets_i$group, true_groups)) {
                                                         # Create design
                                                         design_i <- stats::model.matrix(~ group, data = targets_i)
                                                         # Run DE
                                                         de_results_i <- switch(de_method,
                                                                                edgeR = .runDE.edgeR(pseudobulk = current_pb,
                                                                                                     targets = targets_i,
                                                                                                     design = design_i,
                                                                                                     de_test = de_test),
                                                                                DESeq2 = .runDE.DESeq2(pseudobulk = current_pb,
                                                                                                       targets = targets_i,
                                                                                                       de_test = de_test),
                                                                                limma = .runDE.limma(pseudobulk = current_pb,
                                                                                                     targets = targets_i,
                                                                                                     design = design_i,
                                                                                                     de_test = de_test))
                                                         de_results_i <- de_results_i %>%
                                                           dplyr::mutate(p_adjust = stats::p.adjust(p_value, method = p_adjust_method)) %>%
                                                           dplyr::summarise(n_sig = sum(p_adjust < alpha),
                                                                            min_lfc = min(lfc),
                                                                            max_lfc = max(lfc)) %>%
                                                           data.frame() %>%
                                                           dplyr::mutate(permutation = i,
                                                                         split = current_split) %>%
                                                           dplyr::select(split, permutation, n_sig, min_lfc, max_lfc)
                                                       } else {
                                                         de_results_i <- NULL
                                                       }
                                                       return(de_results_i)
                                                     },
                                                     mc.cores = n_cores,
                                                     mc.set.seed = TRUE)
      permutation_DE_results_s <- do.call(rbind, permutation_DE_results_list) %>% data.frame()
      # If true labels were not removed, remove one random set of labels so total number includes true set
      if (nrow(permutation_DE_results_s) == current_n_iterations) {
        permutation_DE_results_s <- permutation_DE_results_s[-nrow(permutation_DE_results_s),]
      }
      # Conduct permutation test
      true_values <- input$DE_results %>%
        dplyr::filter(split == current_split) %>%
        dplyr::summarise(n_sig = sum(p_adjust < alpha))
      permutation_test_p_value <- 1 - stats::ecdf(c(true_values$true_n_sig, permutation_DE_results_s$n_sig))(true_values$n_sig)
      permutation_test_results_s <- data.frame(split = current_split,
                                               true_n_sig = true_values$n_sig,
                                               p_n_sig = permutation_test_p_value,
                                               n_iterations = current_n_iterations)
      # Add to overall results
      permutation_DE_results <- rbind(permutation_DE_results, permutation_DE_results_s)
      permutation_test_results <- rbind(permutation_test_results, permutation_test_results_s)
    } else {
      # Progress
      if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                           " : Skipping split ", current_split, ", ", dplyr::n_distinct(true_groups), " group(s) present.")
    }
  }

  # ---------------------------------------------------------------------------
  # Wrap up
  # ---------------------------------------------------------------------------

  parameter_list <- list("alpha" = alpha,
                         "n_iterations" = n_iterations,
                         "use_splits" = use_splits,
                         "de_method" = de_method,
                         "de_test" = de_test,
                         "p_adjust_method" = p_adjust_method,
                         "random_seed" = random_seed)

  # Return
  return(list("permutation_test_results" = permutation_test_results,
              "permutation_DE_results" = permutation_DE_results,
              "parameters" = parameter_list))
}
