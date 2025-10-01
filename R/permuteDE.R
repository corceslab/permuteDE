#' Perform permutation test
#'
#' This function performs a permutation test by shuffling the group labels for
#' the given set of replicates and performing pseudobulk differential expression
#' between those permuted groups. The metric used by the permutation test is the
#' number of significantly differentially expressed features for a given
#' significance level \code{alpha} and log fold change threshold
#' \code{lfc_threshold}.
#'
#' As input, this function requires the output from function \code{runDE}
#' containing the pseudobulk (or cell-level) values for each feature and the
#' differential expression results for the true group labels.
#'
#' Permutation is always performed at the replicate level, rather than the cell
#' level, and is performed without replacement, such that each iteration is a
#' unique permutation of the group labels.
#'
#' Notably, this permutation test does not pass judgment on any individual gene,
#' rather, it is intended to assess how many false positive significant
#' differentially expressed features can be expected by chance. In addition, it
#' can be used to characterize the log fold change and significance observed for
#' such false positives to help users prioritize reliable DE results.
#'
#' @param input Output from function \code{runDE}: a list containing
#' differential expression results, pseudobulk values, and parameters used.
#' @param alpha A numeric value indicating the significance level used for
#' permutation test comparisons of the number of differentially expressed
#' features. Defaults to 0.05.
#' @param lfc_threshold A numeric value indicating the minimum absolute value
#' log fold change for a gene to be counted as a "hit". Defaults to 0.5. Set to
#' 0 to disregard log fold change when counting hits.
#' @param n_iterations A numeric value indicating the number of iterations run
#' for the permutation test. Defaults to 1000. Computational time increases
#' approximately linearly with the number of iterations.
#' @param use_splits A vector containing the names of splits to use. Defaults to
#' \code{NULL}, which will try all splits.
#' @param min_DE A numeric value indicating the minimum number of
#' differentially expressed features between the true group labels for a split,
#' below which permutations will not be run. Defaults to 2. Set to 0 to run
#' permutation test for all splits, regardless of the true number of DEGs.
#' @param return_all A Boolean value indicating whether to store and return all
#' DE results (log fold changes and p-values per gene per split) for every
#' single permutation. Defaults = \code{FALSE} will return only high-level
#' permutation test results. Note that setting this to \code{TRUE} will
#' substantially increase the size of the returned output.
#' @param random_seed A numeric value indicating the random seed to be used.
#' Defaults to 1.
#' @param n_cores A numeric value indicating the number of cores to use for
#' parallelization. Default = \code{NULL} will use the number of available cores
#' minus 2.
#' @param verbose A Boolean value indicating whether to use verbose output
#' during the execution of this function. Defaults to \code{TRUE}.
#' Can be set to \code{FALSE} for a cleaner output.
#'
#' @return Returns a list containing the following elements: \describe{
#'   \item{permutation_test_results}{Dataframe containing the permutation test
#'   results by split}
#'   \item{permutation_DE_summary}{Dataframe containing the permutation DE summary
#'   metrics by split}
#'   \item{permutation_DE_results}{If parameter 'return_all' is TRUE,
#'   dataframe DE results for each feature, by split, for each iteration}
#'   \item{metadata}{List recording characteristics of the data and runtime}
#'   \item{parameters}{List recording parameter values used}
#'   }
#'
#' @export
#'
permuteDE <- function(input,
                      alpha = 0.05,
                      lfc_threshold = 0.5,
                      n_iterations = 1000,
                      use_splits = NULL,
                      min_DE = 2,
                      return_all = FALSE,
                      random_seed = 1,
                      n_cores = NULL,
                      verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  time1 <- Sys.time()

  .validInput(input, "input", "permuteDE")
  .validInput(alpha, "alpha")
  .validInput(lfc_threshold, "lfc_threshold")
  .validInput(n_iterations, "n_iterations")
  .validInput(use_splits, "use_splits", input)
  .validInput(min_de, "min_de")
  .validInput(return_all, "return_all")
  .validInput(random_seed, "random_seed")
  .validInput(n_cores, "n_cores")
  .validInput(verbose, "verbose")

  # Set defaults
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 2
  }

  # ---------------------------------------------------------------------------
  # Set up
  # ---------------------------------------------------------------------------

  # Fetch values
  reference_group <- input$parameters$reference_group
  de_method <- input$parameters$de_method
  de_test <- input$parameters$de_test
  de_params <- input$parameters$de_params
  p_adjust_method <- input$parameters$p_adjust_method
  pseudobulk <- input$parameters$pseudobulk
  stored_replicates <- input$parameters$stored_replicates

  if (is.null(use_splits)) {
    if (pseudobulk == "none") {
      use_splits <- names(input$cell_values)
    } else {
      use_splits <- names(input$PB_values)
    }
  }

  # Calculate true number of DE features
  true_DE_values <- input$DE_results |>
    dplyr::group_by(split) |>
    dplyr::summarise(
      n_sig = sum(padj < alpha & abs(lfc) > lfc_threshold, na.rm = TRUE))

  # Remove splits with fewer than the minimum number of DE features
  use_splits <- intersect(use_splits,
                          dplyr::filter(true_DE_values, n_sig >= min_DE)$split)
  remove_splits <- dplyr::filter(true_DE_values, n_sig < min_DE)$split
  if (verbose & length(remove_splits) > 0) {
    message(format(Sys.time(), "%Y-%m-%d %X")," : Will skip ", length(remove_splits), " split label",
            ifelse(length(remove_splits) == 1, "", "s"),
            " due to insufficient differentially expressed features: ",
            paste0(remove_splits,
                   collapse = ", "))
  }

  n_splits <- length(use_splits)

  # Stop if no splits pass min_DE threshold
  if (n_splits == 0) {
    stop("No splits have sufficient differentially expressed features, so no permutation tests were run. Consider setting parameter 'min_DE' to 0.")
  }

  # ---------------------------------------------------------------------------
  # Compute permuted DE results
  # ---------------------------------------------------------------------------

  time2 <- Sys.time()

  # Progress
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                       " : Differential expression results for group ",
                       unique(input$metadata$group_key$group)[1], " vs. ", unique(input$metadata$group_key$group)[2],
                       " across ", n_splits, " pseudobulk matrices..")

  # Set up
  permutation_DE_summary <- data.frame(split = NULL,
                                       permutation = NULL,
                                       n_sig = NULL,
                                       min_lfc = NULL,
                                       max_lfc = NULL)
  permutation_test_results <- data.frame(split = NULL,
                                         true_n_sig = NULL,
                                         pvalue_n_sig = NULL,
                                         n_iterations = NULL)
  if (return_all == TRUE) {
    permutation_DE_results <- data.frame(gene = NULL,
                                         lfc = NULL,
                                         pvalue = NULL,
                                         padj = NULL,
                                         permutation = NULL,
                                         split = NULL)
  }

  # For each split
  for (s in 1:n_splits) {
    # Extract relevant values
    current_split <- use_splits[s]
    if (pseudobulk == "none" & !is.null(stored_replicates)) {
      # If cell-level but permuting biological replicates
      # Grab feature x cell matrix
      current_mat <- input$cell_values[[current_split]]
      current_replicates <- colnames(current_mat)
      # Subset biological replicates to current matrix
      current_stored_replicates <- stored_replicates[current_replicates]
      n_replicates <- length(current_stored_replicates)
      true_groups <- input$metadata$group_key[current_replicates, "group"]
      replicate_set <- unique(current_stored_replicates)
    } else {
      # Grab matrix and get values directly from there
      if (pseudobulk == "none") {
        current_mat <- input$cell_values[[current_split]]
      } else {
        current_mat <- input$PB_values[[current_split]]
      }
      current_replicates <- colnames(current_mat)
      n_replicates <- length(current_replicates)
      true_groups <- input$metadata$group_key[current_replicates, "group"]
    }
    # Proceed if there are two groups present
    if (dplyr::n_distinct(true_groups) == 2) {
      # Set remaining values
      group1 <- sort(unique(true_groups))[1]
      group2 <- sort(unique(true_groups))[2]
      if (pseudobulk == "none" & !is.null(stored_replicates)) {
        # If cell-level but permuting biological replicates
        n_group1 <- dplyr::n_distinct(current_stored_replicates[true_groups == group1])
      } else {
        n_group1 <- sum(true_groups == group1)
      }
      # Compute a set of permuted group labels
      permuted_group_indices <- getCombinations(n_replicates = n_replicates,
                                                n_group1 = n_group1,
                                                n_combinations = n_iterations,
                                                message = paste0(" for split ", current_split, " (", s, "/", n_splits, ")"),
                                                random_seed = random_seed,
                                                verbose = TRUE)
      current_n_iterations <- ncol(permuted_group_indices)
      # If cell-level but permuting biological replicates, convert permuted group indices for each replicate into cell indices
      if (pseudobulk == "none" & !is.null(stored_replicates)) {
        new_permuted_group_indices <- apply(permuted_group_indices, 2,
                                            FUN = function(i) {
                                              # Which replicates are in these indices
                                              replicates_i <- replicate_set[i]
                                              # Which cell indices belong to these replicates
                                              current_cell_indices <- current_stored_replicates %in% replicates_i
                                              return(current_cell_indices)
                                            })
        # Convert to matrix (make lengths match by padding w/ NA)
        max_group_indices <- max(lengths(new_permuted_group_indices))
        new_permuted_group_indices <- lapply(new_permuted_group_indices,
                                             FUN = function(i) {
                                               length(i) <- max_group_indices
                                               return(i)
                                               })
        permuted_group_indices <- do.call(cbind, new_permuted_group_indices)
      }
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
                                                       targets_i[permuted_group_indices[,i][!is.na(permuted_group_indices[,i])], "group"] <- group1
                                                       rownames(targets_i) <- targets_i$replicate

                                                       group_factor <- factor(targets_i$group)
                                                       group_factor <- stats::relevel(group_factor, ref = reference_group)
                                                       targets_i$group <- group_factor

                                                       # Skip if true groups
                                                       if (!identical(targets_i$group, true_groups)) {
                                                         # Create design
                                                         design_i <- stats::model.matrix(~ group, data = targets_i)
                                                         # Run DE
                                                         de_results_i <- switch(de_method,
                                                                                edgeR = .runDE.edgeR(mat = current_mat,
                                                                                                     targets = targets_i,
                                                                                                     design = design_i,
                                                                                                     de_test = de_test,
                                                                                                     de_params = de_params),
                                                                                DESeq2 = .runDE.DESeq2(mat = current_mat,
                                                                                                       targets = targets_i,
                                                                                                       design = design_i,
                                                                                                       de_test = de_test,
                                                                                                       de_params = de_params),
                                                                                limma = .runDE.limma(mat = current_mat,
                                                                                                     targets = targets_i,
                                                                                                     design = design_i,
                                                                                                     de_test = de_test,
                                                                                                     de_params = de_params),
                                                                                presto = .runDE.presto(mat = current_mat,
                                                                                                       targets = targets_i,
                                                                                                       de_test = de_test,
                                                                                                       de_params = de_params))
                                                         de_results_i <- de_results_i |>
                                                           dplyr::mutate(padj = stats::p.adjust(pvalue, method = p_adjust_method),
                                                                         permutation = i,
                                                                         split = current_split) |>
                                                           dplyr::arrange(padj)
                                                         if (return_all != TRUE) {
                                                           de_results_i <- de_results_i |>
                                                             dplyr::group_by(split, permutation) |>
                                                             dplyr::summarise(
                                                               n_sig = sum(padj < alpha & abs(lfc) > lfc_threshold),
                                                               min_lfc_sig = ifelse(n_sig > 0, min(lfc[padj < alpha & abs(lfc) > lfc_threshold]), NA),
                                                               max_lfc_sig = ifelse(n_sig > 0, max(lfc[padj < alpha & abs(lfc) > lfc_threshold]), NA)) |>
                                                             data.frame() |>
                                                             dplyr::select(split, permutation, n_sig, min_lfc_sig, max_lfc_sig)
                                                         }
                                                       } else {
                                                         de_results_i <- NULL
                                                       }
                                                       return(de_results_i)
                                                     },
                                                     mc.cores = n_cores,
                                                     mc.set.seed = TRUE)

      if (return_all == TRUE) {
        permutation_DE_results_s <- do.call(rbind, permutation_DE_results_list) |> data.frame()
        permutation_DE_summary_s <- permutation_DE_results_s |>
          dplyr::group_by(split, permutation) |>
          dplyr::summarise(
            n_sig = sum(padj < alpha & abs(lfc) > lfc_threshold),
            min_lfc_sig = ifelse(n_sig > 0, min(lfc[padj < alpha & abs(lfc) > lfc_threshold]), NA),
            max_lfc_sig = ifelse(n_sig > 0, max(lfc[padj < alpha & abs(lfc) > lfc_threshold]), NA)) |>
          data.frame() |>
          dplyr::select(split, permutation, n_sig, min_lfc_sig, max_lfc_sig)
      } else {
        permutation_DE_summary_s <- do.call(rbind, permutation_DE_results_list) |> data.frame()
      }

      # If true labels were among random permutations, they were then skipped, leaving n_iterations - 1
      # If not, remove one random set of labels, leaving n_iterations - 1
      if (nrow(permutation_DE_summary_s) == current_n_iterations) {
        permutation_DE_summary_s <- permutation_DE_summary_s[-nrow(permutation_DE_summary_s),]
      }

      # Conduct permutation test
      true_n_sig <- dplyr::filter(true_DE_values, split == current_split)$n_sig

      permutation_test_p_value <- sum(c(true_n_sig, permutation_DE_summary_s$n_sig) >= true_n_sig)/current_n_iterations
      permutation_test_results_s <- data.frame(split = current_split,
                                               true_n_sig = true_n_sig,
                                               pvalue_n_sig = permutation_test_p_value,
                                               n_iterations = current_n_iterations)
      # Add to overall results
      permutation_DE_summary <- rbind(permutation_DE_summary, permutation_DE_summary_s)
      permutation_test_results <- rbind(permutation_test_results, permutation_test_results_s)
      if (return_all == TRUE) {
        permutation_DE_results <- rbind(permutation_DE_results, permutation_DE_results_s)
      }

    } else {
      # Progress
      if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                           " : Skipping split ", current_split, ", ", dplyr::n_distinct(true_groups), " group(s) present.")
    }
  }

  # ---------------------------------------------------------------------------
  # Wrap up
  # ---------------------------------------------------------------------------

  time3 <- Sys.time()

  # Metadata
  metadata_list <- list("time" = data.frame(total = difftime(time3, time1, units = "secs"),
                                              step1_setup = difftime(time2, time1, units = "secs"),
                                              step2_permute = difftime(time3, time2, units = "secs")))

  parameter_list <- list("alpha" = alpha,
                         "lfc_threshold" = lfc_threshold,
                         "n_iterations" = n_iterations,
                         "use_splits" = use_splits,
                         "min_DE" = min_DE,
                         "reference_group" = reference_group,
                         "de_method" = de_method,
                         "de_test" = de_test,
                         "de_params" = de_params,
                         "p_adjust_method" = p_adjust_method,
                         "pseudobulk" = pseudobulk,
                         "return_all" = return_all,
                         "random_seed" = random_seed,
                         "n_cores" = n_cores)

  # Return
  if (return_all == TRUE) {
    return(list("permutation_test_results" = permutation_test_results,
                "permutation_DE_summary" = permutation_DE_summary,
                "permutation_DE_results" = permutation_DE_results,
                "metadata" = metadata_list,
                "parameters" = parameter_list))
  } else {
    return(list("permutation_test_results" = permutation_test_results,
                "permutation_DE_summary" = permutation_DE_summary,
                "metadata" = metadata_list,
                "parameters" = parameter_list))
  }
}
