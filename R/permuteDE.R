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
#' Notably, this permutation test does not pass judgment on any individual
#' feature, rather, it is intended to assess how many false positive significant
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
#' log fold change for a feature to be counted as a "hit". Default = 0
#' disregards log fold change when counting hits.
#' @param n_iterations A numeric value indicating the number of iterations run
#' for the permutation test. Defaults to 1000. Computational time increases
#' approximately linearly with the number of iterations.
#' @param use_splits A vector containing the names of splits to use. Defaults to
#' \code{NULL}, which will try all splits.
#' @param permute_by An optional character string indicating the name of the
#' column in \code{input$metadata$group_key} indicating partitions of the data
#' within which group labels should be shuffled together as a unit. Alternately,
#' a vector ordered according to the values in \code{input$metadata$group_key}.
#' Only for specific use cases such as cell-level tests. Default = \code{NULL}
#' will shuffle group labels across all replicates.
#' @param permute_within An optional character string indicating the name of the
#' column in \code{input$metadata$group_key} indicating partitions of the data
#' within which group labels should be shuffled separately. Alternately, a
#' vector ordered according to the values in \code{input$metadata$group_key}.
#' Only for specific use cases such as cell-level tests. Default = \code{NULL}
#' will shuffle group labels across all replicates.
#' @param min_DE A numeric value indicating the minimum number of
#' differentially expressed features between the true group labels for a split,
#' below which permutations will not be run. Defaults to 1. Set to 0 to run
#' permutation test for all splits, regardless of the true number of DEGs.
#' @param return_all A Boolean value indicating whether to store and return all
#' DE results (log fold changes and p-values per feature per split) for every
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
#'   \item{permutation_DE_summary}{Dataframe containing the permutation DE
#'   summary metrics by split}
#'   \item{permutation_DE_results}{If parameter 'return_all' is TRUE,
#'   dataframe DE results for each feature, by split, for each iteration}
#'   \item{metadata}{List recording additional characteristics of the data: the
#'   number of significant DE features from runDE, the group indices for each
#'   iteration of the permutation test, and runtime}
#'   \item{parameters}{List recording parameter values used}
#'   }
#'
#' @export
#'
permuteDE <- function(input,
                      alpha = 0.05,
                      lfc_threshold = 0,
                      n_iterations = 1000,
                      use_splits = NULL,
                      permute_by = NULL,
                      permute_within = NULL,
                      min_DE = 1,
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
  .validInput(use_splits, "use_splits", list(input, "permuteDE"))
  .validInput(min_DE, "min_DE")
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
  non_reference_group <- input$parameters$non_reference_group
  design_formula <- input$parameters$design_formula
  de_method <- input$parameters$de_method
  de_test <- input$parameters$de_test
  de_params <- input$parameters$de_params
  p_adjust_method <- input$parameters$p_adjust_method
  pseudobulk <- input$parameters$pseudobulk

  # Validate retrieved parameters
  .validInput(de_method, "de_method")
  .validInput(de_test, "de_test", de_method)
  .validInput(de_params, "de_params", list(de_method, de_test))
  .validInput(p_adjust_method, "p_adjust_method")
  .validInput(pseudobulk, "pseudobulk", list("permuteDE", input))

  # Additional validation
  .validInput(permute_by, "permute_by", list(input, pseudobulk))
  .validInput(permute_within, "permute_within", list(input, design_formula))

  # Set use_splits
  if (is.null(use_splits)) {
    if (pseudobulk == "none") {
      use_splits <- names(input$cell_values)
    } else {
      use_splits <- names(input$PB_values)
    }
  }

  # Configure group key
  group_key <- input$metadata$group_key
  # If 'permute_by' provided, some replicates will always be kept together
  # (e.g., cells from the same sample)
  if (!is.null(permute_by)) {
    if (length(permute_by) == 1) {
      group_key$permute_by <- .retrieveData(object = NULL,
                                            metadata = group_key,
                                            type = "cell_metadata",
                                            name = permute_by)
    } else {
      group_key$permute_by <- permute_by
    }
  } else {
    group_key$permute_by <- group_key$replicate
  }
  # Check for NA values
  if (any(is.na(group_key$permute_by))) {
    stop("Values provided for 'permute_by' cannot be NA.")
  }
  group_key$permute_by <- as.character(group_key$permute_by)

  # If 'permute_within' provided, replicates will be shuffled separately
  # within each partition
  if (!is.null(permute_within)) {
    if (length(permute_within) == 1) {
      group_key$permute_within <- .retrieveData(object = NULL,
                                            metadata = group_key,
                                            type = "cell_metadata",
                                            name = permute_within)
    } else {
      group_key$permute_within <- permute_within
    }
  } else {
    group_key$permute_within <- "all"
  }
  # Check for NA values
  if (any(is.na(group_key$permute_within))) {
    stop("Values provided for 'permute_within' cannot be NA.")
  }
  group_key$permute_within <- as.character(group_key$permute_within)

  # Calculate true number of DE features
  runDE_values <- input$DE_results |>
    dplyr::group_by(split) |>
    dplyr::summarise(
      runDE_n_sig = sum(padj < alpha & abs(lfc) > lfc_threshold, na.rm = TRUE)) |>
    data.frame()

  # Remove splits with fewer than the minimum number of DE features
  use_splits <- intersect(use_splits,
                          dplyr::filter(runDE_values, runDE_n_sig >= min_DE)$split)
  remove_splits <- dplyr::filter(runDE_values, runDE_n_sig < min_DE)$split
  if (verbose & length(remove_splits) > 0) {
    message(format(Sys.time(), "%Y-%m-%d %X")," : Will skip ", length(remove_splits), " split label",
            ifelse(length(remove_splits) == 1, "", "s"),
            " due to insufficient differentially expressed features: ",
            paste0(remove_splits,
                   collapse = ", "))
  }

  n_splits <- length(use_splits)

  # ---------------------------------------------------------------------------
  # Compute permuted DE results
  # ---------------------------------------------------------------------------

  time2 <- Sys.time()

  # Progress
  proceed <- TRUE
  if (n_splits == 0) {
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : No splits have sufficient differentially expressed features, ",
                         "so permutation tests will not be run. Consider setting parameter 'min_DE' to 0.")
    warning(" No splits have sufficient differentially expressed features, so no permutation tests were run. Consider setting parameter 'min_DE' to 0.")
    proceed <- FALSE
  } else if (verbose) {
    message(format(Sys.time(), "%Y-%m-%d %X"),
            " : Differential expression results for group ",
            unique(input$metadata$group_key$group)[1], " vs. ", unique(input$metadata$group_key$group)[2],
            " across ", n_splits, " pseudobulk matrices..")
  }

  # Set up to collect output
  # Summary of DE results for each iteration
  permutation_DE_summary <- data.frame(split = NULL,
                                       permutation = NULL,
                                       reference_group_overlap = NULL,
                                       non_reference_group_overlap = NULL,
                                       n_sig = NULL,
                                       min_lfc_sig = NULL,
                                       max_lfc_sig = NULL,
                                       min_lfc_all = NULL,
                                       max_lfc_all = NULL)
  # Permutation test results per split
  permutation_test_results <- data.frame(split = NULL,
                                         runDE_n_sig = NULL,
                                         pvalue = NULL,
                                         n_iterations = NULL)
  # Full per-feature DE results for each iteration
  if (return_all == TRUE) {
    permutation_DE_results <- data.frame(feature = NULL,
                                         lfc = NULL,
                                         pvalue = NULL,
                                         padj = NULL,
                                         permutation = NULL,
                                         split = NULL)
  }
  # Permutation test group indices
  permutation_reference_group_indices <- vector(mode = "list", length = n_splits)
  names(permutation_reference_group_indices) <- use_splits

  if (proceed == TRUE) {
    # For each split
    for (s in 1:n_splits) {
      # Extract relevant values
      # Split
      current_split <- use_splits[s]
      # Matrix
      if (pseudobulk == "none") {
        # Grab feature x cell matrix
        current_mat <- input$cell_values[[current_split]]
      } else {
        # Grab PB matrix
        current_mat <- input$PB_values[[current_split]]
      }
      # Our replicates for DE will always be the column names of this matrix
      current_replicates <- colnames(current_mat)
      # Group key for permutation
      current_group_key <- group_key[current_replicates,] |>
        dplyr::select(-replicate) |>
        dplyr::group_by(dplyr::across(dplyr::everything())) |>
        dplyr::summarise(n = dplyr::n()) |>
        dplyr::select(-n) |>
        dplyr::arrange(permute_within) |>
        data.frame()
      # Stop if duplicates
      if (any(duplicated(current_group_key$permute_by))) {
        stop("Each value of parameter 'permute_by' must appear only once in each 'group' and each value of 'permute_within' (if applicable).")
      }
      rownames(current_group_key) <- current_group_key$permute_by

      # Establish permute_within order
      permute_within_ordered <- unique(current_group_key$permute_within)
      # Get permute_by order when arranged by permute_within
      permute_by_ordered <- current_group_key$permute_by
      # Number of replicates (ordered by permute_within order)
      n_replicates <- as.integer(table(current_group_key$permute_within)[permute_within_ordered])
      # Number in reference group (ordered by permute_within order)
      n_reference_group <- as.integer(table(current_group_key$permute_within[current_group_key$group == reference_group])[permute_within_ordered])
      # If input was not supplied to 'permute_by', reorder according to matrix column names
      if (is.null(permute_by)) {
        current_group_key <- current_group_key[current_replicates, ]
      }

      # Proceed if there are two groups present
      if (dplyr::n_distinct(current_group_key$group) != 2) {
        # Progress
        if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                             " : Skipping split ", current_split, ", ", dplyr::n_distinct(current_group_key$group), " group(s) present.")
      } else if (any(n_reference_group == 0) |
                 any((n_replicates - n_reference_group) == 0) |
                 (length(n_replicates) != length(n_reference_group))) {
        # Progress
        if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                             " : Skipping split ", current_split,
                             ", insufficient replicates per group within each value indicated by parameter 'permute_within'.")
      } else {
        # Establish whether there are any potential covariate confounders to check
        confound_check <- NULL
        if (!is.null(design_formula)) {
          # Get terms of design formula
          terms <- attr(terms(design_formula), "term.labels")
          # Remove last term (group)
          terms <- terms[-length(terms)]
          # Create interaction terms if necessary
          interaction_terms <- terms[grepl(":", terms)]
          check_terms <- terms[!grepl(":", terms)]
          if (length(interaction_terms) > 0) {
            for (i in 1:length(interaction_terms)) {
              # Split interaction term
              terms_i <- unlist(strsplit(interaction_terms[i], ":", fixed = TRUE))
              # Add column to current group key
              current_group_key[, paste0(c(terms_i, "interaction"), collapse = "_")] <- paste(current_group_key[,terms_i[1]],
                                                                                              current_group_key[,terms_i[2]],
                                                                                              sep = "_")
              # Add new term
              check_terms <- c(check_terms, paste0(c(terms_i, "interaction"), collapse = "_"))
            }
          }
          # Subset current group key to just columns that are terms of design formula
          confound_check <- current_group_key[, check_terms, drop = FALSE]
          # Subset to categorical only
          categorical_cols <- apply(confound_check, 2, FUN = function(i) {
            intersect(methods::is(i), c("character", "factor", "logical")) > 0
          })
          confound_check <- confound_check[, categorical_cols, drop = FALSE]
          # If there are any columns left, check those for confounds
          if (ncol(confound_check) > 0) {
            # If input was provided to permute_within, reorder
            if (!is.null(permute_within)) {
              confound_check <- confound_check[permute_by_ordered, , drop = FALSE]
            }
          } else {
            # Reset to NULL
            confound_check <- NULL
          }
        }
        # Compute a set of permuted group labels
        permuted_group_indices <- getCombinations(n_replicates = n_replicates,
                                                  n_group1 = n_reference_group,
                                                  n_combinations = n_iterations,
                                                  confound_check = confound_check,
                                                  message = paste0(" for split ", current_split, " (", s, "/", n_splits, ")"),
                                                  random_seed = random_seed,
                                                  verbose = TRUE)
        current_n_iterations <- ncol(permuted_group_indices)

        # Reorder if stratified by permute_within
        if (!is.null(permute_within)) {
          permuted_group_indices <- apply(permuted_group_indices, 2, FUN = function(i) {
            match(permute_by_ordered[i], current_group_key$permute_by)
          })
        }

        # Calculate overlap with true groups
        ref_group_overlap <- apply(permuted_group_indices, 2, FUN = function(i) {
          sum(current_group_key$group[i] == reference_group)/length(current_group_key$group[i])
        })
        non_ref_group_overlap <- apply(permuted_group_indices, 2, FUN = function(i) {
          sum(current_group_key$group[-i] == non_reference_group)/length(current_group_key$group[-i])
        })
        # If true labels were among random permutations, skip, leaving n_iterations - 1
        if (any(ref_group_overlap == 1 & non_ref_group_overlap == 1)) {
          true_index <- which(ref_group_overlap == 1 & non_ref_group_overlap == 1)
          permuted_group_indices <- permuted_group_indices[,-true_index]
          ref_group_overlap <- ref_group_overlap[-true_index]
          non_ref_group_overlap <- non_ref_group_overlap[-true_index]
        } else {
          # If not, remove one random set of labels, leaving n_iterations - 1
          permuted_group_indices <- permuted_group_indices[,-1]
          ref_group_overlap <- ref_group_overlap[-1]
          non_ref_group_overlap <- non_ref_group_overlap[-1]
        }

        # If input supplied to 'permute_by', expand permuted group indices to each replicate within a 'permute_by' set
        # (e.g., cells from the same sample)
        if (!is.null(permute_by)) {
          new_permuted_group_indices <- apply(permuted_group_indices, 2,
                                              FUN = function(i) {
                                                # Which permute_by values are in these indices
                                                permute_by_i <- current_group_key$permute_by[i]
                                                # Which replicate indices belong to these permute_by values
                                                current_replicate_indices <- which(group_key$permute_by %in% permute_by_i)
                                                return(current_replicate_indices)
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
        # Store permuted group indices
        permutation_reference_group_indices[[current_split]] <- permuted_group_indices

        # Progress
        if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                             " : Running ", current_n_iterations, " permutations for split ", current_split,
                             " (", s, "/", n_splits, ")..")
        # Run permutations
        permutation_DE_results_list <- pbmcapply::pbmclapply(seq_len(current_n_iterations-1),
                                                             FUN = function(i) {
                                                               # Create target dataframe by subsetting group_key to current replicates
                                                               targets_i <- group_key[current_replicates,]
                                                               # Switch to permuted group labels
                                                               targets_i$group <- non_reference_group
                                                               targets_i[permuted_group_indices[,i][!is.na(permuted_group_indices[,i])], "group"] <- reference_group

                                                               group_factor <- factor(targets_i$group)
                                                               group_factor <- stats::relevel(group_factor, ref = reference_group)
                                                               targets_i$group <- group_factor

                                                               # Create design
                                                               if (!is.null(design_formula)) {
                                                                 design_i <- stats::model.matrix(design_formula, data = targets_i)
                                                               } else {
                                                                 design_i <- stats::model.matrix(~ group, data = targets_i)
                                                               }

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
                                                                     max_lfc_sig = ifelse(n_sig > 0, max(lfc[padj < alpha & abs(lfc) > lfc_threshold]), NA),
                                                                     min_lfc_all = min(lfc),
                                                                     max_lfc_all = max(lfc)) |>
                                                                   data.frame() |>
                                                                   dplyr::mutate(reference_group_overlap = ref_group_overlap[i],
                                                                                 non_reference_group_overlap = non_ref_group_overlap[i]) |>
                                                                   dplyr::select(split,
                                                                                 permutation,
                                                                                 reference_group_overlap,
                                                                                 non_reference_group_overlap,
                                                                                 n_sig,
                                                                                 min_lfc_sig,
                                                                                 max_lfc_sig,
                                                                                 min_lfc_all,
                                                                                 max_lfc_all)
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
              max_lfc_sig = ifelse(n_sig > 0, max(lfc[padj < alpha & abs(lfc) > lfc_threshold]), NA),
              min_lfc_all = min(lfc),
              max_lfc_all = max(lfc)) |>
            data.frame() |>
            dplyr::mutate(reference_group_overlap = ref_group_overlap,
                          non_reference_group_overlap = non_ref_group_overlap) |>
            dplyr::select(split,
                          permutation,
                          reference_group_overlap,
                          non_reference_group_overlap,
                          n_sig,
                          min_lfc_sig,
                          max_lfc_sig,
                          min_lfc_all,
                          max_lfc_all)
        } else {
          permutation_DE_summary_s <- do.call(rbind, permutation_DE_results_list) |> data.frame()
        }

        # Conduct permutation test
        runDE_n_sig <- dplyr::filter(runDE_values, split == current_split)$runDE_n_sig

        permutation_test_p_value <- sum(c(runDE_n_sig, permutation_DE_summary_s$n_sig) >= runDE_n_sig)/current_n_iterations
        permutation_test_results_s <- data.frame(split = current_split,
                                                 runDE_n_sig = runDE_n_sig,
                                                 pvalue = permutation_test_p_value,
                                                 n_iterations = current_n_iterations)
        # Add to overall results
        permutation_DE_summary <- rbind(permutation_DE_summary, permutation_DE_summary_s)
        permutation_test_results <- rbind(permutation_test_results, permutation_test_results_s)
        if (return_all == TRUE) {
          permutation_DE_results <- rbind(permutation_DE_results, permutation_DE_results_s)
        }
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Wrap up
  # ---------------------------------------------------------------------------

  time3 <- Sys.time()

  # Metadata
  metadata_list <- list("runDE_values" = runDE_values,
                        "permutation_reference_group_indices" = permutation_reference_group_indices,
                        "time" = data.frame(total = difftime(time3, time1, units = "secs"),
                                            step1_setup = difftime(time2, time1, units = "secs"),
                                            step2_permute = difftime(time3, time2, units = "secs")))

  parameter_list <- list("alpha" = alpha,
                         "lfc_threshold" = lfc_threshold,
                         "n_iterations" = n_iterations,
                         "use_splits" = use_splits,
                         "permute_by" = permute_by,
                         "permute_within" = permute_within,
                         "min_DE" = min_DE,
                         "reference_group" = reference_group,
                         "non_reference_group" = non_reference_group,
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
