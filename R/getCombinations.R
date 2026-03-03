#' Generate combinations
#'
#' For a set of replicates that belong to two groups, this function will
#' generate a set of randomized group combinations.
#'
#' Users may provide input in three ways:
#'
#' (1) A vector of replicate labels and a vector of group labels (in order)
#'
#' (2) The total number of replicates (n) and the number of
#' replicates in one group (k) (doesn't matter which group).
#'
#' (3) Column names indicating replicate and group metadata columns in a
#' provided 'Seurat' or 'SingleCellExperiment' object or a provided
#' metadata dataframe.
#'
#' @param object An optional 'Seurat' or 'SingleCellExperiment' object. If
#' \code{NULL}, \code{countCombinations} will expect either (1) vector input to
#' parameters \code{replicate_label} and \code{group_label} or (2) numeric input
#' to parameters \code{n_replicates} and \code{n_group1}.
#' @param metadata An optional dataframe containing relevant metadata columns
#' corresponding to the data provided to parameter \code{object}. Default =
#' \code{NULL} looks for metadata in \code{object} or other provided inputs.
#' @param replicate_labels A string indicating the name of the
#' metadata column containing the biological replicate labels or a character
#' vector containing the biological replicate labels in order.
#' @param group_labels A string indicating the name of the
#' column containing the two comparison group labels or a character vector
#' containing the comparison labels in order.
#' @param use_cells A vector of cell names to subset prior to calculating
#' possible group comibinations. Default = \code{NULL} will use all cells.
#' @param n_replicates A numeric value indicating the total number of
#' replicates. Alternately, a vector can be provided to generate combinations
#' when shuffling separately within multiple sets. Defaults to \code{NULL}.
#' @param n_group1 A numeric value indicating the number of replicates in one
#' group (doesn't matter which). Alternately, a vector can be provided to
#' generate combinations when shuffling separately within multiple sets.
#' Defaults to \code{NULL}.
#' @param n_combinations A numeric value indicating the number of combinations
#' to generate. Defaults to 1000.
#' @param confound_check An optional dataframe of covariates for which to
#' exclude permutations that are perfectly confounded. Defaults to \code{NULL}.
#' @param message A character string indicating additional progress messaging
#' (internal use). Defaults to "".
#' @param random_seed A numeric value indicating the random seed to be used.
#' Defaults to 1.
#' @param verbose A Boolean value indicating whether to use verbose output
#' during the execution of this function. Defaults to \code{TRUE}.
#' Can be set to \code{FALSE} for a cleaner output.
#'
#' @return Returns a matrix where each column contains a combination of index
#' values indicating which replicates to assign to the first group.
#'
#' @export
#'
getCombinations <- function(object = NULL,
                            metadata = NULL,
                            replicate_labels = NULL,
                            group_labels = NULL,
                            use_cells = NULL,
                            n_replicates = NULL,
                            n_group1 = NULL,
                            n_combinations = 1000,
                            confound_check = NULL,
                            message = "",
                            random_seed = 1,
                            verbose = TRUE) {
  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(object, "object", "getCombinations")
  .validInput(metadata, "metadata", object)
  .validInput(replicate_labels, "replicate_labels", list(object, metadata, "not applicable"))
  .validInput(group_labels, "group_labels", list(object, metadata))
  .validInput(use_cells, "use_cells", list(object, "not applicable"))
  .validInput(n_replicates, "n_replicates")
  .validInput(n_group1, "n_group1", n_replicates)
  .validInput(n_combinations, "n_combinations")
  .validInput(message, "message")
  .validInput(random_seed, "random_seed")
  .validInput(verbose, "verbose")

  # Set seed
  set.seed(random_seed)

  # ---------------------------------------------------------------------------
  # Retrieve values
  # ---------------------------------------------------------------------------

  # Either use just 'n_replicates' and 'n_group1', use vectors, or use object and column labels
  if (!is.null(n_replicates) & !is.null(n_group1)) {
    # use 'n_replicates' and 'n_group1'
    if (any(!is.null(object), !is.null(replicate_labels), !is.null(group_labels), !is.null(use_cells))) {
      warning(" When inputs for 'n_replicates' and 'n_group1' are provided, inputs for 'object', 'replicate_labels', 'group_labels', and 'use_cells' are not used.")
    }
  } else {
    # use 'object', 'replicate_labels', 'group_labels'
    if (any(is.null(replicate_labels), is.null(group_labels))) {
      stop("If not using 'n_replicates' and 'n_group1', input must be provided to 'replicate_labels' and 'group_labels'.")
    }
    if (!is.null(n_replicates)) {
      warning(" Input for 'n_replicates' was not used.")
    }
    if (!is.null(n_group1)) {
      warning(" Input for 'n_group1' was not used.")
    }

    # Retrieve metadata
    if (length(replicate_labels) != length(group_labels)) {
      stop("Input to parameters 'replicate_labels' and 'group_labels' must be of the same length. Please supply valid input!")
    }

    if (length(replicate_labels) == 1) {
      if (is.null(object)) {
        stop("When input for 'replicate_labels' and 'group_labels' are single values, input must be provided to parameter 'object'.")
      }
      replicates <- .retrieveData(object = object,
                                  metadata = metadata,
                                  type = "cell_metadata",
                                  name = replicate_labels,
                                  use_cells = use_cells)
      groups <- .retrieveData(object = object,
                              metadata = metadata,
                              type = "cell_metadata",
                              name = group_labels,
                              use_cells = use_cells)
    } else {
      if (!is.null(object)) {
        warning(" When input for 'replicate_labels' and 'group_labels' are vectors, input to parameter 'object' is not used.")
      }
      replicates <- replicate_labels
      groups <- group_labels
    }
    # Check for NA values
    if (any(is.na(replicates))) {
      stop("Values provided for 'replicate_labels' cannot be NA.")
    }
    if (any(is.na(groups))) {
      stop("Values provided for 'group_labels' cannot be NA.")
    }
    # Check that there are exactly two groups present
    if (dplyr::n_distinct(groups) != 2) {
      stop("Input value for parameter 'group_labels' must represent a cell metadata column (or a vector of group labels) that contains exactly 2 groups for the selected data, please supply valid input!")
    }
    # Count unique replicates and unique replicates in one group
    n_replicates <- dplyr::n_distinct(replicates)
    n_group1 <- dplyr::n_distinct(replicates[groups == groups[1]])
  }

  # Additional validation
  .validInput(confound_check, "confound_check", sum(n_replicates))

  # ---------------------------------------------------------------------------
  # If n_replicates and n_group1 contain multiple values, send to sub function
  # ---------------------------------------------------------------------------

  if (length(n_replicates) > 1) {
    output_combinations <- .getStratifiedCombinations(n_replicates = n_replicates,
                                                      n_group1 = n_group1,
                                                      n_combinations = n_combinations,
                                                      confound_check = confound_check,
                                                      message = message,
                                                      verbose = verbose)
  } else {
    # Check whether possible number of combinations is < requested number of combinations
    n_possible_combinations <- choose(n_replicates, n_group1)

    if (n_possible_combinations < n_combinations) {
      if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Input ", n_combinations,
                           " for parameter 'n_combinations' exceeds the number of possible combinations ",
                           n_possible_combinations,". Only ", n_possible_combinations,
                           " combinations will be generated", message, ".")
      n_combinations <- n_possible_combinations
    }

    # ---------------------------------------------------------------------------
    # Generate combinations
    # ---------------------------------------------------------------------------

    # Progress
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Generating ",
                         n_combinations, " of ", n_possible_combinations,
                         " possible combinations", message, "..")

    # When efficient, generate all combinations, then sample
    if (n_possible_combinations < 1000000) {
      all_combinations <- utils::combn(n_replicates, n_group1)
      # If checking for covariate confounds, remove any perfectly confounded combinations before sampling
      if (!is.null(confound_check)) {
        # Nested apply calls to check each combination against each covariate
        keep <- unlist(apply(all_combinations, 2, FUN = function(i) {
          # Convert combination to a 0,1 matrix
          group_i <- rep(1, n_replicates)
          group_i[i] <- 0
          # Test each covariate
          exclude <- apply(confound_check, 2, FUN = function(j) {
            all(colSums(table(group_i, droplevels(factor(j))) > 0) == 1)
          })
          # Retain combination?
          return(!any(exclude))
        }))
        all_combinations <- all_combinations[, keep]
        # Report if any were excluded
        if (ncol(all_combinations) < n_possible_combinations) {
          # Report if number of possible combinations is now lower than what is requested
          if (ncol(all_combinations) < n_combinations) {
            if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Excluding ",
                                 (n_possible_combinations - ncol(all_combinations)), " of ", n_possible_combinations,
                                 " possible combinations due to covariate confounds. Only ", ncol(all_combinations),
                                 " combinations will be generated..")
            n_combinations <- ncol(all_combinations)
          } else if (verbose) {
            message(format(Sys.time(), "%Y-%m-%d %X"), " : Excluding ",
                    (n_possible_combinations - ncol(all_combinations)), " of ", n_possible_combinations,
                    " possible combinations due to covariate confounds..")
          }
        }
      }
      output_combinations <- all_combinations[, sample(1:ncol(all_combinations), n_combinations)]
    } else {
      # Set up collection of combinations
      output_combinations <- vector("list", n_combinations)
      count <- 0
      n_excluded_confound_check <- 0
      # Set up environment to track previously encountered combinations
      seen_combinations <- new.env(hash = TRUE, size = n_combinations)
      # Generate new unique combinations
      while (count < n_combinations) {
        possible_comb <- sample(n_replicates, n_group1)
        # Check if unique
        comb_key <- paste(sort(possible_comb), collapse = "-")
        if (is.null(seen_combinations[[comb_key]])) {
          # Add to seen combinations
          seen_combinations[[comb_key]] <- TRUE

          # If checking for covariate confounds, remove if perfectly confounded
          if (!is.null(confound_check)) {
              # Convert combination to a 0,1 matrix
              group_i <- rep(1, n_replicates)
              group_i[possible_comb] <- 0
              # Test each covariate
              exclude <- apply(confound_check, 2, FUN = function(j) {
                all(colSums(table(group_i, droplevels(factor(j))) > 0) == 1)
              })
              # Retain combination?
              if (any(exclude)) {
                keep <- FALSE
                n_excluded_confound_check <- n_excluded_confound_check + 1
              } else {
                keep <- TRUE
              }
          } else {
            keep <- TRUE
          }

          # Add to output
          if (keep == TRUE) {
            count <- count + 1
            output_combinations[[count]] <- possible_comb
          }

          # Avoid forever loop
          if (length(seen_combinations) > 10*n_combinations) {
            # Stop, even though we have not reached n combinations yet
            count <- n_combinations
            # Report issue
            if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                                 " : >90% of tested combinations have been excluded due to covariate confounds. Could not find ",
                                 n_combinations, " unique combinations. Stopping search at ", length(output_combinations),
                                 " to prevent a potential infinite loop..")
          }
        }
      }
      # Report if any were excluded due to covariate confounds
      if (verbose & (n_excluded_confound_check > 0)) {
        message(format(Sys.time(), "%Y-%m-%d %X"), " : ",
                n_excluded_confound_check, " initially sampled combination(s) excluded due to covariate confounds..")
      }
      # Simplify combinations to a matrix
      output_combinations <- do.call(cbind, output_combinations)
    }
  }
  # Return matrix of combinations
  return(output_combinations)
}

# Generate combinations when shuffling separately within each partition ---------------------------
#
# n_replicates -- Vector of total number of replicates within each partition
# n_group1 -- Vector of number of replicates in group 1 within each partition
# n_combinations -- A numeric value indicating the number of combinations to generate.
# confound_check -- Optional dataframe for checking covariate confounds
# message -- Optional string to report message
# verbose -- Logical value to indicate verbose output

.getStratifiedCombinations <- function(n_replicates,
                                       n_group1,
                                       n_combinations = 1000,
                                       confound_check = NULL,
                                       message = "",
                                       verbose = TRUE) {
  # Number of partitions
  n_partitions <- length(n_replicates)

  # Total possible combinations across partitions (log)
  log_n_possible_combinations <- sum(lchoose(n_replicates, n_group1))
  if (log_n_possible_combinations < log(.Machine$double.xmax)) {
    # Safe to compute real number
    n_possible_combinations <- prod(choose(n_replicates, n_group1))
  } else {
    # Effectively
    n_possible_combinations <- Inf
  }

  # Check whether possible number of combinations is < requested number of combinations
  if (n_possible_combinations < n_combinations) {
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Input ", n_combinations,
                         " for parameter 'n_combinations' exceeds the number of possible combinations ",
                         n_possible_combinations,". Only ", n_possible_combinations,
                         " combinations will be generated", message, ".")
    n_combinations <- n_possible_combinations
  }

  # ---------------------------------------------------------------------------
  # Generate combinations
  # ---------------------------------------------------------------------------

  if (verbose) {
    if (is.finite(n_possible_combinations)) {
      message(format(Sys.time(), "%Y-%m-%d %X"), " : Generating ",
              n_combinations, " of ", n_possible_combinations,
              " possible combinations", message, "..")
    } else {
      message(format(Sys.time(), "%Y-%m-%d %X"), " : Generating ",
              n_combinations, " from a very large set of possible combinations",
              message, "..")
    }
  }

  # Build global index ranges for each partition
  partition_indices <- split(seq_len(sum(n_replicates)),
                             rep.int(seq_along(n_replicates), n_replicates))

  # Decide whether to enumerate all combinations within each partition
  log_n_combinations_per_partition <- lchoose(n_replicates, n_group1)
  generate_all <- log_n_combinations_per_partition <= log(1000000)

  # Precompute per-partition combination sets when feasible
  combinations_by_partition <- vector("list", n_partitions)
  for (p in seq_len(n_partitions)) {
    indices_p <- partition_indices[[p]]
    n_group1_p <- n_group1[[p]]
    if (generate_all[[p]] == TRUE) {
      # Generate a list of all combinations in partition p
      combinations_by_partition[[p]] <- utils::combn(indices_p, n_group1_p, simplify = FALSE)
    } else {
      # Too large to generate all
      combinations_by_partition[[p]] <- NULL
    }
  }

  # Set up collection of combinations
  output_combinations <- vector("list", n_combinations)
  count <- 0
  n_excluded_confound_check <- 0
  # Set up environment to track previously encountered combinations
  seen_combinations <- new.env(hash = TRUE, size = n_combinations)
  # Generate new unique combinations
  while (count < n_combinations) {
    possible_comb <- unlist(lapply(seq_along(partition_indices),
                                   FUN = function(p) {
                                     # Grab pre-generated combinations for partition p (if available)
                                     combinations_p <- combinations_by_partition[[p]]
                                     if (!is.null(combinations_p)) {
                                       # Sample from pre-generated combinations for partition p
                                       combinations_p[[sample.int(length(combinations_p), 1L)]]
                                     } else {
                                       # When not pre-generated, generate random combination for partition p
                                       sample(partition_indices[[p]], n_group1[[p]], replace = FALSE)
                                     }
                                   }))
    # Check if unique
    comb_key <- paste(sort(possible_comb), collapse = "-")
    if (is.null(seen_combinations[[comb_key]])) {
      seen_combinations[[comb_key]] <- TRUE

      # If checking for covariate confounds, remove if perfectly confounded
      if (!is.null(confound_check)) {
        # Convert combination to a 0,1 matrix
        group_i <- rep(1, sum(n_replicates))
        group_i[possible_comb] <- 0
        # Test each covariate
        exclude <- apply(confound_check, 2, FUN = function(j) {
          all(colSums(table(group_i, droplevels(factor(j))) > 0) == 1)
        })
        # Retain combination?
        if (any(exclude)) {
          keep <- FALSE
          n_excluded_confound_check <- n_excluded_confound_check + 1
        } else {
          keep <- TRUE
        }
      } else {
        keep <- TRUE
      }

      # Add to output
      if (keep == TRUE) {
        count <- count + 1
        output_combinations[[count]] <- possible_comb
      }

      # Avoid forever loop
      if (length(seen_combinations) > 10*n_combinations) {
        # Stop, even though we have not reached n combinations yet
        count <- n_combinations
        # Report issue
        if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"),
                             " : Excluded >90% of tested combinations so far due to covariate confounds. Stopping search at ",
                             length(output_combinations), " unique combinations to prevent a potential infinite loop..")
      }
    }
  }

  # Report if any were excluded due to covariate confounds
  if (verbose & (n_excluded_confound_check > 0)) {
    message(format(Sys.time(), "%Y-%m-%d %X"), " : Excluded ",
            n_excluded_confound_check, " initially sampled combination(s) due to covariate confounds..")
  }

  # Simplify combinations to a matrix
  output_combinations <- do.call(cbind, output_combinations)

  # Return matrix of combinations
  return(output_combinations)
}
