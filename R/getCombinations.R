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
#' @param progress_message A character string indicating additional progress messaging
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
                            progress_message = "",
                            random_seed = 1,
                            verbose = TRUE) {
  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(input = object,
              name = "object",
              null_allowed = TRUE,
              class = c("Seurat", "SingleCellExperiment", "matrix", "Matrix"))
  .validInput(input = metadata,
              name = "metadata",
              null_allowed = TRUE,
              class = "data.frame",
              other = object)
  .validInput(input = replicate_labels,
              name = "replicate_labels",
              null_allowed = TRUE,
              class = c("character", "factor", "numeric"),
              other = list(object, metadata, "not applicable"))
  .validInput(input = group_labels,
              name = "group_labels",
              null_allowed = TRUE,
              class = c("character", "factor", "numeric", "logical"),
              other = list(object, metadata))
  .validInput(input = use_cells,
              name = "use_cells",
              null_allowed = TRUE,
              class = "character",
              other = list(object, "not applicable"))
  .validInput(input = n_replicates,
              name = "n_replicates",
              null_allowed = TRUE,
              class = "numeric")
  .validInput(input = n_group1,
              name = "n_group1",
              null_allowed = TRUE,
              class = "numeric",
              other = n_replicates)
  .validInput(input = n_combinations,
              name = "n_combinations",
              class = "numeric",
              len = 1)
  .validInput(input = progress_message,
              name = "progress_message",
              class = "character",
              len = 1)
  .validInput(input = random_seed,
              name = "random_seed",
              class = "numeric",
              len = 1)
  .validInput(input = verbose,
              name = "verbose",
              class = "logical",
              len = 1)

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
    # Check that each replicate belongs to exactly one group
    problem_replicates <- data.frame(replicate = replicates,
                                     group = groups) |>
      dplyr::summarise(n_groups = dplyr::n_distinct(group), .by = replicate) |>
      dplyr::filter(n_groups > 1) |>
      dplyr::pull(replicate)
    if (length(problem_replicates) > 0) {
      stop("Each replicate must belong to exactly one group. ", "The following replicate label(s) appear in multiple groups: ",
           paste(problem_replicates, collapse = ", "), ".")
    }
    # Count unique replicates and unique replicates in one group
    n_replicates <- dplyr::n_distinct(replicates)
    n_group1 <- dplyr::n_distinct(replicates[groups == groups[1]])
  }

  # Additional validation
  .validInput(input = confound_check,
              name = "confound_check",
              null_allowed = TRUE,
              class = "data.frame",
              other = sum(n_replicates))

  # ---------------------------------------------------------------------------
  # If n_replicates and n_group1 contain multiple values, send to sub function
  # ---------------------------------------------------------------------------

  if (length(n_replicates) > 1) {
    output_combinations <- .getStratifiedCombinations(n_replicates = n_replicates,
                                                      n_group1 = n_group1,
                                                      n_combinations = n_combinations,
                                                      confound_check = confound_check,
                                                      progress_message = progress_message,
                                                      verbose = verbose)
  } else {
    # Check whether possible number of combinations is < requested number of combinations
    n_possible_combinations <- choose(n_replicates, n_group1)

    if (n_possible_combinations < n_combinations) {
      if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Input ", n_combinations,
                           " for parameter 'n_combinations' exceeds the number of possible combinations ",
                           n_possible_combinations,". Only ", n_possible_combinations,
                           " combinations will be generated", progress_message, ".")
      n_combinations <- n_possible_combinations
    }

    # ---------------------------------------------------------------------------
    # Generate combinations
    # ---------------------------------------------------------------------------

    # Progress
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Generating ",
                         n_combinations, " of ", n_possible_combinations,
                         " possible combinations", progress_message, "..")

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
        all_combinations <- all_combinations[, keep, drop = FALSE]
        if (ncol(all_combinations) == 0) {
          stop("No valid combinations remain after excluding covariate-confounded combinations.")
        }
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
      output_combinations <- all_combinations[, sample(1:ncol(all_combinations), n_combinations), drop = FALSE]
    } else {
      # Set up collection of combinations
      output_combinations <- vector("list", n_combinations)
      count <- 0
      n_excluded_confound_check <- 0

      # Track how many candidate combinations have been tested
      n_tested <- 0

      # Set up environment to track previously encountered combinations
      seen_combinations <- new.env(hash = TRUE, parent = emptyenv())

      # Generate new unique combinations
      while (count < n_combinations) {
        # Avoid forever loop
        if (n_tested > 10 * n_combinations) {
          if (verbose) {
            message(
              format(Sys.time(), "%Y-%m-%d %X"),
              " : >90% of tested combinations have been excluded due to covariate confounds. Could not find ",
              n_combinations, " unique combinations. Stopping search at ", count,
              " to prevent a potential infinite loop.."
            )
          }
          break
        }

        possible_comb <- sample(n_replicates, n_group1)

        # Compact key for uniqueness check
        # If you have digest installed, this is better than pasting all IDs together:
        comb_key <- digest::digest(sort(possible_comb), algo = "xxhash64")

        # Check if unique
        if (!exists(comb_key, envir = seen_combinations, inherits = FALSE)) {
          assign(comb_key, TRUE, envir = seen_combinations)
          n_tested <- n_tested + 1

          # If checking for covariate confounds, remove if perfectly confounded
          if (!is.null(confound_check)) {
            group_i <- rep(1, n_replicates)
            group_i[possible_comb] <- 0

            exclude <- apply(confound_check, 2, FUN = function(j) {
              all(colSums(table(group_i, droplevels(factor(j))) > 0) == 1)
            })

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
          if (keep) {
            count <- count + 1
            output_combinations[[count]] <- possible_comb
          }
        }
      }

      # Report if any were excluded due to covariate confounds
      if (verbose && n_excluded_confound_check > 0) {
        message(
          format(Sys.time(), "%Y-%m-%d %X"), " : ",
          n_excluded_confound_check, " initially sampled combination(s) excluded due to covariate confounds.."
        )
      }

      # Drop any NULLs if we exited early
      output_combinations <- output_combinations[seq_len(count)]

      # If count = 0, stop
      if (count == 0) {
        stop("No valid combinations remain after excluding covariate-confounded combinations.")
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
# progress_message -- Optional string to report message
# verbose -- Logical value to indicate verbose output

.getStratifiedCombinations <- function(n_replicates,
                                       n_group1,
                                       n_combinations = 1000,
                                       confound_check = NULL,
                                       progress_message = "",
                                       verbose = TRUE) {
  # Number of partitions
  n_partitions <- length(n_replicates)

  # Total possible combinations across partitions (log)
  log_n_possible_combinations <- sum(lchoose(n_replicates, n_group1))
  if (log_n_possible_combinations < log(.Machine$double.xmax)) {
    n_possible_combinations <- prod(choose(n_replicates, n_group1))
  } else {
    n_possible_combinations <- Inf
  }

  # Check whether possible number of combinations is < requested number of combinations
  if (is.finite(n_possible_combinations) && n_possible_combinations < n_combinations) {
    if (verbose) {
      message(
        format(Sys.time(), "%Y-%m-%d %X"), " : Input ", n_combinations,
        " for parameter 'n_combinations' exceeds the number of possible combinations ",
        n_possible_combinations, ". Only ", n_possible_combinations,
        " combinations will be generated", progress_message, "."
      )
    }
    n_combinations <- n_possible_combinations
  }

  # ---------------------------------------------------------------------------
  # Generate combinations
  # ---------------------------------------------------------------------------

  if (verbose) {
    if (is.finite(n_possible_combinations)) {
      message(
        format(Sys.time(), "%Y-%m-%d %X"), " : Generating ",
        n_combinations, " of ", n_possible_combinations,
        " possible combinations", progress_message, ".."
      )
    } else {
      message(
        format(Sys.time(), "%Y-%m-%d %X"), " : Generating ",
        n_combinations, " from a very large set of possible combinations",
        progress_message, ".."
      )
    }
  }

  # Build global index ranges for each partition
  partition_indices <- split(
    seq_len(sum(n_replicates)),
    rep.int(seq_along(n_replicates), n_replicates)
  )

  # Decide whether to enumerate all combinations within each partition
  log_n_combinations_per_partition <- lchoose(n_replicates, n_group1)
  generate_all <- log_n_combinations_per_partition <= log(1000000)

  # Precompute per-partition combination sets when feasible
  combinations_by_partition <- vector("list", n_partitions)
  for (p in seq_len(n_partitions)) {
    indices_p <- partition_indices[[p]]
    n_group1_p <- n_group1[[p]]

    if (isTRUE(generate_all[p])) {
      combinations_by_partition[[p]] <- utils::combn(
        indices_p,
        n_group1_p,
        simplify = FALSE
      )
    } else {
      combinations_by_partition[[p]] <- NULL
    }
  }

  # Set up collection of combinations
  output_combinations <- vector("list", n_combinations)
  count <- 0L
  n_excluded_confound_check <- 0L

  # Set up environment to track previously encountered combinations
  seen_combinations <- new.env(hash = TRUE, parent = emptyenv())

  # Generate new unique combinations
  while (count < n_combinations) {
    possible_comb <- unlist(
      lapply(seq_len(n_partitions), function(p) {
        combinations_p <- combinations_by_partition[[p]]

        if (!is.null(combinations_p)) {
          combinations_p[[sample.int(length(combinations_p), 1L)]]
        } else {
          sample(partition_indices[[p]], n_group1[[p]], replace = FALSE)
        }
      }),
      use.names = FALSE
    )

    # Compact uniqueness key
    comb_key <- digest::digest(sort(possible_comb), algo = "xxhash64")

    # Check if unique
    if (!exists(comb_key, envir = seen_combinations, inherits = FALSE)) {
      assign(comb_key, TRUE, envir = seen_combinations)

      # If checking for covariate confounds, remove if perfectly confounded
      if (!is.null(confound_check)) {
        group_i <- rep(1, sum(n_replicates))
        group_i[possible_comb] <- 0

        exclude <- apply(confound_check, 2, FUN = function(j) {
          all(colSums(table(group_i, droplevels(factor(j))) > 0) == 1)
        })

        if (any(exclude)) {
          keep <- FALSE
          n_excluded_confound_check <- n_excluded_confound_check + 1L
        } else {
          keep <- TRUE
        }
      } else {
        keep <- TRUE
      }

      # Add to output
      if (keep) {
        count <- count + 1L
        output_combinations[[count]] <- possible_comb
      }

      # Avoid forever loop
      if (length(ls(envir = seen_combinations, all.names = TRUE)) > 10 * n_combinations) {
        if (verbose) {
          message(
            format(Sys.time(), "%Y-%m-%d %X"),
            " : Excluded >90% of tested combinations so far due to covariate confounds. Stopping search at ",
            count, " unique combinations to prevent a potential infinite loop.."
          )
        }
        break
      }
    }
  }

  # Report if any were excluded due to covariate confounds
  if (verbose && n_excluded_confound_check > 0) {
    message(
      format(Sys.time(), "%Y-%m-%d %X"), " : Excluded ",
      n_excluded_confound_check, " initially sampled combination(s) due to covariate confounds.."
    )
  }

  # Drop any NULLs if we exited early
  output_combinations <- output_combinations[seq_len(count)]

  # Simplify combinations to a matrix
  output_combinations <- do.call(cbind, output_combinations)

  # Return matrix of combinations
  return(output_combinations)
}
