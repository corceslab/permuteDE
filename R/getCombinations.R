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
#' provided 'Seurat' or 'SingleCellExperiment' object.
#'
#' @param object An optional 'Seurat' or 'SingleCellExperiment' object. If
#' \code{NULL}, \code{countCombinations} will expect either (1) vector input to
#' parameters \code{replicate_label} and \code{group_label} or (2) numeric input
#' to parameters \code{n_replicates} and \code{n_group1}.
#' @param replicate_label A string indicating the name of the
#' metadata column containing the replicate labels or a character vector
#' containing the replicate labels in order.
#' @param group_labels A string indicating the name of the
#' column containing the two comparison group labels or a character vector
#' containing the comparison labels in order.
#' @param use_cells A vector of cell names to subset prior to calculating
#' possible group comibinations. Default = \code{NULL} will use all cells.
#' @param n_replicates A numeric value indicating the total number of
#' replicates. Defaults to \code{NULL}.
#' @param n_group1 A numeric value indicating the number of replicates in one
#' group (doesn't matter which). Defaults to \code{NULL}.
#' @param n_combinations A numeric value indicating the number of combinations
#' to generate. Defaults to 1000.
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
                            replicate_labels = NULL,
                            group_labels = NULL,
                            use_cells = NULL,
                            n_replicates = NULL,
                            n_group1 = NULL,
                            n_combinations = 1000,
                            message = "",
                            random_seed = 1,
                            verbose = TRUE) {
  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(object, "object", "getCombinations")
  .validInput(replicate_labels, "replicate_labels", list(object, "not applicable"))
  .validInput(group_labels, "group_labels", object)
  .validInput(use_cells, "use_cells", list(object, "not applicable"))
  .validInput(n_replicates, "n_replicates")
  .validInput(n_group1, "n_group1", n_replicates)
  .validInput(n_combinations, "n_combinations")
  .validInput(message, "message")
  .validInput(random_seed, "random_seed")
  .validInput(verbose, "verbose")

  # Set seed
  set.seed(random_seed)

  # Either use just 'n_replicates' and 'n_group1', use vectors, or use object and column labels
  if (!is.null(n_replicates) & !is.null(n_group1)) {
    # use 'n_replicates' and 'n_group1'
    if (any(!is.null(object), !is.null(replicate_labels), !is.null(group_labels), !is.null(use_cells))) {
      warning("When inputs for 'n_replicates' and 'n_group1' are provided, inputs for 'object', 'replicate_labels', 'group_labels', and 'use_cells' are not used.")
    }
  } else {
    # use 'object', 'replicate_labels', 'group_labels'
    if (any(is.null(replicate_labels), is.null(group_labels))) {
      stop("If not using 'n_replicates' and 'n_group1', input must be provided to 'replicate_labels' and 'group_labels'.")
    }
    if (!is.null(n_replicates)) {
      warning("Input for 'n_replicates' was not used.")
    }
    if (!is.null(n_group1)) {
      warning("Input for 'n_group1' was not used.")
    }

    # retrieve metadata

    if (length(replicate_labels) != length(group_labels)) {
      stop("Input to parameters 'replicate_labels' and 'group_labels' must be of the same length. Please supply valid input!")
    }

    if (length(replicate_labels) == 1) {
      if (is.null(object)) {
        stop("When input for 'replicate_labels' and 'group_labels' are single values, input must be provided to parameter 'object'.")
      }
      replicates <- .retrieveData(object = object,
                                  type = "cell_metadata",
                                  name = replicate_labels,
                                  use_cells = use_cells)
      groups <- .retrieveData(object = object,
                              type = "cell_metadata",
                              name = group_labels,
                              use_cells = use_cells)
    } else {
      if (!is.null(object)) {
        warning("When input for 'replicate_labels' and 'group_labels' are vectors, input to parameter 'object' is not used.")
      }
      replicates <- replicate_labels
      groups <- group_labels
    }

    # count unique replicates and unique replicates in one group
    n_replicates <- dplyr::n_distinct(replicates)
    n_group1 <- dplyr::n_distinct(replicates[groups == groups[1]])
  }

  # Check whether requested number of combinations is <= possible number of combinations
  n_possible_combinations <- choose(n_replicates, n_group1)

  if (n_possible_combinations < n_combinations) {
    warning("Input ", n_combinations,
            " for parameter 'n_combinations' exceeds the number of possible combinations ",
            n_possible_combinations, ". Only ", n_possible_combinations,
            " combinations will be generated.")
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
    output_combinations <- all_combinations[, sample(1:ncol(all_combinations), n_combinations)]
  } else {
    # Set up collection of combinations
    output_combinations <- vector("list", n_combinations)
    count <- 0
    # Set up environment to track previously encountered combinations
    seen_combinations <- new.env(hash = TRUE, size = n_combinations)
    # Generate new unique combinations
    while (count < n_combinations) {
      possible_comb <- sample(n_replicates, n_group1)
      comb_key <- paste(sort(possible_comb), collapse = "-")
      if (is.null(seen_combinations[[comb_key]])) {
        seen_combinations[[comb_key]] <- TRUE
        count <- count + 1
        output_combinations[[count]] <- possible_comb
      }
    }
    # Simplify combinations to a matrix
    output_combinations <- do.call(cbind, output_combinations)
  }

  # Return matrix of combinations
  return(output_combinations)
}
