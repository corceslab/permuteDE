#' Generate combinations
#'
#' For a set of replicates that belong to two groups, this function will generate
#' a set of randomized group combinations. The input should be either the replicate
#' and group column names from an object's metadata, or simply the total number
#' of replicates (n) and the number of replicates in one of the two groups (k).
#'
#' @param object An object of class 'Seurat' or 'SingleCellExperiment'.
#' @param replicate_labels A character string or vector indicating the name of the
#' column containing the replicate labels. Defaults to \code{NULL}.
#' @param group_labels A character string or vector indicating the name of the
#' column containing the two comparison group labels. Defaults to \code{NULL}.
#' @param use_cells A vector of cell names subset to. Default = \code{NULL} will
#' use all cells.
#' @param n_replicates A numeric value indicating the total number of replicates.
#' Defaults to \code{NULL}.
#' @param n_group1 A numeric value indicating the number of replicates that belong
#' to one of the two groups (doesn't matter which). Defaults to \code{NULL}.
#' @param n_combinations A numeric value indicating the number of combinations
#' to generate. Defaults to 1000.
#' @param message A character string indicating additional progress messaging
#' (internal use). Defaults to "".
#' @param random_seed A numeric value indicating the random seed to be used.
#' Defaults to 1.
#' @param verbose A boolean value indicating whether to use verbose output
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

  .validInput(object, "object", "countCombinations")
  .validInput(replicate_labels, "replicate_labels", object)
  .validInput(group_labels, "group_labels", object)
  .validInput(use_cells, "use_cells", object)
  .validInput(n_replicates, "n_replicates")
  .validInput(n_group1, "n_group1", n_replicates)
  .validInput(n_combinations, "n_combinations")
  .validInput(message, "message")
  .validInput(random_seed, "random_seed")
  .validInput(verbose, "verbose")

  # Set seed
  set.seed(random_seed)

  # Either use just 'n_replicates' and 'n_group1' or use object and column labels
  if (!is.null(n_replicates) & !is.null(n_group1)) {
    if (any(!is.null(object), !is.null(replicate_labels), !is.null(group_labels), !is.null(use_cells))) {
      warning("When inputs for parameters 'n_replicates' and 'n_group1' are provided, inputs for remaining parameters 'object', 'replicate_labels', 'group_labels', and 'use_cells' are not used.")
    }
  } else {
    if (any(is.null(object), is.null(replicate_labels), is.null(group_labels))) {
      stop("When inputs for parameters 'n_replicates' and 'n_group1' are not provided, input must be supplied to parameters 'object', 'replicate_labels', and 'group_labels'. Please supply valid input!")
    }
    if (!is.null(n_replicates)) {
      warning("Input for parameter 'n_replicates' was not used.")
    }
    if (!is.null(n_group1)) {
      warning("Input for parameter 'n_group1' was not used.")
    }
    # Retrieve metadata
    replicates <- .retrieveData(object = object,
                             type = "cell_metadata",
                             name = replicate_labels,
                             use_cells = use_cells)
    groups <- .retrieveData(object = object,
                            type = "cell_metadata",
                            name = group_labels,
                            use_cells = use_cells)
    # Get values
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
