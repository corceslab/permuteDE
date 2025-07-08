#' Generate combinations
#'
#' Given two groups of replicates, this function will generate
#' a set of randomized group combinations. The input takes two forms: (1) replicate
#' and group column names; (2) total number of replicates (n) and the number of
#' replicates in one group (k)(either group is ok).
#'
#' @param object A 'Seurat' or 'SingleCellExperiment' object.
#' @param replicate_labels A string or character vector indicating the column
#' containing replicate labels. Defaults to \code{NULL}.
#' @param group_labels A string or character vector indicating the column
#' containing two comparison group labels. Defaults to \code{NULL}.
#' @param use_cells A vector of cell names to subset. Defaults to \code{NULL},
#' which includes all cells.
#' @param n_replicates A numeric value indicating the total number of replicates.
#' Defaults to \code{NULL}.
#' @param n_group1 A numeric value indicating the number of replicates in one
#' group (doesn't matter which). Defaults to \code{NULL}.
#' @param n_combinations A numeric value indicating the number of combinations
#' to generate. Defaults to 1000.
#' @param message A string indicating additional progress messaging (internal
#' use). Defaults to "".
#' @param random_seed A numeric value indicating the random seed to use.
#' Defaults to 1.
#' @param verbose A boolean value indicating whether to use verbose output.
#' Defaults to \code{TRUE}. Use \code{FALSE} for a cleaner output.
#'
#' @return Returns a matrix of index combinations for assigning replicates to the first group.
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
  .validInput(replicate_labels, "replicate_labels", object)
  .validInput(group_labels, "group_labels", object)
  .validInput(use_cells, "use_cells", object)
  .validInput(n_replicates, "n_replicates")
  .validInput(n_group1, "n_group1", n_replicates)
  .validInput(n_combinations, "n_combinations")
  .validInput(message, "message")
  .validInput(random_seed, "random_seed")
  .validInput(verbose, "verbose")

  # set seed
  set.seed(random_seed)

  # count possible combinations
  n_possible_combinations <- countCombinations(object, replicate_labels,
                                               group_labels, use_cells,
                                               n_replicates, n_group1)

  print(n_possible_combinations)

  # warn if requested combinations exceed the possible threshold
  if (n_possible_combinations < n_combinations) {
    warning("Input ", n_combinations,
            " for parameter 'n_combinations' exceeds the number of possible combinations (",
            n_possible_combinations, "). Only ", n_possible_combinations,
            " combinations will be generated.")
    n_combinations <- n_possible_combinations
  }

  # ---------------------------------------------------------------------------
  # Generate combinations
  # ---------------------------------------------------------------------------

  # print progress message if verbose is enabled
  if (verbose) {
    message(format(Sys.time(), "%Y-%m-%d %X"), " : Generating ",
                       n_combinations, " possible combinations", message, "...")
  }

  if (n_possible_combinations < 1000000) {
    # when 'n_possible_combinations' is not big, generate all combinations, then shuffle
    all_combinations <- utils::combn(n_replicates, n_group1)
    # shuffle 'all_combinations' columns
    output_combinations <- all_combinations[, sample(1:ncol(all_combinations), n_combinations)]
  } else {
    # when 'n_possible_combinations' is extremely big, use a list

    # TODO: this part lacks data to test. expect issues on github.
    # initialize an environments to avoid duplicate combinations
    seen <- new.env(hash = TRUE, size = n_combinations)
    output_combinations <- list()

    # generate unique random combinations until the desired number is reached
    while (length(output_combinations) < n_combinations) {
      comb <- sort(sample(n_replicates, n_group1))

      # create a unique key representing this combination
      key <- paste(comb, collapse = "-")
      if (!exists(key, envir = seen, inherits = FALSE)) {
        # only store the combination if it hasn't been seen before
        seen[[key]] <- TRUE
        output_combinations[[length(output_combinations) + 1]] <- comb
      }
    }

    output_combinations <- do.call(cbind, output_combinations)

  }

  return(output_combinations)
}
