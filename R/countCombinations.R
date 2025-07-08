#' Calculate the number of possible group combinations
#'
#' Given two groups of replicates, this function calculates the number of unique
#' ways to divide the replicates into two groups while preserving group sizes.
#' The formula is: C = n!/(k!(n-k)!). The input takes two forms: (1) replicate
#' and group column names; (2) total number of replicates (n) and the number of
#' replicates in one group (k)(doesn't matter which).
#'
#' (see base function \code{choose})
#'
#' @param object A 'Seurat' or 'SingleCellExperiment' object.
#' @param replicate_label A string or character vector indicating the column
#' containing replicate labels. Defaults to \code{NULL}.
#' @param group_label A string or character vector indicating the column
#' containing two comparison group labels. Defaults to \code{NULL}.
#' @param use_cells A vector of cell names to subset. Defaults to \code{NULL},
#' which includes all cells.
#' @param n_replicates A numeric value indicating the total number of replicates.
#' Defaults to \code{NULL}.
#' @param n_group1 A numeric value indicating the number of replicates in one
#' group (doesn't matter which). Defaults to \code{NULL}.
#'
#' @return Returns a numeric value indicating the number of all possible group
#' combinations.
#'
#' @export
#'
countCombinations <- function(object = NULL,
                              replicate_labels = NULL,
                              group_labels = NULL,
                              use_cells = NULL,
                              n_replicates = NULL,
                              n_group1 = NULL) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(object, "object", "countCombinations")
  .validInput(replicate_labels, "replicate_labels", object)
  .validInput(group_labels, "group_labels", object)
  .validInput(use_cells, "use_cells", object)
  .validInput(n_replicates, "n_replicates")
  .validInput(n_group1, "n_group1", n_replicates)

  # ---------------------------------------------------------------------------
  # Calculate number of combinations
  # ---------------------------------------------------------------------------

  if (!is.null(n_replicates) & !is.null(n_group1)) {
    # use 'n_replicates' and 'n_group1'
    if (any(!is.null(object), !is.null(replicate_labels), !is.null(group_labels), !is.null(use_cells))) {
      warning("When inputs for 'n_replicates' and 'n_group1' are provided, inputs for 'object', 'replicate_labels', 'group_labels', and 'use_cells' are not used.")
    }
  } else {
    # use 'object', 'replicate_labels', 'group_labels'
    if (any(is.null(object), is.null(replicate_labels), is.null(group_labels))) {
      stop("If not using 'n_replicates' and 'n_group1', input must be provided to 'object', 'replicate_labels', and 'group_labels'.")
    }
    if (!is.null(n_replicates)) {
      warning("Input for 'n_replicates' was not used.")
    }
    if (!is.null(n_group1)) {
      warning("Input for 'n_group1' was not used.")
    }

    # retrieve metadata
    replicates <- .retrieveData(object = object,
                             type = "cell_metadata",
                             name = replicate_labels,
                             use_cells = use_cells)
    groups <- .retrieveData(object = object,
                            type = "cell_metadata",
                            name = group_labels,
                            use_cells = use_cells)

    # count unique replicates and unique replicates in one group
    n_replicates <- dplyr::n_distinct(replicates)
    n_group1 <- dplyr::n_distinct(replicates[groups == groups[1]])
  }

  # calculate combination
  n_combinations <- choose(n_replicates, n_group1)

  return(n_combinations)
}


