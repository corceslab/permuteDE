#' Calculate the number of possible group combinations
#'
#' Given two groups of replicates, this function calculates the number of unique
#' ways to divide the replicates into two groups while preserving group sizes.
#' The formula is: C = n!/(k!(n-k)!). (see base function \code{choose})
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
#' @param replicate_labels A string indicating the name of the
#' metadata column containing the biological replicate labels or a character
#' vector containing the biological replicate labels in order.
#' @param group_labels A string indicating the name of the
#' column containing the two comparison group labels or a character vector
#' containing the comparison labels in order.
#' @param use_cells A vector of cell names to subset prior to calculating
#' possible group comibinations. Default = \code{NULL} will use all cells.
#' @param n_replicates A numeric value indicating the total number of
#' replicates. Defaults to \code{NULL}.
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
  .validInput(replicate_labels, "replicate_labels", list(object, "not applicable"))
  .validInput(group_labels, "group_labels", object)
  .validInput(use_cells, "use_cells", list(object, "not applicable"))
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
    if (any(is.null(replicate_labels), is.null(group_labels))) {
      stop("If not using 'n_replicates' and 'n_group1', input must be provided to 'replicate_labels' and 'group_labels'.")
    }
    if (!is.null(n_replicates)) {
      warning("Input for 'n_replicates' was not used.")
    }
    if (!is.null(n_group1)) {
      warning("Input for 'n_group1' was not used.")
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
    # Check that there are exactly two groups present
    if (dplyr::n_distinct(groups) != 2) {
      stop("Input value for parameter 'group_labels' must represent a cell metadata column (or a vector of group labels) that contains exactly 2 groups for the selected data, please supply valid input!")
    }
    # Count unique replicates and unique replicates in one group
    n_replicates <- dplyr::n_distinct(replicates)
    n_group1 <- dplyr::n_distinct(replicates[groups == groups[1]])
  }

  # Calculate number of combinations
  n_combinations <- choose(n_replicates, n_group1)

  return(n_combinations)
}


