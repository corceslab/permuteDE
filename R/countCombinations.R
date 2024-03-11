#' Calculate number of possible group combinations
#' 
#' For a set of replicates that belong to two groups, this function will 
#' calculate how many different randomized group combinations are possible 
#' (i.e., how many ways to divide the replicates into two groups while retaining 
#' the same number of replicates per group). The input should be either the 
#' replicate and group column names from an object's metadata, or simply the 
#' total number of replicates (n) and the number of samples in one of the two 
#' groups (k).
#' 
#' The formula is: C = n!/(k!(n-k)!) 
#' (see base function \code{choose})
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
#'
#' @return Returns a numeric value indicating the number of possible group 
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
  
  # Calculate
  n_combinations <- choose(n_replicates, n_group1)
  
  return(n_combinations)
}


