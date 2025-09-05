#' Calculate pseudobulk values
#'
#' This function will generate a set of pseudobulk expression matrices.
#'
#' This function was inspired by R package neurorestore/Libra
#' (Squair et al. 2021).
#'
#' @param object A 'Seurat' or 'SingleCellExperiment' object.
#' @param replicate_labels A character string or vector indicating the name of the
#' column containing the replicate labels.
#' @param split_labels A character string or vector indicating the name of a
#' column by which to split the cells prior to pseudobulking and performing
#' differential expression (e.g., cell types). Results will be returned for
#' each unique value in the column indicated by 'split_labels'. Default =
#' \code{NULL} will run pseudobulk differential expression on all cells
#' together.
#' @param use_cells A vector of cell names subset to. Default = \code{NULL} will
#' use all cells.
#' @param min_cells_per_split A numeric value indicating the minimum number of
#' cells within one split. Pseudobulk expression will not be calculated for
#' splits with fewer cells. Defaults to 1.
#' @param min_replicates_per_split A numeric value indicating the minimum number
#' of distinct replicates represented within one split. Pseudobulk expression
#' will not be calculated for splits with fewer replicates. Defaults to 1.
#' @param min_cells_per_feature A numeric value indicating the minimum number
#' of cells (within a split) with expression of a gene. Pseudobulk expression
#' will not be calculated for genes with fewer cells. Defaults to 0.
#' @param use_assay A character string indicating the assay to use in the
#' provided object. Default = \code{NULL} will choose the current active assay
#' for Seurat objects and the \code{counts} assay for SingleCellExperiment
#' objects.
#' @param use_layer For Seurat objects, a character string or vector indicating
#' the layer — previously known as slot — to use in the provided object.
#' Default = \code{NULL} will use the \code{counts} layer.
#' @param n_cores A numeric value indicating the number of cores to use for
#' parallelization. Default = \code{NULL} will use the number of available cores
#' minus 2.
#' @param verbose A boolean value indicating whether to use verbose output
#' during the execution of this function. Defaults to \code{TRUE}.
#' Can be set to \code{FALSE} for a cleaner output.
#'
#' @return Returns a list with one pseudobulk (gene x replicate) matrix per
#' split.
#'
#' @export
#'
getPseudobulk <- function(object,
                          replicate_labels,
                          split_labels = NULL,
                          use_cells = NULL,
                          min_cells_per_split = 100,
                          min_replicates_per_split = 6,
                          min_cells_per_feature = 10,
                          use_assay = NULL,
                          use_layer = NULL,
                          n_cores = NULL,
                          verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(object, "object", "countCombinations")
  .validInput(replicate_labels, "replicate_labels", object)
  .validInput(split_labels, "split_labels", object)
  .validInput(use_cells, "use_cells", object)
  .validInput(min_cells_per_split, "min_cells_per_split")
  .validInput(min_replicates_per_split, "min_replicates_per_split")
  .validInput(min_cells_per_feature, "min_cells_per_feature")
  .validInput(use_assay, "use_assay", object)
  .validInput(use_slot, "use_slot", list(object, use_assay))
  .validInput(n_cores, "n_cores")
  .validInput(verbose, "verbose")

  # Set defaults
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 2
  }

  # ---------------------------------------------------------------------------
  # Apply filters
  # ---------------------------------------------------------------------------

  # Fetch cell metadata
  if (is.null(use_cells)) {
    use_cells <- .getCellIDs(object)
  }
  replicates <- .retrieveData(object = object,
                              type = "cell_metadata",
                              name = replicate_labels,
                              use_cells = use_cells)
  if (!is.null(split_labels)) {
    splits <- .retrieveData(object = object,
                            type = "cell_metadata",
                            name = split_labels,
                            use_cells = use_cells)
  } else {
    splits <- rep("all", length(replicates))
  }

  # Filter out splits with too few cells & too few replicates
  filter_table <- table(splits, replicates)
  keep_splits <- rownames(filter_table)[rowSums(filter_table) >= min_cells_per_split &
                                          rowSums(filter_table >= 1) >= min_replicates_per_split]
  keep <- splits %in% keep_splits
  use_cells <- use_cells[keep]
  replicates <- replicates[keep]
  splits <- splits[keep]
  n_splits <- dplyr::n_distinct(splits)

  # ---------------------------------------------------------------------------
  # Pseudobulk
  # ---------------------------------------------------------------------------

  # Extract matrix
  count_matrix <- .getMatrix(object = object,
                             use_assay = use_assay,
                             use_layer = use_layer,
                             use_cells = use_cells,
                             verbose = verbose)

  # Progress
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Generating ", n_splits,
                       " pseudobulk ", ifelse(n_splits == 1, "matrix..", "matrices.."))
  if (verbose & n_splits != nrow(filter_table)) {
    message("Skipped ", nrow(filter_table) - n_splits, " split label",
            ifelse((nrow(filter_table) - n_splits) == 1, "", "s"),
            " due to insufficient cells/replicates: ",
            paste0(rownames(filter_table)[!(rownames(filter_table) %in% unique(splits))],
                   collapse = ", "))
  }

  # Create list of gene x replicate pseudobulk matrices, one per split
  pb_list <- pbmcapply::pbmclapply(keep_splits, FUN = function(s) {
    split_s <- splits == s
    model_mat <- stats::model.matrix(~ 0 + rep_, data = data.frame(rep_ = as.character(replicates[split_s])))
    pb_mat <- count_matrix[, split_s] %*% model_mat
    keep_genes <- rowSums(as.matrix(pb_mat > 0)) >= min_cells_per_feature
    pb_mat <- pb_mat[keep_genes, ]
    return(pb_mat)
  }, mc.cores = n_cores)
  names(pb_list) <- keep_splits

  # Return list
  return(pb_list)
}
