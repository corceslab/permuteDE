#' Calculate pseudobulk values
#'
#' This function will generate pseudobulk expression matri(ces) from a provided
#' Seurat object, SingleCellExperiment object, or feature x cell matrix.
#' Multiple pseudobulk matrices can be generated across different divisions of
#' a dataset, such as cell types.
#'
#' This function was inspired in part by R package neurorestore/Libra and some
#' aspects are adapted therefrom (Squair et al. 2021).
#'
#' @param object An object of class \code{Seurat}, \code{SingleCellExperiment},
#' or \code{matrix}. Data supplied as class \code{matrix} should be a
#' feature x cell matrix.
#' @param replicate_labels A string indicating the name of the
#' metadata column containing the replicate labels or a character vector
#' containing the replicate labels in order.
#' @param split_labels A string indicating the name of a
#' column by which to split the cells prior to pseudobulking and performing
#' differential expression (e.g., cell types). Alternately, a character vector
#' containing the split labels for each cell in order. Results will be returned
#' for each unique value indicated by \code{split_labels}. Default = \code{NULL}
#' will generate a single pseudobulk matrix that includes all cells.
#' @param use_cells A vector of cell names to subset the object to prior to
#' subsequent pseudobulk steps. Default = \code{NULL} will use all cells.
#' @param min_cells_per_split A numeric value indicating the minimum number of
#' cells within one split. Pseudobulk steps will not be performed for splits
#' with fewer cells. Defaults to 100.
#' @param min_replicates_per_split A numeric value indicating the minimum number
#' of distinct replicates represented within one split. Pseudobulk steps
#' will not be performed for splits with fewer replicates. Defaults to 6.
#' @param min_cells_per_feature A numeric value indicating the minimum number
#' of cells (within a split) with expression of a gene. Pseudobulk expression
#' will not be calculated for genes expressed in fewer cells. Defaults to 10.
#' @param min_prop_cells_per_feature A numeric value indicating the minimum
#' proportion of cells (within a split) with expression of a gene. Pseudobulk
#' expression will not be calculated for genes expressed in fewer cells.
#' Defaults to 0.1.
#' @param use_assay A string indicating the assay to use in the
#' provided object. Default = \code{NULL} will choose the current active assay
#' for \code{Seurat} objects and the \code{counts} assay for
#' \code{SingleCellExperiment} objects.
#' @param use_layer For \code{Seurat} objects, a string or vector
#' indicating the layer—previously known as slot—to use in the provided object.
#' Default = \code{NULL} will use the \code{counts} layer.
#' @param n_cores A numeric value indicating the number of cores to use for
#' parallelization. Default = \code{NULL} will use the number of available cores
#' minus 2.
#' @param verbose A Boolean value indicating whether to use verbose output
#' during the execution of this function. Defaults to \code{TRUE}.
#' Can be set to \code{FALSE} for a cleaner output.
#'
#' @return Returns a list with one pseudobulk (feature x replicate) matrix per
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
                          min_prop_cells_per_feature = 0,
                          use_assay = NULL,
                          use_layer = NULL,
                          n_cores = NULL,
                          verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(object, "object", "getPseudobulk")
  .validInput(replicate_labels, "replicate_labels", list(object, "generate"))
  .validInput(split_labels, "split_labels", object)
  .validInput(use_cells, "use_cells", list(object, "generate"))
  .validInput(min_cells_per_split, "min_cells_per_split", "generate")
  .validInput(min_replicates_per_split, "min_replicates_per_split", "generate")
  .validInput(min_cells_per_feature, "min_cells_per_feature", "generate")
  .validInput(min_prop_cells_per_feature, "min_prop_cells_per_feature", "generate")
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
  # Replicate labels
  if (length(replicate_labels) == 1) {
    replicates <- .retrieveData(object = object,
                                type = "cell_metadata",
                                name = replicate_labels,
                                use_cells = use_cells)
  } else {
    replicates <- replicate_labels
    # Check length
    target_length <- length(use_cells)
    if (length(replicates) != target_length) {
      stop("When a vector is provided for 'replicate_labels', it must be the same length and in the same order as the supplied cells.")
    }
  }
  # Split labels
  if (!is.null(split_labels)) {
    if (length(split_labels) == 1) {
      splits <- .retrieveData(object = object,
                              type = "cell_metadata",
                              name = split_labels,
                              use_cells = use_cells)
    } else {
      splits <- split_labels
      # Check length
      if (length(splits) != length(replicates)) {
        stop("When a vector is provided for 'split_labels', it must be the same length and in the same order as the supplied cells.")
      }
    }
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
            paste0(setdiff(rownames(filter_table), unique(splits)),
                   collapse = ", "))
  }

  # Create list of gene x replicate pseudobulk matrices, one per split
  pb_list <- pbmcapply::pbmclapply(keep_splits, FUN = function(s) {
    split_s <- splits == s
    keep_genes_count <- Matrix::rowSums(count_matrix[, split_s, drop = FALSE] > 0) >= min_cells_per_feature
    prop_nonzero <- Matrix::rowMeans((count_matrix[, split_s, drop = FALSE] > 0))
    keep_genes_prop <- prop_nonzero >= min_prop_cells_per_feature
    keep_genes <- which(keep_genes_count & keep_genes_prop)
    model_mat <- stats::model.matrix(~ 0 + rep_, data = data.frame(rep_ = as.character(replicates[split_s])))
    pb_mat <- count_matrix[keep_genes, split_s, drop = FALSE] %*% model_mat

    # Warn if excluded genes are >10% of all genes
    prop_genes_excluded <- 1 - (length(keep_genes)/nrow(count_matrix) )
    if (prop_genes_excluded > 0.1) {
      message("Warning: Excluded ", round(prop_genes_excluded*100, 2),"% of genes in split ", s)
    }

    return(pb_mat)
  }, mc.cores = n_cores)
  names(pb_list) <- keep_splits

  # Return list
  return(pb_list)
}
