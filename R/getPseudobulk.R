#' Calculate pseudobulk values
#'
#' This function will generate pseudobulk expression matri(ces) from a provided
#' Seurat object, SingleCellExperiment object, or feature x cell matrix.
#' Multiple pseudobulk matrices can be generated across different divisions of
#' a dataset, such as cell types.
#'
#' @param object An object of class \code{Seurat}, \code{SingleCellExperiment},
#' or \code{matrix}. Data supplied as class \code{matrix} should be a
#' feature x cell matrix.
#' @param replicate_labels A string indicating the name of the
#' metadata column containing the biological replicate labels or a character
#' vector containing the biological replicate labels in order. The biological
#' replicate labels are used to construct/define the pseudobulks.
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
#' @param min_cells_per_replicate A numeric value indicating the minimum number
#' of cells within one replicate for one split. Pseudobulk steps will not be
#' performed for replicates with fewer cells for that split. Defaults to 10.
#' @param min_replicates_per_split A numeric value indicating the minimum number
#' of distinct replicates represented within one split. Pseudobulk steps
#' will not be performed for splits with fewer replicates. Defaults to 6.
#' @param min_cells_per_feature A numeric value indicating the minimum number
#' of cells (within a split) with expression of a feature. Pseudobulk expression
#' will not be calculated for genes expressed in fewer cells. Defaults to 10.
#' @param min_prop_cells_per_feature A numeric value indicating the minimum
#' proportion of cells (within a split) with expression of a feature. Pseudobulk
#' expression will not be calculated for features expressed in fewer cells.
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
#' @return Returns a list containing the following elements: \describe{
#'   \item{PB_values}{A list of feature x replicate matri(ces) containing
#'   pseudobulk values for each feature, one matrix per split}
#'   \item{metadata}{Dataframe record of quality control metrics for each split}
#'   \item{parameters}{List recording parameter values used}
#'   }
#' @export
#'
getPseudobulk <- function(object,
                          replicate_labels,
                          split_labels = NULL,
                          use_cells = NULL,
                          min_cells_per_split = 100,
                          min_cells_per_replicate = 10,
                          min_replicates_per_split = 6,
                          min_cells_per_feature = 10,
                          min_prop_cells_per_feature = 0.1,
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
  .validInput(min_cells_per_replicate, "min_cells_per_replicate", "generate")
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

  # First, filter based on min cells per replicate per split
  valid_pairs <- filter_table >= min_cells_per_replicate
  valid_combos <- data.frame(which(valid_pairs, arr.ind = TRUE))
  valid_combos$splits <- rownames(filter_table)[valid_combos$splits]
  valid_combos$replicates <- colnames(filter_table)[valid_combos$replicates]
  keep <- mapply(FUN = function(s, r) {
    any(valid_combos$splits == s & valid_combos$replicates == r)},
    splits, replicates)
  use_cells <- use_cells[keep]
  replicates <- replicates[keep]
  splits <- splits[keep]

  # Second, filter based on min (remaining) cells per split
  # And min (remaining) replicates per split
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

  if (n_splits > 0) {
    # Extract matrix
    count_matrix <- .getMatrix(object = object,
                               use_assay = use_assay,
                               use_layer = use_layer,
                               use_cells = use_cells,
                               verbose = verbose)
    # Create list of feature x replicate pseudobulk matrices, one per split
    pb_output <- pbmcapply::pbmclapply(keep_splits, FUN = function(s) {
      split_s <- splits == s
      keep_features_count <- Matrix::rowSums(count_matrix[, split_s, drop = FALSE] > 0) >= min_cells_per_feature
      prop_nonzero <- Matrix::rowMeans((count_matrix[, split_s, drop = FALSE] > 0))
      keep_features_prop <- prop_nonzero >= min_prop_cells_per_feature
      keep_features <- which(keep_features_count & keep_features_prop)
      model_mat <- stats::model.matrix(~ 0 + rep_, data = data.frame(rep_ = as.character(replicates[split_s])))
      pb_mat <- count_matrix[keep_features, split_s, drop = FALSE] %*% model_mat

      # Metadata values
      n_features_excluded <- nrow(count_matrix)-nrow(pb_mat)
      prop_features_excluded <- n_features_excluded/nrow(count_matrix)
      n_reads <- sum(pb_mat)
      if (n_features_excluded > 0) {
        n_reads_pre <- sum(count_matrix[, split_s])
        n_reads_excluded <- n_reads_pre - n_reads
        prop_reads_excluded <- n_reads_excluded/n_reads_pre
      } else {
        n_reads_excluded <- 0
        prop_reads_excluded <- 0
      }
      return(list("pb_mat" = pb_mat,
                  "n_features" = nrow(pb_mat),
                  "n_features_excluded" = n_features_excluded,
                  "prop_features_excluded" = prop_features_excluded,
                  "n_reads" = n_reads,
                  "n_reads_excluded" = n_reads_excluded,
                  "prop_reads_excluded" = prop_reads_excluded))
    }, mc.cores = n_cores)

    pb_list <- do.call(rbind, pb_output)[, "pb_mat"]
    names(pb_list) <- keep_splits

    # Metadata values
    n_features <- unlist(do.call(rbind, pb_output)[, "n_features"])
    n_features_excluded <- unlist(do.call(rbind, pb_output)[, "n_features_excluded"])
    prop_features_excluded <- unlist(do.call(rbind, pb_output)[, "prop_features_excluded"])
    n_reads <- unlist(do.call(rbind, pb_output)[, "n_reads"])
    n_reads_excluded <- unlist(do.call(rbind, pb_output)[, "n_reads_excluded"])
    prop_reads_excluded <- unlist(do.call(rbind, pb_output)[, "prop_reads_excluded"])
    metadata_values <- data.frame(split = keep_splits,
                                  n_features = n_features,
                                  n_features_excluded = n_features_excluded,
                                  prop_features_excluded = prop_features_excluded,
                                  n_reads = n_reads,
                                  n_reads_excluded = n_reads_excluded,
                                  prop_reads_excluded = prop_reads_excluded)

    # Report metadata values
    if (verbose) {
      if (length(pb_list) > 1) {
        # Number of features
        if (min(metadata_values$prop_features_excluded, na.rm = TRUE) == max(metadata_values$prop_features_excluded, na.rm = TRUE)) {
          message("Excluded ", round(min(metadata_values$prop_features_excluded, na.rm = TRUE)*100, 2), "% of features (",
                  metadata_values$n_features_excluded[metadata_values$prop_features_excluded == min(metadata_values$prop_features_excluded, na.rm = TRUE)][1],
                  " features) in each split.")
        } else {
          message("Excluded between ", round(min(metadata_values$prop_features_excluded, na.rm = TRUE)*100, 2), "% (",
                  metadata_values$n_features_excluded[metadata_values$prop_features_excluded == min(metadata_values$prop_features_excluded, na.rm = TRUE)][1],
                  " features) and ", round(max(metadata_values$prop_features_excluded)*100, 2),"% (",
                  metadata_values$n_features_excluded[metadata_values$prop_features_excluded == max(metadata_values$prop_features_excluded, na.rm = TRUE)][1],
                  " features) of features in each split.")
          if (length(metadata_values$split[metadata_values$prop_features_excluded == max(metadata_values$prop_features_excluded, na.rm = TRUE)]) == 1) {
            message("Highest % of features were excluded in split '",
                    metadata_values$split[metadata_values$prop_features_excluded == max(metadata_values$prop_features_excluded, na.rm = TRUE)], "'.")
          } else {
            message("Highest % of features were excluded in splits: ",
                    paste0(metadata_values$split[metadata_values$prop_features_excluded == max(metadata_values$prop_features_excluded, na.rm = TRUE)], collapse = ", "), ".")
          }
        }
        # Number of reads
        if (min(metadata_values$prop_reads_excluded, na.rm = TRUE) == max(metadata_values$prop_reads_excluded, na.rm = TRUE)) {
          message("Excluded ", round(min(metadata_values$prop_reads_excluded, na.rm = TRUE)*100, 2), "% of reads in each split.")
        } else {
          message("Excluded between ", round(min(metadata_values$prop_reads_excluded, na.rm = TRUE)*100, 2), "% (",
                  metadata_values$n_reads_excluded[metadata_values$prop_reads_excluded == min(metadata_values$prop_reads_excluded, na.rm = TRUE)][1],
                  " reads) and ", round(max(metadata_values$prop_reads_excluded)*100, 2),"% (",
                  metadata_values$n_reads_excluded[metadata_values$prop_reads_excluded == max(metadata_values$prop_reads_excluded, na.rm = TRUE)][1],
                  " reads) of reads in each split.")
          if (length(metadata_values$split[metadata_values$prop_reads_excluded == max(metadata_values$prop_reads_excluded, na.rm = TRUE)]) == 1) {
            message("Highest % of reads were excluded in split '",
                    metadata_values$split[metadata_values$prop_reads_excluded == max(metadata_values$prop_reads_excluded, na.rm = TRUE)], "'.")
          } else {
            message("Highest % of reads were excluded in splits: ",
                    paste0(metadata_values$split[metadata_values$prop_reads_excluded == max(metadata_values$prop_reads_excluded, na.rm = TRUE)], collapse = ", "), ".")
          }
        }
      } else {
        # Number of features
        message("Excluded ", round(metadata_values$prop_features_excluded[1]*100, 2),"% of features (",
                metadata_values$n_features_excluded[1],
                " features) in split '", names(pb_list)[1], "'.")
        # Number of reads
        message("Excluded ", round(metadata_values$prop_reads_excluded[1]*100, 2),"% of reads (",
                metadata_values$n_reads_excluded[1],
                " reads) in split '", names(pb_list)[1], "'.")
      }
    }
  } else {
    pb_list <- NULL
    metadata_values <- NULL
  }

  # ---------------------------------------------------------------------------
  # Wrap up
  # ---------------------------------------------------------------------------

  parameter_list <- list("replicate_labels" = replicate_labels,
                         "split_labels" = split_labels,
                         "use_cells" = use_cells,
                         "min_cells_per_split" = min_cells_per_split,
                         "min_replicates_per_split" = min_replicates_per_split,
                         "min_cells_per_feature" = min_cells_per_feature,
                         "min_prop_cells_per_feature" = min_prop_cells_per_feature,
                         "use_assay" = use_assay,
                         "use_layer" = use_layer,
                         "n_cores" = n_cores)

  # Return output
  return(list("PB_values" = pb_list,
              "metadata" = metadata_values,
              "parameters" = parameter_list))
}
