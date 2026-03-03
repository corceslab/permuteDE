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
#' @param metadata An optional dataframe containing relevant metadata columns
#' corresponding to the data provided to parameter \code{object}. Default =
#' \code{NULL} looks for metadata in \code{object} or other provided inputs.
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
#' will not be calculated for features expressed in fewer cells. Defaults to 10.
#' @param min_prop_cells_per_feature A numeric value indicating the minimum
#' proportion of cells (within a split) with expression of a feature. Pseudobulk
#' expression will not be calculated for features expressed in fewer cells.
#' Defaults to 0.1.
#' @param filter A Boolean value indicating whether to remove features from
#' the matrices in the output (\code{TRUE}) or simply list features that do not
#' meet the criteria as part of the metadata (\code{FALSE}). Defaults to
#' \code{TRUE}.
#' @param pseudobulk A string indicating whether to actually pseudobulk
#' ("generate"), or simply apply the filtering thresholds ("none"). Defaults to
#' "generate".
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
#'   \item{exclude_features}{A list of feature x replicate matri(ces) containing
#'   pseudobulk values for each feature, one matrix per split}
#'   \item{metadata}{Dataframe record of quality control metrics for each split}
#'   \item{parameters}{List recording parameter values used}
#'   }
#' @export
#'
getPseudobulk <- function(object,
                          metadata = NULL,
                          replicate_labels,
                          split_labels = NULL,
                          use_cells = NULL,
                          min_cells_per_split = 100,
                          min_cells_per_replicate = 10,
                          min_replicates_per_split = 6,
                          min_cells_per_feature = 10,
                          min_prop_cells_per_feature = 0.1,
                          filter = TRUE,
                          pseudobulk = "generate",
                          use_assay = NULL,
                          use_layer = NULL,
                          n_cores = NULL,
                          verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(object, "object", "getPseudobulk")
  .validInput(metadata, "metadata", object)
  .validInput(replicate_labels, "replicate_labels", list(object, metadata, pseudobulk))
  .validInput(split_labels, "split_labels", list(object, metadata))
  .validInput(use_cells, "use_cells", list(object, pseudobulk))
  .validInput(min_cells_per_split, "min_cells_per_split", list(pseudobulk))
  .validInput(min_cells_per_replicate, "min_cells_per_replicate", list(pseudobulk, "getPseudobulk"))
  .validInput(min_replicates_per_split, "min_replicates_per_split", list(pseudobulk, "getPseudobulk"))
  .validInput(min_cells_per_feature, "min_cells_per_feature", list(pseudobulk))
  .validInput(min_prop_cells_per_feature, "min_prop_cells_per_feature", pseudobulk)
  .validInput(filter, "filter")
  .validInput(pseudobulk, "pseudobulk", list("getPseudobulk", object))
  .validInput(use_assay, "use_assay", object)
  .validInput(use_slot, "use_slot", list(object, use_assay))
  .validInput(n_cores, "n_cores")
  .validInput(verbose, "verbose")

  # Set defaults
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 2
  }

  # Object type
  if (methods::is(object, "Seurat")) {
    object_type <- "Seurat"
  } else if (methods::is(object, "SingleCellExperiment")) {
    object_type <- "SingleCellExperiment"
  } else {
    object_type <- "matrix"
  }

  # ---------------------------------------------------------------------------
  # Apply filters
  # ---------------------------------------------------------------------------

  # Fetch cell metadata
  if (is.null(use_cells)) {
    use_cells <- .getCellIDs(object)
  }

  # Replicate labels
  if (pseudobulk == "generate") {
    if (length(replicate_labels) == 1) {
      replicates <- .retrieveData(object = object,
                                  metadata = metadata,
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
  } else {
    replicates <- use_cells
  }
  # Check for NA values
  if (any(is.na(replicates))) {
    stop("Values provided for 'replicate_labels' cannot be NA.")
  }
  if (!methods::is(replicates, "character")) {
    replicates <- as.character(replicates)
  }

  # Split labels
  if (!is.null(split_labels)) {
    if (length(split_labels) == 1) {
      splits <- .retrieveData(object = object,
                              metadata = metadata,
                              type = "cell_metadata",
                              name = split_labels,
                              use_cells = use_cells)
    } else {
      splits <- split_labels
      # Check length
      if (length(splits) != length(use_cells)) {
        stop("When a vector is provided for 'split_labels', it must be the same length and in the same order as the supplied cells.")
      }
    }
  } else {
    splits <- rep("all", length(use_cells))
  }
  # Check for NA values
  if (any(is.na(splits))) {
    stop("Values provided for 'split_labels' cannot be NA.")
  }
  splits <- as.character(splits)

  if (pseudobulk == "generate") {
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
    keep_splits <- rownames(filter_table)[Matrix::rowSums(filter_table) >= min_cells_per_split &
                                            Matrix::rowSums(filter_table >= 1) >= min_replicates_per_split]
    keep <- splits %in% keep_splits
    use_cells <- use_cells[keep]
    replicates <- replicates[keep]
    splits <- splits[keep]
    n_splits <- dplyr::n_distinct(splits)
  } else {
    # Filter based on min cells per split
    filter_table <- table(splits)
    keep_splits <- names(filter_table)[filter_table >= min_cells_per_split]
    keep <- splits %in% keep_splits
    use_cells <- use_cells[keep]
    splits <- splits[keep]
    n_splits <- dplyr::n_distinct(splits)
  }

  # ---------------------------------------------------------------------------
  # Pseudobulk
  # ---------------------------------------------------------------------------

  # Progress
  if (pseudobulk == "generate") {
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Generating ", n_splits,
                         " pseudobulk ", ifelse(n_splits == 1, "matrix..", "matrices.."))
    if (verbose & n_splits != nrow(filter_table)) {
      message("Skipped ", nrow(filter_table) - n_splits, " split label",
              ifelse((nrow(filter_table) - n_splits) == 1, "", "s"),
              " due to insufficient cells/replicates: ",
              paste0(setdiff(rownames(filter_table), unique(splits)),
                     collapse = ", "))
    }
  } else {
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Filtering ", n_splits,
                         " count ", ifelse(n_splits == 1, "matrix..", "matrices.."))
    if (verbose & n_splits != nrow(filter_table)) {
      message("Skipped ", nrow(filter_table) - n_splits, " split label",
              ifelse((nrow(filter_table) - n_splits) == 1, "", "s"),
              " due to insufficient cells: ",
              paste0(setdiff(rownames(filter_table), unique(splits)),
                     collapse = ", "))
    }
  }

  if (n_splits > 0) {
    # Extract matrix
    count_matrix <- .getMatrix(object = object,
                               use_assay = use_assay,
                               use_layer = use_layer,
                               use_cells = use_cells,
                               verbose = verbose)
    # Create list of feature x replicate pseudobulk matrices, one per split
    output_list <- pbmcapply::pbmclapply(keep_splits, FUN = function(s) {
      split_s <- splits == s
      # Remove features with 0 counts
      count_matrix_s <- count_matrix[Matrix::rowSums(count_matrix[, split_s]) > 0, split_s]

      # Identify features that don't pass input thresholds
      keep_features_count <- Matrix::rowSums(count_matrix_s > 0) >= min_cells_per_feature
      prop_nonzero <- Matrix::rowMeans((count_matrix_s > 0))
      keep_features_prop <- prop_nonzero >= min_prop_cells_per_feature
      exclude_features <- rownames(count_matrix_s)[-which(keep_features_count & keep_features_prop)]

      # Pseudobulk
      if (pseudobulk == "generate") {
        model_mat <- stats::model.matrix(~ 0 + rep_, data = data.frame(rep_ = as.character(replicates[split_s])))
        output_mat <- count_matrix_s %*% model_mat
      } else {
        output_mat <- count_matrix_s
      }

      # Metadata values
      n_all_features <- nrow(count_matrix)
      n_nonzero_features <- nrow(output_mat)
      n_features_exclude <- length(exclude_features)
      n_features_for_DE <- n_nonzero_features-n_features_exclude
      prop_features_exclude <- (n_all_features-n_features_for_DE)/n_all_features
      n_all_reads <- sum(Matrix::rowSums(count_matrix_s))
      n_reads_exclude <- sum(Matrix::rowSums(output_mat[exclude_features,]))
      n_reads_for_DE <- n_all_reads - n_reads_exclude
      prop_reads_exclude <- n_reads_exclude/n_all_reads

      # Exclude features
      if (filter == TRUE) {
        output_mat <- output_mat[!(rownames(output_mat) %in% exclude_features),]
      }

      return(list("output_mat" = output_mat,
                  "exclude_features" = exclude_features,
                  "n_all_features" = n_all_features,
                  "n_nonzero_features" = n_nonzero_features,
                  "n_features_exclude" = n_features_exclude,
                  "n_features_for_DE" = n_features_for_DE,
                  "prop_features_exclude" = prop_features_exclude,
                  "n_all_reads" = n_all_reads,
                  "n_reads_exclude" = n_reads_exclude,
                  "n_reads_for_DE" = n_reads_for_DE,
                  "prop_reads_exclude" = prop_reads_exclude))
    }, mc.cores = n_cores)

    # Combine output
    mat_list <- do.call(rbind, output_list)[, "output_mat"]
    names(mat_list) <- keep_splits
    exclude_features_list <- do.call(rbind, output_list)[, "exclude_features"]
    names(exclude_features_list) <- keep_splits

    # Metadata values
    n_all_features <- unlist(do.call(rbind, output_list)[, "n_all_features"])
    n_nonzero_features <- unlist(do.call(rbind, output_list)[, "n_nonzero_features"])
    n_features_exclude <- unlist(do.call(rbind, output_list)[, "n_features_exclude"])
    n_features_for_DE <- unlist(do.call(rbind, output_list)[, "n_features_for_DE"])
    prop_features_exclude <- unlist(do.call(rbind, output_list)[, "prop_features_exclude"])
    n_all_reads <- unlist(do.call(rbind, output_list)[, "n_all_reads"])
    n_reads_exclude <- unlist(do.call(rbind, output_list)[, "n_reads_exclude"])
    n_reads_for_DE <- unlist(do.call(rbind, output_list)[, "n_reads_for_DE"])
    prop_reads_exclude <- unlist(do.call(rbind, output_list)[, "prop_reads_exclude"])
    metadata_values <- data.frame(split = keep_splits,
                                  n_all_features = n_all_features,
                                  n_nonzero_features = n_nonzero_features,
                                  n_features_exclude = n_features_exclude,
                                  n_features_for_DE = n_features_for_DE,
                                  prop_features_exclude = prop_features_exclude,
                                  n_all_reads = n_all_reads,
                                  n_reads_exclude = n_reads_exclude,
                                  n_reads_for_DE = n_reads_for_DE,
                                  prop_reads_exclude = prop_reads_exclude)
    metadata_list <- list("metrics" = metadata_values,
                          "exclude_features" = exclude_features_list)

    # Report metadata values
    if (verbose) {
      if (filter == TRUE) {
        exclusion_message <- "Excluded"
      } else {
        exclusion_message <- "Will exclude"
      }
      if (length(mat_list) > 1) {
        # Number of features
        if (min(metadata_values$prop_features_exclude, na.rm = TRUE) == max(metadata_values$prop_features_exclude, na.rm = TRUE)) {
          message(exclusion_message, " ", round(min(metadata_values$prop_features_exclude, na.rm = TRUE)*100, 2), "% of features (",
                  metadata_values$n_features_exclude[metadata_values$prop_features_exclude == min(metadata_values$prop_features_exclude, na.rm = TRUE)][1],
                  " features) in each split.")
        } else {
          message(exclusion_message, " between ", round(min(metadata_values$prop_features_exclude, na.rm = TRUE)*100, 2), "% (",
                  metadata_values$n_features_exclude[metadata_values$prop_features_exclude == min(metadata_values$prop_features_exclude, na.rm = TRUE)][1],
                  " features) and ", round(max(metadata_values$prop_features_exclude)*100, 2),"% (",
                  metadata_values$n_features_exclude[metadata_values$prop_features_exclude == max(metadata_values$prop_features_exclude, na.rm = TRUE)][1],
                  " features) of features in each split.")
          if (length(metadata_values$split[metadata_values$prop_features_exclude == max(metadata_values$prop_features_exclude, na.rm = TRUE)]) == 1) {
            message("Highest % of features excluded in split '",
                    metadata_values$split[metadata_values$prop_features_exclude == max(metadata_values$prop_features_exclude, na.rm = TRUE)], "'.")
          } else {
            message("Highest % of features excluded in splits: ",
                    paste0(metadata_values$split[metadata_values$prop_features_exclude == max(metadata_values$prop_features_exclude, na.rm = TRUE)], collapse = ", "), ".")
          }
        }
        # Number of reads
        if (min(metadata_values$prop_reads_exclude, na.rm = TRUE) == max(metadata_values$prop_reads_exclude, na.rm = TRUE)) {
          message(exclusion_message, " ", round(min(metadata_values$prop_reads_exclude, na.rm = TRUE)*100, 2), "% of reads in each split.")
        } else {
          message(exclusion_message, " between ", round(min(metadata_values$prop_reads_exclude, na.rm = TRUE)*100, 2), "% (",
                  metadata_values$n_reads_exclude[metadata_values$prop_reads_exclude == min(metadata_values$prop_reads_exclude, na.rm = TRUE)][1],
                  " reads) and ", round(max(metadata_values$prop_reads_exclude)*100, 2),"% (",
                  metadata_values$n_reads_exclude[metadata_values$prop_reads_exclude == max(metadata_values$prop_reads_exclude, na.rm = TRUE)][1],
                  " reads) of reads in each split.")
          if (length(metadata_values$split[metadata_values$prop_reads_exclude == max(metadata_values$prop_reads_exclude, na.rm = TRUE)]) == 1) {
            message("Highest % of reads excluded in split '",
                    metadata_values$split[metadata_values$prop_reads_exclude == max(metadata_values$prop_reads_exclude, na.rm = TRUE)], "'.")
          } else {
            message("Highest % of reads excluded in splits: ",
                    paste0(metadata_values$split[metadata_values$prop_reads_exclude == max(metadata_values$prop_reads_exclude, na.rm = TRUE)], collapse = ", "), ".")
          }
        }
      } else {
        # Number of features
        message(exclusion_message, " ", round(metadata_values$prop_features_exclude[1]*100, 2),"% of features (",
                metadata_values$n_features_exclude[1],
                " features) in split '", names(mat_list)[1], "'.")
        # Number of reads
        message(exclusion_message, " ", round(metadata_values$prop_reads_exclude[1]*100, 2),"% of reads (",
                metadata_values$n_reads_exclude[1],
                " reads) in split '", names(mat_list)[1], "'.")
      }
    }
  } else {
    mat_list <- NULL
    metadata_list <- NULL
  }

  # ---------------------------------------------------------------------------
  # Wrap up
  # ---------------------------------------------------------------------------

  parameter_list <- list("object_type" = object_type,
                         "metadata_provided" = !is.null(metadata),
                         "replicate_labels" = replicate_labels,
                         "split_labels" = split_labels,
                         "use_cells" = use_cells,
                         "min_cells_per_split" = min_cells_per_split,
                         "min_cells_per_replicate" = min_cells_per_replicate,
                         "min_replicates_per_split" = min_replicates_per_split,
                         "min_cells_per_feature" = min_cells_per_feature,
                         "min_prop_cells_per_feature" = min_prop_cells_per_feature,
                         "use_assay" = use_assay,
                         "use_layer" = use_layer,
                         "n_cores" = n_cores)

  # Return output
  if (pseudobulk == "generate") {
    return(list("PB_values" = mat_list,
                "metadata" = metadata_list,
                "parameters" = parameter_list))
  } else {
    return(list("cell_values" = mat_list,
                "metadata" = metadata_list,
                "parameters" = parameter_list))
  }
}
