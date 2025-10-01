#' Run differential expression analysis to compare two groups
#'
#' This function identifies differentially expressed features between two groups
#' using indicated differential expression analysis methods.
#'
#' By default, pseudobulk matri(ces) are generated or supplied by the user, then
#' used to run pseudobulk differential expression. The following existing tools
#' are supported: \code{edgeR}, \code{DESeq2}, \code{limma}, and the Wilcoxon
#' rank-sum test. Alternately, users may skip pseudobulking and run cell-level
#' differential expression (not recommended in most cases).
#'
#' @param object An object of class \code{Seurat}, \code{SingleCellExperiment},
#' or \code{matrix}. Data supplied as class \code{matrix} may be either a
#' feature x cell matrix or a pre-computed pseudobulk feature x replicate
#' matrix. Note that raw counts are expected, and the normalization method
#' applied during differential expression analysis differs across the methods
#' and tests.
#' @param replicate_labels A string indicating the name of the
#' metadata column containing the biological replicate labels or a character
#' vector containing the biological replicate labels in order. For pseudobulk DE
#' analysis, the biological replicate labels are used to construct/define the
#' pseudobulks. For cell-level DE analysis, the biological replicate labels are
#' not used in this function, but will be passed on to function permuteDE so
#' that the permutation of group labels keeps replicates intact.
#' @param group_labels A string indicating the name of the
#' column containing the two comparison group labels or a character vector
#' containing the comparison labels in order.
#' @param split_labels A string indicating the name of a
#' column by which to split the cells prior to pseudobulking and performing
#' differential expression (e.g., cell types). Alternately, a character vector
#' containing the split labels for each cell in order. Results will be returned
#' for each unique value indicated by \code{split_labels}. Default = \code{NULL}
#' will run pseudobulk differential expression on all cells together.
#' @param reference_group A string specifying the reference group. Defaults to
#' \code{NULL}, in which case the first value alphabetically is used as the
#' reference.
#' @param use_cells A vector of cell names to subset the object to prior to
#' subsequent pseudobulk and differential expression steps. Default =
#' \code{NULL} will use all cells.
#' @param pseudobulk A string indicating pseudobulk handling. Permitted values
#' are: "generate" (pseudobulk matrices will be generated), "supplied"
#' (pseudobulk matrix was supplied by the user to parameter \code{object}), or
#' "none" (pseudobulking will not be used, cell-level differential expression
#' analysis will be run). Defaults to "generate".
#' @param de_method Which tool to use for differential expression analysis.
#' Permitted values are "edgeR", "DESeq2", "limma", and "presto". Defaults to
#' "edgeR".
#' @param de_test Which test to use for differential expression analysis.
#' Available values are dependent on the \code{de_method}: "edgeR" ("LRT",
#' "QLF", "exact"), "DESeq2" ("LRT", "Wald"), "limma" ("trend", "voom"),
#' "presto" ("wilcox_cpm", "wilcox_log_cpm"). Defaults to "LRT".
#' @param de_params A list of lists containing additional parameters to be
#' passed to specific DE functions. The name of each element must be the
#' specific DE function to which those parameters are passed. Defaults to an
#' empty list.
#' @param p_adjust_method A string indicating which multiple comparison
#' adjustment to use. For permitted values, see \code{stats::p.adjust.methods}.
#' Defaults to "fdr" (Benjamini & Hochberg, 1995).
#' @param min_cells_per_split A numeric value indicating the minimum number of
#' cells within one split. Pseudobulk and differential expression steps will not
#' be performed for splits with fewer cells. Defaults to 100.
#' @param min_replicates_per_split A numeric value indicating the minimum number
#' of distinct replicates represented within one split. Pseudobulk expression
#' and differential expression will not be performed for splits with fewer
#' replicates. Defaults to 6.
#' @param min_replicates_per_group A numeric value indicating the minimum number
#' of distinct replicates represented within each of the two comparison groups.
#' Pseudobulk and differential expression steps will not be performed for
#' splits with fewer replicates. Defaults to 3.
#' @param min_cells_per_feature A numeric value indicating the minimum number
#' of cells (within a split) with expression of a gene. Pseudobulk and
#' differential expression will not be calculated for genes expressed in
#' fewer cells. Defaults to 10.
#' @param min_prop_cells_per_feature A numeric value indicating the minimum
#' proportion of cells (within a split) with expression of a gene. Pseudobulk
#' and differential expression will not be calculated for genes
#' expressed in fewer cells. Defaults to 0.1.
#' @param force_balance A boolean indicating whether to force the two comparison
#' groups to have the same sample size. Defaults to \code{FALSE}. If
#' \code{TRUE}, the larger group will be randomly downsampled to the size of the
#' smaller group.
#' @param use_assay A string indicating the assay to use in the
#' provided object. Default = \code{NULL} will choose the current active assay
#' for \code{Seurat} objects and the \code{counts} assay for
#' \code{SingleCellExperiment} objects.
#' @param use_layer For \code{Seurat} objects, a string or vector
#' indicating the layer—previously known as slot—to use in the provided object.
#' Default = \code{NULL} will use the \code{counts} layer.
#' @param random_seed A numerical value indicating the random seed to be used.
#' Defaults to 1. Only relevant in this function when parameter
#' \code{force_balance = TRUE}.
#' @param n_cores A numeric value indicating the number of cores to use for
#' parallelization. Default = \code{NULL} will use the number of available cores
#' minus 2.
#' @param verbose A Boolean value indicating whether to use verbose output
#' during the execution of this function. Defaults to \code{TRUE}.
#' Can be set to \code{FALSE} for a cleaner output.
#'
#' @return Returns a list containing the following elements: \describe{
#'   \item{DE_results}{Dataframe containing DE results for each feature, by
#'   split}
#'   \item{PB_values}{If using pseudobulk data, a list of feature x replicate
#'   matri(ces) containing pseudobulk values for each feature, one matrix per
#'   split}
#'   \item{cell_values}{Alternately, if using cell-level data, a list of feature x
#'   cell matri(ces) containing counts for each feature, one matrix per
#'   split}
#'   \item{group_key}{Dataframe record of group labels corresponding to each
#'   replicate}
#'   \item{metadata}{List recording characteristics of the data and runtime}
#'   \item{parameters}{List recording parameter values used}
#'   }
#'
#' @export
#'
runDE <- function(object,
                  replicate_labels = NULL,
                  group_labels,
                  split_labels = NULL,
                  reference_group = NULL,
                  use_cells = NULL,
                  pseudobulk = "generate",
                  de_method = "edgeR",
                  de_test = "LRT",
                  de_params = list(),
                  p_adjust_method = "fdr",
                  min_cells_per_split = 100,
                  min_replicates_per_split = 6,
                  min_replicates_per_group = 3,
                  min_cells_per_feature = 10,
                  min_prop_cells_per_feature = 0.1,
                  force_balance = FALSE,
                  use_assay = NULL,
                  use_layer = NULL,
                  random_seed = 1,
                  n_cores = NULL,
                  verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  time1 <- Sys.time()

  .validInput(object, "object", "runDE")
  .validInput(pseudobulk, "pseudobulk", object)
  .validInput(replicate_labels, "replicate_labels", list(object, pseudobulk))
  .validInput(group_labels, "group_labels", object)
  .validInput(split_labels, "split_labels", object)
  .validInput(use_cells, "use_cells", list(object, pseudobulk))
  .validInput(reference_group, "reference_group", list(object, group_labels, use_cells))
  .validInput(de_method, "de_method")
  .validInput(de_test, "de_test", de_method)
  .validInput(de_params, "de_params", list(de_method, de_test))
  .validInput(p_adjust_method, "p_adjust_method")
  .validInput(min_cells_per_split, "min_cells_per_split", pseudobulk)
  .validInput(min_replicates_per_split, "min_replicates_per_split", pseudobulk)
  .validInput(min_replicates_per_split, "min_replicates_per_group")
  .validInput(min_cells_per_feature, "min_cells_per_feature", pseudobulk)
  .validInput(min_prop_cells_per_feature, "min_prop_cells_per_feature", pseudobulk)
  .validInput(force_balance, "force_balance", pseudobulk)
  .validInput(use_assay, "use_assay", object)
  .validInput(use_slot, "use_slot", list(object, use_assay))
  .validInput(random_seed, "random_seed")
  .validInput(n_cores, "n_cores")
  .validInput(verbose, "verbose")

  # ---------------------------------------------------------------------------
  # Set up
  # ---------------------------------------------------------------------------

  # Set defaults
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 2
  }

  # Random seed reproducibility
  if (force_balance == TRUE) {
    RNGkind("L'Ecuyer-CMRG")
    set.seed(random_seed)
  }

  # Retrieve metadata
  # Set cells to use
  if (is.null(use_cells) & pseudobulk != "supplied") {
    use_cells <- .getCellIDs(object)
  }

  # Replicate labels
  stored_replicates <- NULL
  if (!is.null(replicate_labels)) {
    if (length(replicate_labels) == 1) {
      replicates <- .retrieveData(object = object,
                                  type = "cell_metadata",
                                  name = replicate_labels,
                                  use_cells = use_cells)
    } else {
      replicates <- replicate_labels
      # Check length
      if (!is.null(use_cells) & pseudobulk != "supplied") {
        target_length <- length(use_cells)
      } else {
        target_length <- ncol(object)
      }
      if (length(replicates) != target_length) {
        if (pseudobulk == "supplied") {
          stop("When a vector is provided for 'replicate_labels', it must be the same length and in the same order as the supplied pseudobulk matrix columns.")
        } else {
          stop("When a vector is provided for 'replicate_labels', it must be the same length and in the same order as the supplied cells.")
        }
      }
    }
    if (pseudobulk == "none") {
      # Store biological replicates for use later if output is passed to permuteDE
      stored_replicates <- replicates
      replicates <- use_cells
      names(stored_replicates) <- replicates
    }
  } else if (pseudobulk == "none") {
    replicates <- use_cells
  }

  # Group labels
  if (length(group_labels) == 1) {
    groups <- .retrieveData(object = object,
                            type = "cell_metadata",
                            name = group_labels,
                            use_cells = use_cells)
  } else {
    groups <- group_labels
    # Check length
    if (length(groups) != length(replicates)) {
      if (pseudobulk == "supplied") {
        stop("When a vector is provided for 'group_labels', it must be the same length and in the same order as the supplied pseudobulk matrix columns.")
      } else {
        stop("When a vector is provided for 'group_labels', it must be the same length and in the same order as the supplied cells.")
      }
    }
  }

  # There must be exactly two comparison groups
  if (dplyr::n_distinct(groups) != 2) {
    stop("Input value for parameter 'group_labels' must represent a cell metadata column (or a vector of group labels) that contains exactly 2 groups for the selected data, please supply valid input!")
  }

  # Reference group
  if (is.null(reference_group)) {
    reference_group <- sort(unique(groups))[1]
  }

  # Group labels for each replicate
  # Returns a dataframe with two columns: replicate, group
  group_key <- data.frame(replicate = paste0("rep_", replicates),
                          group = groups) |>
    dplyr::group_by(replicate, group) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::select(-n) |>
    data.frame()
  rownames(group_key) <- group_key$replicate

  # ---------------------------------------------------------------------------
  # Obtain matrices for DE
  # ---------------------------------------------------------------------------

  time2 <- Sys.time()

  if (pseudobulk == "generate") {
    # Generate pseudobulk matri(ces)
    # Returns a list containing one pseudobulk matrix (gene x replicate) per split
    pseudobulk_output <- getPseudobulk(object = object,
                                       replicate_labels = replicate_labels,
                                       split_labels = split_labels,
                                       use_cells = use_cells,
                                       min_cells_per_split = min_cells_per_split,
                                       min_replicates_per_split = min_replicates_per_split,
                                       min_cells_per_feature = min_cells_per_feature,
                                       min_prop_cells_per_feature = min_prop_cells_per_feature,
                                       n_cores = n_cores,
                                       verbose = verbose)
    matrix_list <- pseudobulk_output[["PB_values"]]
    pseudobulk_metadata <- pseudobulk_output[["metadata"]]
  } else {
    # If necessary, separate the supplied pseudobulk matrix by split
    # Returns a list containing one pseudobulk matrix (gene x replicate) per split
    if (is.null(split_labels)) {
      split_labels <- rep("all", length(replicates))
    } else if (length(split_labels) != length(replicates)) { # Check length
      if (pseudobulk == "supplied") {
        stop("When a vector is provided for 'split_labels', it must be the same length and in the same order as the supplied pseudobulk matrix columns.")
      } else {
        stop("When a vector is provided for 'split_labels', it must be the same length and in the same order as the supplied cells.")
      }
    }
    split_indices <- split(seq_along(split_labels), split_labels)
    # Filter
    if (pseudobulk == "supplied") {
      keep_indices <- split_indices[lengths(split_indices) >= min_replicates_per_split]
      n_splits <- length(keep_indices)
      # Progress
      if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Extracting ", n_splits,
                           " pseudobulk ", ifelse(n_splits == 1, "matrix..", "matrices.."))
      # Create matrix list
      matrix_list <- lapply(keep_indices, function(i) object[, i, drop = FALSE])
      # Progress
      if (verbose & n_splits != dplyr::n_distinct(split_labels)) {
        message("Skipped ", dplyr::n_distinct(split_labels) - n_splits, " split label",
                ifelse((dplyr::n_distinct(split_labels) - n_splits) == 1, "", "s"),
                " due to insufficient cells/replicates: ",
                paste0(setdiff(unique(split_labels), names(matrix_list)),
                       collapse = ", "))
      }
    } else if (pseudobulk == "none") {
      keep_indices <- split_indices[lengths(split_indices) >= min_cells_per_split]
      n_splits <- length(keep_indices)
      # Extract matrix
      count_matrix <- .getMatrix(object = object,
                                 use_assay = use_assay,
                                 use_layer = use_layer,
                                 use_cells = use_cells,
                                 verbose = verbose)
      # Progress
      if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Extracting ", n_splits,
                           " pseudobulk ", ifelse(n_splits == 1, "matrix..", "matrices.."))
      # Create matrix list
      matrix_list <- lapply(keep_indices, function(i) {
        keep_genes_count <- Matrix::rowSums(count_matrix[, i, drop = FALSE] > 0) >= min_cells_per_feature
        prop_nonzero <- Matrix::rowMeans((count_matrix[, i, drop = FALSE] > 0))
        keep_genes_prop <- prop_nonzero >= min_prop_cells_per_feature
        keep_genes <- which(keep_genes_count & keep_genes_prop)
        filtered_mat <- count_matrix[keep_genes, i, drop = FALSE]

        # Warn if excluded genes are >10% of all genes
        prop_genes_excluded <- 1 - (length(keep_genes)/nrow(count_matrix) )
        if (prop_genes_excluded > 0.1) {
          message("Warning: Excluded ", round(prop_genes_excluded*100, 2),"% of genes in split ",
                  unique(split_labels[i]))
        }

        return(filtered_mat)
      })

      if (verbose & n_splits != dplyr::n_distinct(split_labels)) {
        message("Skipped ", dplyr::n_distinct(split_labels) - n_splits, " split label",
                ifelse((dplyr::n_distinct(split_labels) - n_splits) == 1, "", "s"),
                " due to insufficient cells/replicates: ",
                paste0(setdiff(unique(split_labels), names(matrix_list)),
                       collapse = ", "))
      }
    }
  }

  # Run through each matrix
  # Downsample when appropriate if parameter force_balance = TRUE
  # Check group composition against min_replicates_per_group
  matrix_list <- lapply(matrix_list, function(m) {
    current_replicates <- colnames(m)
    current_groups <- group_key$group[match(current_replicates, group_key$replicate)]
    replicates_per_group <- table(current_groups)
    # Check whether either group is < min_replicates_per_group or only 1 group is present
    if (any(replicates_per_group < min_replicates_per_group) | (length(replicates_per_group) < 2)) {
      m <- NULL
    } else if (force_balance == TRUE) {
      # If force balancing, downsample larger group to match size of smaller group
      if (any(replicates_per_group > min(replicates_per_group))) {
        n_exclude <- max(replicates_per_group) - min(replicates_per_group)
        downsample_group <- names(replicates_per_group[replicates_per_group == max(replicates_per_group)])
        exclude_indices <- sample(which(current_groups == downsample_group), n_exclude)
        m <- m[, -exclude_replicates, drop = FALSE]
        current_groups <- current_groups[-exclude_replicates]
      }
    }
    return(m)
  })
  # Remove NULL elements
  pre_labels <- names(matrix_list)
  matrix_list <- matrix_list[lengths(matrix_list) > 0]
  if (verbose & any(lengths(matrix_list) == 0)) {
    message("Skipped ", sum(lengths(matrix_list) == 0), " split label",
            ifelse((lengths(matrix_list) == 0) == 1, "", "s"),
            " due to insufficient replicates per group: ",
            paste0(setdiff(pre_labels, names(matrix_list)),
                   collapse = ", "))
  }

  # Create corresponding list of replicates/groups
  target_list <- lapply(matrix_list, FUN = function(i) {
    replicates_i <- colnames(i)
    groups_i <- group_key[replicates_i, c("replicate", "group")]
    return(groups_i)
  })
  names(target_list) <- names(matrix_list)

  # ---------------------------------------------------------------------------
  # Perform pseudobulk differential expression
  # ---------------------------------------------------------------------------

  time3 <- Sys.time()

  # Progress
  if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Running DE on ", length(matrix_list),
                       ifelse(length(matrix_list) == 1, " matrix..", " matrices.."))

  # for each item in matrix_list
  de_results_list <- pbmcapply::pbmclapply(seq_len(length(matrix_list)),
                                           FUN = function(i) {
                                             n_groups <- dplyr::n_distinct(target_list[[i]]$group)
                                             if (n_groups == 2) {
                                               group_factor <- factor(target_list[[i]]$group)
                                               group_factor <- stats::relevel(group_factor, ref = reference_group)
                                               target_list[[i]]$group <- group_factor

                                               design_i <- stats::model.matrix(~ group, data = target_list[[i]])
                                               de_results_i <- switch(de_method,
                                                                      edgeR = .runDE.edgeR(mat = matrix_list[[i]],
                                                                                           targets = target_list[[i]],
                                                                                           design = design_i,
                                                                                           de_test = de_test,
                                                                                           de_params = de_params),
                                                                      DESeq2 = .runDE.DESeq2(mat = matrix_list[[i]],
                                                                                             targets = target_list[[i]],
                                                                                             design = design_i,
                                                                                             de_test = de_test,
                                                                                             de_params = de_params),
                                                                      limma = .runDE.limma(mat = matrix_list[[i]],
                                                                                           targets = target_list[[i]],
                                                                                           design = design_i,
                                                                                           de_test = de_test,
                                                                                           de_params = de_params),
                                                                      presto = .runDE.presto(mat = matrix_list[[i]],
                                                                                             targets = target_list[[i]],
                                                                                             de_test = de_test,
                                                                                             de_params = de_params))
                                               de_results_i <- de_results_i |>
                                                 dplyr::mutate(padj = stats::p.adjust(pvalue, method = p_adjust_method),
                                                               split = names(matrix_list)[i]) |>
                                                 dplyr::arrange(padj)
                                             } else {
                                               de_results_i <- NULL
                                               if (verbose) message("Skipped split label ", names(matrix_list)[i],
                                                                    ", only ", n_groups,
                                                                    " group (", unique(target_list[[i]]$group),
                                                                    ") present.")
                                             }
                                             return(de_results_i)
                                           },
                                           mc.cores = n_cores,
                                           mc.set.seed = TRUE)
  de_results <- do.call(rbind, de_results_list)
  de_results <- de_results |>
    data.frame()

  # ---------------------------------------------------------------------------
  # Wrap up
  # ---------------------------------------------------------------------------

  time4 <- Sys.time()

  # Metadata
  if (pseudobulk == "generate") {
    metadata_list <- list("PB_metadata" = pseudobulk_metadata,
                          "time" = data.frame(total = difftime(time4, time1, units = "secs"),
                                              step1_setup = difftime(time2, time1, units = "secs"),
                                              step2_get_matrices = difftime(time3, time2, units = "secs"),
                                              step3_DE = difftime(time4, time3, units = "secs")))
  } else {
    metadata_list <- list("time" = data.frame(total = difftime(time4, time1, units = "secs"),
                                              step1_setup = difftime(time2, time1, units = "secs"),
                                              step2_get_matrices = difftime(time3, time2, units = "secs"),
                                              step3_DE = difftime(time4, time3, units = "secs")))
  }

  # Parameters
  parameter_list <- list("replicate_labels" = replicate_labels,
                         "group_labels" = group_labels,
                         "split_labels" = split_labels,
                         "use_cells" = use_cells,
                         "pseudobulk" = pseudobulk,
                         "de_method" = de_method,
                         "de_test" = de_test,
                         "de_params" = de_params,
                         "p_adjust_method" = p_adjust_method,
                         "min_cells_per_split" = min_cells_per_split,
                         "min_replicates_per_split" = min_replicates_per_split,
                         "min_replicates_per_group" = min_replicates_per_group,
                         "min_cells_per_feature" = min_cells_per_feature,
                         "min_prop_cells_per_feature" = min_prop_cells_per_feature,
                         "force_balance" = force_balance,
                         "use_assay" = use_assay,
                         "use_layer" = use_layer,
                         "stored_replicates" = stored_replicates,
                         "random_seed" = random_seed,
                         "n_cores" = n_cores)

  # Return output
  if (pseudobulk == "none") {
    return(list("DE_results" = de_results,
                "cell_values" = matrix_list,
                "group_key" = group_key,
                "metadata" = metadata_list,
                "parameters" = parameter_list))
  } else {
    return(list("DE_results" = de_results,
                "PB_values" = matrix_list,
                "group_key" = group_key,
                "metadata" = metadata_list,
                "parameters" = parameter_list))
  }
}


# Run edgeR differential expression ---------------------------
#
# mat -- A feature x replicate pseudobulk matrix or a feature x cell matrix
# targets -- A dataframe containing sample to group key
# design -- A model.matrix design object
# de_test -- Which test to use for differential expression
# de_params -- Additional parameters to pass

.runDE.edgeR <- function(mat,
                         targets,
                         design,
                         de_test = "LRT",
                         de_params = list()) {
  tryCatch({
    # Normalization
    dge <- do.call(edgeR::DGEList, c(list("counts" = mat,
                                          "group" = targets$group),
                                     de_params[["DGEList"]]))
    dge <- do.call(edgeR::calcNormFactors, c(list("object" = dge),
                                             de_params[["calcNormFactors"]]))
    y <- do.call(edgeR::estimateDisp, c(list("y" = dge,
                                             "design" = design),
                                        de_params[["estimateDisp"]]))
    # Run test
    fit <- switch(de_test,
                  QLF = do.call(edgeR::glmQLFit, c(list("y" = y,
                                                        "design" = design),
                                                   de_params[["glmQLFit"]])),
                  LRT = do.call(edgeR::glmFit, c(list("y" = y,
                                                      "design" = design),
                                                 de_params[["glmFit"]])),
                  exact = do.call(edgeR::exactTest, c(list("object" = y),
                                                      de_params[["exactTest"]])))
    test <- switch(de_test,
                   QLF = do.call(edgeR::glmQLFTest, c(list("glmfit" = fit,
                                                           "coef" = 2),
                                                    de_params[["glmQLFTest"]])),
                   LRT = do.call(edgeR::glmLRT, c(list("glmfit" = fit,
                                                       "coef" = 2),
                                                  de_params[["glmLRT"]])),
                   exact = fit)
    # Compile results
    edgeR_results <- edgeR::topTags(object = test,
                                    n = Inf,
                                    adjust.method = "none") |>
      data.frame()
    edgeR_results <- edgeR_results |>
      dplyr::transmute(gene = rownames(edgeR_results),
                       lfc = logFC,
                       pvalue = PValue)
    rownames(edgeR_results) <- NULL
  }, error = function(e) message(e))

  if (!exists("edgeR_results")) {
    edgeR_results <- NULL
  }
  return(edgeR_results)
}

# Run DESeq2 differential expression ---------------------------
#
# mat -- A feature x replicate pseudobulk matrix or a feature x cell matrix
# targets -- A dataframe containing sample to group key (splits to keep)
# design -- A model.matrix design object
# de_test -- Which test to use for differential expression
# de_params -- Additional parameters to pass

.runDE.DESeq2 <- function(mat,
                          targets,
                          design,
                          de_test = "LRT",
                          de_params = list()) {

  .requirePackage("DESeq2", source = "bioc")

  # Construct DESeq2 dataset
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = mat,
    colData = targets,
    design = design
  )

  # Run DESeq
  if (de_test == "LRT") {
    reduced_design <- stats::model.matrix(~ 1, data = targets)
    dds <- do.call(DESeq2::DESeq, c(list("object" = dds,
                                         "test" = "LRT",
                                         "reduced" = reduced_design),
                                    de_params[["DESeq"]]))
  } else if (de_test == "Wald") {
    dds <- do.call(DESeq2::DESeq, c(list("object" = dds,
                                         "test" = "Wald"),
                                    de_params[["DESeq"]]))
  }

  DESeq2_results <- DESeq2::results(dds) |>
    data.frame()
  DESeq2_results <- DESeq2_results |>
    dplyr::transmute(gene = rownames(DESeq2_results),
                     lfc = log2FoldChange,
                     pvalue = pvalue)
  rownames(DESeq2_results) <- NULL

  return(DESeq2_results)
}

# Run limma differential expression ---------------------------
#
# mat -- A feature x replicate pseudobulk matrix or a feature x cell matrix
# targets -- A dataframe containing sample to group key
# design -- A model.matrix design object
# de_test -- Which test to use for differential expression
# de_params -- Additional parameters to pass

.runDE.limma <- function(mat,
                         targets,
                         design,
                         de_test = "voom",
                         de_params = list()) {

  .requirePackage("limma", source = "bioc")

  if (de_test %in% c("trend", "voom")) {
    # Create a DGE list
    dge <- do.call(edgeR::DGEList, c(list("counts" = mat),
                                     de_params[["DGEList"]]))
    # Apply TMM normalization
    dge <- do.call(edgeR::calcNormFactors, c(list("object" = dge),
                                             de_params[["calcNormFactors"]]))
    # Transform
    if (de_test == "trend") {
      # Apply logCPM
      transformed_dge <- do.call(edgeR::cpm, c(list("y" = dge,
                                                    "log" = TRUE),
                                               de_params[["cpm"]]))
    } else if (de_test == "voom") {
      # Apply voom transformation
      transformed_dge <- do.call(limma::voom, c(list("counts" = dge,
                                                     "design" = design),
                                                de_params[["voom"]]))
    }
    # limma DE
    fit <- do.call(limma::lmFit, c(list("object" = transformed_dge,
                                        "design" = design),
                                   de_params[["lmFit"]]))
    fit <- do.call(limma::eBayes, c(list("fit" = fit),
                                    de_params[["eBayes"]]))

    limma_results <- limma::topTable(fit, coef = ncol(design), number = Inf) |>
      data.frame()
    limma_results <- limma_results |>
      dplyr::transmute(gene = rownames(limma_results),
                       lfc = logFC,
                       pvalue = P.Value)
    rownames(limma_results) <- NULL
  } else if (de_test %in% c("wilcox_cpm", "wilcox_log_cpm")) {
    # Normalization
    if (de_test == "wilcox_cpm") {
      cpm_mat <- do.call(edgeR::cpm, c(list("y" = mat),
                                       de_params[["cpm"]]))
    } else if (de_test == "wilcox_log_cpm") {
      cpm_mat <- do.call(edgeR::cpm, c(list("y" = mat,
                                            "log" = TRUE),
                                       de_params[["cpm"]]))
    }
    # Group indices
    group1_indices <- which(targets$group == levels(targets$group)[1])
    group2_indices <- which(targets$group == levels(targets$group)[2])

    # Wilcoxon rank sum test p-values
    pvalues <- apply(cpm_mat, 1,
      FUN = function(x) {
        return(min(2 * min(do.call(limma::rankSumTestWithCorrelation, c(list("index" = group1_indices,
                                                                             "statistics" = x),
                                                                        de_params[["rankSumTestWithCorrelation"]]))), 1))
      }
    )
    # LFC
    if ("pseudocount" %in% de_params[["lfc"]]) {
      pseudocount <- de_params$lfc$pseudocount
    } else {
      pseudocount <- 1
    }
    lfcs <- log2(rowMeans(mat[, group1_indices, drop = FALSE] + pseudocount)/rowMeans(mat[, group2_indices, drop = FALSE] + pseudocount))

    wilcox_results <- data.frame(gene = rownames(cpm_mat),
                                 lfc = lfcs,
                                 pvalue = pvalues)
    rownames(wilcox_results) <- NULL
  }
  return(limma_results)
}

# Run Wilcoxon rank sum test differential expression using presto ---------------------------
#
# mat -- A feature x replicate pseudobulk matrix or a feature x cell matrix
# targets -- A dataframe containing sample to group key
# de_test -- Which test to use for differential expression
# de_params -- Additional parameters to pass to cpm and/or wilcoxauc

.runDE.presto <- function(mat,
                          targets,
                          de_test = "wilcox_cpm",
                          de_params = list()) {

  .requirePackage("presto", installInfo = 'devtools::install_github("immunogenomics/presto")')

  # Normalization
  if (de_test == "wilcox_cpm") {
    cpm_mat <- do.call(edgeR::cpm, c(list("y" = mat),
                                     de_params[["cpm"]]))
  } else if (de_test == "wilcox_log_cpm") {
    cpm_mat <- do.call(edgeR::cpm, c(list("y" = mat,
                                          "log" = TRUE),
                                     de_params[["cpm"]]))
  }
  # Run Wilcoxon rank sum test
  wilcox_results <- do.call(presto::wilcoxauc, c(list("X" = cpm_mat,
                                                      "y" = targets$group),
                                                 de_params[["wilcoxauc"]]))

  wilcox_results <- wilcox_results |>
    dplyr::transmute(gene = feature,
                     lfc = logFC,
                     pvalue = pval)
  rownames(wilcox_results) <- NULL

  return(wilcox_results)
}
