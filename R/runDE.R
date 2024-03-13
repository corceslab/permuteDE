#' Run pseudobulk differential expression
#'
#' Generate pseudobulk matri(ces) and run pseudobulk differential expression
#' between two groups for each using existing tools: \code{edgeR},
#' \code{DESeq2}, or \code{limma}.
#'
#' This function was inspired by R package neurorestore/Libra and some aspects
#' are adapted therefrom (Squair et al. 2021).
#'
#' @param object An object of class 'Seurat' or 'SingleCellExperiment'.
#' @param replicate_labels A character string or vector indicating the name of the
#' column containing the replicate labels.
#' @param group_labels A character string or vector indicating the name of the
#' column containing the two comparison group labels.
#' @param split_labels A character string or vector indicating the name of a
#' column by which to split the cells prior to pseudobulking and performing
#' differential expression (e.g., cell types). Results will be returned for
#' each unique value in the column indicated by 'split_labels'. Default =
#' \code{NULL} will run pseudobulk differential expression on all cells
#' together.
#' @param use_cells A vector of cell names subset to. Default = \code{NULL} will
#' use all cells.
#' @param min_cells_per_split A numeric value indicating the minimum number of
#' cells within one split. Pseudobulk expression and differential expression
#' will not be performed for splits with fewer cells. Defaults to 4.
#' @param min_replicates_per_split A numeric value indicating the minimum number
#' of distinct replicates represented within one split. Pseudobulk expression
#' and differential expression will not be performed for splits with fewer
#' replicates. Defaults to 4.
#' @param min_replicates_per_group A numeric value indicating the minimum number
#' of distinct replicates represented within each of the two comparison groups.
#' Pseudobulk expression and differential expression will not be performed for
#' splits with fewer replicates. Defaults to 2.
#' @param min_cells_per_feature A numeric value indicating the minimum number
#' of cells (within a split) with expression of a gene. Pseudobulk expression
#' and differential expression will not be calculated for genes with fewer
#' cells. Defaults to 1.
#' @param de_method Which tool to use for differential expression. Permitted
#' values are 'edgeR', 'DESeq2', and 'limma'. Defaults to 'edgeR'.
#' @param de_test Which test to use for differential expression. Defaults to
#' 'LRT'.
#' @param p_adjust_method A string indicating which multiple comparison
#' adjustment to use. For permitted values, see \code{stats::p.adjust}. Defaults
#' to 'fdr' (Benjamini & Hochberg, 1995).
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
#' @return Returns a list containing the following elements: \describe{
#'   \item{DE_results}{Dataframe containing DE results for each feature, by
#'   split}
#'   \item{PB_values}{List of feature x replicate matri(ces) containing
#'   pseudobulk values for each feature, one per split}
#'   \item{group_key}{Dataframe record of group labels corresponding to each replicate}
#'   \item{parameters}{Dataframe record of parameter values used}
#'   }
#'
#' @export
#'
runDE <- function(object,
                  replicate_labels,
                  group_labels,
                  split_labels = NULL,
                  use_cells = NULL,
                  min_cells_per_split = 4,
                  min_replicates_per_split = 4,
                  min_replicates_per_group = 2,
                  min_cells_per_feature = 1,
                  de_method = "edgeR",
                  de_test = "LRT",
                  p_adjust_method = "fdr",
                  use_assay = NULL,
                  use_layer = NULL,
                  n_cores = NULL,
                  verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Check input validity
  # ---------------------------------------------------------------------------

  .validInput(object, "object", "countCombinations")
  .validInput(replicate_labels, "replicate_labels", object)
  .validInput(group_labels, "group_labels", object)
  .validInput(split_labels, "split_labels", object)
  .validInput(use_cells, "use_cells", object)
  .validInput(min_cells_per_split, "min_cells_per_split")
  .validInput(min_replicates_per_split, "min_replicates_per_split")
  .validInput(min_replicates_per_split, "min_replicates_per_group")
  .validInput(min_cells_per_feature, "min_cells_per_feature")
  .validInput(de_method, "de_method")
  .validInput(de_test, "de_test", de_method)
  .validInput(p_adjust_method, "p_adjust_method")
  .validInput(use_assay, "use_assay", object)
  .validInput(use_slot, "use_slot", list(object, use_assay))
  .validInput(n_cores, "n_cores")
  .validInput(verbose, "verbose")

  # Set defaults
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 2
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
  if (dplyr::n_distinct(groups) != 2) {
    stop("Input value '", group_labels,
         "' for parameter 'group_labels' must represent a cell metadata column that contains exactly 2 groups for the selected cells, please supply valid input!")
  }

  # ---------------------------------------------------------------------------
  # Calculate pseudobulk values
  # ---------------------------------------------------------------------------

  # Pseudobulk matri(ces)
  pb_list <- getPseudobulk(object = object,
                           replicate_labels = replicate_labels,
                           split_labels = split_labels,
                           use_cells = use_cells,
                           min_cells_per_split = min_cells_per_split,
                           min_replicates_per_split = min_replicates_per_split,
                           min_cells_per_feature = min_cells_per_feature,
                           n_cores = n_cores,
                           verbose = verbose)

  # Group labels for each replicate, by split
  group_key <- data.frame(replicate = paste0("rep_", replicates),
                          group = groups) %>%
    dplyr::group_by(replicate, group) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::select(-n) %>%
    data.frame()
  rownames(group_key) <- group_key$replicate
  remove_splits <- c()
  target_list <- lapply(pb_list, FUN = function(i) {
    replicates_i <- colnames(i)
    groups_i <- group_key[replicates_i, c("replicate", "group")]
    # Remove splits with fewer than required number of replicates per group
    if (min(table_groups_i) < min_replicates_per_group) {
      remove_splits <- c(remove_splits, i)
    }
  })
  # Remove splits with fewer than required number of replicates per group
  if (length(remove_splits) > 0) {
    pb_list[[remove_splits]] <- NULL
    target_list[[remove_splits]] <- NULL
    message("Skipped ", length(remove_splits), " split label",
            ifelse(length(remove_splits) == 1, "", "s"),
            " due to insufficient replicates per group: ",
            paste0(remove_splits,
                   collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # Perform pseudobulk differential expression
  # ---------------------------------------------------------------------------

  de_results_list <- pbmcapply::pbmclapply(seq_len(length(pb_list)),
                                           FUN = function(i) {
                                             n_groups <- dplyr::n_distinct(target_list[[i]]$group)
                                             if (n_groups == 2) {
                                               design_i <- stats::model.matrix(~ group, data = target_list[[i]])
                                               de_results_i <- switch(de_method,
                                                                      edgeR = .runDE.edgeR(pseudobulk = pb_list[[i]],
                                                                                          targets = target_list[[i]],
                                                                                          design = design_i,
                                                                                          de_test = de_test),
                                                                      DESeq2 = .runDE.DESeq2(pseudobulk = pb_list[[i]],
                                                                                            targets = target_list[[i]],
                                                                                            de_test = de_test),
                                                                      limma = .runDE.limma(pseudobulk = pb_list[[i]],
                                                                                          targets = target_list[[i]],
                                                                                          design = design_i,
                                                                                          de_test = de_test))
                                               de_results_i <- de_results_i %>%
                                                 dplyr::mutate(p_adjust = stats::p.adjust(p_value, method = p_adjust_method),
                                                               split = names(pb_list)[i])
                                             } else {
                                               de_results_i <- NULL
                                               if (verbose) message("Skipped split label ", names(pb_list)[i],
                                                                    ", only ", n_groups,
                                                                    " group (", unique(target_list[[i]]$group),
                                                                    ") present.")
                                             }
                                             return(de_results_i)
                                           },
                                           mc.cores = n_cores,
                                           mc.set.seed = TRUE)
  de_results <- do.call(rbind, de_results_list)
  de_results <- de_results %>% data.frame() %>% dplyr::arrange(split, p_adjust)

  # ---------------------------------------------------------------------------
  # Wrap up
  # ---------------------------------------------------------------------------

  parameter_list <- list("replicate_labels" = replicate_labels,
                         "group_labels" = group_labels,
                         "split_labels" = split_labels,
                         "use_cells" = use_cells,
                         "min_cells_per_split" = min_cells_per_split,
                         "min_replicates_per_split" = min_replicates_per_split,
                         "min_replicates_per_group" = min_replicates_per_group,
                         "min_cells_per_feature" = min_cells_per_feature,
                         "de_method" = de_method,
                         "de_test" = de_test,
                         "p_adjust_method" = p_adjust_method)

  # Return output
  return(list("DE_results" = de_results,
              "PB_values" = pb_list,
              "group_key" = group_key,
              "parameters" = parameter_list))
}


# Run edgeR differential expression ---------------------------
#
# pseudobulk -- A feature x replicate pseudobulk matrix
# targets -- A dataframe containing sample to group key
# design -- A model.matrix design object
# de_test -- Which test to use for differential expression

.runDE.edgeR <- function(pseudobulk,
                        targets,
                        design,
                        de_test = "LRT") {
  tryCatch({
    y <- edgeR::DGEList(counts = pseudobulk, group = targets$group) %>%
      edgeR::calcNormFactors(method = 'TMM') %>%
      edgeR::estimateDisp(design)
    fit <- switch(de_test,
                  QLF = edgeR::glmQLFit(y, design),
                  LRT = edgeR::glmFit(y, design = design))
    test <- switch(de_test,
                   QLF = edgeR::glmQLFTest(fit, coef = 2),
                   LRT = edgeR::glmLRT(fit))
    edgeR_results <- edgeR::topTags(object = test,
                                    n = Inf,
                                    adjust.method = "none") %>%
      data.frame()
    edgeR_results <- edgeR_results %>%
      dplyr::transmute(gene = rownames(edgeR_results),
                       lfc = logFC,
                       p_value = PValue)
  }, error = function(e) message(e))

  if (!exists("edgeR_results")) {
    edgeR_results <- NULL
  }
  return(edgeR_results)
}

# Run DESeq2 differential expression ---------------------------
#
# pseudobulk -- A feature x replicate pseudobulk matrix
# targets -- A dataframe containing sample to group key
# design -- A model.matrix design object
# de_test -- Which test to use for differential expression

runDE.DESeq2 <- function(pseudobulk,
                         targets,
                         design,
                         de_test = "LRT") {
  NULL
}

# Run limma differential expression ---------------------------
#
# pseudobulk -- A feature x replicate pseudobulk matrix
# targets -- A dataframe containing sample to group key
# design -- A model.matrix design object
# de_test -- Which test to use for differential expression

runDE.limma <- function(pseudobulk,
                        targets,
                        design,
                        de_test = "voom") {
  NULL
}

