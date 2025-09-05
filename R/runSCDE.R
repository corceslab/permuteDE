#' Run single-cell differential expression between two groups using `Seurat` defaults.
#'
#' This function performs differential expression analysis between two 
#' groups of cells in a Seurat object using Seurat's `FindMarkers()`. 
#'
#' @param object An object of class 'Seurat'.
#' @param group_labels A character string or vector indicating the name of the
#' column containing the two comparison group labels.
#' @param reference_group A string specifying the reference group. Default to 
#' \code{NULL}, in which case the first value in the group column is used as the reference.
#' @param force_balance A boolean indicating if two groups have equal sample size.
#' Default to \code{FALSE}. If TRUE, and the two groups have unequal sample sizes, 
#' the larger group will be randomly downsampled to match the size of the smaller group.
#' @param use_cells A vector of cell names subset to. Default = \code{NULL} will
#' use all cells.
#' @param min_replicates_per_group A numeric value indicating the minimum number
#' of distinct replicates represented within each of the two comparison groups.
#' Pseudobulk expression and differential expression will not be performed for
#' splits with fewer replicates. Defaults to 2.
#' @param min_cells_per_feature A numeric value indicating the minimum number
#' of cells (within a split) with expression of a gene. Pseudobulk expression
#' and differential expression will not be calculated for genes with fewer
#' cells. Defaults to 1.
#' @param use_assay A character string indicating the assay to use in the
#' provided object. Default = \code{NULL} will choose the `RNA` assay.
#' @param use_layer  a character string or vector indicating
#' the layer — previously known as slot — to use in the provided object
#' @param n_cores A numeric value indicating the number of cores to use for
#' parallelization. Default = \code{NULL} will use the number of available cores
#' minus 2.
#' @param verbose A boolean value indicating whether to use verbose output
#' during the execution of this function. Defaults to \code{TRUE}.
#' Can be set to \code{FALSE} for a cleaner output.
#'
#' @return A data.frame of differentially expressed genes from `FindMarkers()`.
#' @export
#'
#' @examples
runSCDE <- function(object,
                    group_labels,
                    reference_group,
                    force_balance = FALSE,
                    use_cells = NULL,
                    min_replicates_per_group = 2,
                    min_cells_per_feature = 1,
                    use_assay = NULL,
                    use_layer = NULL,
                    n_cores = NULL,
                    verbose = TRUE) {
  library(Seurat)
  
  # normalize count matrix if haven't
  if (!'data' %in% Assays(object)) {
    message("The Seurat object has not been normalized yet.")
    object <- NormalizeData(object)
  }
  
  # sanity check
  Idents(object) <- group_labels
  if (length(table(Idents(object))) != 2){
    stop('Group labels are not binary.')
  }
  
  if (any(table(Idents(object)) < min_replicates_per_group)){
    stop('Replicates per group are not sufficient.')
  }
  
  x <- GetAssayData(object, layer = "data") 
  keep_genes <- Matrix::rowSums(x > 0) >= min_cells_per_feature
  object <- subset(object, features = rownames(object)[keep_genes])
  
  # get case group value
  case_group <- setdiff(levels(Idents(object)), reference_group)
  
  # perform DE
  de.markers <- FindMarkers(object, layer = 'data', logfc.threshold = 0,
                            ident.1 = case_group, ident.2 = reference_group)
  
  return(de.markers)
  
}