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
#' @param metadata An optional dataframe containing relevant metadata columns
#' corresponding to the data provided to parameter \code{object}. Default =
#' \code{NULL} looks for metadata in \code{object} or other provided inputs.
#' @param replicate_labels A string indicating the name of the
#' metadata column containing the biological replicate labels or a vector
#' containing the biological replicate labels in order. For pseudobulk DE
#' analysis, the biological replicate labels are used to construct/define the
#' pseudobulks. Input is not required for cell-level DE analysis.
#' @param group_labels A string indicating the name of the metadata column
#' containing the two comparison group labels or a vector containing the
#' comparison labels in order.
#' @param split_labels A string indicating the name of a metadata column by
#' which to split the cells prior to pseudobulking and performing differential
#' expression (e.g., cell types). Alternately, a vector containing the split
#' labels for each cell in order. Results will be returned for each unique value
#' indicated by \code{split_labels}. Default = \code{NULL} will run pseudobulk
#' differential expression on all cells together.
#' @param reference_group A string specifying the reference group. Defaults to
#' \code{NULL}, in which case the first value alphabetically is used as the
#' reference.
#' @param design An optional string specifying a model formula for more complex
#' designs. Last term in formula must correspond to group labels. Default =
#' \code{NULL} will run a pairwise group comparison (~ group) based on the input
#' provided to parameter \code{group_labels}.
#' @param use_cells A vector of cell names to subset the object to prior to
#' subsequent pseudobulk and differential expression steps. Default =
#' \code{NULL} will use all cells.
#' @param pseudobulk A string indicating pseudobulk handling. Permitted values
#' are: "generate" (pseudobulk matrices will be generated), "supplied"
#' (pseudobulk matrix was supplied by the user to parameter \code{object}), or
#' "none" (pseudobulking will not be used, cell-level differential expression
#' analysis will be run). Defaults to "generate".
#' @param de_method Which tool to use for differential expression analysis.
#' Permitted values are "edgeR", "DESeq2", "limma", "presto", and "BPCells".
#' Defaults to "edgeR".
#' @param de_test Which test to use for differential expression analysis.
#' Available values are dependent on the \code{de_method}: "edgeR" ("LRT",
#' "QLF", "exact"), "DESeq2" ("LRT", "Wald"), "limma" ("trend", "voom",
#' "wilcox_cpm", "wilcox_log_cpm"), "presto" ("wilcox_cpm", "wilcox_log_cpm"),
#' and "BPCells" ("wilcox_cpm", "wilcox_log_cpm"). Defaults to "LRT".
#' @param de_params A list of lists containing additional parameters to be
#' passed to specific DE functions. The name of each element must be the
#' specific DE function to which those parameters are passed. Defaults to an
#' empty list.
#' @param return_raw_de A Boolean value indicating whether to also return the
#' raw output from the selected DE method/test. Defaults to \code{FALSE}.
#' @param normalize_prefilter A Boolean value indicating whether
#' normalization should be applied before (\code{TRUE}) or after (\code{FALSE})
#' filtering out features with low counts. Defaults to \code{FALSE}.
#' @param p_adjust_method A string indicating which multiple comparison
#' adjustment to use. For permitted values, see \code{stats::p.adjust.methods}.
#' Defaults to "fdr" (Benjamini & Hochberg, 1995).
#' @param min_cells_per_split A numeric value indicating the minimum number of
#' cells within one split. Pseudobulk and differential expression steps will not
#' be performed for splits with fewer cells. Defaults to 100.
#' @param min_cells_per_replicate A numeric value indicating the minimum number
#' of cells within one replicate for one split. Pseudobulk steps will not be
#' performed for replicates with fewer cells for that split. Defaults to 10.
#' @param min_replicates_per_split A numeric value indicating the minimum number
#' of distinct replicates represented within one split. Pseudobulk expression
#' and differential expression will not be performed for splits with fewer
#' replicates. Defaults to 6.
#' @param min_replicates_per_group A numeric value indicating the minimum number
#' of distinct replicates represented within each of the two comparison groups.
#' Pseudobulk and differential expression steps will not be performed for
#' splits with fewer replicates. Defaults to 3.
#' @param min_cells_per_feature A numeric value indicating the minimum number
#' of cells (within a split) with expression of a feature. Pseudobulk and
#' differential expression will not be calculated for features expressed in
#' fewer cells. Defaults to 10.
#' @param min_prop_cells_per_feature A numeric value indicating the minimum
#' proportion of cells (within a split) with expression of a feature. Pseudobulk
#' and differential expression will not be calculated for features expressed in
#' fewer cells. Defaults to 0.1.
#' @param force_balance A boolean indicating whether to force the two comparison
#' groups to have the same sample size. Defaults to \code{FALSE}. If
#' \code{TRUE}, the larger group will be randomly downsampled to the size of the
#' smaller group.
#' @param use_assay A string indicating the assay to use in the
#' provided object. Default = \code{NULL} will choose the current active assay
#' for \code{Seurat} objects and the \code{counts} assay for
#' \code{SingleCellExperiment} objects.
#' @param use_layer For \code{Seurat} objects, a string or vector
#' indicating the layer (previously known as slot) to use in the provided
#' object. Default = \code{NULL} will use the \code{counts} layer.
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
#'   \item{metadata}{List recording characteristics of the data and runtime}
#'   \item{parameters}{List recording parameter values used}
#'   }
#'
#' @export
#'
runDE <- function(object,
                  metadata = NULL,
                  replicate_labels = NULL,
                  group_labels,
                  split_labels = NULL,
                  reference_group = NULL,
                  design = NULL,
                  use_cells = NULL,
                  pseudobulk = "generate",
                  de_method = "edgeR",
                  de_test = "LRT",
                  de_params = list(),
                  return_raw_de = FALSE,
                  normalize_prefilter = FALSE,
                  p_adjust_method = "fdr",
                  min_cells_per_split = 100,
                  min_cells_per_replicate = 10,
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

  .validInput(input = object,
              name = "object",
              class = c("Seurat", "SingleCellExperiment", "matrix", "Matrix"))
  .validInput(input = metadata,
              name = "metadata",
              null_allowed = TRUE,
              class = "data.frame",
              other = object)
  .validInput(input = pseudobulk,
              name = "pseudobulk",
              class = "character",
              len = 1,
              caller = "runDE",
              other = object)
  .validInput(input = replicate_labels,
              name = "replicate_labels",
              null_allowed = pseudobulk == "none",
              class = c("character", "factor", "numeric"),
              other = list(object, metadata, pseudobulk))
  .validInput(input = group_labels,
              name = "group_labels",
              class = c("character", "factor", "numeric", "logical"),
              other = list(object, metadata))
  .validInput(input = split_labels,
              name = "split_labels",
              null_allowed = TRUE,
              class = c("character", "factor", "numeric", "logical"),
              other = list(object, metadata, "runDE"))
  .validInput(input = use_cells,
              name = "use_cells",
              null_allowed = TRUE,
              class = "character",
              other = list(object, pseudobulk))
  .validInput(input = return_raw_de,
              name = "return_raw_de",
              class = "logical",
              len = 1)
  .validInput(input = reference_group,
              name = "reference_group",
              null_allowed = TRUE,
              len = 1,
              other = list(object, metadata, group_labels, use_cells))
  .validInput(input = normalize_prefilter,
              name = "normalize_prefilter",
              class = "logical",
              len = 1)
  .validInput(input = p_adjust_method,
              name = "p_adjust_method",
              class = "character",
              len = 1)
  .validInput(input = min_cells_per_split,
              name = "min_cells_per_split",
              class = "numeric",
              len = 1,
              other = pseudobulk)
  .validInput(input = min_cells_per_replicate,
              name = "min_cells_per_replicate",
              class = "numeric",
              len = 1,
              caller = "runDE",
              other = pseudobulk)
  .validInput(input = min_replicates_per_split,
              name = "min_replicates_per_split",
              class = "numeric",
              len = 1,
              caller = "runDE",
              other = pseudobulk)
  .validInput(input = min_replicates_per_group,
              name = "min_replicates_per_group",
              class = "numeric",
              len = 1)
  .validInput(input = min_cells_per_feature,
              name = "min_cells_per_feature",
              class = "numeric",
              len = 1,
              other = pseudobulk)
  .validInput(input = min_prop_cells_per_feature,
              name = "min_prop_cells_per_feature",
              class = "numeric",
              len = 1,
              other = pseudobulk)
  .validInput(input = force_balance,
              name = "force_balance",
              class = "logical",
              len = 1,
              other = pseudobulk)
  .validInput(input = use_assay,
              name = "use_assay",
              null_allowed = TRUE,
              class = "character",
              len = 1,
              other = object)
  .validInput(input = use_layer,
              name = "use_layer",
              null_allowed = TRUE,
              class = "character",
              len = 1,
              other = list(object, use_assay))
  .validInput(input = de_method,
              name = "de_method",
              class = "character",
              len = 1,
              caller = "runDE",
              other = list(pseudobulk, object, use_assay, use_layer, use_cells))
  .validInput(input = de_test,
              name = "de_test",
              class = "character",
              len = 1,
              other = de_method)
  .validInput(input = de_params,
              name = "de_params",
              class = "list",
              other = list(de_method, de_test))
  .validInput(input = design,
              name = "design",
              null_allowed = TRUE,
              class = "character",
              len = 1,
              other = list(object, metadata, group_labels, de_test))
  .validInput(input = random_seed,
              name = "random_seed",
              class = "numeric",
              len = 1)
  .validInput(input = n_cores,
              name = "n_cores",
              null_allowed = TRUE,
              class = "numeric",
              len = 1)
  .validInput(input = verbose,
              name = "verbose",
              class = "logical",
              len = 1)

  # ---------------------------------------------------------------------------
  # Set up
  # ---------------------------------------------------------------------------

  # Object type
  if (methods::is(object, "Seurat")) {
    object_type <- "Seurat"
  } else if (methods::is(object, "SingleCellExperiment")) {
    object_type <- "SingleCellExperiment"
  } else {
    object_type <- "matrix"
  }

  # Set defaults
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 2)
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
  if (pseudobulk == "none") {
    replicates <- use_cells
  } else if (!is.null(replicate_labels)) {
    if (length(replicate_labels) == 1) {
      replicates <- .retrieveData(object = object,
                                  metadata = metadata,
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
  }
  # Check for NA values
  if (any(is.na(replicates))) {
    stop("Values provided for 'replicate_labels' cannot be NA.")
  }
  replicates <- as.character(replicates)

  # Group labels
  if (length(group_labels) == 1) {
    groups <- .retrieveData(object = object,
                            metadata = metadata,
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
  # Check for NA values
  if (any(is.na(groups))) {
    stop("Values provided for 'group_labels' cannot be NA.")
  }
  if (!methods::is(groups, "character") & !methods::is(groups, "factor")) {
    groups <- as.character(groups)
  }

  # There must be exactly two comparison groups
  if (dplyr::n_distinct(groups) != 2) {
    stop("Input value for parameter 'group_labels' must represent a cell metadata column (or a vector of group labels) that contains exactly 2 groups for the selected data, please supply valid input!")
  }

  # Reference group
  if (is.null(reference_group)) {
    reference_group <- sort(unique(groups))[1]
  }
  non_reference_group <- unique(groups)[unique(groups) != reference_group]

  # Group labels for each replicate
  # Returns a dataframe with two columns: replicate, group
  if (pseudobulk == "generate") {
    replicate_prefix <- "rep_"
  } else {
    replicate_prefix <- ""
  }
  group_key <- data.frame(replicate = paste0(replicate_prefix, replicates),
                          group = groups)

  # Additional variables when design is provided
  if (!is.null(design)) {
    # Replace last term of design formula string with "group"
    # Convert to formula and extract terms
    design_formula <- stats::as.formula(sub(" [^ ]+$", " group", design))
    terms <- attr(terms(design_formula), "term.labels")
    # Remove last term (because that's refers to groups, so it's already in group_key)
    terms <- terms[-length(terms)]
    # Break up interaction terms
    terms <- unique(unlist(strsplit(terms, ":", fixed = TRUE)))
    # Remove term "replicate" if present (already in group_key)
    terms <- terms[terms != "replicate"]
    # Add data for each term to group_key
    if (length(terms) > 0) {
      for (t in 1:length(terms)) {
        term_t <- .retrieveData(object = object,
                                metadata = metadata,
                                type = "cell_metadata",
                                name = terms[t],
                                use_cells = use_cells)
        # Check for NA values
        if (any(is.na(term_t))) {
          stop("Values provided for covariate '", terms[t], "' cannot be NA.")
        }
        if (length(term_t) != length(groups)) {
          stop("Input provided to parameter 'design' requires metadata column '", terms[t], "', but it is the wrong length. Please supply valid input!")
        } else {
          group_key <- cbind(group_key, tmp = term_t)
          colnames(group_key)[ncol(group_key)] <- terms[t]
        }
      }
    }
  } else {
    design_formula <- stats::as.formula("~ group")
  }

  group_key <- group_key |>
    dplyr::group_by(dplyr::across(dplyr::everything())) |>
    dplyr::summarise(n = dplyr::n()) |>
    dplyr::select(-n) |>
    data.frame()

  # Check for duplicated replicates across groupings
  if (any(duplicated(group_key$replicate))) {
    if (!is.null(design)) {
      stop("Input to parameter 'replicate_labels' must be the atomic grouping unit, no replicate value may appear in multiple groups/partitions.")
    } else {
      stop("Input to parameter 'replicate_labels' must be the atomic grouping unit, no replicate value may appear in multiple groups.")
    }
  }

  rownames(group_key) <- group_key$replicate

  # ---------------------------------------------------------------------------
  # Obtain matrices for DE
  # ---------------------------------------------------------------------------

  time2 <- Sys.time()

  exclude_features <- NULL

  if (pseudobulk %in% c("generate", "none")) {
    # Generate pseudobulk matri(ces) or filter count matrix
    # Returns a list containing one matrix per split
    output_list <- getPseudobulk(object = object,
                                 metadata = metadata,
                                 replicate_labels = replicate_labels,
                                 split_labels = split_labels,
                                 use_cells = use_cells,
                                 min_cells_per_split = min_cells_per_split,
                                 min_cells_per_replicate = min_cells_per_replicate,
                                 min_replicates_per_split = min_replicates_per_split,
                                 min_cells_per_feature = min_cells_per_feature,
                                 min_prop_cells_per_feature = min_prop_cells_per_feature,
                                 filter = !normalize_prefilter,
                                 pseudobulk = pseudobulk,
                                 use_assay = use_assay,
                                 use_layer = use_layer,
                                 n_cores = n_cores,
                                 verbose = verbose)
    if (pseudobulk == "generate") {
      matrix_list <- output_list[["PB_values"]]
    } else {
      matrix_list <- output_list[["cell_values"]]
    }
    feature_metrics <- output_list[["metadata"]][["metrics"]]
    if (normalize_prefilter == TRUE) {
      exclude_features <- output_list[["metadata"]][["exclude_features"]]
    }
    rm(output_list)
    gc(verbose = FALSE)
  } else {
    # If necessary, separate the supplied pseudobulk matrix by split
    # Returns a list containing one pseudobulk matrix (feature x replicate) per split
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
    keep_indices <- split_indices[lengths(split_indices) >= min_replicates_per_split]
    n_splits <- length(keep_indices)
    # Progress
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Extracting ", n_splits,
                         " pseudobulk ", ifelse(n_splits == 1, "matrix..", "matrices.."))
    # Create matrix list & remove features with 0 counts in a split
    matrix_list <- lapply(keep_indices, FUN = function(i) {
      object[Matrix::rowSums(object[,i]) > 0, i, drop = FALSE]
    })
    # Progress
    if (verbose & n_splits != dplyr::n_distinct(split_labels)) {
      message("Skipped ", dplyr::n_distinct(split_labels) - n_splits, " split label",
              ifelse((dplyr::n_distinct(split_labels) - n_splits) == 1, "", "s"),
              " due to insufficient cells/replicates: ",
              paste0(setdiff(unique(split_labels), names(matrix_list)),
                     collapse = ", "))
    }
  }

  # Run through each matrix
  # Downsample when appropriate if parameter force_balance = TRUE
  # Check group composition against min_replicates_per_group
  pre_labels <- names(matrix_list)
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
        m <- m[, -exclude_indices, drop = FALSE]
        current_groups <- current_groups[-exclude_indices]
      }
    }
    return(m)
  })
  # Remove NULL elements
  keep_matrix <- lengths(matrix_list) > 0
  skipped_labels <- names(matrix_list)[!keep_matrix]
  matrix_list <- matrix_list[keep_matrix]
  if (verbose && length(skipped_labels) > 0) {
    message("Skipped ", length(skipped_labels), " split label",
            ifelse(length(skipped_labels) == 1, "", "s"),
            " due to insufficient replicates per group: ",
            paste0(skipped_labels, collapse = ", "))
  }

  # Create corresponding list of replicates/groups
  target_list <- lapply(matrix_list, FUN = function(i) {
    replicates_i <- colnames(i)
    groups_i <- group_key[replicates_i, ]
    return(groups_i)
  })
  names(target_list) <- names(matrix_list)

  # ---------------------------------------------------------------------------
  # Perform pseudobulk differential expression
  # ---------------------------------------------------------------------------

  time3 <- Sys.time()

  # Progress
  proceed <- TRUE
  if (length(matrix_list) == 0) {
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : No matrices had sufficient cells/replicates, DE will not be run.")
    warning(" No matrices had sufficient cells/replicates, DE was not run.")
    proceed <- FALSE
  } else if (verbose) {
    message(format(Sys.time(), "%Y-%m-%d %X"), " : Running DE on ", length(matrix_list),
            ifelse(length(matrix_list) == 1, " matrix..", " matrices.."))
  }

  if (proceed == TRUE) {
    # Define function to call correct DE method
    .callDE <- function(i) {
      worker_warnings <- character()
      raw_results_i <- NULL

      de_results_i <- withCallingHandlers(
        {
          n_groups <- dplyr::n_distinct(target_list[[i]]$group)
          if (n_groups == 2) {
            group_factor <- factor(target_list[[i]]$group)
            group_factor <- stats::relevel(group_factor, ref = reference_group)
            target_list[[i]]$group <- group_factor

            design_i <- stats::model.matrix(design_formula, data = target_list[[i]])
            exclude_features_i <- if (is.null(exclude_features)) {
              NULL
            } else {
              exclude_features[[i]]
            }

            de_output_i <- switch(de_method,
                                   edgeR = .runDE.edgeR(mat = matrix_list[[i]],
                                                        targets = target_list[[i]],
                                                        design = design_i,
                                                        de_test = de_test,
                                                        de_params = de_params,
                                                        normalize_prefilter = normalize_prefilter,
                                                        exclude_features = exclude_features_i),
                                   DESeq2 = .runDE.DESeq2(mat = matrix_list[[i]],
                                                          targets = target_list[[i]],
                                                          design = design_i,
                                                          design_formula = design_formula,
                                                          de_test = de_test,
                                                          de_params = de_params,
                                                          normalize_prefilter = normalize_prefilter,
                                                          exclude_features = exclude_features_i),
                                   limma = .runDE.limma(mat = matrix_list[[i]],
                                                        targets = target_list[[i]],
                                                        design = design_i,
                                                        de_test = de_test,
                                                        de_params = de_params,
                                                        normalize_prefilter = normalize_prefilter,
                                                        exclude_features = exclude_features_i),
                                   presto = .runDE.presto(mat = matrix_list[[i]],
                                                          targets = target_list[[i]],
                                                          de_test = de_test,
                                                          de_params = de_params,
                                                          normalize_prefilter = normalize_prefilter,
                                                          exclude_features = exclude_features_i,
                                                          non_reference_group = non_reference_group),
                                   BPCells = .runDE.BPCells(mat = matrix_list[[i]],
                                                            targets = target_list[[i]],
                                                            de_test = de_test,
                                                            de_params = de_params,
                                                            normalize_prefilter = normalize_prefilter,
                                                            exclude_features = exclude_features_i,
                                                            non_reference_group = non_reference_group))
            de_results_i <- de_output_i$results |>
              dplyr::mutate(padj = stats::p.adjust(pvalue, method = p_adjust_method),
                            split = names(matrix_list)[i]) |>
              dplyr::arrange(padj)
            raw_results_i <- de_output_i$raw_results
          } else {
            de_results_i <- NULL
            raw_results_i <- NULL

            if (verbose) {
              message("Skipped split label ", names(matrix_list)[i],
                      ", only ", n_groups,
                      " group (", unique(target_list[[i]]$group),
                      ") present.")
            }
          }
          de_results_i
        },
        warning = function(w) {
          worker_warnings <<- c(worker_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      return(list(results = de_results_i,
                  raw_results = raw_results_i,
                  warnings = worker_warnings))
    }

    # Use pblapply vs pbmcapply depending on # of cores
    if (n_cores == 1) {
      de_output_list <- lapply(
        X = seq_len(length(matrix_list)),
        FUN = .callDE
      )
    } else {
      de_output_list <- pbmcapply::pbmclapply(
        X = seq_len(length(matrix_list)),
        FUN = .callDE,
        mc.cores = n_cores,
        mc.set.seed = TRUE
      )
    }

    # Check that every worker returned the expected wrapper structure
    bad_outputs <- !vapply(de_output_list,
                           FUN = function(x) {
                             is.list(x) &&
                               all(c("results", "raw_results", "warnings") %in% names(x))
                           },
                           FUN.VALUE = logical(1))
    if (any(bad_outputs)) {
      stop("DE failed for one or more splits. ",
           "At least one worker returned an unexpected result. ",
           "Try rerunning with n_cores = 1 for a clearer error message.")
    }

    # Extract the DE results
    de_results_list <- lapply(de_output_list, `[[`, "results")

    # If requested, extract raw DE results
    if (return_raw_de) {
      raw_de_results_list <- lapply(de_output_list, `[[`, "raw_results")
      names(raw_de_results_list) <- names(matrix_list)
      raw_de_results_list <- raw_de_results_list[!vapply(raw_de_results_list, is.null, logical(1))]
    }

    # Extract and combine warnings captured inside each worker
    de_warnings <- unlist(lapply(de_output_list, `[[`, "warnings"), use.names = FALSE)
    # Report unique warnings if verbose output is enabled
    if (verbose && length(de_warnings) > 0) {
      warning("Warnings occurred during DE:\n",
              paste(unique(de_warnings), collapse = "\n"))
    }

    # Check that each split returned either NULL or a dataframe
    bad_results <- !vapply(de_results_list,
                           FUN = function(x) is.null(x) || is.data.frame(x),
                           FUN.VALUE = logical(1))
    if (any(bad_results)) {
      stop("DE failed for one or more splits. ",
           "At least one worker returned a non-data.frame result. ",
           "Try rerunning with n_cores = 1 for a clearer error message.")
    }

    # Remove skipped results
    de_results_list <- Filter(Negate(is.null), de_results_list)

    # Combine split-specific DE tables into one result dataframe
    if (length(de_results_list) == 0) {
      de_results <- NULL
    } else {
      de_results <- do.call(rbind, de_results_list)
      de_results <- data.frame(de_results)
    }
  } else {
    de_results <- NULL
    raw_de_results_list <- NULL
  }

  # ---------------------------------------------------------------------------
  # Wrap up
  # ---------------------------------------------------------------------------

  time4 <- Sys.time()

  # Metadata
  if (pseudobulk == "generate") {
    metadata_list <- list("group_key" = group_key,
                          "feature_metrics" = feature_metrics,
                          "time" = data.frame(total = difftime(time4, time1, units = "secs"),
                                              step1_setup = difftime(time2, time1, units = "secs"),
                                              step2_get_matrices = difftime(time3, time2, units = "secs"),
                                              step3_DE = difftime(time4, time3, units = "secs")))
  } else {
    metadata_list <- list("group_key" = group_key,
                          "time" = data.frame(total = difftime(time4, time1, units = "secs"),
                                              step1_setup = difftime(time2, time1, units = "secs"),
                                              step2_get_matrices = difftime(time3, time2, units = "secs"),
                                              step3_DE = difftime(time4, time3, units = "secs")))
  }

  # Parameters
  parameter_list <- list("object_type" = object_type,
                         "metadata_provided" = !is.null(metadata),
                         "replicate_labels" = replicate_labels,
                         "group_labels" = group_labels,
                         "split_labels" = split_labels,
                         "reference_group" = reference_group,
                         "non_reference_group" = non_reference_group,
                         "design" = design,
                         "design_formula" = design_formula,
                         "use_cells" = use_cells,
                         "pseudobulk" = pseudobulk,
                         "de_method" = de_method,
                         "de_test" = de_test,
                         "de_params" = de_params,
                         "return_raw_de" = return_raw_de,
                         "normalize_prefilter" = normalize_prefilter,
                         "p_adjust_method" = p_adjust_method,
                         "min_cells_per_split" = min_cells_per_split,
                         "min_cells_per_replicate" = min_cells_per_replicate,
                         "min_replicates_per_split" = min_replicates_per_split,
                         "min_replicates_per_group" = min_replicates_per_group,
                         "min_cells_per_feature" = min_cells_per_feature,
                         "min_prop_cells_per_feature" = min_prop_cells_per_feature,
                         "force_balance" = force_balance,
                         "use_assay" = use_assay,
                         "use_layer" = use_layer,
                         "random_seed" = random_seed,
                         "n_cores" = n_cores)

  # Return output
  if (return_raw_de) {
    if (pseudobulk == "none") {
      return(list("DE_results" = de_results,
                  "raw_DE_results" = raw_de_results_list,
                  "cell_values" = matrix_list,
                  "metadata" = metadata_list,
                  "parameters" = parameter_list))
    } else {
      return(list("DE_results" = de_results,
                  "raw_DE_results" = raw_de_results_list,
                  "PB_values" = matrix_list,
                  "metadata" = metadata_list,
                  "parameters" = parameter_list))
    }
  } else {
    if (pseudobulk == "none") {
      return(list("DE_results" = de_results,
                  "cell_values" = matrix_list,
                  "metadata" = metadata_list,
                  "parameters" = parameter_list))
    } else {
      return(list("DE_results" = de_results,
                  "PB_values" = matrix_list,
                  "metadata" = metadata_list,
                  "parameters" = parameter_list))
    }
  }
}


# Run edgeR differential expression ---------------------------
#
# mat                 -- A feature x replicate pseudobulk matrix or a feature x cell matrix
# targets             -- A dataframe containing sample to group key
# design              -- A model.matrix design object
# de_test             -- Which test to use for differential expression
# de_params           -- Additional parameters to pass
# normalize_prefilter -- Whether to get normalization/size factors before filtering out features
# exclude_features    -- A vector of feature names to filter out if normalize_prefilter is TRUE

.runDE.edgeR <- function(mat,
                         targets,
                         design,
                         de_test = "LRT",
                         de_params = list(),
                         normalize_prefilter = FALSE,
                         exclude_features = NULL) {

  # Create edgeR object
  dge <- do.call(edgeR::DGEList, c(list("counts" = mat,
                                        "group" = targets$group),
                                   de_params[["DGEList"]]))
  # Get size factors
  dge <- do.call(edgeR::calcNormFactors, c(list("object" = dge),
                                           de_params[["calcNormFactors"]]))
  # If filtering, do so
  if (normalize_prefilter & !is.null(exclude_features)) {
    mat <- mat[!(rownames(mat) %in% exclude_features),]
    dge_filtered <- do.call(edgeR::DGEList, c(list("counts" = mat,
                                                   "group" = targets$group),
                                              de_params[["DGEList"]]))
    dge_filtered$samples <- dge$samples
    dge <- dge_filtered
  }
  # Estimate dispersion
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
                 QLF = do.call(edgeR::glmQLFTest, c(list("glmfit" = fit),
                                                    de_params[["glmQLFTest"]])),
                 LRT = do.call(edgeR::glmLRT, c(list("glmfit" = fit),
                                                de_params[["glmLRT"]])),
                 exact = fit)
  # Compile results
  raw_results <- edgeR::topTags(object = test,
                                n = Inf,
                                adjust.method = "none") |>
    data.frame()
  edgeR_results <- raw_results |>
    dplyr::transmute(feature = rownames(raw_results),
                     lfc = logFC,
                     pvalue = PValue)
  rownames(edgeR_results) <- NULL

  return(list("results" = edgeR_results,
              "raw_results" = raw_results))
}

# Run DESeq2 differential expression ---------------------------
#
# mat                 -- A feature x replicate pseudobulk matrix or a feature x cell matrix
# targets             -- A dataframe containing sample to group key (splits to keep)
# design              -- A model.matrix design object
# design_formula      -- The design formula
# de_test             -- Which test to use for differential expression
# de_params           -- Additional parameters to pass
# normalize_prefilter -- Whether to get normalization/size factors before filtering out features
# exclude_features    -- A vector of feature names to filter out if normalize_prefilter is TRUE

.runDE.DESeq2 <- function(mat,
                          targets,
                          design,
                          design_formula,
                          de_test = "LRT",
                          de_params = list(),
                          normalize_prefilter = FALSE,
                          exclude_features = NULL) {

  .requirePackage("DESeq2", source = "bioc")

  # Default DESeq2 to quiet output
  if (is.null(de_params[["DESeq"]])) {
    de_params[["DESeq"]] <- list()
  }
  if (is.null(de_params[["DESeq"]][["quiet"]])) {
    de_params[["DESeq"]][["quiet"]] <- TRUE
  }

  # Construct DESeq2 dataset
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = mat,
                                        colData = targets,
                                        design = design)

  # Estimate size factors
  dds <- do.call(DESeq2::estimateSizeFactors, c(list(object = dds),
                                                de_params[["estimateSizeFactors"]]))

  # If filtering, do so after size factor estimation so the size factors are
  # calculated from the original matrix.
  if (normalize_prefilter && !is.null(exclude_features)) {
    mat <- mat[!(rownames(mat) %in% exclude_features), , drop = FALSE]

    dds_filtered <- DESeq2::DESeqDataSetFromMatrix(countData = mat,
                                                   colData = targets,
                                                   design = design)

    DESeq2::sizeFactors(dds_filtered) <- DESeq2::sizeFactors(dds)
    dds <- dds_filtered
  }

  # Run DESeq
  if (de_test == "LRT") {
    # Reduce formula
    term_labels <- attr(stats::terms(design_formula), "term.labels")
    if (length(term_labels) <= 1) {
      reduced_formula <- stats::as.formula("~ 1")
    } else {
      reduced_terms <- term_labels[-length(term_labels)]
      reduced_formula <- stats::as.formula(paste("~", paste(reduced_terms, collapse = " + ")))
    }
    reduced_design <- stats::model.matrix(reduced_formula, data = targets)

    dds <- do.call(DESeq2::DESeq, c(list(object = dds,
                                         test = "LRT",
                                         reduced = reduced_design), de_params[["DESeq"]]))
  } else if (de_test == "Wald") {
    dds <- do.call(DESeq2::DESeq, c(list(object = dds,
                                         test = "Wald"),
                                    de_params[["DESeq"]]))
  }

  raw_results <- DESeq2::results(dds) |>
    data.frame()

  DESeq2_results <- raw_results |>
    dplyr::transmute(feature = rownames(raw_results),
                     lfc = log2FoldChange,
                     pvalue = pvalue)

  rownames(DESeq2_results) <- NULL

  return(list("results" = DESeq2_results,
              "raw_results" = raw_results))
}

# Run limma differential expression ---------------------------
#
# mat                 -- A feature x replicate pseudobulk matrix or a feature x cell matrix
# targets             -- A dataframe containing sample to group key
# design              -- A model.matrix design object
# de_test             -- Which test to use for differential expression
# de_params           -- Additional parameters to pass
# normalize_prefilter -- Whether to get normalization/size factors before filtering out features
# exclude_features    -- A vector of feature names to filter out if normalize_prefilter is TRUE

.runDE.limma <- function(mat,
                         targets,
                         design,
                         de_test = "voom",
                         de_params = list(),
                         normalize_prefilter = FALSE,
                         exclude_features = NULL) {

  .requirePackage("limma", source = "bioc")

  # Create edgeR object
  dge <- do.call(edgeR::DGEList, c(list("counts" = mat,
                                        "group" = targets$group),
                                   de_params[["DGEList"]]))

  if (de_test %in% c("trend", "voom")) {
    # Get size factors
    dge <- do.call(edgeR::calcNormFactors, c(list("object" = dge),
                                             de_params[["calcNormFactors"]]))
    # If filtering, do so
    if (normalize_prefilter & !is.null(exclude_features)) {
      mat <- mat[!(rownames(mat) %in% exclude_features),]
      dge_filtered <- do.call(edgeR::DGEList, c(list("counts" = mat,
                                                     "group" = targets$group),
                                                de_params[["DGEList"]]))
      dge_filtered$samples <- dge$samples
      dge <- dge_filtered
    }

    # Transform
    if (de_test == "trend") {
      # Apply logCPM (uses stored $samples$lib.size)
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

    raw_results <- limma::topTable(fit, coef = ncol(design), number = Inf) |>
      data.frame()
    limma_results <- raw_results |>
      dplyr::transmute(feature = rownames(raw_results),
                       lfc = logFC,
                       pvalue = P.Value)
    rownames(limma_results) <- NULL
  } else if (de_test %in% c("wilcox_cpm", "wilcox_log_cpm")) {
    # Get library sizes
    lib_sizes <- dge$samples[colnames(mat), "lib.size"]
    # If filtering, do so
    if (normalize_prefilter & !is.null(exclude_features)) {
      mat <- mat[!(rownames(mat) %in% exclude_features),]
    }
    # Normalization
    if (de_test == "wilcox_cpm") {
      cpm_mat <- do.call(edgeR::cpm, c(list("y" = mat,
                                            "lib.size" = lib_sizes),
                                       de_params[["cpm"]]))
    } else if (de_test == "wilcox_log_cpm") {
      cpm_mat <- do.call(edgeR::cpm, c(list("y" = mat,
                                            "lib.size" = lib_sizes,
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
    lfcs <- log2(rowMeans(cpm_mat[, group1_indices, drop = FALSE] + pseudocount)/rowMeans(cpm_mat[, group2_indices, drop = FALSE] + pseudocount))

    # Results
    limma_results <- data.frame(feature = rownames(cpm_mat),
                                lfc = lfcs,
                                pvalue = pvalues)
    rownames(limma_results) <- NULL

    # Raw results are the same here
    raw_results <- limma_results
  }

  return(list("results" = limma_results,
              "raw_results" = raw_results))
}

# Run Wilcoxon rank sum test differential expression using presto ---------------------------
#
# mat                 -- A feature x replicate pseudobulk matrix or a feature x cell matrix
# targets             -- A dataframe containing sample to group key
# de_test             -- Which test to use for differential expression
# de_params           -- Additional parameters to pass to cpm and/or wilcoxauc
# normalize_prefilter -- Whether to get normalization/size factors before filtering out features
# exclude_features    -- A vector of feature names to filter out if normalize_prefilter is TRUE
# non_reference_group -- String indicating the non-reference group
.runDE.presto <- function(mat,
                          targets,
                          de_test = "wilcox_cpm",
                          de_params = list(),
                          normalize_prefilter = FALSE,
                          exclude_features = NULL,
                          non_reference_group) {

  .requirePackage("presto", installInfo = 'devtools::install_github("immunogenomics/presto")')

  mat <- methods::as(mat, "dgCMatrix")

  # Get library sizes
  dge <- do.call(edgeR::DGEList, c(list("counts" = mat,
                                        "group" = targets$group),
                                   de_params[["DGEList"]]))
  lib_sizes <- dge$samples[colnames(mat), "lib.size"]
  # If filtering, do so
  if (normalize_prefilter & !is.null(exclude_features)) {
    mat <- mat[!(rownames(mat) %in% exclude_features),]
  }
  # Normalization
  if (de_test == "wilcox_cpm") {
    cpm_mat <- do.call(edgeR::cpm, c(list("y" = mat,
                                          "lib.size" = lib_sizes),
                                     de_params[["cpm"]]))
  } else if (de_test == "wilcox_log_cpm") {
    cpm_mat <- do.call(edgeR::cpm, c(list("y" = mat,
                                          "lib.size" = lib_sizes,
                                          "log" = TRUE),
                                     de_params[["cpm"]]))
  }
  # Run Wilcoxon rank sum test
  raw_results <- do.call(presto::wilcoxauc, c(list("X" = cpm_mat,
                                                   "y" = targets$group),
                                              de_params[["wilcoxauc"]]))
  presto_results <- raw_results |> dplyr::filter(group == non_reference_group)

  presto_results <- presto_results |>
    dplyr::transmute(feature,
                     lfc = logFC,
                     pvalue = pval)
  rownames(presto_results) <- NULL

  return(list("results" = presto_results,
              "raw_results" = raw_results))
}

# Run Wilcoxon rank sum test differential expression using BPCells ---------------------------
#
# mat                 -- A feature x replicate pseudobulk matrix or a feature x cell matrix (must be IterableMatrix)
# targets             -- A dataframe containing sample to group key
# de_test             -- Which test to use for differential expression
# de_params           -- Additional parameters to pass to cpm and/or wilcoxauc
# normalize_prefilter -- Whether to get normalization/size factors before filtering out features
# exclude_features    -- A vector of feature names to filter out if normalize_prefilter is TRUE
# non_reference_group -- String indicating the non-reference group
.runDE.BPCells <- function(mat,
                           targets,
                           de_test = "wilcox_cpm",
                           de_params = list(),
                           normalize_prefilter = FALSE,
                           exclude_features = NULL,
                           non_reference_group) {
  .requirePackage("BPCells", installInfo = "devtools::install_github('bnprks/BPCells/r')")
  if (!methods::is(mat, "IterableMatrix")) {
    stop("Input matrix to '.runDE.BPCells()' must be class 'IterableMatrix'.")
  }
  # Library sizes are column sums of the raw count matrix.
  lib_sizes <- Matrix::colSums(mat)
  if (any(lib_sizes == 0)) {
    stop("Cannot calculate CPM: at least one cell/sample has zero total counts.")
  }
  # If requested, filter features after calculating library sizes.
  if (normalize_prefilter && !is.null(exclude_features)) {
    mat <- mat[!(rownames(mat) %in% exclude_features), , drop = FALSE]
  }

  # CPM
  cpm_mat <- BPCells::multiply_cols(mat,
                                    1e6 / lib_sizes[colnames(mat)])

  # Log transformation
  if (de_test == "wilcox_log_cpm") {
    cpm_mat <- log1p(cpm_mat)
  }

  # Run Wilcoxon rank sum test using BPCells.
  raw_results <- do.call(BPCells::marker_features,
                             c(list(mat = cpm_mat,
                                    groups = targets$group),
                               de_params[["marker_features"]]))

  BPCells_results <- raw_results |>
    dplyr::filter(foreground == non_reference_group)

  # LFC
  if ("pseudocount" %in% de_params[["lfc"]]) {
    pseudocount <- de_params$lfc$pseudocount
  } else {
    pseudocount <- 1
  }

  if (de_test == "wilcox_log_cpm") {
    # cpm_mat is log1p-transformed using the natural log.
    # Difference of means is in natural-log units, so divide by log(2)
    # to report log2-scale fold change (see marker_features docs)
    BPCells_results <- BPCells_results |>
      dplyr::mutate(lfc = (foreground_mean - background_mean) / log(2))
  } else {
    # cpm_mat is CPM-scale, so calculate log2 fold-change from group means
    # (see marker_features docs)
    BPCells_results <- BPCells_results |>
      dplyr::mutate(lfc = log2((foreground_mean + pseudocount) /(background_mean + pseudocount)))
  }

  # Results
  BPCells_results <- BPCells_results |>
    dplyr::transmute(feature,
                     lfc,
                     pvalue = p_val_raw)

  rownames(BPCells_results) <- NULL

  return(list("results" = BPCells_results,
              "raw_results" = raw_results))
}
