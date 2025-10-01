# ---------------------------------------------------------------------------
# Validation methods
# ---------------------------------------------------------------------------

# Check parameter input validity ---------------------------
#
# Checks validity of provided values for input parameters
#
# input -- input value
# name -- name of parameter
# other -- other inputs
.validInput <- function(input = NULL,
                        name = NULL,
                        other = NULL) {

  # object
  if (name == "object") {
    # Only allowed to be NULL for functions countCombinations and getCombinations
    if (any(!(other %in% c("countCombinations", "getCombinations")), !is.null(input))) {
      # Otherwise, must be of type Seurat, SingleCellExperiment, or matrix
      if (length(intersect(methods::is(input), c("Seurat", "SingleCellExperiment", "matrix"))) < 1) {
        stop("Input value for '", name, "' is not one of classes 'Seurat', 'SingleCellExperiment', or 'matrix'. Please supply valid input!")
      }
      # If object is of type matrix, must have row names
      if ("matrix" %in% methods::is(input) & is.null(rownames(input))) {
        stop("When input value for '", name, "' is of class 'matrix', row names cannot be NULL. Please set row names to feature names.")
      }
    }
  }

  # replicate_labels, group_labels, split_labels
  if (name %in% c("replicate_labels", "group_labels", "split_labels")) {
    # If not NULL
    if (!is.null(input)) {
      # Should be of class "character"
      if (!methods::is(input, "character")) {
        stop("Input value for '", name, "' must be of class 'character', please supply valid input!")
      }
      # replicate_labels are not used for cell-level tests
      if (name == "replicate_labels") {
        if (other[[2]] == "none") {
          stop("Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'none'. Please supply valid input!")
        }
      }
      # If single value, must be a column name
      if (length(input) == 1) {
        # Must be present in metadata of provided object
        if (methods::is(other[[1]], "Seurat")) {
          if (!(input %in% colnames(other[[1]]@meta.data))) {
            stop("When a single input value is provided for '", name, "', it must indicate a column present in the 'meta.data' of the provided object, please supply valid input!")
          }
        } else if (methods::is(other[[1]], "SingleCellExperiment")) {
          if (!(input %in% colnames(other[[1]]@colData))) {
            stop("When a single input value is provided for '", name, "', it must indicate a column present in the 'colData' of the provided object, please supply valid input!")
          }
        } else if (methods::is(other[[1]], "matrix")) {
          stop("When input to parameter 'object' is of class 'matrix', input value for '", name, "' cannot be a single value, it must be a vector.")
        }
      }
    } else {
      # replicate_labels are required for pseudobulk tests
      if (name == "replicate_labels") {
        if (other[[2]] %in% c("generate", "supplied")) {
          stop("Input value for '", name, "' cannot be NULL when parameter 'pseudobulk' is set to '", other[[2]], "'. Please supply valid input!")
        } else if (other[[2]] == "none") {
          warning("If the output of this function is intended to be passed to function 'permuteDE', input value for '", name, "' cannot be NULL when parameter 'pseudobulk' is set to '", other[[2]], "'. Please supply valid input!")
        }
      }
    }
  }

  # use_cells
  if (name == "use_cells") {
    # If not NULL & object is not NULL
    if (!is.null(input) & !is.null(other[[1]])) {
      if (other[[2]] == "supplied") {
        warning("Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'supplied' (when a pre-computed pseudobulk matrix is supplied by the user).")
      } else {
        cell_ids <- colnames(other[[1]])
        if (length(intersect(input, cell_ids)) != length(input)) {
          stop("Not all provided cells are present in the provided object, please supply valid input!")
        }
      }
    }
  }

  # Single positive integer
  # n_replicates, n_group1, n_cores, min_cells_per_split, min_replicates_per_split, min_replicates_per_group, min_cells_per_feature, random_seed
  if (name %in% c("min_cells_per_split", "min_replicates_per_split", "min_replicates_per_group", "min_cells_per_feature",
                  "n_replicates", "n_group1",
                  "n_combinations", "n_iterations",
                  "random_seed", "n_cores")) {
    # n_cores, n_replicates, n_group1, can be NULL
    if (!(name %in% c("n_replicates", "n_group1", "n_cores") & is.null(input))) {
      # Should be of class 'numeric', must be a single value
      if (!methods::is(input, "numeric") | length(input) != 1) {
        stop("Input value for '", name, "' must be a single value of class 'numeric'. Please supply valid input!")
      }
      # Must be positive integer
      if (input %% 1 != 0 | input < 1) {
        stop("Input value for '", name, "' must be a positive integer. Please supply valid input!")
      }
      # n_group1 must be < n_replicates
      if (name == "n_group1") {
        if (input >= other) {
          stop("Input value for '", name, "' must be less than input value for 'n_replicates'. Please supply valid input!")
        }
      }
      # min_cells_per_split, min_cells_per_feature are not applicable when pre-computed pseudobulk matrix is supplied by the user
      if (name %in% c("min_cells_per_split", "min_cells_per_feature")) {
        if (other == "supplied") {
          warning("Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'supplied' (when a pre-computed pseudobulk matrix is supplied by the user).")
        }
      }
      # min_replicates_per_split is not applicable when doing cell-level tests
      if (name %in% c("min_replicates_per_split")) {
        if (other == "none") {
          warning("Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'none'.")
        }
      }
    }
  }

  # Single non-negative integer
  # min_DE
  if (name %in% c("min_DE")) {
    # Should be of class 'numeric', must be a single value
    if (!methods::is(input, "numeric") | length(input) != 1) {
      stop("Input value for '", name, "' must be a single value of class 'numeric'. Please supply valid input!")
    }
    # Must be positive integer
    if (input %% 1 != 0 | input < 0) {
      stop("Input value for '", name, "' must be a non-negative integer. Please supply valid input!")
    }
  }

  # Single Boolean value
  # force_balance, return_all, verbose
  if (name %in% c("force_balance", "return_all", "verbose")) {
    # Must be T/F
    if (!methods::is(input, "logical") | length(input) != 1) {
      stop("Input value for '", name, "' is not a single value of class 'logical', please supply valid input!")
    }
    # force_balance is not applicable for cell-level tests
    if (name %in% c("force_balance")) {
      if (other == "none" & input == TRUE) {
        stop("Parameter '", name, "' cannot be set to '", input, "' when parameter 'pseudobulk' is set to 'none'. Please supply valid input!")
      }
    }
  }

  # use_assay
  if (name == "use_assay") {
    # If not NULL
    if (!is.null(input)) {
      # Only relevant for Seurat and SingleCellExperiment objects
      if (is.null(other)) {
        warning("Input value for '", name, "' is not used when 'object' is NULL.")
      } else if (methods::is(other, "matrix")) {
        warning("Input value for '", name, "' is not used when 'object' is of class 'matrix'.")
      } else if (methods::is(other, "Seurat") | methods::is(other, "SingleCellExperiment")) {
        # Should be of class 'character'
        if (!methods::is(input, "character") | length(input) != 1) {
          stop("Input for '", name, "' must be a single value of class 'character', please supply valid input!")
        }
        # Must exist in provided object
        if (!(input %in% names(other@assays))) {
          stop("Assay '", input, "' provided for parameter '", name, "' is not present in provided object, please supply valid input!")
        }
      }
    }
  }

  # use_layer
  if (name == "use_layer") {
    # If not NULL
    if (!is.null(input)) {
      # Only relevant for Seurat objects
      if (methods::is(other[[1]], "Seurat")) {
        # Should be of class 'character'
        if (!methods::is(input, "character") | length(input) != 1) {
          stop("Input for '", name, "' must be a single value of class 'character', please supply valid input!")
        }
        # Check correspondence with input for 'use_assay'
        if (is.null(other[[2]])) {
          other[[2]] <- Seurat::DefaultAssay(other[[1]])
        }
        # Layer must be present under specified assay in provided object
        if (length(other[[2]]) == 1) {
          # Check if object is Seurat v5
          if ("Assay5" %in% methods::is(other[[1]][[other[[2]]]])) {
            if (!(input %in% names(other[[1]][[other[[2]]]]@layers))) {
              stop("Layer '", input, "' is not present in assay '", other[[2]], "' of provided Seurat v5 object, please supply valid input!")
            }
          } else {
            try(slot_exists <- methods::validObject(methods::slot(other[[1]][[other[[2]]]], input)))
            if (!exists("slot_exists")) {
              stop("Slot '", input, "' is not present in assay '", other[[2]], "' of provided Seurat object, please supply valid input!")
            } else if (slot_exists == FALSE) {
              stop("Slot '", input, "' is not present in assay '", other[[2]], "' of provided Seurat object, please supply valid input!")
            }
          }
        }
      } else {
        warning("Input value(s) for '", name, "' are not used when provided object is not of class 'Seurat'.")
      }
    }
  }

  # de_method
  if (name == "de_method") {
    # Should be of class 'character'
    if (!methods::is(input, "character") | length(input) != 1) {
      stop("Input for '", name, "' must be a single value of class 'character', please supply valid input!")
    }
    # Must be among permitted values
    if (!(input %in% c("edgeR", "DESeq2", "limma", "presto"))) {
      stop("Input for '", name, "' must be among permitted values (", paste0(c("edgeR", "DESeq2", "limma", "presto"), collapse = ", "), "), please supply valid input!")
    }
  }

  # de_test
  if (name == "de_test") {
    # Should be of class 'character'
    if (!methods::is(input, "character") | length(input) != 1) {
      stop("Input for '", name, "' must be a single value of class 'character', please supply valid input!")
    }
    # Must be among permitted values
    if (other == "edgeR") {
      if (!(input %in% c("LRT", "QLF"))) {
        stop("When input for 'de_method' is '", other, "', input for '", name, "' must be among permitted values (", paste0(c("LRT", "QLF"), collapse = ", "), "), please supply valid input!")
      }
    } else if (other == "DESeq2") {
      if (!(input %in% c("LRT", "Wald"))) {
        stop("When input for 'de_method' is '", other, "', input for '", name, "' must be among permitted values (", paste0(c("LRT", "Wald"), collapse = ", "), "), please supply valid input!")
      }
    } else if (other == "limma") {
      if (!(input %in% c("trend", "voom"))) {
        stop("When input for 'de_method' is '", other, "', input for '", name, "' must be among permitted values (", paste0(c("trend", "voom"), collapse = ", "), "), please supply valid input!")
      }
    } else if (other == "presto") {
      if (!(input %in% c("wilcox_cpm", "wilcox_log_cpm"))) {
        stop("When input for 'de_method' is '", other, "', input for '", name, "' must be among permitted values (", paste0(c("wilcox_cpm", "wilcox_log_cpm"), collapse = ", "), "), please supply valid input!")
      }
    }
  }

  # p_adjust_method
  if (name == "p_adjust_method") {
    # Should be of class 'character'
    if (!methods::is(input, "character") | length(input) != 1) {
      stop("Input for '", name, "' must be a single value of class 'character', please supply valid input!")
    }
    # Must be among permitted values
    if (!(input %in% stats::p.adjust.methods)) {
      stop("Input for '", name, "' must be among permitted values (", paste0(stats::p.adjust.methods, collapse = ", "), "), please supply valid input!")
    }
  }

  # input
  if (name == "input") {
    # Must be of type list
    if (!methods::is(input, "list")) {
      stop("Parameter '", name, "' must be a list, please supply valid input!")
    }
    # Must have expected elements with set names
    if (!identical(names(input), c("DE_results", "PB_values",  "group_key", "metadata", "parameters")) &
        !identical(names(input), c("DE_results", "cell_values",  "group_key", "metadata", "parameters"))) {
      stop("Structure of list provided for parameter 'input' is unexpected. It should be a list with five named elements ",
           "('DE_results', 'PB_values' (or 'cell_values'), 'group_key', 'metadata', and 'parameters'). Please supply valid input!")
    }
  }

  # Single number from 0-1
  # alpha, min_prop_cells_per_feature
  if (name %in% c("alpha", "min_prop_cells_per_feature")) {
    if (!methods::is(input, "numeric") | length(input) != 1) {
      stop("Input value for '", name, "' must be a single value of class 'numeric'. Please supply valid input!")
    }
    if (input < 0) {
      stop("Input value for '", name, "' cannot be negative. Please supply valid input!")
    }
    if (input > 1) {
      stop("Input value for '", name, "' cannot be greater than 1. Please supply valid input!")
    }
    # min_prop_cells_per_feature is not applicable when pre-computed pseudobulk matrix is supplied by the user
    if (name %in% c("min_prop_cells_per_feature")) {
      if (other == "supplied") {
        warning("Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'supplied' (when a pre-computed pseudobulk matrix is supplied by the user).")
      }
    }
  }

  # use_splits
  if (name == "use_splits") {
    # If not NULL
    if (!is.null(input)) {
      # Values must all be among those in input
      if (!all(input %in% names(other$PB_values))) {
        stop("Input value(s) for '", name, "' must all be present among provided pseudobulk matrices. Please supply valid input!")
      }
    }
  }

  # reference_group
  if (name == "reference_group") {
    # If not NULL
    if (!is.null(input)) {
      # Value must be among those indicated by group_labels
      groups <- .retrieveData(object = other[[1]],
                              type = "cell_metadata",
                              name = other[[2]],
                              use_cells = other[[3]])
      if (!input %in% groups) {
        stop("Input value for '", name, "' must be present among provided group labels. Please supply valid input!")
      }
    }
  }

  # pseudobulk
  if (name == "pseudobulk") {
    # Should be of class 'character'
    if (!methods::is(input, "character") | length(input) != 1) {
      stop("Input for '", name, "' must be a single value of class 'character', please supply valid input!")
    }
    # Must be among permitted values
    if (!(input %in% c("generate", "supplied", "none"))) {
      stop("Input for '", name, "' must be among permitted values (", paste0(c("generate", "supplied", "none"), collapse = ", "), "), please supply valid input!")
    }
    # Issue warning for cell-level tests
    if (input == "none") {
      warning("Cell-level tests are not recommended in most cases, proceed with caution.")
    }
    # If supplied, object cannot be Seurat or SingleCellExperiment
    if (input == "supplied" & length(intersect(methods::is(other), c("Seurat", "SingleCellExperiment"))) > 0) {
      stop("When input for '", name, "' is '", input, "', parameter 'object' must be of class 'matrix'. Please supply valid input!")
    }
  }

  # message
  if (name == "message") {
    # Should be of class 'character'
    if (!methods::is(input, "character") | length(input) != 1) {
      stop("Input for '", name, "' must be a single value of class 'character'.",
           " This parameter is intended for internal use, we recommend leaving it as the default value.")
    }
  }

  # lfc_threshold
  if (name == "lfc_threshold") {
    # Single non-negative number
    if (!methods::is(input, "numeric") | length(input) != 1) {
      stop("Input value for '", name, "' must be a single value of class 'numeric'. Please supply valid input!")
    }
    if (input < 0) {
      stop("Input value for '", name, "' cannot be negative. Please supply valid input!")
    }
  }

  # de_params
  if (name == "de_params") {
    # Check that it is a list
    if (!methods::is(input, "list")) {
      stop("Input value for '", name, "' is not of class 'list', please supply valid input!")
    }
    # Check whether it is empty
    if (length(input) > 0) {
      # Allowed functions to pass parameters to depends on selected DE method/test
      if (other[[1]] == "edgeR") {
        if (other[[2]] == "LRT") {
          allowed_functions <- c("DGEList", "calcNormFactors", "estimateDisp", "glmQLFit", "glmQLFTest")
        } else if (other[[2]] == "QLF") {
          allowed_functions <- c("DGEList", "calcNormFactors", "estimateDisp", "glmFit", "glmLRT")
        } else if (other[[2]] == "exact") {
          allowed_functions <- c("DGEList", "calcNormFactors", "estimateDisp", "exactTest")
        }
      } else if (other[[1]] == "DESeq2") {
        allowed_functions <- c("DESeq")
      } else if (other[[1]] == "limma") {
        if (other[[2]] == "trend") {
          allowed_functions <- c("DGEList", "calcNormFactors", "cpm", "lmFit", "eBayes")
        } else if (other[[2]] == "voom") {
          allowed_functions <- c("DGEList", "calcNormFactors", "voom", "lmFit", "eBayes")
        } else if (other[[2]] %in% c("wilcox_cpm", "wilcox_log_cpm")) {
          allowed_functions <- c("cpm", "rankSumTestWithCorrelation", "lfc")
        }
      } else if (other[[1]] == "presto") {
        allowed_functions <- c("cpm", "wilcoxauc")
      }
      if (!all(names(input) %in% allowed_functions)) {
        stop("When supplying additional parameters to '", name, "' for use",
             " with '", de_method, ": ", de_test, "',",
             " please provide a list of lists, where each secondary list is named according to",
             " the allowed functions (", paste0(allowed_functions, collapse = ", "), ").")
      }
    }
  }
}
