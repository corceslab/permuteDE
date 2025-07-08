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
    # Only allowed to be NULL for function countCombinations() or getCombinations()
    # Otherwise, must be either of type Seurat or SingleCellExperiment
    if (any(other != "countCombinations", !is.null(input)) &
        length(intersect(methods::is(input), c("Seurat", "SingleCellExperiment"))) < 1) {
      stop("Input value for '", name, "' is not one of classes Seurat or SingleCellExperiment. Please supply valid input!")
    }
  }

  # replicate_labels, group_labels, split_labels, message
  if (name %in% c("replicate_labels", "group_labels", "split_labels", "message")) {
    # If not NULL
    if (!is.null(input)) {
      # Should be of class "character"
      if (!methods::is(input, "character") | length(input) != 1) {
        stop("Input value for '", name, "' must be a single value of class 'character', please supply valid input!")
      }
      # Must be present in metadata of provided object
      if (methods::is(other, "Seurat")) {
        if (!(input %in% colnames(other@meta.data))) {
          stop("Input value for '", name, "' must indicate a column present in the 'meta.data' of the provided object, please supply valid input!")
        }
      } else if (methods::is(other, "SingleCellExperiment")) {
        if (!(input %in% colnames(other@colData))) {
          stop("Input value for '", name, "' must indicate a column present in the 'colData' of the provided object, please supply valid input!")
        }
      }
    }
  }

  # use_cells
  if (name == "use_cells") {
    # If not NULL & object is not NULL
    if (!is.null(input) & (methods::is(other, "Seurat") | methods::is(other, "SingleCellExperiment"))) {
      cell_ids <- colnames(other)
      if (length(intersect(input, cell_ids)) != length(input)) {
        stop("Not all provided cells are present in the provided object, please supply valid input!")
      }
    }
  }

  # Single positive integer
  # n_replicates, n_group1, n_cores, min_cells_per_split, min_replicates_per_split, min_replicates_per_group, random_seed
  if (name %in% c("min_cells_per_split", "min_replicates_per_split", "min_replicates_per_group",
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
    }
  }

  # Single non-negative integer
  # min_cells_per_feature, min_DE
  if (name %in% c("min_cells_per_feature", "min_DE")) {
    # Should be of class 'numeric', must be a single value
    if (!methods::is(input, "numeric") | length(input) != 1) {
      stop("Input value for '", name, "' must be a single value of class 'numeric'. Please supply valid input!")
    }
    # Must be positive integer
    if (input %% 1 != 0 | input < 0) {
      stop("Input value for '", name, "' must be a non-negative integer. Please supply valid input!")
    }
  }

  # Single logical value
  # verbose, return_all
  if (name %in% c("verbose")) {
    # Must be T/F
    if (!methods::is(input, "logical") | length(input) != 1) {
      stop("Input value for '", name, "' is not a single value of class 'logical', please supply valid input!")
    }
  }

  # use_assay
  if (name == "use_assay") {
    # If not NULL
    if (!is.null(input)) {
      # Only relevant for Seurat and SingleCellExperiment objects
      if (is.null(other)) {
        warning("Input value for '", name, "' is not used when 'object' is NULL.")
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
    if (!(input %in% c("edgeR", "DESeq2", "limma"))) {
      stop("Input for '", name, "' must be among permitted values (", paste0(c("edgeR", "DESeq2", "limma"), collapse = ", "), "), please supply valid input!")
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
    if (!identical(names(input), c("DE_results", "PB_values",  "group_key", "parameters"))) {
      stop("Structure of list provided for parameter 'input' is unexpected. It should be a list with four named elements ",
           "('DE_results', 'PB_values', 'group_key', and 'parameters'). Please supply valid input!")
    }
  }

  # Single non-negative number
  # alpha
  if (name %in% c("min_cells_per_feature")) {
    if (!methods::is(input, "numeric") | length(input) != 1) {
      stop("Input value for '", name, "' must be a single value of class 'numeric'. Please supply valid input!")
    }
    if (input < 0) {
      stop("Input value for '", name, "' cannot be negative. Please supply valid input!")
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
}
