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
      if (length(intersect(methods::is(input), c("Seurat", "SingleCellExperiment", "matrix", "Matrix", "dgCMatrix"))) < 1) {
        stop("Input value for '", name, "' is not one of classes 'Seurat', 'SingleCellExperiment', or 'matrix'. Please supply valid input!")
      }
      # If object is of type matrix, must have row names
      if ((length(intersect(methods::is(input), c("matrix", "Matrix", "dgCMatrix"))) >= 1) & is.null(rownames(input))) {
        stop("When input value for '", name, "' is of class 'matrix', row names cannot be NULL. Please set row names to feature names.")
      }
    }
  }

  # replicate_labels, group_labels, split_labels
  if (name %in% c("replicate_labels", "group_labels", "split_labels")) {
    # If not NULL
    if (!is.null(input)) {
      # Should be of class "character" or "factor"
      if (!methods::is(input, "character") & !methods::is(input, "factor")) {
        if (methods::is(input, "numeric") | methods::is(input, "logical")) {
          warning(" Input value for '", name, "' will be converted to class 'character'.")
        } else {
          stop("Input value for '", name, "' must be a single value of class 'character' or a vector of labels, please supply valid input!")
        }
      }
      # replicate_labels are not used for cell-level tests
      if (name == "replicate_labels") {
        if (other[[3]] == "none") {
          warning(" Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'none'.")
        }
      }
      # If single value, must be a column name
      if (length(input) == 1) {
        # If metadata is provided, must be found there
        if (!is.null(other[[2]])) {
          if (!(input %in% colnames(other[[2]]))) {
            stop("When a single input value is provided for '", name, "', it must indicate a column present in the provided 'metadata', please supply valid input!")
          }
        } else {
          # Otherwise must be present in metadata of provided object
          if (methods::is(other[[1]], "Seurat")) {
            if (!(input %in% colnames(other[[1]]@meta.data))) {
              stop("When a single input value is provided for '", name, "', it must indicate a column present in the 'meta.data' of the provided object, please supply valid input!")
            }
          } else if (methods::is(other[[1]], "SingleCellExperiment")) {
            if (!(input %in% colnames(other[[1]]@colData))) {
              stop("When a single input value is provided for '", name, "', it must indicate a column present in the 'colData' of the provided object, please supply valid input!")
            }
          } else if (methods::is(other[[1]], "matrix")) {
            stop("When input to parameter 'object' is of class 'matrix' and input to parameter 'metadata' is NULL, input value for '",
                 name, "' cannot be a single value, it must be a vector.")
          }
        }
      }
    } else {
      # replicate_labels are required for pseudobulk tests
      if (name == "replicate_labels") {
        if (other[[3]] %in% c("generate", "supplied")) {
          stop("Input value for '", name, "' cannot be NULL when parameter 'pseudobulk' is set to '", other[[3]], "'. Please supply valid input!")
        }
      }
    }
  }

  # use_cells
  if (name == "use_cells") {
    # If not NULL & object is not NULL
    if (!is.null(input) & !is.null(other[[1]])) {
      if (other[[2]] == "supplied") {
        warning(" Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'supplied' (when a pre-computed pseudobulk matrix is supplied by the user).")
      } else {
        cell_ids <- colnames(other[[1]])
        if (length(intersect(input, cell_ids)) != length(input)) {
          stop("Not all provided cells are present in the provided object, please supply valid input!")
        }
      }
    }
  }

  # Single positive integer
  # n, n_cores, min_cells_per_split, min_cells_per_replicate, min_replicates_per_split,
  # min_replicates_per_group, min_cells_per_feature, random_seed, n_combinations
  if (name %in% c("min_cells_per_split", "min_cells_per_replicate", "min_replicates_per_split", "min_replicates_per_group", "min_cells_per_feature",
                  "n_combinations", "n_iterations",
                  "random_seed", "n_cores", "n")) {
    # n_cores & n can be NULL
    if (!(name %in% c("n_cores", "n") & is.null(input))) {
      # Should be of class 'numeric', must be a single value
      if (!methods::is(input, "numeric") | length(input) != 1) {
        stop("Input value for '", name, "' must be a single value of class 'numeric'. Please supply valid input!")
      }
      # Must be positive integer
      if (input %% 1 != 0 | input < 1) {
        stop("Input value for '", name, "' must be a positive integer. Please supply valid input!")
      }
      # min_cells_per_split, min_cells_per_replicate, min_cells_per_feature
      # are not applicable when pre-computed pseudobulk matrix is supplied by the user
      if (name %in% c("min_cells_per_split", "min_cells_per_replicate", "min_cells_per_feature")) {
        if (other[[1]] == "supplied") {
          warning(" Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'supplied' (when a pre-computed pseudobulk matrix is supplied by the user).")
        }
      }
      # min_cells_per_replicate, min_replicates_per_split are not applicable when doing cell-level tests
      if (name %in% c("min_cells_per_replicate", "min_replicates_per_split")) {
        if (other[[1]] == "none" & other[[2]] != "runDE") {
          warning(" Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'none'.")
        }
      }
    }
  }

  # n_replicates
  if (name %in% "n_replicates") {
    # If not NULL
    if (!is.null(input)) {
      # Should be of class 'numeric'
      if (!methods::is(input, "numeric")) {
        stop("Input value for '", name, "' must be of class 'numeric'. Please supply valid input!")
      }
      # Each value must be positive integer
      for (i in 1:length(input)) {
        if (i %% 1 != 0 | i < 1) {
          stop("Input value(s) for '", name, "' must be positive integer(s). Please supply valid input!")
        }
      }
    }
  }

  # n_group1
  if (name %in% "n_group1") {
    # If not NULL
    if (!is.null(input)) {
      # Should be of class 'numeric'
      if (!methods::is(input, "numeric")) {
        stop("Input value for '", name, "' must be of class 'numeric'. Please supply valid input!")
      }
      # Each value must be positive integer
      for (i in 1:length(input)) {
        if (input[i] %% 1 != 0 | input[i] < 1) {
          stop("Input value(s) for '", name, "' must be positive integer(s). Please supply valid input!")
        }
      }
      # n_replicates cannot be NULL
      if (is.null(other)) {
        warning(" Input for '", name, "' is not used when parameter 'n_replicates' is NULL.")
      } else {
        # Must have same length as n_replicates
        if (length(input) != length(other)) {
          stop("Input for '", name, "' must be a vector of the same length as input to 'n_replicates'. Please supply valid input!")
        }
        # Each value < corresponding value of n_replicates
        for (i in 1:length(input)) {
          if (input[i] >= other[i]) {
            stop("Each input value for '", name,
                 "' must be less than the corresponding input value for 'n_replicates'. Please supply valid input!")
          }
        }
      }
    } else {
      # n_replicates must also be NULL
      if (!is.null(other)) {
        warning(" Input for 'n_replicates' is not used when parameter 'n_group1' is NULL.")
      }
    }
  }

  # Single non-negative integer
  # min_DE, n_max_label
  if (name %in% c("min_DE", "n_max_label")) {
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
  # center, force_balance, label_pvalue, label_replicates, label_splits, label_statistics,
  # return_all, verbose, swatch, fix_coords, normalize_prefilter, filter
  if (name %in% c("center", "force_balance", "label_pvalue", "label_replicates", "label_splits", "label_statistics",
                  "return_all", "verbose", "swatch", "fix_coords", "normalize_prefilter", "filter")) {
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
        warning(" Input value for '", name, "' is not used when 'object' is NULL.")
      } else if (methods::is(other, "matrix")) {
        warning(" Input value for '", name, "' is not used when 'object' is of class 'matrix'.")
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
        warning(" Input value(s) for '", name, "' are not used when provided object is not of class 'Seurat'.")
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
      if (!(input %in% c("trend", "voom", "wilcox_cpm", "wilcox_log_cpm"))) {
        stop("When input for 'de_method' is '", other, "', input for '", name, "' must be among permitted values (", paste0(c("trend", "voom", "wilcox_cpm", "wilcox_log_cpm"), collapse = ", "), "), please supply valid input!")
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
    # If not NULL
    if (!is.null(input)) {
      # Must be of type list
      if (!methods::is(input, "list")) {
        stop("Parameter '", name, "' must be a list, please supply valid input!")
      }
      if (other == "permuteDE") {
        # Must have expected elements with set names
        if (!identical(names(input), c("DE_results", "PB_values",  "metadata", "parameters")) &
            !identical(names(input), c("DE_results", "cell_values",  "metadata", "parameters"))) {
          stop("Structure of list provided for parameter 'input' is unexpected. ",
               "It should be the output returned by function 'runDE()' or a list following the same structure with four named elements ",
               "('DE_results', 'PB_values' (or 'cell_values'), 'metadata', and 'parameters'). Please supply valid input!")
        }
        # Metadata must contain group key
        if (!("group_key" %in% names(input$metadata))) {
          stop("Structure of list provided for parameter 'input' is unexpected. ",
               "It should be the output returned by function 'runDE()' or a list following the same structure. ",
               "Element 'metadata' must contain a dataframe under name 'group_key'. Please supply valid input!")
        }
        # Parameters must contain certain values
        if (!all(c("reference_group", "non_reference_group", "design_formula", "de_method", "de_test", "de_params",
                   "p_adjust_method", "pseudobulk") %in% names(input$parameters))) {
          stop("Structure of list provided for parameter 'input' is unexpected. ",
               "It should be the output returned by function 'runDE()' or a list following the same structure. ",
               "Element 'parameters' must be a list containing a minimal set of named elements ",
               "('reference_group', 'non_reference_group', 'design_formula', 'de_method', 'de_test', 'de_params', ",
               "'p_adjust_method', 'pseudobulk'). Please supply valid input!")
        }
        # Is the content consistent?
        if (length(input$parameters$reference_group) != 1 | length(input$parameters$non_reference_group) != 1) {
          stop("Content of list provided for parameter 'input' is unexpected. ",
               "It should be the output returned by function 'runDE()' or a list following the same structure. ",
               "The reference/non-reference groups provided under element 'parameters' must be single values. Please supply valid input!")
        } else {
          if (!(input$parameters$reference_group %in% input$metadata$group_key[, "group"])) {
            stop("Content of list provided for parameter 'input' is unexpected. ",
                 "It should be the output returned by function 'runDE()' or a list following the same structure. ",
                 "The reference group provided under element 'parameters' was not found within the 'group_key' provided under element 'metadata'. Please supply valid input!")
          }
          if (!(input$parameters$non_reference_group %in% input$metadata$group_key[, "group"])) {
            stop("Content of list provided for parameter 'input' is unexpected. ",
                 "It should be the output returned by function 'runDE()' or a list following the same structure. ",
                 "The non-reference group provided under element 'parameters' was not found within the 'group_key' provided under element 'metadata'. Please supply valid input!")
          }
        }
        if (!is.null(input$parameters$design_formula)) {
          # Must be formula
          if(!methods::is(input$parameters$design_formula, "formula")) {
            stop("Content of list provided for parameter 'input' is unexpected. ",
                 "It should be the output returned by function 'runDE()' or a list following the same structure. ",
                 "The 'design_formula' provided under element 'parameters' must be of class 'formula'. Please supply valid input!")
          }
          # Terms must be within column names of group key
          terms <- attr(terms(input$parameters$design_formula), "term.labels")
          if (any(c(grepl("group:", terms), grepl(":group", terms)))) {
            stop("Content of list provided for parameter 'input' is unexpected. ",
                 "The terms within the 'design_formula' provided under element 'parameters' include interaction term(s) that involve ",
                 "the primary comparison groups to be permuted. The 'permuteDE()' function is not compatible with these interaction terms. ",
                 "Please supply valid input!")
          }
          terms <- unique(unlist(strsplit(terms, ":", fixed = TRUE)))
          if (!all(terms %in% colnames(input$metadata$group_key))) {
            stop("Content of list provided for parameter 'input' is unexpected. ",
                 "It should be the output returned by function 'runDE()' or a list following the same structure. ",
                 "The terms within the 'design_formula' provided under element 'parameters' must correspond to ",
                 "the column names of dataframe 'group_key' provided under element 'metadata'. Please supply valid input!")
          }
          # For each term, check for NA values
          for (t in terms) {
            if (any(is.na(input$metadata$group_key[, t]))) {
              stop("Content of list provided for parameter 'input' is unexpected. ",
                   "It should be the output returned by function 'runDE()' or a list following the same structure. ",
                   "Dataframe 'group_key' provided under element 'metadata' cannot contain NA values. Please supply valid input!")
            }
          }
        }
      } else if (other == "plotVolcano") {
        # Must have expected elements with set names
        if (!("DE_results" %in% names(input))) {
          stop("Structure of list provided for parameter 'input' is unexpected, ",
               "it should be the output returned by function 'runDE()' or a list containing (at minimum) a dataframe named 'DE_results'. ",
               "Please supply valid input!")
        }
      } else if (other == "plotHistogram") {
        # Must have expected elements with set names
        if (!("permutation_test_results" %in% names(input)) | !("permutation_DE_summary" %in% names(input))) {
          stop("Structure of list provided for parameter 'input' is unexpected, it should be the output returned by function 'permuteDE()' ",
               "or a list containing (at minimum) dataframes named 'permutation_test_results' and 'permutation_DE_summary'. ",
               "Please supply valid input!")
        }
      } else if (other == "plotFeature") {
        # Must have expected elements with set names
        if (!("metadata" %in% names(input)) | !("parameters" %in% names(input))) {
          stop("Structure of list provided for parameter 'input' is unexpected, it should be the output returned by function 'runDE()' ",
               "or a list containing (at minimum) elements named 'PB_values' (or 'cell_values'), 'metadata', and 'parameters'. Please supply valid input!")
        }
        # Metadata must contain group key
        if (!("group_key" %in% names(input$metadata))) {
          stop("Structure of list provided for parameter 'input' is unexpected. ",
               "It should be the output returned by function 'runDE()' or a list following the same structure. ",
               "Element 'metadata' must contain a dataframe under name 'group_key'. Please supply valid input!")
        }
      } else if (other == "plotDimReduction") {
        # Must have expected elements with set names
        if (!("permutation_test_results" %in% names(input))) {
          stop("Structure of list provided for parameter 'input' is unexpected, it should be the output returned by function 'permuteDE()' ",
               "or a list containing (at minimum) an element named 'permutation_test_results'. Please supply valid input!")
        }
      }
    } else {
      if (other != "plotDimReduction") {
        stop("Parameter '", name, "' cannot be NULL, please supply valid input!")
      }
    }
  }

  # Single number from 0-1
  # alpha, min_prop_cells_per_feature, permutation_test_alpha
  if (name %in% c("alpha", "min_prop_cells_per_feature", "permutation_test_alpha")) {
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
        warning(" Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'supplied' (when a pre-computed pseudobulk matrix is supplied by the user).")
      }
    } else if ((name == "permutation_test_alpha") & (input < 1)) {
      if (!(other %in% c("pvalue", "split"))) {
        warning(" Input value for '", name, "' is not used when parameter 'color_by' is set to '", other, "'.")
      }
    }
  }

  # use_splits
  if (name == "use_splits") {
    # If not NULL
    if (!is.null(input)) {
      # Values must all be among those in input
      if (other[[2]] %in% c("permuteDE", "plotFeature")) {
        if ("PB_values" %in% names(other[[1]])) {
          if (!all(input %in% names(other[[1]]$PB_values))) {
            stop("Input value(s) for '", name, "' must all be present among provided pseudobulk matrices. Please supply valid input!")
          }
        } else if ("cell_values" %in% names(other[[1]])) {
          if (!all(input %in% names(other[[1]]$cell_values))) {
            stop("Input value(s) for '", name, "' must all be present among provided matrices. Please supply valid input!")
          }
        }
      } else if (other[[2]] %in% c("plotVolcano")) {
        if (!all(input %in% other[[1]]$DE_results$split)) {
          stop("Input value(s) for '", name, "' must all be present in DE results. Please supply valid input!")
        }
      } else if (other[[2]] %in% c("plotHistogram")) {
        if (!all(input %in% other[[1]]$permutation_DE_summary$split)) {
          stop("Input value(s) for '", name, "' must all be present in permutation test summary. Please supply valid input!")
        }
      }
    }
  }

  # reference_group
  if (name == "reference_group") {
    # If not NULL
    if (!is.null(input)) {
      if (length(input) != 1) {
        stop("Input for '", name, "' must be a single value, please supply valid input!")
      } else {
        if (length(other[[3]]) == 1) {
          # Value must be among those indicated by group_labels
          groups <- .retrieveData(object = other[[1]],
                                  metadata = other[[2]],
                                  type = "cell_metadata",
                                  name = other[[3]],
                                  use_cells = other[[4]])
        } else {
          groups <- other[[3]]
        }
        if (!input %in% groups) {
          stop("Input value for '", name, "' must be present among provided group labels. Please supply valid input!")
        }
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
    # Issue warning for cell-level tests (when function is runDE, will warn during internal getPseudobulk call)
    if (other[[1]] != "runDE" & input == "none") {
      warning(" Cell-level tests are not recommended in most cases, proceed with caution.")
    }
    if (other[[1]] == "permuteDE") {
      if (input == "none") {
        if (!("cell_values" %in% names(other[[2]]))) {
          stop("Structure of list provided for parameter 'input' is unexpected. When conducting cell-level tests, ",
               "it should be the output returned by function 'runDE()' or a list following the same structure with four named elements ",
               "('DE_results', 'cell_values', 'metadata', and 'parameters'). Please supply valid input!")
        }
      } else {
        if (!("PB_values" %in% names(other[[2]]))) {
          stop("Structure of list provided for parameter 'input' is unexpected. When conducting pseudobulk tests, ",
               "it should be the output returned by function 'runDE()' or a list following the same structure with four named elements ",
               "('DE_results', 'PB_values', 'metadata', and 'parameters'). Please supply valid input!")
        }
      }
    } else {
      # If supplied, object cannot be Seurat or SingleCellExperiment
      if (input == "supplied" & length(intersect(methods::is(other[[2]]), c("Seurat", "SingleCellExperiment"))) > 0) {
        stop("When input for '", name, "' is '", input, "', parameter 'object' must be of class 'matrix'. Please supply valid input!")
      }
    }
  }

  # message, feature
  if (name %in% c("message", "feature")) {
    # Should be of class 'character'
    if (!methods::is(input, "character") | length(input) != 1) {
      stop("Input for '", name, "' must be a single value of class 'character'.")
    }
  }

  # title, subtitle, feature_name
  if (name %in% c("title", "subtitle", "feature_name")) {
    # If not NULL
    if (!is.null(input)) {
      # Should be of class 'character'
      if (!methods::is(input, "character") | length(input) != 1) {
        stop("Input for '", name, "' must be a single value of class 'character'. Please supply valid input!")
      }
      if (name == "feature_name") {
        if (other != "feature") {
          warning(" Input for '", name, "' is not used when parameter 'color_by' is set to '", other, "'.")
        }
      }
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
        allowed_functions <- c("estimateSizeFactors", "DESeq")
      } else if (other[[1]] == "limma") {
        if (other[[2]] == "trend") {
          allowed_functions <- c("DGEList", "calcNormFactors", "cpm", "lmFit", "eBayes")
        } else if (other[[2]] == "voom") {
          allowed_functions <- c("DGEList", "calcNormFactors", "voom", "lmFit", "eBayes")
        } else if (other[[2]] %in% c("wilcox_cpm", "wilcox_log_cpm")) {
          allowed_functions <- c("DGEList", "cpm", "rankSumTestWithCorrelation", "lfc")
        }
      } else if (other[[1]] == "presto") {
        allowed_functions <- c("DGEList", "cpm", "wilcoxauc")
      }
      if (!all(names(input) %in% allowed_functions)) {
        stop("When supplying additional parameters to '", name, "' for use",
             " with '", other[[1]], ": ", other[[2]], "',",
             " please provide a list of lists, where each secondary list is named according to",
             " the allowed functions (", paste0(allowed_functions, collapse = ", "), ").")
      }
    }
  }

  # type
  if (name == "type") {
    # Should be of class 'character'
    if (!methods::is(input, "character") | length(input) != 1) {
      stop("Input for '", name, "' must be a single value of class 'character', please supply valid input!")
    }
    # Must be among permitted values
    if (!(input %in% c("discrete", "gradient"))) {
      stop("Input for '", name, "' must be among permitted values (", paste0(c("discrete", "gradient"), collapse = ", "), "), please supply valid input!")
    }
  }

  # reduction
  if (name == "reduction") {
    # Should be of class "matrix"
    if (!methods::is(input, "matrix")) {
      stop("Input value for '", name, "' is not of class 'matrix', please supply valid input!")
    }
    # Must have at least 2 columns
    if (ncol(input) < 2) {
      stop("Input value for '", name, "' must have at least 2 columns, please supply valid input!")
    }
  }

  # color_by
  if (name == "color_by") {
    # If not NULL
    if (!is.null(input)) {
      # Should be of class 'character'
      if (!methods::is(input, "character") | length(input) != 1) {
        stop("Input for '", name, "' must be a single value of class 'character', please supply valid input!")
      }
      # Must be among permitted values
      if (!(input %in% c("split", "n_sig", "pvalue", "feature"))) {
        stop("Input for '", name, "' must be among permitted values (", paste0(c("split", "n_sig", "pvalue", "feature"), collapse = ", "), "), please supply valid input!")
      }
      # Other input can't be NULL
      if (input != "feature" & is.null(other[[1]])) {
        stop("Input for 'split_labels' cannot be NULL when input for '", name, "' is '", input, "', please supply valid input!")
      }
      if (input %in% c("n_sig", "pvalue") & is.null(other[[2]])) {
        stop("Parameter 'input' cannot be NULL when parameter '", name, "' is '", input, "', please supply valid input!")
      }
    }
  }

  # name
  if (name == "palette_name") {
    # If not NULL
    if (!is.null(input)) {
      # Should be of class 'character'
      if (!methods::is(input, "character") | length(input) != 1) {
        stop("Input for '", name, "' must be a single value of class 'character', please supply valid input!")
      }
      # Must be among permitted values for palette type
      if (other == "discrete") {
        if (!(input %in% c("choir", "archr"))) {
          stop("When input for 'type' is 'discrete', input for '", name, "' must be among permitted values (",
               paste0(c("choir", "archr"), collapse = ", "), "), please supply valid input!")
        }
      } else if (other == "gradient") {
        if (!(input %in% c("frozen", "inferno"))) {
          stop("When input for 'type' is 'gradient', input for '", name, "' must be among permitted values (",
               paste0(c("frozen", "inferno"), collapse = ", "), "), please supply valid input!")
        }
      }
    }
  }

  # normalization_method
  if (name == "normalization_method") {
    # Should be of class 'character'
    if (!methods::is(input, "character") | length(input) != 1) {
      stop("Input for '", name, "' must be a single value of class 'character', please supply valid input!")
    }
    # Must be among permitted values
    if (!(input %in% c("cpm", "log_cpm", "none"))) {
      stop("Input for '", name, "' must be among permitted values (", paste0(c("cpm", "log_cpm", "none"), collapse = ", "), "), please supply valid input!")
    }
  }

  # color_limits
  if (name == "color_limits") {
    # If not NULL
    if (!is.null(input)) {
      if (!(other %in% c("pvalue", "n_sig", "feature"))) {
        warning(" Input values for '", name, "' are not used when parameter 'color_by' is set to '", other, "'.")
      } else {
        # Should be of class 'numeric'
        if (!methods::is(input, "numeric") | length(input) != 2) {
          stop("Input for '", name, "' must be a vector containing two values of class 'numeric', please supply valid input!")
        }
      }
    }
  }

  # feature_values
  if (name %in% c("feature_values")) {
    # If not NULL
    if (!is.null(input)) {
      if (other[[3]] != "feature") {
        warning(" Input values for '", name, "' are not used when parameter 'color_by' is set to '", other[[3]], "'.")
      } else {
        # Should be of class "numeric"
        if (!methods::is(input, "numeric")) {
          stop("Input value for '", name, "' must be of class 'numeric', please supply valid input!")
        }
        # Must be of same length as either use cells or reduction matrix
        if (!is.null(other[[2]])) {
          if (length(input) != length(other[[2]])) {
            stop("Input value for '", name, "' must be the same length as input to parameter 'use_cells', please supply valid input!")
          }
        } else {
          if (length(input) != nrow(other[[1]])) {
            stop("Input value for '", name, "' must be the same length as there are rows in the input to parameter 'reduction', please supply valid input!")
          }
        }
      }
    } else {
      # feature_values are required when color_by = "feature"
      if (other[[3]] == "feature") {
        stop("Input value for '", name, "' cannot be NULL when parameter 'color_by' is set to '", other[[3]], "'. Please supply valid input!")
      }
    }
  }

  # label_features
  if (name %in% c("label_features")) {
    # If not NULL
    if (!is.null(input)) {
      # Should be of class 'character'
      if (!methods::is(input, "character")) {
        stop("Input for '", name, "' must be of class 'character', please supply valid input!")
      }
      # If n_max_label < length(input), issue warning
      if (other < length(input)) {
        warning(" When input to parameter '", name, "' is provided, input to parameter 'n_max_label' is disregarded.")
      }
    }
  }

  # design
  if (name %in% c("design")) {
    # If not NULL
    if (!is.null(input)) {
      # Should be a single value of class 'character'
      if (methods::is(input, "character") & length(input) == 1) {
        # Check formula syntax
        if (!grepl("^\\s*~", input)) {
          stop("Input value for '", name, "' must be a one-sided formula starting with '~'. Please supply valid input!")
        }
        try(input_formula <- stats::as.formula(input), silent = TRUE)
        if (!exists("input_formula")) {
          stop("Input value for '", name, "' must be a one-sided formula of proper syntax. Please supply valid input!")
        } else {
          # Parse terms
          terms <- attr(terms(input_formula), "term.labels")
          # Last term must be same as 'group_labels'
          if (terms[length(terms)] != other[[3]]) {
            stop("When input is provided to parameter '", name,
                 "', the last term in the formula must be 'group' or the same as input provided to parameter 'group_labels'. Please supply valid input!")
          }
          terms <- terms[-length(terms)]
          if (length(terms) > 0) {
            # Break up interaction terms
            terms <- unique(unlist(strsplit(terms, ":", fixed = TRUE)))
            # Warn if term "replicate" is in formula, because it will be used to refer to the replicate_labels
            if ("replicate" %in% terms) {
              terms <- terms[terms != "replicate"]
              warning(" Formula provided for parameter '", name,
                      "' includes the term 'replicate'. This will be used to refer to the input provided to parameter 'replicate_labels'. If that is not your intention, please rename the term.")
            }
            # Check for presence of terms in metadata of provided object
            if (!is.null(other[[2]])) {
              if (!all(terms %in% colnames(other[[2]]))) {
                stop("When input for '", name, "' is a character string, the terms must indicate column(s) present in the provided 'metadata', please supply valid input!")
              }
            } else if (methods::is(other[[1]], "Seurat")) {
              if (!all(terms %in% colnames(other[[1]]@meta.data))) {
                stop("When input for '", name, "' is a character string, the terms must indicate column(s) present in the 'meta.data' of the provided object, please supply valid input!")
              }
            } else if (methods::is(other[[1]], "SingleCellExperiment")) {
              if (!all(terms %in% colnames(other[[1]]@colData))) {
                stop("When input for '", name, "' is a character string, the terms must indicate column(s) present in the 'colData' of the provided object, please supply valid input!")
              }
            }
          }
          # Clean up
          rm(input_formula)
        }
      } else {
        stop("Input for '", name, "' must be a single value of class 'character', please supply valid input!")
      }
    }
  }

  # metadata
  if (name %in% c("metadata")) {
    # If not NULL
    if (!is.null(input)) {
      # Must be a dataframe
      if (!methods::is(input, "data.frame")) {
        stop("Input for '", name, "' must be of class 'data.frame', please supply valid input!")
      }
      # Must have rownames corresponding to columns in the provided object
      if (!identical(rownames(input), colnames(other))) {
        stop("Input for '", name, "' must have row names corresponding to each column name in the provided 'object', please supply valid input!")
      }
    }
  }

  # permute_by, permute_within
  if (name %in% c("permute_by", "permute_within")) {
    # If not NULL
    if (!is.null(input)) {
      # Should be of class "character" or "factor"
      if (!methods::is(input, "character") & !methods::is(input, "factor")) {
        if (methods::is(input, "numeric") | methods::is(input, "logical")) {
          warning(" Input value for '", name, "' will be converted to class 'character'.")
        } else {
          stop("Input value for '", name, "' must be a single value of class 'character' or a vector of labels, please supply valid input!")
        }
      }
      # If single value, must be a column name
      if (length(input) == 1) {
        # Must not be 'group'
        if (input == "group") {
          stop("Input for '", name, "' cannot be 'group', please supply valid input!")
        }
        # Must correspond to a column in group_key
        if (!(input %in% colnames(other[[1]]$metadata$group_key))) {
          stop("When a single input value is provided for '", name, "', it ",
               "' must be among the column names of dataframe 'group_key' provided under element ",
               "'metadata' of list provided to parameter 'input', please supply valid input!")
        }
      } else {
        # Length must correspond to group_key rows
        if (length(input) != nrow(other[[1]]$metadata$group_key))
          stop("When a vector is provided for '", name, "', it must contain values corresponding ",
               "to each row of dataframe 'group_key' provided under element ",
               "'metadata' of list provided to parameter 'input', please supply valid input!")
      }
      if (name == "permute_by") {
        # If not a cell-level test, issue warning
        if (other[[2]] != "none") {
          warning(" Input for '", name, "' is intended for use with cell-level tests, where parameter 'pseudobulk' is 'none'.")
        }
      } else if (name == "permute_within") {
        if (is.null(other[[2]])) {
          warning(" Input for '", name, "' is intended for use with complex design formulas, such as paired tests.")
        }
      }
    } else {
      # If a cell-level test, issue warning
      if (name == "permute_by") {
        if (other[[2]] == "none") {
          warning(" When running a cell-level test (parameter 'pseudobulk' is 'none'), consider providing biological replicates labels to ",
                  "parameter '", name, "', such that group labels for each cell in a biological replicate are shuffled together as a unit.")
        }
      }
    }
  }

  # confound_check
  if (name %in% c("confound_check")) {
    # If not NULL
    if (!is.null(input)) {
      # Should be a dataframe
      if (!methods::is(input, "data.frame")) {
        stop("Input value for '", name, "' must be of class 'data.frame', please supply valid input!")
      }
      # Check for correct number of rows
      if (nrow(input) != other) {
        stop("Input value for '", name, "' must a 'data.frame' with rows corresponding to the number of replicates, please supply valid input!")
      }
      # Check class of each column, should not be numeric
      for (i in 1:ncol(input)) {
        if (intersect(methods::is(input[,i]), c("character", "factor", "logical")) < 1) {
          stop("Input value for '", name, "' must a 'data.frame' with columns of classes 'character', 'factor', or 'logical', please supply valid input!")
        }
      }
    }
  }
}
