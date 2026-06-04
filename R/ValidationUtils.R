# -----------------------------------------------------------------------------
# Validation methods
# -----------------------------------------------------------------------------

# Check parameter input validity ---------------------------
#
# Checks validity of provided values for input parameters
#
# input        -- input value
# name         -- name of parameter
# null_allowed -- logical, whether input may be NULL
# class        -- allowed class or classes
# len          -- required length; use "single value" for length 1
# caller       -- function calling parameter (when relevant)
# other        -- other inputs, using the same structure as the original function

.validInput <- function(input,
                        name,
                        null_allowed = FALSE,
                        class = NULL,
                        len = NULL,
                        caller = NULL,
                        other = NULL) {

  # ---------------------------------------------------------------------------
  # General NULL, length, and class checks
  # ---------------------------------------------------------------------------

  if (is.null(input)) {
    if (null_allowed == TRUE) {
      # Warnings
      if (name == "permute_by" && other[[2]] == "none") {
        warning(" When running a cell-level test (parameter 'pseudobulk' is 'none'), consider providing biological replicates labels to parameter '",
                name, "', such that group labels for each cell in a biological replicate are shuffled together as a unit.")
      }
      if (name == "n_group1" && !is.null(other)) {
        warning(" Input for 'n_replicates' is not used when parameter 'n_group1' is NULL.")
      }
      return(invisible(NULL))
    } else {
      stop("Input value for '", name, "' cannot be NULL. Please supply valid input!")
    }
  }

  if (!is.null(class)) {
    if (length(intersect(methods::is(input), class)) < 1) {
      stop("Input value for '", name, "' must be among allowed classes (",
           paste(class, collapse = ", "), "). Please supply valid input!")
    }
  }

  if (!is.null(len)) {
    if (length(input) != len) {
      if (len == 1) {
        stop("Input value for '", name, "' must be a single value. Please supply valid input!")
      } else {
        stop("Input value for '", name, "' must be ", len, " given other provided values. Please supply valid input!")
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Parameters for which no further checks are needed
  # ---------------------------------------------------------------------------

  if (name %in% c("center",
                  "filter",
                  "fix_coords",
                  "label_pvalue",
                  "label_replicates",
                  "label_splits",
                  "label_statistics",
                  "normalize_prefilter",
                  "progress_message",
                  "return_all",
                  "return_raw_de",
                  "swatch",
                  "verbose")) {
    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # Large objects
  # ---------------------------------------------------------------------------

  if (name == "object") {
    if (length(intersect(methods::is(input), c("matrix", "Matrix", "dgCMatrix"))) >= 1 && is.null(rownames(input))) {
      stop("When input value for '", name, "' is of class 'matrix', row names cannot be NULL. Please set row names to feature names.")
    }
    return(invisible(NULL))
  }

  if (name == "metadata") {
    if (!is.null(other) && !identical(rownames(input), colnames(other))) {
      stop("Input for '", name, "' must have row names corresponding to each column name in the provided 'object', please supply valid input!")
    }
    return(invisible(NULL))
  }

  if (name == "reduction") {
    if (ncol(input) < 2) {
      stop("Input value for '", name, "' must have at least 2 columns, please supply valid input!")
    }
    return(invisible(NULL))
  }


  if (name == "input") {
    # For permuteDE
    if (caller == "permuteDE") {
      valid_names <- identical(names(input), c("DE_results", "PB_values", "metadata", "parameters")) ||
        identical(names(input), c("DE_results", "cell_values", "metadata", "parameters"))
      if (!valid_names) {
        stop("Structure of list provided for parameter 'input' is unexpected. It should be the output returned by function 'runDE()' or a list following the same structure with four named elements ('DE_results', 'PB_values' (or 'cell_values'), 'metadata', and 'parameters'). Please supply valid input!")
      }
      if (!("group_key" %in% names(input$metadata))) {
        stop("Structure of list provided for parameter 'input' is unexpected. It should be the output returned by function 'runDE()' or a list following the same structure. Element 'metadata' must contain a dataframe under name 'group_key'. Please supply valid input!")
      }
      required_parameters <- c("reference_group", "non_reference_group", "design_formula", "de_method", "de_test", "de_params", "p_adjust_method", "pseudobulk")
      if (!all(required_parameters %in% names(input$parameters))) {
        stop("Structure of list provided for parameter 'input' is unexpected. It should be the output returned by function 'runDE()' or a list following the same structure. Element 'parameters' must be a list containing a minimal set of named elements ('reference_group', 'non_reference_group', 'design_formula', 'de_method', 'de_test', 'de_params', 'p_adjust_method', 'pseudobulk'). Please supply valid input!")
      }
      if (length(input$parameters$reference_group) != 1 || length(input$parameters$non_reference_group) != 1) {
        stop("Content of list provided for parameter 'input' is unexpected. It should be the output returned by function 'runDE()' or a list following the same structure. The reference/non-reference groups provided under element 'parameters' must be single values. Please supply valid input!")
      }
      if (!(input$parameters$reference_group %in% input$metadata$group_key[, "group"])) {
        stop("Content of list provided for parameter 'input' is unexpected. It should be the output returned by function 'runDE()' or a list following the same structure. The reference group provided under element 'parameters' was not found within the 'group_key' provided under element 'metadata'. Please supply valid input!")
      }
      if (!(input$parameters$non_reference_group %in% input$metadata$group_key[, "group"])) {
        stop("Content of list provided for parameter 'input' is unexpected. It should be the output returned by function 'runDE()' or a list following the same structure. The non-reference group provided under element 'parameters' was not found within the 'group_key' provided under element 'metadata'. Please supply valid input!")
      }
      if (!is.null(input$parameters$design_formula)) {
        if (!methods::is(input$parameters$design_formula, "formula")) {
          stop("Content of list provided for parameter 'input' is unexpected. It should be the output returned by function 'runDE()' or a list following the same structure. The 'design_formula' provided under element 'parameters' must be of class 'formula'. Please supply valid input!")
        }
        terms <- attr(terms(input$parameters$design_formula), "term.labels")
        if (any(grepl("group:", terms) | grepl(":group", terms))) {
          stop("Content of list provided for parameter 'input' is unexpected. The terms within the 'design_formula' provided under element 'parameters' include interaction term(s) that involve the primary comparison groups to be permuted. The 'permuteDE()' function is not compatible with these interaction terms. Please supply valid input!")
        }
        terms <- unique(unlist(strsplit(terms, ":", fixed = TRUE)))
        if (!all(terms %in% colnames(input$metadata$group_key))) {
          stop("Content of list provided for parameter 'input' is unexpected. It should be the output returned by function 'runDE()' or a list following the same structure. The terms within the 'design_formula' provided under element 'parameters' must correspond to the column names of dataframe 'group_key' provided under element 'metadata'. Please supply valid input!")
        }
        for (t in terms) {
          if (any(is.na(input$metadata$group_key[, t]))) {
            stop("Content of list provided for parameter 'input' is unexpected. It should be the output returned by function 'runDE()' or a list following the same structure. Dataframe 'group_key' provided under element 'metadata' cannot contain NA values. Please supply valid input!")
          }
        }
      }
    }
    # For plots
    if (caller == "plotVolcano" && !("DE_results" %in% names(input))) {
      stop("Structure of list provided for parameter 'input' is unexpected, it should be the output returned by function 'runDE()' or a list containing (at minimum) a dataframe named 'DE_results'. Please supply valid input!")
    }
    if (caller == "plotHistogram" && (!("permutation_test_results" %in% names(input)) || !("permutation_DE_summary" %in% names(input)))) {
      stop("Structure of list provided for parameter 'input' is unexpected, it should be the output returned by function 'permuteDE()' or a list containing (at minimum) dataframes named 'permutation_test_results' and 'permutation_DE_summary'. Please supply valid input!")
    }
    if (caller == "plotFeature") {
      if (!all(c("metadata", "parameters") %in% names(input))) {
        stop("Structure of list provided for parameter 'input' is unexpected, it should be the output returned by function 'runDE()' or a list containing (at minimum) elements named 'PB_values' (or 'cell_values'), 'metadata', and 'parameters'. Please supply valid input!")
      }
      if (!("group_key" %in% names(input$metadata))) {
        stop("Structure of list provided for parameter 'input' is unexpected. It should be the output returned by function 'runDE()' or a list following the same structure. Element 'metadata' must contain a dataframe under name 'group_key'. Please supply valid input!")
      }
    }
    if (caller == "plotDimReduction" && !("permutation_test_results" %in% names(input))) {
      stop("Structure of list provided for parameter 'input' is unexpected, it should be the output returned by function 'permuteDE()' or a list containing (at minimum) an element named 'permutation_test_results'. Please supply valid input!")
    }
    return(invisible(NULL))
  }

  if (name == "confound_check") {
    if (nrow(input) != other) {
      stop("Input value for '", name, "' must a 'data.frame' with rows corresponding to the number of replicates, please supply valid input!")
    }
    valid_cols <- vapply(input, function(x) {
      length(intersect(methods::is(x), c("character", "factor", "logical"))) >= 1
    }, logical(1))
    if (!all(valid_cols)) {
      stop("Input value for '", name, "' must a 'data.frame' with columns of classes 'character', 'factor', or 'logical', please supply valid input!")
    }
    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # Labels / metadata column parameters
  # ---------------------------------------------------------------------------

  if (name %in% c("replicate_labels",
                  "group_labels",
                  "split_labels")) {
    # split_label has particular requirements for function plotDimReduction
    if (name == "split_labels" && !is.null(caller) && caller == "plotDimReduction") {
      if (is.null(other[[2]]) || !(other[[2]] %in% c("split", "n_sig", "pvalue"))) {
        warning(" Input for '", name, "' is not used for value '", other[[2]], "' of parameter 'color_by'.")
        return(invisible(NULL))
      }
      if (length(input) != nrow(other[[1]])) {
        stop("Input value for '", name,
             "' must be a vector of labels corresponding to each row of the input provided to parameter 'reduction', please supply valid input!")
      }
      return(invisible(NULL))
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
    # replicate_labels are not used for cell-level tests
    if (name == "replicate_labels") {
      if (other[[3]] == "none") {
        warning(" Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'none'.")
      }
    }
    # Numeric/logical values will be converted to character
    if (methods::is(input, "numeric") | methods::is(input, "logical")) {
      warning(" Input value for '", name, "' will be converted to class 'character'.")
    }
    return(invisible(NULL))
  }

  if (name %in% c("permute_by",
                  "permute_within")) {
    # If single value, must be a column name
    if (length(input) == 1) {
      # Must not be 'group'
      if (input == "group") {
        stop("Input for '", name, "' cannot be 'group', please supply valid input!")
      }
      # Must correspond to a column in group_key
      if (!(input %in% colnames(other[[1]]$metadata$group_key))) {
        stop("When a single input value is provided for '", name,
             "', it must be among the column names of dataframe 'group_key' provided under element 'metadata' of list provided to parameter 'input', please supply valid input!")
      }
    } else {
      # Length must correspond to group_key rows
      if (length(input) != nrow(other[[1]]$metadata$group_key)) {
        stop("When a vector is provided for '", name,
             "', it must contain values corresponding to each row of dataframe 'group_key' provided under element 'metadata' of list provided to parameter 'input', please supply valid input!")
      }
    }
    if (name == "permute_by" && other[[2]] != "none") {
      warning(" Input for '", name, "' is intended for use with cell-level tests, where parameter 'pseudobulk' is 'none'.")
    }
    if (name == "permute_within" && is.null(other[[2]])) {
      warning(" Input for '", name, "' is intended for use with complex design formulas, such as paired tests.")
    }
    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # Cell, assay, and layer selection
  # ---------------------------------------------------------------------------

  if (name == "use_cells") {
    if (other[[2]] == "supplied") {
      warning(" Input value for '", name,
              "' is not used when parameter 'pseudobulk' is set to 'supplied' (when a pre-computed pseudobulk matrix is supplied by the user).")
      return(invisible(NULL))
    }
    if (is.null(other[[1]])) {
      warning(" Input value for '", name,
              "' is not used when parameter 'object' is set to 'NULL'.")
      return(invisible(NULL))
    }
    if (length(intersect(input, colnames(other[[1]]))) != length(input)) {
      stop("Not all provided cells are present in the provided object, please supply valid input!")
    }
    return(invisible(NULL))
  }

  if (name == "use_assay") {
    if (is.null(other)) {
      warning(" Input value for '", name, "' is not used when 'object' is NULL.")
      return(invisible(NULL))
    }
    if (length(intersect(methods::is(other), c("matrix", "Matrix", "dgCMatrix"))) >= 1) {
      warning(" Input value for '", name, "' is not used when 'object' is of class 'matrix'.")
      return(invisible(NULL))
    }
    if (methods::is(other, "Seurat") || methods::is(other, "SingleCellExperiment")) {
      if (!(input %in% names(other@assays))) {
        stop("Assay '", input, "' provided for parameter '", name, "' is not present in provided object, please supply valid input!")
      }
    }
    return(invisible(NULL))
  }

  if (name == "use_layer") {
    if (!methods::is(other[[1]], "Seurat")) {
      warning(" Input value for '", name, "' are not used when provided object is not of class 'Seurat'.")
      return(invisible(NULL))
    }
    if (is.null(other[[2]])) {
      other[[2]] <- Seurat::DefaultAssay(other[[1]])
    }
    if ("Assay5" %in% methods::is(other[[1]][[other[[2]]]])) {
      if (!(input %in% names(other[[1]][[other[[2]]]]@layers))) {
        stop("Layer '", input, "' is not present in assay '", other[[2]], "' of provided Seurat v5 object, please supply valid input!")
      }
    } else {
      slot_exists <- tryCatch({
        methods::slot(other[[1]][[other[[2]]]], input)
        TRUE
      }, error = function(e) FALSE)
      if (!slot_exists) {
        stop("Slot '", input, "' is not present in assay '", other[[2]], "' of provided Seurat object, please supply valid input!")
      }
    }
    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # Numeric and logical parameters
  # ---------------------------------------------------------------------------

  # Positive integer
  if (name %in% c("min_cells_per_feature",
                  "min_cells_per_split",
                  "min_cells_per_replicate",
                  "min_replicates_per_split",
                  "min_replicates_per_group",
                  "n",
                  "n_combinations",
                  "n_cores",
                  "n_iterations",
                  "n_replicates",
                  "random_seed")) {
    # Must be positive integer
    for (i in input) {
      if (i %% 1 != 0 | i < 1) {
        if (name %in% c("n_replicates")) {
          stop("Input value for '", name, "' must be positive integer(s). Please supply valid input!")
        } else {
          stop("Input value for '", name, "' must be a positive integer. Please supply valid input!")
        }
      }
    }
    # n_iterations >= 2
    if (name == "n_iterations" && input < 2) {
        stop("Input value for '", name, "' must be at least 2. Please supply valid input!")
    }
    # Warnings
    if (name %in% c("min_cells_per_split", "min_cells_per_replicate", "min_cells_per_feature") && other == "supplied") {
      warning(" Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'supplied' (when a pre-computed pseudobulk matrix is supplied by the user).")
    }
    if (name %in% c("min_cells_per_replicate", "min_replicates_per_split") && other == "none" && !is.null(caller) && caller != "runDE") {
      warning(" Input value for '", name, "' is not used when parameter 'pseudobulk' is set to 'none'.")
    }
    return(invisible(NULL))
  }

  if (name == "n_group1") {
    for (i in input) {
      if (i %% 1 != 0 | i < 1) {
        stop("Input value for '", name, "' must be positive integer(s). Please supply valid input!")
      }
    }
    if (is.null(other)) {
      warning(" Input for '", name, "' is not used when parameter 'n_replicates' is NULL.")
      return(invisible(NULL))
    }
    if (length(input) != length(other)) {
      stop("Input for '", name, "' must be a vector of the same length as input to 'n_replicates'. Please supply valid input!")
    }
    if (any(input >= other)) {
      stop("Each input value for '", name,
           "' must be less than the corresponding input value for 'n_replicates'. Please supply valid input!")
    }
    return(invisible(NULL))
  }

  # Non-negative integer
  if (name %in% c("min_DE", "n_max_label")) {
    if (input %% 1 != 0 | input < 0) {
      stop("Input value for '", name, "' must be a non-negative integer. Please supply valid input!")
    }
    return(invisible(NULL))
  }

  # Number from 0-1
  if (name %in% c("alpha", "min_prop_cells_per_feature", "permutation_test_alpha")) {
    if (input < 0) {
      stop("Input value for '", name, "' cannot be negative. Please supply valid input!")
    }
    if (input > 1) {
      stop("Input value for '", name, "' cannot be greater than 1. Please supply valid input!")
    }
    if (name == "min_prop_cells_per_feature" && other == "supplied") {
      warning(" Input value for '", name,
              "' is not used when parameter 'pseudobulk' is set to 'supplied' (when a pre-computed pseudobulk matrix is supplied by the user).")
    }
    if (name == "permutation_test_alpha" && input < 1) {
      if (is.null(other)) {
        warning(" Input value for '", name, "' is not used when parameter 'color_by' is set to NULL.")
      } else if (!(other %in% c("pvalue", "split"))) {
        warning(" Input value for '", name, "' is not used when parameter 'color_by' is set to '", other, "'.")
      }
    }
    return(invisible(NULL))
  }

  # Non-negative number
  if (name == "lfc_threshold") {
    if (input < 0) {
      stop("Input value for '", name, "' cannot be negative. Please supply valid input!")
    }
    return(invisible(NULL))
  }

  if (name == "color_limits") {
    if (is.null(other)) {
      warning(" Input values for '", name, "' are not used when parameter 'color_by' is set to NULL.")
      return(invisible(NULL))
    }
    if (!(other %in% c("pvalue", "n_sig", "feature"))) {
      warning(" Input values for '", name, "' are not used when parameter 'color_by' is set to '", color_by, "'.")
      return(invisible(NULL))
    }
    return(invisible(NULL))
  }

  if (name == "feature_values") {
    if (is.null(other[[3]])) {
      warning(" Input values for '", name, "' are not used when parameter 'color_by' is set to NULL.")
      return(invisible(NULL))
    }
    if (other[[3]] != "feature") {
      warning(" Input values for '", name, "' are not used when parameter 'color_by' is set to '", other[[3]], "'.")
      return(invisible(NULL))
    }
    if (!is.null(other[[2]])) {
      if (length(input) != length(other[[2]])) {
        stop("Input value for '", name, "' must be the same length as input to parameter 'use_cells', please supply valid input!")
      }
    } else {
      if (length(input) != nrow(other[[1]])) {
        stop("Input value for '", name, "' must be the same length as there are rows in the input to parameter 'reduction', please supply valid input!")
      }
    }
    return(invisible(NULL))
  }

  # Logical (w/ additional requirements)
  if (name %in% c("force_balance")) {
    if (other == "none" && isTRUE(input)) {
      stop("Parameter '", name, "' cannot be set to '", input, "' when parameter 'pseudobulk' is set to 'none'. Please supply valid input!")
    }
    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # Character parameters with allowed values
  # ---------------------------------------------------------------------------

  # Simple cases
  if (name %in% c("normalization_method",
                  "p_adjust_method",
                  "plot_type",
                  "type")) {
    allowed_values <- switch(name,
                             normalization_method = c("cpm", "log_cpm", "none"),
                             p_adjust_method = stats::p.adjust.methods,
                             plot_type = c("boxplot", "bar_se", "bar_sd", "beeswarm"),
                             type = c("discrete", "gradient"))
    if (!(input %in% allowed_values)) {
      stop("Input for '", name, "' must be among allowed values (",
           paste0(allowed_values, collapse = ", "), "), please supply valid input!")
    }
    return(invisible(NULL))
  }

  if (name == "de_method") {
    allowed_values <- c("edgeR", "DESeq2", "limma", "presto", "BPCells")
    if (!(input %in% allowed_values)) {
      stop("Input for '", name, "' must be among allowed values (",
           paste0(allowed_values, collapse = ", "), "), please supply valid input!")
    }
    input_matrix <- NULL
    if (caller == "runDE") {
      input_matrix <- .getMatrix(object = other[[2]],
                                 use_assay = other[[3]],
                                 use_layer = other[[4]],
                                 use_cells = other[[5]],
                                 verbose = FALSE)
    } else if (caller == "permuteDE") {
      input_matrix <- if (other[[1]] == "none") {
        other[[2]]$cell_values[[1]]
      } else {
        other[[2]]$PB_values[[1]]
      }
    }
    # BPCells checks
    if (input == "BPCells" && !(other[[1]] %in% c("supplied", "none"))) {
      stop("Input for '", name, "' can only be 'BPCells' when parameter ",
           "'pseudobulk' is 'supplied' or 'none' and counts are provided as ",
           "class 'IterableMatrix'. Please supply valid input!")
    }
    if (other[[1]] %in% c("supplied", "none") && input == "BPCells" && !(methods::is(input_matrix, "IterableMatrix"))) {
      stop("When input for '", name, "' is 'BPCells' and parameter 'pseudobulk' ",
           "is 'supplied' or 'none', counts must be provided as class ",
           "'IterableMatrix'. Please supply valid input!")
    }

    if (other[[1]] %in% c("supplied", "none") && input != "BPCells" && methods::is(input_matrix, "IterableMatrix")) {
      stop("When counts are provided as class 'IterableMatrix' and input to ",
           "parameter 'pseudobulk' is '", other[[1]], "', input for '",
           name, "' must be 'BPCells'. Please supply valid input!")
    }
    return(invisible(NULL))
  }

  if (name == "de_test") {
    allowed_values <- switch(other,
                             edgeR = c("LRT", "QLF", "exact"),
                             DESeq2 = c("LRT", "Wald"),
                             limma = c("trend", "voom", "wilcox_cpm", "wilcox_log_cpm"),
                             presto = c("wilcox_cpm", "wilcox_log_cpm"),
                             BPCells = c("wilcox_cpm", "wilcox_log_cpm"))
    if (!(input %in% allowed_values)) {
      stop("When input for 'de_method' is '", other, "', input for '", name, "' must be among allowed values (",
           paste0(allowed_values, collapse = ", "), "), please supply valid input!")
    }
    return(invisible(NULL))
  }

  if (name == "de_params") {
    if (length(input) == 0) return(invisible(NULL))

    de_method <- other[[1]]
    de_test <- other[[2]]
    return_raw_de <- other[[3]]

    function_values <- switch(de_method,
      edgeR = switch(de_test,
        LRT = c("DGEList", "calcNormFactors", "estimateDisp", "glmFit", "glmLRT"),
        QLF = c("DGEList", "calcNormFactors", "estimateDisp", "glmQLFit", "glmQLFTest"),
        exact = c("DGEList", "calcNormFactors", "estimateDisp", "exactTest")),
      DESeq2 = c("estimateSizeFactors", "DESeq"),
      limma = switch(de_test,
        trend = c("DGEList", "calcNormFactors", "cpm", "lmFit", "eBayes"),
        voom = c("DGEList", "calcNormFactors", "voom", "lmFit", "eBayes"),
        wilcox_cpm = c("DGEList", "cpm", "rankSumTestWithCorrelation", "lfc"),
        wilcox_log_cpm = c("DGEList", "cpm", "rankSumTestWithCorrelation", "lfc")),
      presto = c("DGEList", "cpm", "wilcoxauc"),
      BPCells = c("DGEList", "cpm", "marker_features", "lfc"))

    allowed_values <- c(function_values, "return_all_coefficients")

    if (!all(names(input) %in% allowed_values)) {
      stop("When supplying additional parameters to '", name, "' for use with '",
        de_method, ": ", de_test,
        "', please provide a list where each element is named according to the allowed functions/options (",
        paste0(allowed_values, collapse = ", "), ").")
    }

    if ("return_all_coefficients" %in% names(input) &&
        (!is.logical(input[["return_all_coefficients"]]) ||
         length(input[["return_all_coefficients"]]) != 1)) {
      stop("Input 'de_params$return_all_coefficients' must be a single value of class logical.")
    }

    supports_all_coefficients <- ((de_method == "edgeR" && de_test %in% c("LRT", "QLF")) ||
        (de_method == "limma" && de_test %in% c("trend", "voom")) ||
        (de_method == "DESeq2" && de_test == "Wald"))

    if (isTRUE(input[["return_all_coefficients"]]) &&
        !supports_all_coefficients) {
      stop("Input 'de_params$return_all_coefficients = TRUE' is only supported for ",
        "coefficient-based model tests: edgeR LRT/QLF, limma trend/voom, and DESeq2 Wald.")
    }

    if ("return_all_coefficients" %in% names(input) &&
        isFALSE(return_raw_de)) {
      warning("Value provided for 'return_all_coefficients' in 'de_params' will not be used ",
        "when parameter 'return_raw_de' is FALSE.")
    }

    # All backend function parameters should still be supplied as lists.
    # Special top-level options such as return_all_coefficients are checked separately.
    function_param_names <- intersect(names(input), function_values)

    if (!all(vapply(input[function_param_names], is.list, logical(1)))) {
      stop("When supplying function-specific parameters to '", name, "', each ",
        "function-specific element must be a list. The only non-list option currently ",
        "allowed is 'return_all_coefficients'.")
    }

    # The tested coefficient/contrast is controlled internally
    # so standardized DE_results always correspond to the group effect
    # Cannot override this through de_params
    if (de_method == "edgeR" && de_test == "LRT" && "glmLRT" %in% names(input) &&
        any(c("coef", "contrast") %in% names(input[["glmLRT"]]))) {
      stop("Do not supply 'coef' or 'contrast' in 'de_params$glmLRT'. ",
        "The tested coefficient is controlled internally so 'DE_results' always ",
        "corresponds to the group effect. Use 'de_params$return_all_coefficients = TRUE' ",
        "with 'return_raw_de = TRUE' to inspect all coefficient-level results.")
    }
    if (de_method == "edgeR" && de_test == "QLF" && "glmQLFTest" %in% names(input) &&
        any(c("coef", "contrast") %in% names(input[["glmQLFTest"]]))) {
      stop("Do not supply 'coef' or 'contrast' in 'de_params$glmQLFTest'. ",
        "The tested coefficient is controlled internally so 'DE_results' always ",
        "corresponds to the group effect. Use 'de_params$return_all_coefficients = TRUE' ",
        "with 'return_raw_de = TRUE' to inspect all coefficient-level results.")
    }
    if (de_method == "DESeq2" && "DESeq" %in% names(input) &&
        any(c("test", "reduced", "full") %in% names(input[["DESeq"]]))) {
      stop( "Do not supply 'test', 'reduced', or 'full' in 'de_params$DESeq'. ",
        "These are controlled internally so 'DE_results' always corresponds to ",
        "the requested DESeq2 test and group effect.")
    }

    return(invisible(NULL))
  }

  if (name == "pseudobulk") {
    allowed_values <- c("generate", "supplied", "none")
    if (!(input %in% allowed_values)) {
      stop("Input for '", name, "' must be among allowed values (",
           paste0(allowed_values, collapse = ", "), "), please supply valid input!")
    }
    if (caller != "runDE" && input == "none") {
      warning(" Cell-level tests are not recommended in most cases, proceed with caution.")
    }
    if (caller == "permuteDE") {
      if (input == "none" && !("cell_values" %in% names(other))) {
        stop("Structure of list provided for parameter 'input' is unexpected. When conducting cell-level tests, it should be the output returned by function 'runDE()' or a list following the same structure with four named elements ('DE_results', 'cell_values', 'metadata', and 'parameters'). Please supply valid input!")
      }
      if (input != "none" && !("PB_values" %in% names(other))) {
        stop("Structure of list provided for parameter 'input' is unexpected. When conducting pseudobulk tests, it should be the output returned by function 'runDE()' or a list following the same structure with four named elements ('DE_results', 'PB_values', 'metadata', and 'parameters'). Please supply valid input!")
      }
    } else {
      if (input == "supplied" && length(intersect(methods::is(other), c("Seurat", "SingleCellExperiment"))) > 0) {
        stop("When input for '", name, "' is '", input, "', parameter 'object' must be of class 'matrix'. Please supply valid input!")
      }
    }
    return(invisible(NULL))
  }

  if (name == "use_splits") {
    if (caller %in% c("permuteDE", "plotFeature")) {
      if ("PB_values" %in% names(other)) {
        if (!all(input %in% names(other$PB_values))) {
          stop("Input value(s) for '", name, "' must all be present among provided pseudobulk matrices. Please supply valid input!")
        }
      } else if ("cell_values" %in% names(other)) {
        if (!all(input %in% names(other$cell_values))) {
          stop("Input value(s) for '", name, "' must all be present among provided matrices. Please supply valid input!")
        }
      }
    } else if (caller == "plotVolcano") {
      if (!all(input %in% other$DE_results$split)) {
        stop("Input value(s) for '", name, "' must all be present in DE results. Please supply valid input!")
      }
    } else if (caller == "plotHistogram") {
      if (!all(input %in% other$permutation_DE_summary$split)) {
        stop("Input value(s) for '", name, "' must all be present in permutation test summary. Please supply valid input!")
      }
    }
    return(invisible(NULL))
  }

  if (name == "reference_group") {
    if (length(input) != 1) {
      stop("Input for '", name, "' must be a single value, please supply valid input!")
    }
    if (length(other[[3]]) == 1) {
      groups <- .retrieveData(object = other[[1]],
                              metadata = other[[2]],
                              type = "cell_metadata",
                              name = other[[3]],
                              use_cells = other[[4]])
    } else {
      groups <- other[[3]]
    }
    if (!(input %in% groups)) {
      stop("Input value for '", name, "' must be present among provided group labels. Please supply valid input!")
    }
    return(invisible(NULL))
  }

  if (name == "color_by") {
    allowed_values <- c("split", "n_sig", "pvalue", "feature")
    if (!(input %in% allowed_values)) {
      stop("Input for '", name, "' must be among allowed values (",
           paste0(allowed_values, collapse = ", "), "), please supply valid input!")
    }
    if (input != "feature" && is.null(other[[1]])) {
      stop("Input for 'split_labels' cannot be NULL when input for '", name, "' is '", input, "', please supply valid input!")
    }
    if (input %in% c("n_sig", "pvalue") && is.null(other[[2]])) {
      stop("Parameter 'input' cannot be NULL when parameter '", name, "' is '", input, "', please supply valid input!")
    }
    return(invisible(NULL))
  }

  if (name == "palette") {
    is_color <- vapply(input, function(i) {
      tryCatch(is.matrix(grDevices::col2rgb(i)), error = function(e) FALSE)
    }, logical(1))
    if (all(is_color)) return(invisible(NULL))
    if (length(input) != 1) {
      stop("Input for '", name, "' must be either a palette name or a vector of color values, please supply valid input!")
    }
    if (other == "discrete" && !(input %in% c("choir", "archr"))) {
      stop("When input for 'type' is 'discrete', input for '", name, "' must be among allowed values (choir, archr), please supply valid input!")
    }
    if (other == "gradient" && !(input %in% c("frozen", "inferno"))) {
      stop("When input for 'type' is 'gradient', input for '", name, "' must be among allowed values (frozen, inferno), please supply valid input!")
    }
    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # Other character parameters
  # ---------------------------------------------------------------------------

  if (name == "feature_name") {
    if (!is.null(caller) && caller == "plotDimReduction") {
      if (is.null(other)) {
        warning(" Input for '", name, "' is not used when parameter 'color_by' is set to NULL.")
      } else if (other != "feature") {
        warning(" Input for '", name, "' is not used when parameter 'color_by' is set to '", color_by, "'.")
      }
    }
    return(invisible(NULL))
  }

  if (name == "label_features") {
    if (other < length(input)) {
      warning(" When input to parameter '", name, "' is provided, input to parameter 'n_max_label' is disregarded.")
    }
    return(invisible(NULL))
  }

  if (name == "design") {
    # Stop if DE test is incompatible
    if (other[[4]] %in% c("wilcox_cpm", "wilcox_log_cpm", "exact")) {
      stop("Complex formulas for parameter '", name,
           "' are not supported when input to parameter 'de_test' is '", other[[4]],
           "'. This test does not model covariates.")
    }

    if (!grepl("^\\s*~", input)) {
      stop("Input value for '", name, "' must be a one-sided formula starting with '~'. Please supply valid input!")
    }
    input_formula <- tryCatch(stats::as.formula(input), error = function(e) NULL)
    if (is.null(input_formula)) {
      stop("Input value for '", name, "' must be a one-sided formula of proper syntax. Please supply valid input!")
    }
    terms <- attr(terms(input_formula), "term.labels")
    if (terms[length(terms)] != other[[3]]) {
      stop("When input is provided to parameter '", name,
           "', the last term in the formula must be 'group' or the same as input provided to parameter 'group_labels'. Please supply valid input!")
    }
    terms <- terms[-length(terms)]
    if (any(grepl("\\(", terms))) {
      stop("permuteDE is not currently compatible with non-standard terms (or terms based on function calls) within input provided for '",
           name, "'. Please supply valid input!")
    }
    if (length(terms) > 0) {
      terms <- unique(unlist(strsplit(terms, ":", fixed = TRUE)))
      if ("replicate" %in% terms) {
        terms <- terms[terms != "replicate"]
        warning(" Formula provided for parameter '", name,
                "' includes the term 'replicate'. This will be used to refer to the input provided to parameter 'replicate_labels'. If that is not your intention, please rename the term.")
      }
      if (!is.null(other[[2]])) {
        if (!all(terms %in% colnames(other[[2]]))) {
          stop("When input is provided for '", name, "', the terms must indicate column(s) present in the provided 'metadata', please supply valid input!")
        }
      } else if (methods::is(other[[1]], "Seurat")) {
        if (!all(terms %in% colnames(other[[1]]@meta.data))) {
          stop("When input is provided for '", name, "', the terms must indicate column(s) present in the 'meta.data' of the provided object, please supply valid input!")
        }
      } else if (methods::is(other[[1]], "SingleCellExperiment")) {
        if (!all(terms %in% colnames(other[[1]]@colData))) {
          stop("When input is provided for '", name, "', the terms must indicate column(s) present in the 'colData' of the provided object, please supply valid input!")
        }
      }
    }
    return(invisible(NULL))
  }

  # ---------------------------------------------------------------------------
  # Unrecognized name
  # ---------------------------------------------------------------------------

  stop("No validation rule has been implemented for parameter '", name, "'.")
}

