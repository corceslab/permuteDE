# ---------------------------------------------------------------------------
# Tests for function runDE
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------

# Make a test count matrix
#
# n_genes     -- Number of genes
# random_seed -- Random seed
makeCounts <- function(n_genes = 100,
                       random_seed = 1) {
  set.seed(random_seed)
  n_cells <- 12
  group <- rep(c("A", "A", "A", "B", "B", "B"), each = 2)
  # Baseline gene means vary across genes
  base_mean <- stats::rgamma(n_genes, shape = 2, rate = 0.08)
  # Add gene-specific dispersion; larger size = less dispersion
  size <- stats::runif(n_genes, min = 5, max = 30)
  # Some genes are up in A, some up in B, most unchanged
  lfc <- rep(0, n_genes)
  lfc[seq_len(10)] <- 1
  lfc[11:20] <- -1
  mat <- matrix(0L, nrow = n_genes, ncol = n_cells)
  for (g in seq_len(n_genes)) {
    mu <- base_mean[g] * ifelse(group == "A", 2^lfc[g], 1)
    # Add mild cell-level library size variation
    lib_factor <- stats::runif(n_cells, min = 0.8, max = 1.2)
    mat[g, ] <- stats::rnbinom(
      n = n_cells,
      mu = mu * lib_factor,
      size = size[g]
    )
  }
  rownames(mat) <- paste0("gene", seq_len(n_genes))
  colnames(mat) <- paste0("cell", seq_len(n_cells))
  # Ensure no all-zero genes in this small test matrix
  zero_genes <- Matrix::rowSums(mat) == 0
  if (any(zero_genes)) {
    mat[zero_genes, ] <- 1L
  }
  return(mat)
}

# Make a test cell metadata dataframe
makeMetaData <- function() {
  data.frame(replicate = rep(paste0("rep", 1:6), each = 2),
             group = rep(c("A", "A", "A", "B", "B", "B"), each = 2),
             split = rep("all", 12),
             row.names = paste0("cell", 1:12))
}

# Make test pseudobulk matrix
makePB <- function() {
  counts <- makeCounts()
  metadata <- makeMetaData()
  reps <- unique(metadata$replicate)
  pb <- vapply(
    reps,
    FUN = function(r) {
      Matrix::rowSums(counts[, metadata$replicate == r, drop = FALSE])
    },
    FUN.VALUE = numeric(nrow(counts)))
  rownames(pb) <- rownames(counts)
  colnames(pb) <- reps
  return(pb)
}

# Make test group labels
makeGroups <- function() {
  c("A", "A", "A", "B", "B", "B")
}

# Get expected output value name
#
# pseudobulk -- Pseudobulk setting
expectedValueName <- function(pseudobulk) {
  if (identical(pseudobulk, "none")) {
    return("cell_values")
  }
  return("PB_values")
}

# Make set of expectations
#
# output -- runDE output
expectValidOutput <- function(output) {
  expected_name <- expectedValueName(output$parameters$pseudobulk)
  expect_named(output, c("DE_results", expected_name, "metadata", "parameters"))

  # Check DE_results structure
  expect_s3_class(output$DE_results, "data.frame")
  expect_true(all(c("feature", "lfc", "pvalue", "padj", "split") %in% colnames(output$DE_results)))
  expect_gt(nrow(output$DE_results), 0)

  # Check DE_results content
  expect_type(output$DE_results$feature, "character")
  expect_false(anyNA(output$DE_results$feature))
  expect_equal(anyDuplicated(output$DE_results[, c("feature", "split")]), 0L)
  expect_true(is.numeric(output$DE_results$lfc))
  expect_true(is.numeric(output$DE_results$pvalue))
  expect_true(is.numeric(output$DE_results$padj))
  expect_true(all(is.finite(output$DE_results$lfc) | is.na(output$DE_results$lfc)))
  expect_true(all(output$DE_results$pvalue >= 0 | is.na(output$DE_results$pvalue)))
  expect_true(all(output$DE_results$pvalue <= 1 | is.na(output$DE_results$pvalue)))
  expect_true(all(output$DE_results$padj >= 0 | is.na(output$DE_results$padj)))
  expect_true(all(output$DE_results$padj <= 1 | is.na(output$DE_results$padj)))
  expect_type(output$DE_results$split, "character")
  expect_false(anyNA(output$DE_results$split))

  # Check that split labels match output matrix/cell list names
  expect_setequal(unique(output$DE_results$split), names(output[[expected_name]]))

  # Check adjusted p-values are consistent with p-values within each split
  for (split_i in unique(output$DE_results$split)) {
    rows_i <- output$DE_results$split == split_i

    expect_equal(
      output$DE_results$padj[rows_i],
      stats::p.adjust(
        output$DE_results$pvalue[rows_i],
        method = output$parameters$p_adjust_method
      ),
      tolerance = 1e-12
    )
  }

  # Check output matrices/cell values
  expect_type(output[[expected_name]], "list")
  expect_gt(length(output[[expected_name]]), 0)
  if (identical(expected_name, "PB_values")) {
    expect_false("cell_values" %in% names(output))
  } else {
    expect_false("PB_values" %in% names(output))
  }
  expect_true("group_key" %in% names(output$metadata))
  expect_true("time" %in% names(output$metadata))
}

# DE method/test grid
runDE_method_grid <- data.frame(
  de_method = c("edgeR", "edgeR", "edgeR",
                "DESeq2", "DESeq2",
                "limma", "limma", "limma", "limma",
                "presto", "presto",
                "BPCells", "BPCells"),
  de_test = c("LRT", "QLF", "exact",
              "LRT", "Wald",
              "trend", "voom", "wilcox_cpm", "wilcox_log_cpm",
              "wilcox_cpm", "wilcox_log_cpm",
              "wilcox_cpm", "wilcox_log_cpm"),
  stringsAsFactors = FALSE)

# Skip tests for uninstalled packages
#
# de_method -- DE method
skipUninstalled <- function(de_method) {
  if (de_method %in% c("edgeR", "limma", "presto", "BPCells")) {
    testthat::skip_if_not_installed("edgeR")
  }
  if (de_method == "DESeq2") {
    testthat::skip_if_not_installed("DESeq2")
  }
  if (de_method == "limma") {
    testthat::skip_if_not_installed("limma")
  }
  if (de_method == "presto") {
    testthat::skip_if_not_installed("presto")
  }
  if (de_method == "BPCells") {
    testthat::skip_if_not_installed("BPCells")
  }
  testthat::skip_if_not_installed("pbmcapply")
}

# ---------------------------------------------------------------------------
# Set up input types
# ---------------------------------------------------------------------------

# Set up cell-level matrix input
setInput.CellMatrix <- function() {
  list(object = makeCounts(),
       metadata = makeMetaData(),
       replicate_labels = "replicate",
       group_labels = "group",
       split_labels = NULL,
       use_assay = NULL,
       use_layer = NULL)
}

# Set up PB matrix input
setInput.PBMatrix <- function() {
  pb <- makePB()
  list(object = pb,
       metadata = NULL,
       replicate_labels = colnames(pb),
       group_labels = makeGroups(),
       split_labels = NULL,
       use_assay = NULL,
       use_layer = NULL)
}

# Set up Seurat input with dgCMatrix
setInput.Seurat <- function() {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("Matrix")
  counts <- Matrix::Matrix(makeCounts(), sparse = TRUE)
  metadata <- makeMetaData()
  object <- Seurat::CreateSeuratObject(counts = counts,
                                       meta.data = metadata)
  list(object = object,
       metadata = NULL,
       replicate_labels = "replicate",
       group_labels = "group",
       split_labels = NULL,
       use_assay = Seurat::DefaultAssay(object),
       use_layer = "counts")
}

# Set up Seurat input with BPCells
setInput.BPCells <- function() {
  testthat::skip_if_not_installed("Seurat")
  testthat::skip_if_not_installed("BPCells")
  testthat::skip_if_not_installed("Matrix")
  counts <- Matrix::Matrix(makeCounts(), sparse = TRUE)
  metadata <- makeMetaData()
  tmp_dir <- tempfile("bpcells_counts_")
  BPCells::write_matrix_dir(mat = counts,
                            dir = tmp_dir)
  counts_bp <- BPCells::open_matrix_dir(tmp_dir)
  object <- Seurat::CreateSeuratObject(counts = counts_bp,
                                       meta.data = metadata)
  list(object = object,
       metadata = NULL,
       replicate_labels = "replicate",
       group_labels = "group",
       split_labels = NULL,
       use_assay = Seurat::DefaultAssay(object),
       use_layer = "counts")
}

# Set up SingleCellExperiment input
setInput.SCE <- function() {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  counts <- makeCounts()
  metadata <- makeMetaData()
  object <- suppressWarnings(SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts),
                                                                        colData = S4Vectors::DataFrame(metadata)))
  list(object = object,
       metadata = NULL,
       replicate_labels = "replicate",
       group_labels = "group",
       split_labels = NULL,
       use_assay = "counts",
       use_layer = NULL)
}

# All input setups
runDE_input_list <- list(cell_matrix = list(fun = setInput.CellMatrix,
                                            is_BPCells_input = FALSE,
                                            pseudobulk_values = c("generate", "none")),
                         Seurat_dgCMatrix = list(fun = setInput.Seurat,
                                                 is_BPCells_input = FALSE,
                                                 pseudobulk_values = c("generate", "none")),
                         Seurat_BPCells = list(fun = setInput.BPCells,
                                               is_BPCells_input = TRUE,
                                               pseudobulk_values = c("generate", "none")),
                         SCE = list(fun = setInput.SCE,
                                    is_BPCells_input = FALSE,
                                    pseudobulk_values = c("generate", "none")),
                         PB_matrix = list(fun = setInput.PBMatrix,
                                          is_BPCells_input = FALSE,
                                          pseudobulk_values = "supplied"))

# ---------------------------------------------------------------------------
# More test helpers
# ---------------------------------------------------------------------------

# Check whether this combination should error due to BPCells validation
#
# input_info  -- Entry from runDE_input_list
# de_method   -- DE method
# pseudobulk  -- Pseudobulk setting
expectBPCellsError <- function(input_info,
                               de_method,
                               pseudobulk) {
  # BPCells DE should only be used on existing IterableMatrix input.
  if (de_method == "BPCells") {
    return(!(input_info$is_BPCells_input && pseudobulk %in% c("supplied", "none")))
  }
  # Existing IterableMatrix input with supplied/none must use BPCells.
  if (input_info$is_BPCells_input && pseudobulk %in% c("supplied", "none") && de_method != "BPCells") {
    return(TRUE)
  }
  return(FALSE)
}

# Get test grid for input
#
# input_info -- Entry from runDE_input_list
getInputGrid <- function(input_info) {
  grid <- do.call(rbind, lapply(input_info$pseudobulk_values, function(p) {
    transform(runDE_method_grid, pseudobulk = p)
  }))
  grid$expect_error <- vapply(seq_len(nrow(grid)),
                              FUN = function(i) {
                                expectBPCellsError(input_info = input_info,
                                                   de_method = grid$de_method[[i]],
                                                   pseudobulk = grid$pseudobulk[[i]])
                              },
                              FUN.VALUE = logical(1))
  return(grid)
}

# Run one test case
#
# input      -- Input to runDE
# de_method  -- DE method
# de_test    -- DE test
# pseudobulk -- Pseudobulk setting
runTestCase <- function(input,
                        de_method,
                        de_test,
                        pseudobulk) {
  skipUninstalled(de_method)

  replicate_labels <- input$replicate_labels
  if (identical(pseudobulk, "none")) {
    replicate_labels <- NULL
  }

  warnings <- character()

  output <- withCallingHandlers(
    runDE(object = input$object,
          metadata = input$metadata,
          replicate_labels = replicate_labels,
          group_labels = input$group_labels,
          split_labels = input$split_labels,
          pseudobulk = pseudobulk,
          de_method = de_method,
          de_test = de_test,
          min_cells_per_split = 1,
          min_cells_per_replicate = 1,
          min_replicates_per_split = 1,
          min_replicates_per_group = 1,
          min_cells_per_feature = 1,
          min_prop_cells_per_feature = 0,
          use_assay = input$use_assay,
          use_layer = input$use_layer,
          n_cores = 2,
          verbose = TRUE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (identical(pseudobulk, "none")) {
    expect_true(any(grepl("Cell-level tests are not recommended", warnings)))
    expect_true(any(grepl("min_cells_per_replicate.*not used", warnings)))
    expect_true(any(grepl("min_replicates_per_split.*not used", warnings)))
  } else {
    expect_false(any(grepl("Cell-level tests are not recommended", warnings)))
  }

  expectValidOutput(output)
  expect_equal(output$parameters$de_method, de_method)
  expect_equal(output$parameters$de_test, de_test)
  expect_equal(output$parameters$pseudobulk, pseudobulk)

  return(output)
}

# Standardize DE_results before comparison
#
# de_results -- DE_results dataframe from runDE
standardizeDEResults <- function(de_results) {
  de_results |>
    dplyr::select(feature, split, lfc, pvalue, padj) |>
    dplyr::arrange(split, feature) |>
    data.frame()
}

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

# Set up counter
total_tests <- sum(vapply(
  runDE_input_list,
  FUN = function(x) nrow(getInputGrid(x)),
  FUN.VALUE = numeric(1)
))

test_counter <- 0L

# Store DE results for cross-input comparison
de_results_by_case <- new.env(parent = emptyenv())

# Loop through all input, DE method, DE test, and pseudobulk combinations
for (input_name in names(runDE_input_list)) {
  input_grid <- getInputGrid(runDE_input_list[[input_name]])
  for (i in seq_len(nrow(input_grid))) {
    local({
      test_counter <<- test_counter + 1L
      test_index_i <- test_counter
      total_tests_i <- total_tests
      input_name_i <- input_name
      input_info_i <- runDE_input_list[[input_name]]
      input_i <- input_info_i$fun
      de_method_i <- input_grid$de_method[[i]]
      de_test_i <- input_grid$de_test[[i]]
      pseudobulk_i <- input_grid$pseudobulk[[i]]
      expect_error_i <- input_grid$expect_error[[i]]
      test_that(paste("runDE works for",
                      input_name_i,
                      "with",
                      de_method_i,
                      de_test_i,
                      "and pseudobulk",
                      pseudobulk_i),
                {
                  message("[", test_index_i, "/", total_tests_i, "] ",
                          "runDE test: input = ", input_name_i,
                          ", method = ", de_method_i,
                          ", test = ", de_test_i,
                          ", pseudobulk = ", pseudobulk_i,
                          ", expect_error = ", expect_error_i)
                  input <- input_i()
                  if (expect_error_i) {
                    expect_error(runTestCase(input = input,
                                             de_method = de_method_i,
                                             de_test = de_test_i,
                                             pseudobulk = pseudobulk_i),
                                 "BPCells|IterableMatrix|de_method")
                  } else {
                    output <- runTestCase(input = input,
                                          de_method = de_method_i,
                                          de_test = de_test_i,
                                          pseudobulk = pseudobulk_i)
                    case_key <- paste(de_method_i, de_test_i, pseudobulk_i, sep = "__")
                    if (is.null(de_results_by_case[[case_key]])) {
                      de_results_by_case[[case_key]] <- list()
                    }
                    case_key <- paste(de_method_i, de_test_i, pseudobulk_i, sep = "__")
                    case_results <- de_results_by_case[[case_key]]
                    if (is.null(case_results)) {
                      case_results <- list()
                    }
                    case_results[[input_name_i]] <- standardizeDEResults(output$DE_results)
                    de_results_by_case[[case_key]] <- case_results
                    expect_equal(output$parameters$de_method, de_method_i)
                    expect_equal(output$parameters$de_test, de_test_i)
                    expect_equal(output$parameters$pseudobulk, pseudobulk_i)
                    expect_gt(nrow(output$DE_results), 0)
                  }
                })
    })
  }
}

# Compare DE results across input types for matching method/test/pseudobulk cases
message("Stored cross-input DE result counts:")

for (case_key in ls(de_results_by_case)) {
  message(
    case_key,
    ": ",
    paste(names(de_results_by_case[[case_key]]), collapse = ", ")
  )
}

test_that("runDE DE_results match across comparable input types", {
  testthat::skip_if(length(ls(de_results_by_case)) == 0)
  for (case_key in ls(de_results_by_case)) {
    results_i <- de_results_by_case[[case_key]]
    if (length(results_i) < 2) {
      next
    }
    reference_name <- names(results_i)[[1]]
    reference_results <- results_i[[reference_name]]
    for (input_name in names(results_i)[-1]) {
      expect_equal(results_i[[input_name]],
                   reference_results,
                   tolerance = 1e-8,
                   info = paste("Case:", case_key,
                                "| Compared", input_name, "to", reference_name))
    }
  }
})
