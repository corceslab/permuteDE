# ---------------------------------------------------------------------------
# Shared test helpers for testthat scripts
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Small deterministic data for getPseudobulk arithmetic tests
# ---------------------------------------------------------------------------

# Make a small deterministic test count matrix ---------------------------
makeSmallCounts <- function() {
  matrix(
    c(
      # c1 c2 c3 c4 c5 c6 c7 c8
      1,  2,  0,  0,  5,  5,  0,  0,  # gene1
      0,  0,  3,  4,  0,  0,  6,  6,  # gene2
      1,  1,  1,  1,  1,  1,  1,  1,  # gene3
      0,  0,  0,  0,  0,  0,  0,  0   # gene4
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("gene", 1:4),
                    paste0("cell", 1:8)))
}

# Make small deterministic replicate labels ---------------------------
makeSmallReplicates <- function() {
  c("1", "1", "2", "2", "3", "3", "4", "4")
}

# Make small deterministic split labels ---------------------------
makeSmallSplits <- function() {
  c("split1", "split1", "split1", "split1",
    "split2", "split2", "split2", "split2")
}

# ---------------------------------------------------------------------------
# Larger data for runDE / permuteDE tests
# ---------------------------------------------------------------------------

# Make a larger test count matrix ---------------------------
#
# n_genes     -- Number of genes
# random_seed -- Random seed
makeCounts <- function(n_genes = 100,
                       random_seed = 1) {
  set.seed(random_seed)

  n_reps <- 12
  cells_per_rep <- 2
  n_cells <- n_reps * cells_per_rep

  group_by_rep <- rep(c("A", "B"), each = 6)
  group <- rep(group_by_rep, each = cells_per_rep)

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

  zero_genes <- Matrix::rowSums(mat) == 0
  if (any(zero_genes)) {
    mat[zero_genes, ] <- 1L
  }

  return(mat)
}

# Make a test cell metadata dataframe ---------------------------
makeMetaData <- function() {
  data.frame(
    replicate = rep(paste0("rep", 1:12),
                    each = 2),
    group = rep(c(rep("A", 6),
                  rep("B", 6)),
                each = 2),
    split = rep("all", 24),
    split_multi = rep(c(rep("split1", 3),
                        rep("split2", 3),
                        rep("split1", 3),
                        rep("split2", 3)),
                      each = 2),
    batch = rep(c("batch1", "batch2", "batch1",
                  "batch1", "batch2", "batch1",
                  "batch1", "batch2", "batch1",
                  "batch1", "batch2", "batch1"),
                each = 2),
    age = rep(c(50, 55, 60,
                52, 57, 62,
                51, 56, 61,
                53, 58, 63),
              each = 2),
    row.names = paste0("cell", 1:24))
}

# Make test pseudobulk matrix ---------------------------
makePB <- function() {
  counts <- makeCounts()
  metadata <- makeMetaData()
  reps <- unique(metadata$replicate)

  pb <- vapply(reps,
               FUN = function(r) {
                 Matrix::rowSums(counts[, metadata$replicate == r, drop = FALSE])
               },
               FUN.VALUE = numeric(nrow(counts)))

  rownames(pb) <- rownames(counts)
  colnames(pb) <- reps

  return(pb)
}

# ---------------------------------------------------------------------------
# Output expectations
# ---------------------------------------------------------------------------

# Get expected output value name ---------------------------
#
# pseudobulk -- Pseudobulk setting
expectedValueName <- function(pseudobulk) {
  if (identical(pseudobulk, "none")) {
    return("cell_values")
  }
  return("PB_values")
}

# Make set of expectations for runDE ---------------------------
#
# output -- runDE output
expectValidOutput.runDE <- function(output) {
  # Check names
  expected_name <- expectedValueName(output$parameters$pseudobulk)

  expected_names <- c("DE_results", expected_name, "metadata", "parameters")
  if (isTRUE(output$parameters$return_raw_de)) {
    expected_names <- c("DE_results", "raw_DE_results", expected_name, "metadata", "parameters")
  }

  expect_named(output, expected_names)

  # Check DE_results structure
  expect_s3_class(output$DE_results, "data.frame")
  expect_true(all(c("feature", "lfc", "pvalue", "padj", "split") %in%
                    colnames(output$DE_results)))
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

  # Check raw_DE_results
  if (isTRUE(output$parameters$return_raw_de)) {
    expect_type(output$raw_DE_results, "list")
    expect_gt(length(output$raw_DE_results), 0)
    expect_setequal(names(output$raw_DE_results), names(output[[expected_name]]))
  }

  # Check that split labels match output matrix/cell list names
  expect_setequal(unique(output$DE_results$split), names(output[[expected_name]]))

  # Check adjusted p-values are consistent with p-values within each split
  for (split_i in unique(output$DE_results$split)) {
    rows_i <- output$DE_results$split == split_i

    if (identical(output$parameters$p_adjust_method, "fdrtool")) {
      expect_true(all(is.na(output$DE_results$padj[rows_i]) |
          (output$DE_results$padj[rows_i] >= 0 &
              output$DE_results$padj[rows_i] <= 1)))
    } else {
      expect_equal(output$DE_results$padj[rows_i],
                   stats::p.adjust(
                     output$DE_results$pvalue[rows_i],
                     method = output$parameters$p_adjust_method),
                   tolerance = 1e-12)
    }

  }

  # Check output matrices/cell values
  expect_type(output[[expected_name]], "list")
  expect_gt(length(output[[expected_name]]), 0)

  if (identical(expected_name, "PB_values")) {
    expect_false("cell_values" %in% names(output))
  } else {
    expect_false("PB_values" %in% names(output))
  }

  # Check metadata names
  expect_true("group_key" %in% names(output$metadata))
  expect_true("replicates_by_split" %in% names(output$metadata))
  expect_true("n_replicates_by_split" %in% names(output$metadata))
  expect_true("time" %in% names(output$metadata))

  # Check group_key
  expect_s3_class(output$metadata$group_key, "data.frame")
  expect_true(all(c("replicate", "group") %in%
                    colnames(output$metadata$group_key)))
  expect_type(output$metadata$group_key$replicate, "character")
  expect_false(anyNA(output$metadata$group_key$replicate))

  # Check replicates_by_split
  expect_s3_class(output$metadata$replicates_by_split, "data.frame")
  expect_true(all(c("split", "replicate", "group") %in%
                    colnames(output$metadata$replicates_by_split)))

  expect_type(output$metadata$replicates_by_split$split, "character")
  expect_type(output$metadata$replicates_by_split$replicate, "character")
  expect_false(anyNA(output$metadata$replicates_by_split$split))
  expect_false(anyNA(output$metadata$replicates_by_split$replicate))

  expect_setequal(output$metadata$replicates_by_split$split,
                  names(output[[expected_name]]))

  for (split_i in names(output[[expected_name]])) {
    rows_i <- output$metadata$replicates_by_split$split == split_i
    expect_setequal(output$metadata$replicates_by_split$replicate[rows_i],
                    colnames(output[[expected_name]][[split_i]]))
  }

  # Check n_replicates_by_split
  expect_s3_class(output$metadata$n_replicates_by_split, "data.frame")
  expect_true(all(c("split", "n_reference_group", "n_non_reference_group") %in%
                    colnames(output$metadata$n_replicates_by_split)))

  expect_type(output$metadata$n_replicates_by_split$split, "character")
  expect_true(is.numeric(output$metadata$n_replicates_by_split$n_reference_group))
  expect_true(is.numeric(output$metadata$n_replicates_by_split$n_non_reference_group))

  expect_false(anyNA(output$metadata$n_replicates_by_split$split))
  expect_false(anyNA(output$metadata$n_replicates_by_split$n_reference_group))
  expect_false(anyNA(output$metadata$n_replicates_by_split$n_non_reference_group))

  expect_setequal(output$metadata$n_replicates_by_split$split,
                  names(output[[expected_name]]))

  expect_true(all(output$metadata$n_replicates_by_split$n_reference_group >= 1))
  expect_true(all(output$metadata$n_replicates_by_split$n_non_reference_group >= 1))

  for (split_i in names(output[[expected_name]])) {
    replicate_rows_i <- output$metadata$replicates_by_split$split == split_i
    count_rows_i <- output$metadata$n_replicates_by_split$split == split_i

    expect_equal(sum(count_rows_i), 1L)

    split_groups_i <- output$metadata$replicates_by_split$group[replicate_rows_i]

    expect_equal(output$metadata$n_replicates_by_split$n_reference_group[count_rows_i],
                 sum(split_groups_i == output$parameters$reference_group))

    expect_equal(output$metadata$n_replicates_by_split$n_non_reference_group[count_rows_i],
                 sum(split_groups_i == output$parameters$non_reference_group))
  }

  # Check time metadata
  expect_s3_class(output$metadata$time, "data.frame")
}

# Make set of expectations for permuteDE output ---------------------------
#
# output     -- permuteDE output
# return_all -- Whether permutation_DE_results should be present
expectValidOutput.permuteDE <- function(output,
                                        return_all = FALSE) {
  expected_names <- c("permutation_test_results",
                      "permutation_DE_summary",
                      "metadata",
                      "parameters")

  if (isTRUE(return_all)) {
    expected_names <- c("permutation_test_results",
                        "permutation_DE_summary",
                        "permutation_DE_results",
                        "metadata",
                        "parameters")
  }

  expect_named(output, expected_names)

  # Check permutation_test_results
  expect_s3_class(output$permutation_test_results, "data.frame")
  expect_true(all(c("split",
                    "runDE_n_sig",
                    "pvalue",
                    "n_iterations") %in%
                    colnames(output$permutation_test_results)))

  expect_type(output$permutation_test_results$split, "character")
  expect_true(is.numeric(output$permutation_test_results$runDE_n_sig))
  expect_true(is.numeric(output$permutation_test_results$pvalue))
  expect_true(is.numeric(output$permutation_test_results$n_iterations))

  expect_true(all(output$permutation_test_results$pvalue >= 0 |
                    is.na(output$permutation_test_results$pvalue)))
  expect_true(all(output$permutation_test_results$pvalue <= 1 |
                    is.na(output$permutation_test_results$pvalue)))
  expect_true(all(output$permutation_test_results$n_iterations >= 1))

  # Check permutation_DE_summary
  expect_s3_class(output$permutation_DE_summary, "data.frame")
  expect_true(all(c("split",
                    "permutation",
                    "reference_group_overlap",
                    "non_reference_group_overlap",
                    "n_sig",
                    "min_lfc_sig",
                    "max_lfc_sig",
                    "min_lfc_all",
                    "max_lfc_all") %in%
                    colnames(output$permutation_DE_summary)))

  expect_type(output$permutation_DE_summary$split, "character")
  expect_true(is.numeric(output$permutation_DE_summary$permutation))
  expect_true(is.numeric(output$permutation_DE_summary$reference_group_overlap))
  expect_true(is.numeric(output$permutation_DE_summary$non_reference_group_overlap))
  expect_true(is.numeric(output$permutation_DE_summary$n_sig))

  expect_true(all(output$permutation_DE_summary$reference_group_overlap >= 0 |
                    is.na(output$permutation_DE_summary$reference_group_overlap)))
  expect_true(all(output$permutation_DE_summary$reference_group_overlap <= 1 |
                    is.na(output$permutation_DE_summary$reference_group_overlap)))
  expect_true(all(output$permutation_DE_summary$non_reference_group_overlap >= 0 |
                    is.na(output$permutation_DE_summary$non_reference_group_overlap)))
  expect_true(all(output$permutation_DE_summary$non_reference_group_overlap <= 1 |
                    is.na(output$permutation_DE_summary$non_reference_group_overlap)))

  # Check metadata
  expect_true("runDE_values" %in% names(output$metadata))
  expect_true("permutation_reference_group_indices" %in% names(output$metadata))
  expect_true("time" %in% names(output$metadata))

  expect_s3_class(output$metadata$runDE_values, "data.frame")
  expect_true(all(c("split", "runDE_n_sig") %in%
                    colnames(output$metadata$runDE_values)))
  expect_type(output$metadata$permutation_reference_group_indices, "list")

  # Check parameters
  expect_true(all(c("alpha",
                    "lfc_threshold",
                    "n_iterations",
                    "use_splits",
                    "permute_by",
                    "permute_within",
                    "min_DE",
                    "reference_group",
                    "non_reference_group",
                    "de_method",
                    "de_test",
                    "de_params",
                    "p_adjust_method",
                    "pseudobulk",
                    "return_all",
                    "random_seed",
                    "n_cores") %in%
                    names(output$parameters)))

  expect_equal(output$parameters$return_all, return_all)

  # Check full permutation DE results if requested
  if (isTRUE(return_all)) {
    expect_s3_class(output$permutation_DE_results, "data.frame")
    expect_true(all(c("feature",
                      "lfc",
                      "pvalue",
                      "padj",
                      "permutation",
                      "split") %in%
                      colnames(output$permutation_DE_results)))

    expect_type(output$permutation_DE_results$feature, "character")
    expect_true(is.numeric(output$permutation_DE_results$lfc))
    expect_true(is.numeric(output$permutation_DE_results$pvalue))
    expect_true(is.numeric(output$permutation_DE_results$padj))
    expect_true(is.numeric(output$permutation_DE_results$permutation))
    expect_type(output$permutation_DE_results$split, "character")
  } else {
    expect_false("permutation_DE_results" %in% names(output))
  }
}

# ---------------------------------------------------------------------------
# DE method grids and skips
# ---------------------------------------------------------------------------

# DE method/test grid
runDE_method_grid <- data.frame(de_method = c("edgeR", "edgeR", "edgeR",
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

# Split grid
runDE_split_grid <- data.frame(split_condition = c("single_split", "multiple_splits"),
                               split_labels = c(NA_character_, "split_multi"),
                               stringsAsFactors = FALSE)

# Design grid
runDE_design_grid <- data.frame(design_condition = c("default",
                                                     "categorical",
                                                     "numeric",
                                                     "categorical_numeric"),
                                design = c(NA_character_,
                                           "~ batch + group",
                                           "~ age + group",
                                           "~ batch + age + group"),
                                stringsAsFactors = FALSE)

# Skip tests for uninstalled packages ---------------------------
#
# de_method -- DE method
skipUninstalled <- function(de_method) {
  if (de_method %in% c("edgeR", "limma", "presto")) {
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
# Input for runDE / permuteDE integration tests
# ---------------------------------------------------------------------------

# Set up cell-level matrix input ---------------------------
setInput.CellMatrix <- function() {
  list(object = makeCounts(),
       metadata = makeMetaData(),
       replicate_labels = "replicate",
       group_labels = "group",
       split_labels = NULL,
       use_assay = NULL,
       use_layer = NULL)
}

# Set up PB matrix input ---------------------------
setInput.PBMatrix <- function() {
  pb <- makePB()

  metadata <- data.frame(replicate = colnames(pb),
                         group = c(rep("A", 6), rep("B", 6)),
                         split = rep("all", 12),
                         split_multi = c(rep("split1", 3),
                                         rep("split2", 3),
                                         rep("split1", 3),
                                         rep("split2", 3)),
                         batch = c("batch1", "batch2", "batch1",
                                   "batch1", "batch2", "batch1",
                                   "batch1", "batch2", "batch1",
                                   "batch1", "batch2", "batch1"),
                         age = c(50, 55, 60,
                                 52, 57, 62,
                                 51, 56, 61,
                                 53, 58, 63),
                         row.names = colnames(pb))

  list(object = pb,
       metadata = metadata,
       replicate_labels = "replicate",
       group_labels = "group",
       split_labels = NULL,
       use_assay = NULL,
       use_layer = NULL)
}

# Set up Seurat input with dgCMatrix ---------------------------
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

# Set up Seurat input with BPCells ---------------------------
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

# Set up SingleCellExperiment input ---------------------------
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

# Check whether this combination should error due to BPCells validation ---------------------------
#
# input_info -- Entry from runDE_input_list
# de_method  -- DE method
# pseudobulk -- Pseudobulk setting
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

# Check whether this combination should error because the DE test does not
# support covariate-adjusted designs
#
# de_test          -- DE test
# design_condition -- Design condition from runDE_design_grid
expectDesignError <- function(de_test,
                              design_condition) {
  if (identical(design_condition, "default")) {
    return(FALSE)
  }

  de_test %in% c("exact", "wilcox_cpm", "wilcox_log_cpm")
}

# Get test grid for input ---------------------------
#
# input_info -- Entry from runDE_input_list
# Get test grid for input
#
# input_info -- Entry from runDE_input_list
getInputGrid <- function(input_info) {
  grid <- do.call(rbind, lapply(input_info$pseudobulk_values, function(p) {
    transform(runDE_method_grid, pseudobulk = p)
  }))

  grid <- merge(grid, runDE_split_grid, by = NULL)
  grid <- merge(grid, runDE_design_grid, by = NULL)

  grid$expect_error <- vapply(seq_len(nrow(grid)),
                              FUN = function(i) {
                                expectBPCellsError(input_info = input_info,
                                                   de_method = grid$de_method[[i]],
                                                   pseudobulk = grid$pseudobulk[[i]]) ||
                                  expectDesignError(de_test = grid$de_test[[i]],
                                                    design_condition = grid$design_condition[[i]])
                              },
                              FUN.VALUE = logical(1))

  return(grid)
}

# ---------------------------------------------------------------------------
# Test case helpers
# ---------------------------------------------------------------------------

# Run one runDE test case ---------------------------
#
# input                     -- Input to runDE
# split_labels              -- runDE parameter
# pseudobulk                -- runDE parameter
# de_method                 -- runDE parameter
# de_test                   -- runDE parameter
# de_params                 -- runDE parameter
# return_raw_de             -- runDE parameter
# min_replicates_per_group  -- runDE parameter
# n_cores                   -- runDE parameter
# verbose                   -- runDE parameter
testCase.runDE <- function(input,
                           split_labels = NULL,
                           design = NULL,
                           pseudobulk,
                           de_method,
                           de_test,
                           de_params = list(),
                           return_raw_de = FALSE,
                           min_replicates_per_group = 1,
                           n_cores = 1,
                           verbose = FALSE) {
  skipUninstalled(de_method)

  replicate_labels <- input$replicate_labels
  if (identical(pseudobulk, "none")) {
    replicate_labels <- NULL
  }
  if (!is.null(split_labels) &&
      identical(pseudobulk, "supplied") &&
      length(split_labels) == 1 &&
      !is.null(input$metadata) &&
      split_labels %in% colnames(input$metadata)) {
    split_labels <- input$metadata[[split_labels]]
  }

  warnings <- character()

  output <- withCallingHandlers(runDE(object = input$object,
                                      metadata = input$metadata,
                                      replicate_labels = replicate_labels,
                                      group_labels = input$group_labels,
                                      split_labels = split_labels,
                                      design = design,
                                      pseudobulk = pseudobulk,
                                      de_method = de_method,
                                      de_test = de_test,
                                      de_params = de_params,
                                      return_raw_de = return_raw_de,
                                      min_cells_per_split = 1,
                                      min_cells_per_replicate = 1,
                                      min_replicates_per_split = 1,
                                      min_replicates_per_group = min_replicates_per_group,
                                      min_cells_per_feature = 1,
                                      min_prop_cells_per_feature = 0,
                                      use_assay = input$use_assay,
                                      use_layer = input$use_layer,
                                      n_cores = n_cores,
                                      verbose = verbose),
                                warning = function(w) {
                                  warnings <<- c(warnings, conditionMessage(w))
                                  invokeRestart("muffleWarning")
                                })

  if (identical(pseudobulk, "none")) {
    expect_true(any(grepl("Cell-level tests are not recommended", warnings)))
    expect_true(any(grepl("min_cells_per_replicate.*not used", warnings)))
    expect_true(any(grepl("min_replicates_per_split.*not used", warnings)))
  } else {
    expect_false(any(grepl("Cell-level tests are not recommended", warnings)))
  }

  expectValidOutput.runDE(output)
  expect_equal(output$parameters$de_method, de_method)
  expect_equal(output$parameters$de_test, de_test)
  expect_equal(output$parameters$pseudobulk, pseudobulk)

  expect_s3_class(output$parameters$design_formula, "formula")
  expected_design_formula <- if (is.null(design)) {
    stats::as.formula("~ group")
  } else {
    stats::as.formula(sub(" [^ ]+$", " group", design))
  }
  expect_equal(paste(deparse(output$parameters$design_formula), collapse = ""),
               paste(deparse(expected_design_formula), collapse = ""))

  return(output)
}

# Run one permuteDE test case ---------------------------
#
# input        -- Input to runDE
# de_method    -- DE method
# de_test      -- DE test
# pseudobulk   -- Pseudobulk setting
# split_labels -- Split labels for runDE
# design       -- Design formula for runDE
# return_all   -- Whether permuteDE should return all per-feature DE results
testCase.permuteDE <- function(input,
                               de_method,
                               de_test,
                               pseudobulk,
                               split_labels = NULL,
                               design = NULL,
                               return_all = FALSE,
                               n_iterations = 5,
                               min_DE = 0,
                               n_cores = 1,
                               verbose = FALSE) {
  runDE_output <- testCase.runDE(input = input,
                                 de_method = de_method,
                                 de_test = de_test,
                                 pseudobulk = pseudobulk,
                                 split_labels = split_labels,
                                 design = design,
                                 n_cores = 1,
                                 verbose = verbose)

  warnings <- character()

  output <- withCallingHandlers(
    permuteDE(input = runDE_output,
              alpha = 0.05,
              lfc_threshold = 0,
              n_iterations = n_iterations,
              min_DE = min_DE,
              return_all = return_all,
              random_seed = 1,
              n_cores = n_cores,
              verbose = verbose),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (identical(pseudobulk, "none")) {
    expect_true(any(grepl("Cell-level tests are not recommended", warnings)))
    expect_true(any(grepl("consider providing biological replicates labels to parameter 'permute_by'", warnings)))
  } else {
    expect_false(any(grepl("Cell-level tests are not recommended", warnings)))
    expect_false(any(grepl("consider providing biological replicates labels to parameter 'permute_by'", warnings)))
  }

  expectValidOutput.permuteDE(output, return_all = return_all)

  expect_equal(output$parameters$de_method, de_method)
  expect_equal(output$parameters$de_test, de_test)
  expect_equal(output$parameters$pseudobulk, pseudobulk)
  expect_equal(output$parameters$return_all, return_all)
  expect_equal(output$parameters$n_iterations, n_iterations)
  expect_equal(output$parameters$min_DE, min_DE)

  return(output)
}

# ---------------------------------------------------------------------------
# Inter-input type comparison helpers
# ---------------------------------------------------------------------------

# Standardize DE_results before comparison ---------------------------
#
# de_results -- DE_results dataframe from runDE
standardizeDEResults <- function(de_results) {
  de_results |>
    dplyr::select(feature, split, lfc, pvalue, padj) |>
    dplyr::arrange(split, feature) |>
    data.frame()
}

# Standardize permuteDE results before comparison ---------------------------
#
# output -- permuteDE output
standardizePermutationResults <- function(output) {
  list(permutation_test_results = output$permutation_test_results |>
         dplyr::select(split, runDE_n_sig, pvalue, n_iterations) |>
         dplyr::arrange(split) |>
         data.frame(),
       permutation_DE_summary = output$permutation_DE_summary |>
         dplyr::select(split,
                       permutation,
                       reference_group_overlap,
                       non_reference_group_overlap,
                       n_sig,
                       min_lfc_sig,
                       max_lfc_sig,
                       min_lfc_all,
                       max_lfc_all) |>
         dplyr::arrange(split, permutation) |>
         data.frame())
}

# ---------------------------------------------------------------------------
# Other helpers
# ---------------------------------------------------------------------------

# Capture warnings ---------------------------
#
# expr -- expression
capture_warnings <- function(expr) {
  warnings <- character()

  value <- withCallingHandlers(expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = value, warnings = warnings)
}
