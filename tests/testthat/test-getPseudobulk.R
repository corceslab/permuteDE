# ---------------------------------------------------------------------------
# Tests for function getPseudobulk
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------

# Make a test count matrix
makeCounts <- function() {
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

# Make test replicate labels
makeReplicates <- function() {
  c("1", "1", "2", "2", "3", "3", "4", "4")
}

# Make test split labels
makeSplits <- function() {
  c("split1", "split1", "split1", "split1", "split2", "split2", "split2", "split2")
}

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

# Matrix input

test_that("getPseudobulk generates one pseudobulk matrix from matrix input", {
  counts <- makeCounts()
  output <- getPseudobulk(object = counts,
                          replicate_labels = makeReplicates(),
                          min_cells_per_split = 1,
                          min_cells_per_replicate = 1,
                          min_replicates_per_split = 1,
                          min_cells_per_feature = 1,
                          min_prop_cells_per_feature = 0,
                          n_cores = 1)
  expect_named(output, c("PB_values", "metadata", "parameters"))
  expect_type(output$PB_values, "list")
  expect_named(output$PB_values, "all")
  expect_true(is.matrix(output$PB_values$all) || inherits(output$PB_values$all, "Matrix"))
  expect_equal(colnames(output$PB_values$all), paste0("rep_", 1:4))
  expected <- cbind(rep_1 = Matrix::rowSums(counts[, c("cell1", "cell2"), drop = FALSE]),
                    rep_2 = Matrix::rowSums(counts[, c("cell3", "cell4"), drop = FALSE]),
                    rep_3 = Matrix::rowSums(counts[, c("cell5", "cell6"), drop = FALSE]),
                    rep_4 = Matrix::rowSums(counts[, c("cell7", "cell8"), drop = FALSE]))
  expected <- expected[Matrix::rowSums(expected) > 0, , drop = FALSE]
  expect_equal(as.matrix(output$PB_values$all), as.matrix(expected))
})

test_that("getPseudobulk generates one pseudobulk matrix per split with proper output structure", {
  counts <- makeCounts()
  output <- getPseudobulk(object = counts,
                          replicate_labels = makeReplicates(),
                          split_labels = makeSplits(),
                          min_cells_per_split = 1,
                          min_cells_per_replicate = 1,
                          min_replicates_per_split = 1,
                          min_cells_per_feature = 1,
                          min_prop_cells_per_feature = 0,
                          n_cores = 1)

  # Check that PB values has split structure
  expect_named(output$PB_values, c("split1", "split2"))
  expect_equal(colnames(output$PB_values$split1), c("rep_1", "rep_2"))
  expect_equal(colnames(output$PB_values$split2), c("rep_3", "rep_4"))
  expect_equal(as.matrix(output$PB_values$split1),
               as.matrix(cbind(rep_1 = Matrix::rowSums(counts[, c("cell1", "cell2"), drop = FALSE]),
                               rep_2 = Matrix::rowSums(counts[, c("cell3", "cell4"), drop = FALSE]))[c("gene1", "gene2", "gene3"), , drop = FALSE]))
  expect_equal(as.matrix(output$PB_values$split2),
               as.matrix(cbind(rep_3 = Matrix::rowSums(counts[, c("cell5", "cell6"), drop = FALSE]),
                               rep_4 = Matrix::rowSums(counts[, c("cell7", "cell8"), drop = FALSE]))[c("gene1", "gene2", "gene3"), , drop = FALSE]))

  # Check metadata structure
  expect_named(output$metadata, c("metrics", "exclude_features"))
  expect_true(all(c("split",
                    "n_all_features",
                    "n_nonzero_features",
                    "n_features_exclude",
                    "n_features_for_DE",
                    "prop_features_exclude",
                    "n_all_reads",
                    "n_reads_exclude",
                    "n_reads_for_DE",
                    "prop_reads_exclude") %in% colnames(output$metadata$metrics)))
  expect_equal(output$metadata$metrics$split, c("split1", "split2"))
  expect_type(output$metadata$exclude_features, "list")

  # Check that features with zero counts are removed
  expect_false("gene4" %in% rownames(output$PB_values$all))
  expect_equal(output$metadata$metrics$n_all_features, c(4, 4))
  expect_equal(output$metadata$metrics$n_nonzero_features, c(3, 3))
})

test_that("getPseudobulk filters features using min_cells_per_feature", {
  counts <- makeCounts()
  output <- getPseudobulk(object = counts,
                          replicate_labels = makeReplicates(),
                          min_cells_per_split = 1,
                          min_cells_per_replicate = 1,
                          min_replicates_per_split = 1,
                          min_cells_per_feature = 8,
                          min_prop_cells_per_feature = 0,
                          n_cores = 1)
  expect_equal(rownames(output$PB_values$all), "gene3")
  expect_equal(output$metadata$exclude_features$all, c("gene1", "gene2"))
})

test_that("getPseudobulk reports excluded features without removing them when filter is FALSE", {
  counts <- makeCounts()
  output <- getPseudobulk(object = counts,
                          replicate_labels = makeReplicates(),
                          min_cells_per_split = 1,
                          min_cells_per_replicate = 1,
                          min_replicates_per_split = 1,
                          min_cells_per_feature = 8,
                          min_prop_cells_per_feature = 0,
                          filter = FALSE,
                          n_cores = 1)
  expect_equal(rownames(output$PB_values$all), c("gene1", "gene2", "gene3"))
  expect_equal(output$metadata$exclude_features$all, c("gene1", "gene2"))
})

test_that("getPseudobulk with pseudobulk none returns cell-level matrices", {
  counts <- makeCounts()
  expect_warning(
    expect_warning(
      expect_warning(
        output <- getPseudobulk(object = counts,
                                split_labels = makeSplits(),
                                min_cells_per_split = 1,
                                min_cells_per_replicate = 1,
                                min_replicates_per_split = 1,
                                min_cells_per_feature = 1,
                                min_prop_cells_per_feature = 0,
                                pseudobulk = "none",
                                n_cores = 1),
        "Cell-level tests are not recommended"),
      "min_cells_per_replicate.*not used"),
    "min_replicates_per_split.*not used")
  expect_named(output, c("cell_values", "metadata", "parameters"))
  expect_named(output$cell_values, c("split1", "split2"))
  expect_equal(colnames(output$cell_values$split1), paste0("cell", 1:4))
  expect_equal(colnames(output$cell_values$split2), paste0("cell", 5:8))
})

test_that("getPseudobulk with pseudobulk none reports skipped splits correctly", {
  counts <- makeCounts()
  split_labels <- c("split1", "split1", "split1", "split1", "split1", "split1",
                    "split2", "split2")
  expect_warning(
    expect_warning(
      expect_warning(
        expect_message(
          output <- getPseudobulk(object = counts,
                                  replicate_labels = NULL,
                                  split_labels = split_labels,
                                  min_cells_per_split = 5,
                                  min_cells_per_replicate = 1,
                                  min_replicates_per_split = 1,
                                  min_cells_per_feature = 1,
                                  min_prop_cells_per_feature = 0,
                                  filter = TRUE,
                                  pseudobulk = "none",
                                  n_cores = 1,
                                  verbose = TRUE),
          "Skipped 1 split label due to insufficient cells: split2"),
        "Cell-level tests are not recommended"),
      "min_cells_per_replicate.*not used"),
    "min_replicates_per_split.*not used")
  expect_named(output$cell_values, "split1")
  expect_false("split2" %in% names(output$cell_values))
  expect_equal(colnames(output$cell_values$split1), paste0("cell", 1:6))
})

test_that("getPseudobulk respects use_cells", {
  counts <- makeCounts()
  output <- getPseudobulk(object = counts,
                          replicate_labels = makeReplicates()[1:4],
                          use_cells = paste0("cell", 1:4),
                          min_cells_per_split = 1,
                          min_cells_per_replicate = 1,
                          min_replicates_per_split = 1,
                          min_cells_per_feature = 1,
                          min_prop_cells_per_feature = 0,
                          n_cores = 1)
  expect_equal(colnames(output$PB_values$all), c("rep_1", "rep_2"))
  expect_equal(output$parameters$use_cells, paste0("cell", 1:4))
})

test_that("getPseudobulk errors when label vectors have the wrong length", {
  counts <- makeCounts()
  expect_error(getPseudobulk(object = counts,
                             replicate_labels = c("rep1", "rep2"),
                             min_cells_per_split = 1,
                             min_cells_per_replicate = 1,
                             min_replicates_per_split = 1,
                             n_cores = 1),
               "same length and in the same order as the supplied cells")
  expect_error(getPseudobulk(object = counts,
                             replicate_labels = makeReplicates(),
                             split_labels = c("split1", "split2"),
                             min_cells_per_split = 1,
                             min_cells_per_replicate = 1,
                             min_replicates_per_split = 1,
                             n_cores = 1),
               "same length and in the same order as the supplied cells")
})

test_that("getPseudobulk errors when labels contain NA", {
  counts <- makeCounts()
  replicates <- makeReplicates()
  replicates[1] <- NA
  expect_error(getPseudobulk(object = counts,
                             replicate_labels = replicates,
                             min_cells_per_split = 1,
                             min_cells_per_replicate = 1,
                             min_replicates_per_split = 1,
                             n_cores = 1),
               "replicate_labels.*cannot be NA")
  replicates <- makeReplicates()
  splits <- makeSplits()
  splits[1] <- NA
  expect_error(getPseudobulk(object = counts,
                             replicate_labels = makeReplicates(),
                             split_labels = splits,
                             min_cells_per_split = 1,
                             min_cells_per_replicate = 1,
                             min_replicates_per_split = 1,
                             n_cores = 1),
               "split_labels.*cannot be NA")
})

test_that("getPseudobulk removes replicate-split pairs with too few cells", {
  counts <- makeCounts()
  output <- getPseudobulk(object = counts,
                          replicate_labels = makeReplicates(),
                          split_labels = makeSplits(),
                          min_cells_per_split = 1,
                          min_cells_per_replicate = 3,
                          min_replicates_per_split = 1,
                          min_cells_per_feature = 1,
                          min_prop_cells_per_feature = 0,
                          n_cores = 1)
  expect_null(output$PB_values)
  expect_null(output$metadata)
})

test_that("getPseudobulk retrieves replicate and split labels from metadata", {
  counts <- makeCounts()
  metadata <- data.frame(replicate = makeReplicates(),
                         split = makeSplits(),
                         row.names = colnames(counts))
  output <- getPseudobulk(object = counts,
                          metadata = metadata,
                          replicate_labels = "replicate",
                          split_labels = "split",
                          min_cells_per_split = 1,
                          min_cells_per_replicate = 1,
                          min_replicates_per_split = 1,
                          min_cells_per_feature = 1,
                          min_prop_cells_per_feature = 0,
                          n_cores = 1)
  expect_named(output$PB_values, c("split1", "split2"))
  expect_equal(colnames(output$PB_values$split1), c("rep_1", "rep_2"))
  expect_equal(colnames(output$PB_values$split2), c("rep_3", "rep_4"))
})

# SingleCellExperiment input
test_that("getPseudobulk works with SingleCellExperiment input", {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("S4Vectors")
  counts <- makeCounts()
  col_data <- S4Vectors::DataFrame(replicate = makeReplicates(),
                                   split = makeSplits(),
                                   row.names = colnames(counts))
  object <- suppressWarnings(SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts),
                                                                        colData = col_data))
  output <- getPseudobulk(object = object,
                          replicate_labels = "replicate",
                          split_labels = "split",
                          min_cells_per_split = 1,
                          min_cells_per_replicate = 1,
                          min_replicates_per_split = 1,
                          min_cells_per_feature = 1,
                          min_prop_cells_per_feature = 0,
                          n_cores = 1)
  expect_named(output, c("PB_values", "metadata", "parameters"))
  expect_named(output$PB_values, c("split1", "split2"))
  expect_equal(output$parameters$object_type, "SingleCellExperiment")
  expect_equal(colnames(output$PB_values$split1), c("rep_1", "rep_2"))
  expect_equal(colnames(output$PB_values$split2), c("rep_3", "rep_4"))
  expected_split1 <- cbind(rep_1 = Matrix::rowSums(counts[, c("cell1", "cell2"), drop = FALSE]),
                           rep_2 = Matrix::rowSums(counts[, c("cell3", "cell4"), drop = FALSE]))
  expected_split2 <- cbind(rep_3 = Matrix::rowSums(counts[, c("cell5", "cell6"), drop = FALSE]),
                           rep_4 = Matrix::rowSums(counts[, c("cell7", "cell8"), drop = FALSE]))
  expected_split1 <- expected_split1[Matrix::rowSums(expected_split1) > 0, , drop = FALSE]
  expected_split2 <- expected_split2[Matrix::rowSums(expected_split2) > 0, , drop = FALSE]
  expect_equal(as.matrix(output$PB_values$split1), as.matrix(expected_split1))
  expect_equal(as.matrix(output$PB_values$split2), as.matrix(expected_split2))
})

# Seurat input
test_that("getPseudobulk works with Seurat input", {
  testthat::skip_if_not_installed("Seurat")
  counts <- makeCounts()
  object <- suppressWarnings(Seurat::CreateSeuratObject(counts = counts,
                                                        meta.data = data.frame(
                                                          replicate = makeReplicates(),
                                                          split = makeSplits(),
                                                          row.names = colnames(counts))))
  output <- getPseudobulk(object = object,
                          replicate_labels = "replicate",
                          split_labels = "split",
                          min_cells_per_split = 1,
                          min_cells_per_replicate = 1,
                          min_replicates_per_split = 1,
                          min_cells_per_feature = 1,
                          min_prop_cells_per_feature = 0,
                          n_cores = 1)
  expect_named(output, c("PB_values", "metadata", "parameters"))
  expect_named(output$PB_values, c("split1", "split2"))
  expect_equal(output$parameters$object_type, "Seurat")
  expect_equal(colnames(output$PB_values$split1), c("rep_1", "rep_2"))
  expect_equal(colnames(output$PB_values$split2), c("rep_3", "rep_4"))
  expected_split1 <- cbind(rep_1 = Matrix::rowSums(counts[, c("cell1", "cell2"), drop = FALSE]),
                           rep_2 = Matrix::rowSums(counts[, c("cell3", "cell4"), drop = FALSE]))
  expected_split2 <- cbind(rep_3 = Matrix::rowSums(counts[, c("cell5", "cell6"), drop = FALSE]),
                           rep_4 = Matrix::rowSums(counts[, c("cell7", "cell8"), drop = FALSE]))
  expected_split1 <- expected_split1[Matrix::rowSums(expected_split1) > 0, , drop = FALSE]
  expected_split2 <- expected_split2[Matrix::rowSums(expected_split2) > 0, , drop = FALSE]
  expect_equal(as.matrix(output$PB_values$split1), as.matrix(expected_split1))
  expect_equal(as.matrix(output$PB_values$split2), as.matrix(expected_split2))
})


