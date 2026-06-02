# ---------------------------------------------------------------------------
# Tests for plotting functions
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# permuteDEtheme
# ---------------------------------------------------------------------------

test_that("permuteDEtheme returns a ggplot theme", {
  theme <- permuteDEtheme()

  expect_s3_class(theme, "theme")
  expect_s3_class(theme, "gg")
})

# ---------------------------------------------------------------------------
# permuteDEpalette
# ---------------------------------------------------------------------------

test_that("permuteDEpalette returns discrete palettes", {
  pal <- permuteDEpalette(type = "discrete", n = 5, palette = "choir")

  expect_type(pal, "character")
  expect_length(pal, 5)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal)))
})

test_that("permuteDEpalette returns gradient palettes", {
  pal <- permuteDEpalette(type = "gradient", n = 5, palette = "inferno")

  expect_type(pal, "character")
  expect_length(pal, 5)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal)))
})

test_that("permuteDEpalette accepts custom palettes", {
  pal <- permuteDEpalette(type = "gradient",
                          n = 4,
                          palette = c("black", "white"))

  expect_type(pal, "character")
  expect_length(pal, 4)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", pal)))
})

# ---------------------------------------------------------------------------
# plotVolcano
# ---------------------------------------------------------------------------

test_that("plotVolcano returns a ggplot for one split", {
  input <- setInput.CellMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 pseudobulk = "generate",
                                 split_labels = "split_multi",
                                 design = NULL,
                                 n_cores = 1,
                                 verbose = FALSE)

  p <- plotVolcano(input = runDE_output,
                   use_splits = "split1",
                   n_max_label = 3)

  expect_s3_class(p, "ggplot")
})

test_that("plotVolcano returns a named list for multiple splits", {
  input <- setInput.CellMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 pseudobulk = "generate",
                                 split_labels = "split_multi",
                                 design = NULL,
                                 n_cores = 1,
                                 verbose = FALSE)

  p <- plotVolcano(input = runDE_output,
                   use_splits = c("split1", "split2"),
                   n_max_label = 3)

  expect_type(p, "list")
  expect_named(p, c("split1", "split2"))
  expect_true(all(vapply(p, inherits, logical(1), what = "ggplot")))
})

test_that("plotVolcano works with explicit labels and uncentered axis", {
  input <- setInput.CellMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 pseudobulk = "generate",
                                 split_labels = "split_multi",
                                 design = NULL,
                                 n_cores = 1,
                                 verbose = FALSE)

  label_features <- head(runDE_output$DE_results$feature, 2)

  p <- plotVolcano(input = runDE_output,
                   use_splits = "split1",
                   label_features = label_features,
                   center = FALSE)

  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------------
# plotHistogram
# ---------------------------------------------------------------------------

test_that("plotHistogram returns a ggplot for one split", {
  input <- setInput.CellMatrix()

  permute_output <- testCase.permuteDE(input = input,
                                       de_method = "edgeR",
                                       de_test = "LRT",
                                       pseudobulk = "generate",
                                       split_labels = "split_multi",
                                       design = NULL,
                                       return_all = FALSE,
                                       n_iterations = 5,
                                       min_DE = 0,
                                       n_cores = 1,
                                       verbose = FALSE)

  p <- plotHistogram(input = permute_output,
                     use_splits = "split1")

  expect_s3_class(p, "ggplot")
})

test_that("plotHistogram returns a named list for multiple splits", {
  input <- setInput.CellMatrix()

  permute_output <- testCase.permuteDE(input = input,
                                       de_method = "edgeR",
                                       de_test = "LRT",
                                       pseudobulk = "generate",
                                       split_labels = "split_multi",
                                       design = NULL,
                                       return_all = FALSE,
                                       n_iterations = 5,
                                       min_DE = 0,
                                       n_cores = 1,
                                       verbose = FALSE)

  p <- plotHistogram(input = permute_output,
                     use_splits = c("split1", "split2"),
                     label_pvalue = FALSE)

  expect_type(p, "list")
  expect_named(p, c("split1", "split2"))
  expect_true(all(vapply(p, inherits, logical(1), what = "ggplot")))
})

# ---------------------------------------------------------------------------
# plotFeature
# ---------------------------------------------------------------------------

test_that("plotFeature returns a ggplot for one split", {
  input <- setInput.CellMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 pseudobulk = "generate",
                                 split_labels = "split_multi",
                                 design = NULL,
                                 n_cores = 1,
                                 verbose = FALSE)

  feature_name <- runDE_output$DE_results$feature[[1]]

  p <- plotFeature(input = runDE_output,
                   feature_name = feature_name,
                   use_splits = "split1",
                   normalization_method = "cpm",
                   plot_type = "boxplot")

  expect_s3_class(p, "ggplot")
})

test_that("plotFeature returns a named list for multiple splits", {
  input <- setInput.CellMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 pseudobulk = "generate",
                                 split_labels = "split_multi",
                                 design = NULL,
                                 n_cores = 1,
                                 verbose = FALSE)

  feature_name <- runDE_output$DE_results$feature[[1]]

  p <- plotFeature(input = runDE_output,
                   feature_name = feature_name,
                   use_splits = c("split1", "split2"),
                   normalization_method = "none",
                   plot_type = "bar_se")

  expect_type(p, "list")
  expect_named(p, c("split1", "split2"))
  expect_true(all(vapply(p, inherits, logical(1), what = "ggplot")))
})

test_that("plotFeature supports plot variants", {
  input <- setInput.CellMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 pseudobulk = "generate",
                                 split_labels = "split_multi",
                                 design = NULL,
                                 n_cores = 1,
                                 verbose = FALSE)

  feature_name <- runDE_output$DE_results$feature[[1]]

  p_bar_sd <- plotFeature(input = runDE_output,
                          feature_name = feature_name,
                          use_splits = "split1",
                          normalization_method = "log_cpm",
                          plot_type = "bar_sd")

  p_beeswarm <- plotFeature(input = runDE_output,
                            feature_name = feature_name,
                            use_splits = "split1",
                            normalization_method = "none",
                            plot_type = "beeswarm",
                            label_statistics = FALSE)

  expect_s3_class(p_bar_sd, "ggplot")
  expect_s3_class(p_beeswarm, "ggplot")
})

test_that("plotFeature warns and returns NULL for missing feature", {
  input <- setInput.CellMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 pseudobulk = "generate",
                                 split_labels = "split_multi",
                                 design = NULL,
                                 n_cores = 1,
                                 verbose = FALSE)

  expect_warning(
    p <- plotFeature(input = runDE_output,
                     feature_name = "not_a_real_gene",
                     use_splits = "split1"),
    "not_a_real_gene")

  expect_null(p)
})

# ---------------------------------------------------------------------------
# plotDimReduction
# ---------------------------------------------------------------------------

test_that("plotDimReduction returns ggplot with no color_by", {
  reduction <- matrix(rnorm(24 * 2),
                      ncol = 2,
                      dimnames = list(paste0("cell", 1:24), c("UMAP_1", "UMAP_2")))

  p <- plotDimReduction(reduction = reduction,
                        color_by = NULL)

  expect_s3_class(p, "ggplot")
})

test_that("plotDimReduction returns ggplot colored by split", {
  reduction <- matrix(rnorm(24 * 2),
                      ncol = 2,
                      dimnames = list(paste0("cell", 1:24), c("UMAP_1", "UMAP_2")))

  split_labels <- makeMetaData()$split_multi

  p <- plotDimReduction(reduction = reduction,
                        split_labels = split_labels,
                        color_by = "split",
                        label_splits = TRUE)

  expect_s3_class(p, "ggplot")
})

test_that("plotDimReduction returns ggplot colored by permutation metrics", {
  input <- setInput.CellMatrix()

  permute_output <- testCase.permuteDE(input = input,
                                       de_method = "edgeR",
                                       de_test = "LRT",
                                       pseudobulk = "generate",
                                       split_labels = "split_multi",
                                       design = NULL,
                                       return_all = FALSE,
                                       n_iterations = 5,
                                       min_DE = 0,
                                       n_cores = 1,
                                       verbose = FALSE)

  reduction <- matrix(rnorm(24 * 2),
                      ncol = 2,
                      dimnames = list(paste0("cell", 1:24), c("UMAP_1", "UMAP_2")))

  split_labels <- makeMetaData()$split_multi

  # Avoid Seurat warning when all plotted p-values are identical.
  permute_output_for_plot <- permute_output
  permute_output_for_plot$permutation_test_results$pvalue <- seq(from = 0.2, to = 0.8,
                                                                 length.out = nrow(permute_output_for_plot$permutation_test_results))

  p_nsig <- plotDimReduction(reduction = reduction,
                             input = permute_output,
                             split_labels = split_labels,
                             color_by = "n_sig",
                             label_statistics = TRUE)

  p_pvalue <- plotDimReduction(reduction = reduction,
                               input = permute_output_for_plot,
                               split_labels = split_labels,
                               color_by = "pvalue",
                               permutation_test_alpha = 1)

  expect_s3_class(p_nsig, "ggplot")
  expect_s3_class(p_pvalue, "ggplot")
})

test_that("plotDimReduction returns ggplot colored by feature values", {
  reduction <- matrix(rnorm(24 * 2),
                      ncol = 2,
                      dimnames = list(paste0("cell", 1:24), c("UMAP_1", "UMAP_2")))

  feature_values <- stats::rpois(24, lambda = 5)

  p <- plotDimReduction(reduction = reduction,
                        feature_values = feature_values,
                        feature_name = "gene1",
                        color_by = "feature")

  expect_s3_class(p, "ggplot")
})
