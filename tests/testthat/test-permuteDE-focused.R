# ---------------------------------------------------------------------------
# Tests for function permuteDE: Focused
# ---------------------------------------------------------------------------

test_that("permuteDE supports return_all TRUE", {
  testthat::skip_on_os("windows")
  testthat::skip_if(parallel::detectCores() < 2, "Need at least 2 cores")

  input <- setInput.PBMatrix()

  output <- testCase.permuteDE(input = input,
                               de_method = "edgeR",
                               de_test = "LRT",
                               pseudobulk = "supplied",
                               split_labels = NULL,
                               design = NULL,
                               return_all = TRUE,
                               n_iterations = 5,
                               min_DE = 0,
                               n_cores = 2,
                               verbose = FALSE)
  expect_true("permutation_DE_results" %in% names(output))
  expect_gt(nrow(output$permutation_DE_results), 0)
})

test_that("permuteDE respects use_splits", {
  input <- setInput.CellMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 pseudobulk = "generate",
                                 split_labels = "split_multi",
                                 design = NULL,
                                 n_cores = 2,
                                 verbose = FALSE)

  output <- permuteDE(input = runDE_output,
                      n_iterations = 5,
                      use_splits = "split1",
                      min_DE = 0,
                      n_cores = 2,
                      verbose = FALSE)

  expectValidOutput.permuteDE(output, return_all = FALSE)
  expect_equal(output$parameters$use_splits, "split1")
  expect_named(output$metadata$permutation_reference_group_indices, "split1")
  expect_equal(unique(output$permutation_test_results$split), "split1")
  expect_equal(unique(output$permutation_DE_summary$split), "split1")
})

test_that("permuteDE supports permute_within", {
  testthat::skip_on_os("windows")
  testthat::skip_if(parallel::detectCores() < 2, "Need at least 2 cores")

  input <- setInput.CellMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 pseudobulk = "generate",
                                 split_labels = NULL,
                                 design = "~ batch + group",
                                 n_cores = 2,
                                 verbose = FALSE)

  output <- permuteDE(input = runDE_output,
                      n_iterations = 5,
                      permute_within = "batch",
                      min_DE = 0,
                      n_cores = 2,
                      verbose = FALSE)

  expectValidOutput.permuteDE(output, return_all = FALSE)
  expect_equal(output$parameters$permute_within, "batch")
  expect_gt(nrow(output$permutation_test_results), 0)
})

test_that("permuteDE supports permute_by for cell-level tests", {
  testthat::skip_on_os("windows")
  testthat::skip_if(parallel::detectCores() < 2, "Need at least 2 cores")

  input <- setInput.CellMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 pseudobulk = "none",
                                 split_labels = NULL,
                                 design = NULL,
                                 n_cores = 2,
                                 verbose = FALSE)

  warnings <- character()

  output <- withCallingHandlers(permuteDE(input = runDE_output,
                                          n_iterations = 5,
                                          permute_by = "replicate",
                                          min_DE = 0,
                                          n_cores = 2,
                                          verbose = FALSE),
                                warning = function(w) {
                                  warnings <<- c(warnings, conditionMessage(w))
                                  invokeRestart("muffleWarning")
                                })

  expect_true(any(grepl("Cell-level tests are not recommended", warnings)))
  expect_false(any(grepl("consider providing biological replicates labels to parameter 'permute_by'", warnings)))

  expectValidOutput.permuteDE(output, return_all = FALSE)
  expect_equal(output$parameters$permute_by, "replicate")
  expect_gt(nrow(output$permutation_test_results), 0)
})

test_that("permuteDE gives identical output with the same random seed", {
  testthat::skip_on_os("windows")
  testthat::skip_if(parallel::detectCores() < 2, "Need at least 2 cores")

  input <- setInput.PBMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 pseudobulk = "supplied",
                                 split_labels = NULL,
                                 design = NULL,
                                 n_cores = 2,
                                 verbose = FALSE)

  output_1 <- permuteDE(input = runDE_output,
                        n_iterations = 5,
                        min_DE = 0,
                        random_seed = 123,
                        n_cores = 2,
                        verbose = FALSE)

  output_2 <- permuteDE(input = runDE_output,
                        n_iterations = 5,
                        min_DE = 0,
                        random_seed = 123,
                        n_cores = 2,
                        verbose = FALSE)

  expectValidOutput.permuteDE(output_1, return_all = FALSE)
  expectValidOutput.permuteDE(output_2, return_all = FALSE)

  expect_equal(output_1$permutation_test_results,
               output_2$permutation_test_results)

  expect_equal(output_1$permutation_DE_summary,
               output_2$permutation_DE_summary)

  expect_equal(output_1$metadata$runDE_values,
               output_2$metadata$runDE_values)

  expect_equal(output_1$metadata$permutation_reference_group_indices,
               output_2$metadata$permutation_reference_group_indices)

  expect_equal(output_1$parameters,
               output_2$parameters)
})

test_that("permuteDE uses covariate design during DESeq2 LRT permutations", {
  testthat::skip_on_os("windows")
  testthat::skip_if(parallel::detectCores() < 2, "Need at least 2 cores")

  input <- setInput.CellMatrix()

  runDE_output <- testCase.runDE(input = input,
                                 de_method = "DESeq2",
                                 de_test = "LRT",
                                 pseudobulk = "generate",
                                 split_labels = NULL,
                                 design = "~ batch + group",
                                 n_cores = 2,
                                 verbose = FALSE)

  output <- permuteDE(input = runDE_output,
                      n_iterations = 5,
                      min_DE = 0,
                      n_cores = 2,
                      verbose = FALSE)

  expectValidOutput.permuteDE(output, return_all = FALSE)
  expect_equal(output$parameters$de_method, "DESeq2")
  expect_equal(output$parameters$de_test, "LRT")
})

test_that("permuteDE handles small permutation spaces without dropping permutation matrix dimensions", {
  input <- setInput.PBMatrix()

  output_runDE <- testCase.runDE(input = input,
                                 split_labels = NULL,
                                 pseudobulk = "supplied",
                                 de_method = "edgeR",
                                 de_test = "LRT",
                                 return_raw_de = FALSE,
                                 min_replicates_per_group = 2,
                                 n_cores = 1,
                                 verbose = FALSE)

  output_perm <- permuteDE(input = output_runDE,
                           alpha = 1,
                           lfc_threshold = 0,
                           n_iterations = 2,
                           min_DE = 0,
                           return_all = FALSE,
                           n_cores = 1,
                           verbose = FALSE)

  expectValidOutput.permuteDE(output_perm)

  expect_s3_class(output_perm$permutation_test_results, "data.frame")
  expect_s3_class(output_perm$permutation_DE_summary, "data.frame")

  expect_true(all(output_perm$permutation_test_results$n_iterations >= 1))
  expect_true(all(output_perm$permutation_DE_summary$permutation >= 1))

  expect_true(is.list(output_perm$metadata$permutation_reference_group_indices))

  for (split_i in names(output_perm$metadata$permutation_reference_group_indices)) {
    indices_i <- output_perm$metadata$permutation_reference_group_indices[[split_i]]

    expect_true(is.matrix(indices_i))
    expect_gt(ncol(indices_i), 0)
  }
})

test_that("permuteDE uses runDE excluded features when normalize_prefilter is TRUE", {
  input <- setInput.CellMatrix()

  forced_excluded_feature <- rownames(input$object)[1]

  input$object[forced_excluded_feature, ] <- 0
  input$object[forced_excluded_feature, seq_len(3)] <- 1

  output_runDE <- runDE(object = input$object,
                        metadata = input$metadata,
                        replicate_labels = "replicate",
                        group_labels = "group",
                        split_labels = NULL,
                        pseudobulk = "generate",
                        de_method = "edgeR",
                        de_test = "LRT",
                        normalize_prefilter = TRUE,
                        min_cells_per_split = 1,
                        min_cells_per_replicate = 1,
                        min_replicates_per_split = 2,
                        min_replicates_per_group = 1,
                        min_cells_per_feature = 4,
                        min_prop_cells_per_feature = 0,
                        n_cores = 1,
                        verbose = FALSE)

  expectValidOutput.runDE(output_runDE)

  expect_true(isTRUE(output_runDE$parameters$normalize_prefilter))
  expect_true("exclude_features" %in% names(output_runDE$metadata))

  excluded <- output_runDE$metadata$exclude_features

  expect_type(excluded, "list")
  expect_true(any(vapply(excluded,
                         FUN = function(x) forced_excluded_feature %in% x,
                         FUN.VALUE = logical(1))))

  expect_false(forced_excluded_feature %in% output_runDE$DE_results$feature)

  output_perm <- permuteDE(input = output_runDE,
                           alpha = 1,
                           lfc_threshold = 0,
                           n_iterations = 3,
                           min_DE = 0,
                           return_all = TRUE,
                           n_cores = 1,
                           verbose = FALSE)

  expectValidOutput.permuteDE(output_perm, return_all = TRUE)

  expect_false(
    forced_excluded_feature %in% output_perm$permutation_DE_results$feature
  )
})

test_that("permuteDE uses fdrtool p-value adjustment by default", {
  skip_if_not_installed("fdrtool")

  input <- setInput.PBMatrix()

  captured_runDE <- capture_warnings(runDE(object = input$object,
                                           metadata = input$metadata,
                                           replicate_labels = "replicate",
                                           group_labels = "group",
                                           split_labels = NULL,
                                           pseudobulk = "supplied",
                                           de_method = "edgeR",
                                           de_test = "LRT",
                                           p_adjust_method = "fdrtool",
                                           min_replicates_per_split = 2,
                                           min_replicates_per_group = 1,
                                           n_cores = 1,
                                           verbose = FALSE))

  output_runDE <- captured_runDE$value

  expect_true(any(grepl("fdrtool.*discretion", captured_runDE$warnings)))
  expectValidOutput.runDE(output_runDE)
  expect_equal(output_runDE$parameters$p_adjust_method, "fdrtool")

  captured_perm <- capture_warnings(permuteDE(input = output_runDE,
                                              alpha = 1,
                                              lfc_threshold = 0,
                                              n_iterations = 3,
                                              min_DE = 0,
                                              return_all = TRUE,
                                              n_cores = 1,
                                              verbose = FALSE))

  output_perm <- captured_perm$value

  expectValidOutput.permuteDE(output_perm, return_all = TRUE)
  expect_equal(output_perm$parameters$p_adjust_method, "fdrtool")
  expect_equal(output_perm$parameters$de_method, "edgeR")
  expect_equal(output_perm$parameters$de_test, "LRT")
  expect_true(all(is.na(output_perm$permutation_DE_results$padj) |
                    (output_perm$permutation_DE_results$padj >= 0 &
                       output_perm$permutation_DE_results$padj <= 1)))
  expect_true(all(is.na(output_perm$permutation_DE_summary$n_sig) |
                    output_perm$permutation_DE_summary$n_sig >= 0))
})

test_that("permuteDE uses fdrtool zscore adjustment for DESeq2 Wald", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("fdrtool")

  input <- setInput.PBMatrix()

  captured_runDE <- capture_warnings(runDE(object = input$object,
                                           metadata = input$metadata,
                                           replicate_labels = "replicate",
                                           group_labels = "group",
                                           split_labels = NULL,
                                           pseudobulk = "supplied",
                                           de_method = "DESeq2",
                                           de_test = "Wald",
                                           p_adjust_method = "fdrtool",
                                           de_params = list(fdrtool = list(statistic = "zscore")),
                                           min_replicates_per_split = 2,
                                           min_replicates_per_group = 1,
                                           n_cores = 1,
                                           verbose = FALSE))

  output_runDE <- captured_runDE$value

  expect_true(any(grepl("fdrtool.*discretion", captured_runDE$warnings)))
  expectValidOutput.runDE(output_runDE)

  expect_equal(output_runDE$parameters$p_adjust_method, "fdrtool")
  expect_equal(output_runDE$parameters$de_params$fdrtool$statistic, "zscore")

  captured_perm <- capture_warnings(permuteDE(input = output_runDE,
                                              alpha = 1,
                                              lfc_threshold = 0,
                                              n_iterations = 5,
                                              min_DE = 0,
                                              return_all = TRUE,
                                              n_cores = 1,
                                              verbose = FALSE))

  output_perm <- captured_perm$value

  expectValidOutput.permuteDE(output_perm, return_all = TRUE)
  expect_equal(output_perm$parameters$p_adjust_method, "fdrtool")
  expect_equal(output_perm$parameters$de_method, "DESeq2")
  expect_equal(output_perm$parameters$de_test, "Wald")
  expect_equal(output_perm$parameters$de_params$fdrtool$statistic, "zscore")
  expect_true(all(is.na(output_perm$permutation_DE_results$padj) |
                    (output_perm$permutation_DE_results$padj >= 0 &
                       output_perm$permutation_DE_results$padj <= 1)))
  expect_true(all(is.na(output_perm$permutation_DE_summary$n_sig) |
                    output_perm$permutation_DE_summary$n_sig >= 0))
})
