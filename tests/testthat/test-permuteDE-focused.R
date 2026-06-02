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
