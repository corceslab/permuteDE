# ---------------------------------------------------------------------------
# Tests for function permuteDE: Grid
# ---------------------------------------------------------------------------

# Set up counter
total_tests <- sum(vapply(runDE_input_list,
                          FUN = function(x) nrow(getInputGrid(x)),
                          FUN.VALUE = numeric(1)))

test_counter <- 0L

# Store permutation results for cross-input comparison
permute_results_by_case <- new.env(parent = emptyenv())

# Loop through all input, DE method, DE test, pseudobulk, split, and design combinations
for (input_name in names(runDE_input_list)) {
  input_grid <- getInputGrid(runDE_input_list[[input_name]])

  for (i in seq_len(nrow(input_grid))) {
    local({
      test_counter <<- test_counter + 1L
      test_index_i <- test_counter
      total_tests_i <- total_tests

      input_name_i <- input_name
      input_i <- runDE_input_list[[input_name]]$fun

      de_method_i <- input_grid$de_method[[i]]
      de_test_i <- input_grid$de_test[[i]]
      pseudobulk_i <- input_grid$pseudobulk[[i]]

      split_condition_i <- input_grid$split_condition[[i]]
      split_labels_i <- input_grid$split_labels[[i]]
      if (is.na(split_labels_i)) {
        split_labels_i <- NULL
      }

      design_condition_i <- input_grid$design_condition[[i]]
      design_i <- input_grid$design[[i]]
      if (is.na(design_i)) {
        design_i <- NULL
      }

      expect_error_i <- input_grid$expect_error[[i]]

      test_that(paste("permuteDE works for", input_name_i,
                      "with", de_method_i, de_test_i,
                      "pseudobulk", pseudobulk_i,
                      "split", split_condition_i,
                      "design", design_condition_i),
                {
                  testthat::skip_on_os("windows")
                  testthat::skip_if(parallel::detectCores() < 2, "Need at least 2 cores")

                  message("[", test_index_i, "/", total_tests_i, "] ",
                          "permuteDE test: input = ", input_name_i,
                          ", method = ", de_method_i,
                          ", test = ", de_test_i,
                          ", pseudobulk = ", pseudobulk_i,
                          ", split = ", split_condition_i,
                          ", design = ", design_condition_i,
                          ", expect_error = ", expect_error_i)

                  input <- input_i()

                  if (expect_error_i) {
                    expect_error(testCase.permuteDE(input = input,
                                                    de_method = de_method_i,
                                                    de_test = de_test_i,
                                                    pseudobulk = pseudobulk_i,
                                                    split_labels = split_labels_i,
                                                    design = design_i,
                                                    return_all = FALSE,
                                                    n_iterations = 5,
                                                    min_DE = 0,
                                                    n_cores = 2,
                                                    verbose = FALSE),
                                 "BPCells|IterableMatrix|de_method|design|covariate|Wilcoxon|exact")
                  } else {
                    output <- testCase.permuteDE(input = input,
                                                 de_method = de_method_i,
                                                 de_test = de_test_i,
                                                 pseudobulk = pseudobulk_i,
                                                 split_labels = split_labels_i,
                                                 design = design_i,
                                                 return_all = FALSE,
                                                 n_iterations = 5,
                                                 min_DE = 0,
                                                 n_cores = 2,
                                                 verbose = FALSE)

                    case_key <- paste(de_method_i,
                                      de_test_i,
                                      pseudobulk_i,
                                      split_condition_i,
                                      design_condition_i,
                                      sep = "__")

                    if (is.null(permute_results_by_case[[case_key]])) {
                      permute_results_by_case[[case_key]] <- list()
                    }

                    permute_results_by_case[[case_key]][[input_name_i]] <-
                      standardizePermutationResults(output)

                    expect_equal(output$parameters$de_method, de_method_i)
                    expect_equal(output$parameters$de_test, de_test_i)
                    expect_equal(output$parameters$pseudobulk, pseudobulk_i)
                    expect_equal(output$parameters$return_all, FALSE)
                    expect_equal(output$parameters$n_iterations, 5)
                    expect_equal(output$parameters$min_DE, 0)

                    expect_gt(nrow(output$permutation_test_results), 0)
                    expect_gt(nrow(output$permutation_DE_summary), 0)
                  }
                }
      )
    })
  }
}

# Compare permuteDE results across input types for matching method/test/pseudobulk/split/design cases
message("Stored cross-input permuteDE result counts:")

for (case_key in ls(permute_results_by_case)) {
  message(case_key, ": ",
          paste(names(permute_results_by_case[[case_key]]), collapse = ", "))
}

test_that("permuteDE results match across comparable input types", {
  testthat::skip_if(length(ls(permute_results_by_case)) == 0)

  for (case_key in ls(permute_results_by_case)) {
    results_i <- permute_results_by_case[[case_key]]

    if (length(results_i) < 2) {
      next
    }

    reference_name <- names(results_i)[[1]]
    reference_results <- results_i[[reference_name]]

    for (input_name in names(results_i)[-1]) {
      expect_equal(
        results_i[[input_name]],
        reference_results,
        tolerance = 1e-8,
        info = paste(
          "Case:", case_key,
          "| Compared", input_name, "to", reference_name
        )
      )
    }
  }
})
