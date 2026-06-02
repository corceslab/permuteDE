# ---------------------------------------------------------------------------
# Tests for function runDE: Grid
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

      test_that(
        paste("runDE works for", input_name_i,
              "with", de_method_i, de_test_i,
              "pseudobulk", pseudobulk_i,
              "split", split_condition_i,
              "design", design_condition_i),
        {
          testthat::skip_on_os("windows")
          testthat::skip_if(parallel::detectCores() < 2, "Need at least 2 cores")

          message("[", test_index_i, "/", total_tests_i, "] ",
                  "runDE test: input = ", input_name_i,
                  ", method = ", de_method_i,
                  ", test = ", de_test_i,
                  ", pseudobulk = ", pseudobulk_i,
                  ", split = ", split_condition_i,
                  ", design = ", design_condition_i,
                  ", expect_error = ", expect_error_i)

          input <- input_i()

          if (expect_error_i) {
            expect_error(testCase.runDE(input = input,
                                        de_method = de_method_i,
                                        de_test = de_test_i,
                                        pseudobulk = pseudobulk_i,
                                        split_labels = split_labels_i,
                                        design = design_i,
                                        n_cores = 2),
                         "BPCells|IterableMatrix|de_method|design|covariate|Wilcoxon|exact")
          } else {
            output <- testCase.runDE(input = input,
                                     de_method = de_method_i,
                                     de_test = de_test_i,
                                     pseudobulk = pseudobulk_i,
                                     split_labels = split_labels_i,
                                     design = design_i,
                                     n_cores = 2)

            case_key <- paste(de_method_i,
                              de_test_i,
                              pseudobulk_i,
                              split_condition_i,
                              design_condition_i,
                              sep = "__")

            if (is.null(de_results_by_case[[case_key]])) {
              de_results_by_case[[case_key]] <- list()
            }

            de_results_by_case[[case_key]][[input_name_i]] <- standardizeDEResults(
              output$DE_results
            )

            expect_equal(output$parameters$de_method, de_method_i)
            expect_equal(output$parameters$de_test, de_test_i)
            expect_equal(output$parameters$pseudobulk, pseudobulk_i)
            expect_gt(nrow(output$DE_results), 0)
          }
        }
      )
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

test_that("runDE requires replicate_labels to match supplied pseudobulk column names", {
  input <- setInput.PBMatrix()

  # Positive control: exact column-name labels should work.
  metadata_ok <- input$metadata
  metadata_ok$replicate_id <- rownames(metadata_ok)

  warnings <- character()

  output <- withCallingHandlers(runDE(object = input$object,
                                      metadata = metadata_ok,
                                      replicate_labels = "replicate_id",
                                      group_labels = input$group_labels,
                                      split_labels = NULL,
                                      pseudobulk = "supplied",
                                      de_method = "edgeR",
                                      de_test = "LRT",
                                      min_replicates_per_split = 1,
                                      min_replicates_per_group = 1,
                                      n_cores = 1,
                                      verbose = FALSE),
                                warning = function(w) {
                                  warnings <<- c(warnings, conditionMessage(w))
                                  invokeRestart("muffleWarning")
                                })

  expect_true(any(grepl("min_cells_per_split.*not used", warnings)))
  expect_true(any(grepl("min_cells_per_replicate.*not used", warnings)))
  expect_true(any(grepl("min_cells_per_feature.*not used", warnings)))
  expect_true(any(grepl("min_prop_cells_per_feature.*not used", warnings)))

  expectValidOutput.runDE(output)
  expect_equal(output$parameters$pseudobulk, "supplied")
  expect_true(identical(output$metadata$group_key$replicate, colnames(input$object)))

  # Negative control: sample labels are same length/order, but do not match
  # the supplied pseudobulk matrix column names.
  expect_error(suppressWarnings(runDE(object = input$object,
                                      metadata = input$metadata,
                                      replicate_labels = "replicate",
                                      group_labels = input$group_labels,
                                      split_labels = NULL,
                                      pseudobulk = "supplied",
                                      de_method = "edgeR",
                                      de_test = "LRT",
                                      min_replicates_per_split = 1,
                                      min_replicates_per_group = 1,
                                      n_cores = 1,
                                      verbose = FALSE)),
               "replicate_labels.*match.*column names|column names.*replicate_labels")

  # Negative control: same labels as column names, but wrong order.
  metadata_wrong_order <- metadata_ok
  metadata_wrong_order$replicate_id <- rev(metadata_wrong_order$replicate_id)

  expect_error(suppressWarnings(runDE(object = input$object,
                                      metadata = metadata_wrong_order,
                                      replicate_labels = "replicate_id",
                                      group_labels = input$group_labels,
                                      split_labels = NULL,
                                      pseudobulk = "supplied",
                                      de_method = "edgeR",
                                      de_test = "LRT",
                                      min_replicates_per_split = 1,
                                      min_replicates_per_group = 1,
                                      n_cores = 1,
                                      verbose = FALSE)),
               "replicate_labels.*match.*column names|column names.*replicate_labels")
})
