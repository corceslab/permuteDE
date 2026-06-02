# ---------------------------------------------------------------------------
# Tests for function runDE: Focused
# ---------------------------------------------------------------------------

test_that("runDE errors when split_labels contains NA values", {
  input <- setInput.CellMatrix()

  split_labels <- input$metadata$split_multi
  split_labels[split_labels == "split2"] <- NA

  expect_error(testCase.runDE(input = input,
                              de_method = "edgeR",
                              de_test = "LRT",
                              pseudobulk = "generate",
                              split_labels = split_labels,
                              design = NULL,
                              n_cores = 1,
                              verbose = FALSE),
               "split_labels.*cannot be NA")
})

test_that("runDE errors when unsupported tests are used with covariate design", {
  input <- setInput.CellMatrix()

  expect_error(testCase.runDE(input = input,
                              de_method = "edgeR",
                              de_test = "exact",
                              pseudobulk = "generate",
                              split_labels = NULL,
                              design = "~ batch + group"),
               "design|covariate|exact")

  expect_error(testCase.runDE(input = input,
                              de_method = "limma",
                              de_test = "wilcox_cpm",
                              pseudobulk = "generate",
                              split_labels = NULL,
                              design = "~ batch + group"),
               "design|covariate|Wilcoxon|wilcox")
})

test_that("runDE returns raw DE output for all method/test combinations with Seurat BPCells input", {
  input <- setInput.BPCells()

  for (i in seq_len(nrow(runDE_method_grid))) {
    de_method_i <- runDE_method_grid$de_method[[i]]
    de_test_i <- runDE_method_grid$de_test[[i]]

    # BPCells DE requires existing IterableMatrix input
    # so test it with pseudobulk = "none"
    # Other methods use generated pseudobulk matrices
    pseudobulk_i <- if (identical(de_method_i, "BPCells")) {
      "none"
    } else {
      "generate"
    }

    message("[raw DE output BPCells input] method = ", de_method_i,
            ", test = ", de_test_i,
            ", pseudobulk = ", pseudobulk_i)

    output <- testCase.runDE(input = input,
                             de_method = de_method_i,
                             de_test = de_test_i,
                             pseudobulk = pseudobulk_i,
                             split_labels = NULL,
                             design = NULL,
                             return_raw_de = TRUE,
                             n_cores = 1,
                             verbose = FALSE)

    expected_name <- expectedValueName(pseudobulk_i)

    expect_true("raw_DE_results" %in% names(output))
    expect_true(isTRUE(output$parameters$return_raw_de))

    expect_type(output$raw_DE_results, "list")
    expect_gt(length(output$raw_DE_results), 0)
    expect_named(output$raw_DE_results, names(output[[expected_name]]))

    expect_true(all(vapply(output$raw_DE_results,
                           FUN = is.data.frame,
                           FUN.VALUE = logical(1))))

    expect_true(all(vapply(output$raw_DE_results,
                           FUN = nrow,
                           FUN.VALUE = integer(1)) > 0))

    expect_equal(output$parameters$de_method, de_method_i)
    expect_equal(output$parameters$de_test, de_test_i)
    expect_equal(output$parameters$pseudobulk, pseudobulk_i)

    expectValidOutput.runDE(output)
  }
})

test_that("runDE resolves split_labels from metadata for supplied pseudobulk matrices", {
  input <- setInput.PBMatrix()

  warnings <- character()

  output <- withCallingHandlers(
    runDE(object = input$object,
      metadata = input$metadata,
      replicate_labels = input$replicate_labels,
      group_labels = input$group_labels,
      split_labels = "split_multi",
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
  expect_equal(output$parameters$split_labels, input$metadata$split_multi)
  expect_setequal(names(output$PB_values), c("split1", "split2"))
  expect_setequal(unique(output$DE_results$split), c("split1", "split2"))
})

test_that("runDE accepts split_labels vector for supplied pseudobulk matrices", {
  input <- setInput.PBMatrix()

  output <- testCase.runDE(input = input,
    de_method = "edgeR",
    de_test = "LRT",
    pseudobulk = "supplied",
    split_labels = input$metadata$split_multi,
    design = NULL,
    return_raw_de = FALSE,
    n_cores = 1,
    verbose = FALSE)

  expectValidOutput.runDE(output)

  expect_equal(output$parameters$pseudobulk, "supplied")
  expect_setequal(names(output$PB_values), c("split1", "split2"))
  expect_setequal(unique(output$DE_results$split), c("split1", "split2"))
})
