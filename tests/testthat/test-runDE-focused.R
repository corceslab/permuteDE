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

test_that("runDE requires replicate_labels to match supplied pseudobulk column names", {
  input <- setInput.PBMatrix()

  # Positive control: exact column-name labels should work.
  metadata_ok <- input$metadata
  metadata_ok$replicate_id <- colnames(input$object)

  warnings <- character()

  output <- withCallingHandlers(
    runDE(object = input$object,
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
  expect_setequal(output$metadata$group_key$replicate, colnames(input$object))

  # Negative control: same length/order, but values do not match column names.
  metadata_bad_names <- metadata_ok
  metadata_bad_names$bad_replicate_id <- paste0("wrong_", seq_len(ncol(input$object)))

  expect_error(suppressWarnings(
    runDE(object = input$object,
          metadata = metadata_bad_names,
          replicate_labels = "bad_replicate_id",
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
  metadata_wrong_order$replicate_id <- rev(colnames(input$object))

  expect_error(suppressWarnings(runDE(
    object = input$object,
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

test_that("runDE returns raw coefficient results for coefficient-based model tests", {
  cases <- data.frame(de_method = c("edgeR", "limma", "DESeq2"),
                      de_test = c("LRT", "voom", "Wald"),
                      stringsAsFactors = FALSE)

  for (i in seq_len(nrow(cases))) {
    de_method_i <- cases$de_method[[i]]
    de_test_i <- cases$de_test[[i]]

    skipUninstalled(de_method_i)

    input <- setInput.PBMatrix()

    output <- testCase.runDE(input = input,
                             de_method = de_method_i,
                             de_test = de_test_i,
                             pseudobulk = "supplied",
                             split_labels = NULL,
                             design = "~ age + batch + group",
                             return_raw_de = TRUE,
                             n_cores = 1,
                             verbose = FALSE,
                             de_params = list(return_all_coefficients = TRUE))

    expectValidOutput.runDE(output)

    expect_true(isTRUE(output$parameters$return_raw_de))
    expect_true(isTRUE(output$parameters$de_params$return_all_coefficients))

    expect_type(output$raw_DE_results, "list")
    expect_setequal(names(output$raw_DE_results), names(output$PB_values))

    for (split_i in names(output$raw_DE_results)) {
      raw_i <- output$raw_DE_results[[split_i]]

      expect_type(raw_i, "list")
      expect_true(all(c("group", "coefficients") %in% names(raw_i)))

      expect_s3_class(raw_i$group, "data.frame")
      expect_gt(nrow(raw_i$group), 0)

      expect_type(raw_i$coefficients, "list")
      expect_gt(length(raw_i$coefficients), 0)

      # All three of these model-based cases should include age, at least one
      # batch coefficient, and the group coefficient.
      coefficient_names <- names(raw_i$coefficients)

      expect_true(any(grepl("^age$", coefficient_names)))

      if (de_method_i == "DESeq2") {
        expect_true(any(grepl("batch", coefficient_names)))
        expect_true(any(grepl("group", coefficient_names)))
      } else {
        expect_true(any(grepl("^batch", coefficient_names)))
        expect_true(any(grepl("^group", coefficient_names)))
      }

      expect_true(all(vapply(raw_i$coefficients,
                             FUN = is.data.frame,
                             FUN.VALUE = logical(1))))

      expect_true(all(vapply(raw_i$coefficients,
                             FUN = nrow,
                             FUN.VALUE = integer(1)) > 0))
    }
  }
})

test_that("runDE rejects return_all_coefficients for unsupported tests", {
  input <- setInput.PBMatrix()

  expect_error(testCase.runDE(input = input,
                              de_method = "DESeq2",
                              de_test = "LRT",
                              pseudobulk = "supplied",
                              split_labels = NULL,
                              design = "~ age + batch + group",
                              return_raw_de = TRUE,
                              n_cores = 1,
                              verbose = FALSE,
                              de_params = list(return_all_coefficients = TRUE)),
               "return_all_coefficients.*only supported|coefficient-based model tests")

  expect_error(testCase.runDE(input = input,
                              de_method = "limma",
                              de_test = "wilcox_cpm",
                              pseudobulk = "supplied",
                              split_labels = NULL,
                              design = NULL,
                              return_raw_de = TRUE,
                              n_cores = 1,
                              verbose = FALSE,
                              de_params = list(return_all_coefficients = TRUE)),
               "return_all_coefficients.*only supported|coefficient-based model tests")
})

test_that("runDE stores excluded features when normalize_prefilter is TRUE", {
  input <- setInput.CellMatrix()

  output_runDE <- runDE(object = input$object,
                        metadata = input$metadata,
                        replicate_labels = "replicate",
                        group_labels = "group",
                        split_labels = "split",
                        pseudobulk = "generate",
                        de_method = "edgeR",
                        de_test = "LRT",
                        normalize_prefilter = TRUE,
                        min_cells_per_split = 1,
                        min_cells_per_replicate = 1,
                        min_replicates_per_split = 2,
                        min_replicates_per_group = 1,
                        min_cells_per_feature = 2,
                        min_prop_cells_per_feature = 0,
                        n_cores = 1,
                        verbose = FALSE)

  expectValidOutput.runDE(output_runDE)

  expect_true(isTRUE(output_runDE$parameters$normalize_prefilter))
  expect_true("exclude_features" %in% names(output_runDE$metadata))

  expect_true(
    is.null(output_runDE$metadata$exclude_features) ||
      is.list(output_runDE$metadata$exclude_features)
  )
})

test_that("fdrtool p-value adjustment works by default", {
  skip_if_not_installed("fdrtool")

  input <- setInput.PBMatrix()

  captured <- capture_warnings(runDE(object = input$object,
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

  output <- captured$value

  expect_true(any(grepl("fdrtool.*discretion", captured$warnings)))
  expect_true(any(grepl("min_cells_per_split.*not used", captured$warnings)))
  expect_true(any(grepl("min_cells_per_replicate.*not used", captured$warnings)))
  expect_true(any(grepl("min_cells_per_feature.*not used", captured$warnings)))
  expect_true(any(grepl("min_prop_cells_per_feature.*not used", captured$warnings)))

  expectValidOutput.runDE(output)
  expect_equal(output$parameters$p_adjust_method, "fdrtool")
  expect_true(all(is.na(output$DE_results$padj) |
                    (output$DE_results$padj >= 0 & output$DE_results$padj <= 1)))
})

test_that("fdrtool zscore adjustment is supported for DESeq2 Wald", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("fdrtool")

  input <- setInput.PBMatrix()

  captured <- capture_warnings(runDE(object = input$object,
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

  output <- captured$value

  expect_true(any(grepl("fdrtool.*discretion", captured$warnings)))
  expect_true(any(grepl("min_cells_per_split.*not used", captured$warnings)))
  expect_true(any(grepl("min_cells_per_replicate.*not used", captured$warnings)))
  expect_true(any(grepl("min_cells_per_feature.*not used", captured$warnings)))
  expect_true(any(grepl("min_prop_cells_per_feature.*not used", captured$warnings)))

  expectValidOutput.runDE(output)
  expect_equal(output$parameters$p_adjust_method, "fdrtool")
  expect_equal(output$parameters$de_params$fdrtool$statistic, "zscore")
  expect_true(all(is.na(output$DE_results$padj) |
                    (output$DE_results$padj >= 0 & output$DE_results$padj <= 1)))
})

test_that("fdrtool zscore adjustment is rejected outside DESeq2 Wald", {
  input <- setInput.PBMatrix()

  captured <- capture_warnings(expect_error(runDE(object = input$object,
                                                  metadata = input$metadata,
                                                  replicate_labels = "replicate",
                                                  group_labels = "group",
                                                  split_labels = NULL,
                                                  pseudobulk = "supplied",
                                                  de_method = "edgeR",
                                                  de_test = "LRT",
                                                  p_adjust_method = "fdrtool",
                                                  de_params = list(fdrtool = list(statistic = "zscore")),
                                                  min_replicates_per_split = 2,
                                                  min_replicates_per_group = 1,
                                                  n_cores = 1,
                                                  verbose = FALSE),
                                            "zscore.*DESeq2.*Wald"))

  expect_true(any(grepl("fdrtool.*discretion", captured$warnings)))
  expect_true(any(grepl("min_cells_per_split.*not used", captured$warnings)))
  expect_true(any(grepl("min_cells_per_replicate.*not used", captured$warnings)))
  expect_true(any(grepl("min_cells_per_feature.*not used", captured$warnings)))
  expect_true(any(grepl("min_prop_cells_per_feature.*not used", captured$warnings)))
})
