# ---------------------------------------------------------------------------
# Tests for function getCombinations
# ---------------------------------------------------------------------------

# Numeric input
test_that("getCombinations works with numeric input", {
  output <- getCombinations(n_replicates = 4,
                            n_group1 = 2,
                            n_combinations = 3)
  expect_true(is.matrix(output))
  expect_equal(nrow(output), 2)
  expect_equal(ncol(output), 3)
  expect_true(all(output %in% seq_len(4)))
  # Columns should be unique combinations.
  unique_combinations <- unique(t(apply(output, 2, sort)))
  expect_equal(nrow(unique_combinations), ncol(output))
})

# Vector input
test_that("getCombinations works with replicate and group label vectors", {
  replicate_labels <- rep(paste0("rep", 1:4), each = 3)
  group_labels <- rep(c("A", "A", "B", "B"), each = 3)
  output <- getCombinations(replicate_labels = replicate_labels,
                            group_labels = group_labels,
                            n_combinations = 3)
  expect_true(is.matrix(output))
  expect_equal(nrow(output), 2)
  expect_equal(ncol(output), 3)
  expect_true(all(output %in% seq_len(4)))
})

# Different length inputs
test_that("getCombinations errors when replicate and group label vectors have different lengths", {
  expect_error(getCombinations(replicate_labels = c("rep1", "rep2", "rep3"),
                               group_labels = c("A", "B")),
               "same length")
})

# Number of groups
test_that("getCombinations errors when group labels do not contain exactly two groups", {
  expect_error(getCombinations(replicate_labels = paste0("rep", 1:4),
                               group_labels = c("A", "A", "A", "A")),
               "exactly 2 groups")
  expect_error(getCombinations(replicate_labels = paste0("rep", 1:4),
                               group_labels = c("A", "B", "C", "C")),
               "exactly 2 groups")
})

# NA inputs
test_that("getCombinations rejects NA replicate or group labels", {
  expect_error(getCombinations(replicate_labels = c("rep1", NA, "rep3", "rep4"),
                               group_labels = c("A", "A", "B", "B")),
               "cannot be NA"
  )
  expect_error(getCombinations(replicate_labels = paste0("rep", 1:4),
                               group_labels = c("A", "A", NA, "B")),
               "cannot be NA")
})

# Random seed reproducible
test_that("getCombinations is reproducible with the same random seed", {
  output_1 <- getCombinations(n_replicates = 6,
                              n_group1 = 3,
                              n_combinations = 5,
                              random_seed = 42)
  output_2 <- getCombinations(n_replicates = 6,
                              n_group1 = 3,
                              n_combinations = 5,
                              random_seed = 42)
  expect_identical(output_1, output_2)
})

# Request more than possible combinations
test_that("getCombinations returns all possible combinations when requested number is too large", {
  expect_message(output <- getCombinations(n_replicates = 4,
                                           n_group1 = 2,
                                           n_combinations = 100),
                 "exceeds the number of possible combinations")
  expect_true(is.matrix(output))
  expect_equal(nrow(output), 2)
  expect_equal(ncol(output), choose(4, 2))
})

# One combination
test_that("getCombinations preserves matrix output for one requested combination", {
  output <- getCombinations(n_replicates = 5,
                            n_group1 = 2,
                            n_combinations = 1)
  expect_true(is.matrix(output))
  expect_equal(nrow(output), 2)
  expect_equal(ncol(output), 1)
})

# Warning when both numeric and object inputs provided
test_that("getCombinations ignores object input when numeric counts are supplied", {
  data("sample_data")
  expect_warning(output <- getCombinations(object = sample_data,
                                           n_replicates = 4,
                                           n_group1 = 2,
                                           n_combinations = 3),
                 "not used")
  expect_true(is.matrix(output))
  expect_equal(nrow(output), 2)
  expect_equal(ncol(output), 3)
})

# Covariate confounds
test_that("getCombinations excludes perfectly confounded combinations", {
  # With 4 replicates and choosing 2, these covariate levels would be perfectly
  # separated for combinations c(1, 2) and c(3, 4).
  confound_check <- data.frame(batch = c("B1", "B1", "B2", "B2"))
  output <- getCombinations(n_replicates = 4,
                            n_group1 = 2,
                            n_combinations = 6,
                            confound_check = confound_check)
  expect_true(is.matrix(output))
  sorted_cols <- apply(output, 2, function(x) paste(sort(x), collapse = ","))
  expect_false("1,2" %in% sorted_cols)
  expect_false("3,4" %in% sorted_cols)
})

# Stratified input
test_that("getCombinations handles stratified numeric input", {
  output <- getCombinations(n_replicates = c(4, 4),
                            n_group1 = c(2, 2),
                            n_combinations = 5,
                            random_seed = 1,
                            verbose = FALSE)
  expect_true(is.matrix(output))
  expect_equal(nrow(output), 4)
  expect_equal(ncol(output), 5)
  # Each combination should choose 2 from partition 1 and 2 from partition 2.
  expect_true(all(apply(output, 2, function(x) sum(x %in% 1:4) == 2)))
  expect_true(all(apply(output, 2, function(x) sum(x %in% 5:8) == 2)))
})

# Request more than possible combinations -- stratified
test_that("getCombinations returns all possible stratified combinations when requested number is too large", {
  expect_message(output <- getCombinations(n_replicates = c(3, 3),
                                           n_group1 = c(1, 1),
                                           n_combinations = 100),
                 "exceeds the number of possible combinations")
  expect_true(is.matrix(output))
  expect_equal(nrow(output), 2)
  expect_equal(ncol(output), choose(3, 1) * choose(3, 1))
})
