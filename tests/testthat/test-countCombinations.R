# ---------------------------------------------------------------------------
# Tests for function countCombinations
# ---------------------------------------------------------------------------

# Numeric input
test_that("countCombinations works with numeric n_replicates and n_group1", {
  expect_equal(countCombinations(n_replicates = 4, n_group1 = 2), 6)
  expect_equal(countCombinations(n_replicates = 5, n_group1 = 2), 10)
  expect_equal(countCombinations(n_replicates = 6, n_group1 = 3), 20)
})

# Vector input
test_that("countCombinations works with replicate and group label vectors", {
  replicate_labels <- c("rep1", "rep2", "rep3", "rep4")
  group_labels <- c("A", "A", "B", "B")
  expect_equal(countCombinations(replicate_labels = replicate_labels,
                                 group_labels = group_labels),
               6)
})

# No input
test_that("countCombinations requires either numeric inputs or label inputs", {
  expect_error(countCombinations(),
               "input must be provided to 'replicate_labels' and 'group_labels'")
})

# Different length inputs
test_that("countCombinations requires replicate_labels and group_labels to have the same length", {
  expect_error(countCombinations(replicate_labels = c("rep1", "rep2", "rep3"),
                                 group_labels = c("A", "B")),
               "must be of the same length")})

# Number of groups
test_that("countCombinations requires exactly two groups", {
  expect_error(countCombinations(replicate_labels = c("rep1", "rep2", "rep3"),
                                 group_labels = c("A", "A", "A")),
               "contains exactly 2 groups")

  expect_error(countCombinations(replicate_labels = c("rep1", "rep2", "rep3"),
                                 group_labels = c("A", "B", "C")),
               "contains exactly 2 groups")
})

# NA inputs
test_that("countCombinations rejects NA replicate or group labels", {
  expect_error(countCombinations(replicate_labels = c("rep1", NA, "rep3", "rep4"),
                                 group_labels = c("A", "A", "B", "B")),
               "replicate_labels.*cannot be NA")

  expect_error(countCombinations(replicate_labels = c("rep1", "rep2", "rep3", "rep4"),
                                 group_labels = c("A", NA, "B", "B")),
               "group_labels.*cannot be NA")
})

# Warning when both numeric and label inputs provided
test_that("countCombinations warns when numeric inputs override label inputs", {
  expect_warning(
    output <- countCombinations(
      replicate_labels = c("rep1", "rep2", "rep3", "rep4"),
      group_labels = c("A", "A", "B", "B"),
      n_replicates = 5,
      n_group1 = 2
    ), "not used")
  expect_equal(output, 10)
})
