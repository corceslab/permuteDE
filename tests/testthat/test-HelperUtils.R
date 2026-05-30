# ---------------------------------------------------------------------------
# Tests for helper functions
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# .retrieveData
# ---------------------------------------------------------------------------

test_that(".retrieveData retrieves metadata column from data.frame", {
  metadata <- data.frame(replicate = c("rep1", "rep2", "rep3"),
                         group = c("A", "A", "B"),
                         row.names = c("cell1", "cell2", "cell3"))
  output <- .retrieveData(object = NULL,
                          metadata = metadata,
                          type = "cell_metadata",
                          name = "replicate")
  expect_equal(output, c("rep1", "rep2", "rep3"))
})

test_that(".retrieveData subsets metadata by use_cells", {
  metadata <- data.frame(replicate = c("rep1", "rep2", "rep3"),
                         group = c("A", "A", "B"),
                         row.names = c("cell1", "cell2", "cell3"))
  output <- .retrieveData(object = NULL,
                          metadata = metadata,
                          type = "cell_metadata",
                          name = "group",
                          use_cells = c("cell1", "cell3"))
  expect_equal(output, c("A", "B"))
})

test_that(".retrieveData errors if object and metadata are both NULL", {
  expect_error(.retrieveData(object = NULL,
                             metadata = NULL,
                             type = "cell_metadata",
                             name = "group"),
               "both 'object' and 'metadata' are NULL")
})

# ---------------------------------------------------------------------------
# .getCellIDs
# ---------------------------------------------------------------------------

test_that(".getCellIDs retrieves column names from matrix input", {
  mat <- matrix(1:6,
                nrow = 2,
                dimnames = list(c("gene1", "gene2"),
                                c("cell1", "cell2", "cell3")))
  output <- .getCellIDs(mat)
  expect_equal(output, c("cell1", "cell2", "cell3"))
})

test_that(".getCellIDs works with SingleCellExperiment objects", {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  counts <- matrix(1:6,
                   nrow = 2,
                   dimnames = list(c("gene1", "gene2"),
                                   c("cell1", "cell2", "cell3")))
  suppressWarnings(object <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts)))
  output <- .getCellIDs(object)
  expect_equal(output, c("cell1", "cell2", "cell3"))
})

test_that(".getCellIDs works with Seurat objects", {
  testthat::skip_if_not_installed("Seurat")
  counts <- matrix(1:6,
                   nrow = 2,
                   dimnames = list(c("gene1", "gene2"),
                                   c("cell1", "cell2", "cell3")))
  suppressWarnings(object <- Seurat::CreateSeuratObject(counts = counts))
  output <- .getCellIDs(object)
  expect_equal(output, c("cell1", "cell2", "cell3"))
})

# ---------------------------------------------------------------------------
# .getMatrix
# ---------------------------------------------------------------------------

# Matrix objects
test_that(".getMatrix.matrix returns the full matrix when no subsets are supplied", {
  mat <- matrix(1:9,
                nrow = 3,
                dimnames = list(c("gene1", "gene2", "gene3"),
                                c("cell1", "cell2", "cell3")))
  output <- .getMatrix(object = mat,
                       verbose = TRUE)
  expect_equal(output, mat)
})

test_that(".getMatrix.matrix subsets features and cells", {
  mat <- matrix(1:12,
                nrow = 3,
                dimnames = list(c("gene1", "gene2", "gene3"),
                                c("cell1", "cell2", "cell3", "cell4")))
  output <- .getMatrix(object = mat,
                       use_matrix = NULL,
                       use_features = c("gene1", "gene3"),
                       exclude_features = NULL,
                       use_cells = c("cell2", "cell4"),
                       verbose = TRUE)
  expect_equal(rownames(output), c("gene1", "gene3"))
  expect_equal(colnames(output), c("cell2", "cell4"))
  expect_equal(output, mat[c("gene1", "gene3"), c("cell2", "cell4")])
})

test_that(".getMatrix.matrix excludes requested features", {
  mat <- matrix(1:12,
                nrow = 3,
                dimnames = list(c("gene1", "gene2", "gene3"),
                                c("cell1", "cell2", "cell3", "cell4")))
  output <- .getMatrix(object = mat,
                       exclude_features = "gene2",
                       verbose = TRUE)
  expect_equal(rownames(output), c("gene1", "gene3"))
  expect_false("gene2" %in% rownames(output))
})

test_that(".getMatrix.matrix applies use_features before exclude_features", {
  mat <- matrix(1:12,
                nrow = 3,
                dimnames = list(c("gene1", "gene2", "gene3"),
                                c("cell1", "cell2", "cell3", "cell4")))
  output <- .getMatrix(object = mat,
                       use_features = c("gene1", "gene2"),
                       exclude_features = "gene2",
                       verbose = TRUE)
  expect_equal(rownames(output), "gene1")
})

test_that(".getMatrix.matrix errors when no features remain", {
  mat <- matrix(1:6,
                nrow = 2,
                dimnames = list(c("gene1", "gene2"),
                                c("cell1", "cell2", "cell3")))
  expect_error(.getMatrix(object = mat,
                          exclude_features = c("gene1", "gene2"),
                          verbose = TRUE),
               "No remaining features")
})

test_that(".getMatrix.matrix warns about missing requested features", {
  mat <- matrix(1:6,
                nrow = 2,
                dimnames = list(c("gene1", "gene2"),
                                c("cell1", "cell2", "cell3")))
  expect_warning(output <- .getMatrix(object = mat,
                                      use_features = c("gene1", "missing_gene"),
                                      verbose = TRUE),
                 "Could not find the following")
  expect_equal(rownames(output), "gene1")
})

test_that(".getMatrix.matrix warns about missing requested cells", {
  mat <- matrix(1:6,nrow = 2,
                dimnames = list(c("gene1", "gene2"),
                                c("cell1", "cell2", "cell3")))
  expect_warning(output <- .getMatrix(object = mat,
                                      use_cells = c("cell1", "missing_cell"),
                                      verbose = TRUE),
                 "Could not find the following")
  expect_equal(colnames(output), "cell1")
})

test_that(".getMatrix.matrix errors when matrix has no rownames and features are requested", {
  mat <- matrix(1:6, nrow = 2)
  colnames(mat) <- c("cell1", "cell2", "cell3")
  expect_error(.getMatrix(object = mat,
                          use_features = "gene1",
                          verbose = TRUE),
               "has no row names")
})

test_that(".getMatrix.matrix errors when matrix has no colnames and cells are requested", {
  mat <- matrix(1:6, nrow = 2)
  rownames(mat) <- c("gene1", "gene2")
  expect_error(.getMatrix(object = mat,
                          use_cells = "cell1",
                          verbose = TRUE),
               "has no column names")
})

test_that(".getMatrix.matrix adds row and column names when missing and no subsetting is requested", {
  mat <- matrix(1:6, nrow = 2)
  output <- .getMatrix(object = mat,
                       verbose = TRUE)
  expect_equal(rownames(output), c("1", "2"))
  expect_equal(colnames(output), c("1", "2", "3"))
})

test_that(".getMatrix.matrix dispatches to matrix helper for matrix input", {
  mat <- matrix(1:6,
                nrow = 2,
                dimnames = list(c("gene1", "gene2"),
                                c("cell1", "cell2", "cell3")))
  output <- .getMatrix(object = mat,
                       use_features = "gene2",
                       use_cells = "cell3",
                       verbose = TRUE)
  expect_equal(output, mat["gene2", "cell3", drop = FALSE])
})

# SingleCellExperiment object
test_that(".getMatrix.SingleCellExperiment retrieves and subsets assay matrix", {
  testthat::skip_if_not_installed("SingleCellExperiment")
  testthat::skip_if_not_installed("SummarizedExperiment")
  counts <- matrix(1:12,
                   nrow = 3,
                   dimnames = list(c("gene1", "gene2", "gene3"),
                                   c("cell1", "cell2", "cell3", "cell4")))
  suppressWarnings(object <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts)))
  output <- .getMatrix(object = object,
                       use_assay = "counts",
                       use_features = c("gene1", "gene3"),
                       use_cells = c("cell2", "cell4"),
                       verbose = TRUE)
  expect_equal(output, counts[c("gene1", "gene3"), c("cell2", "cell4")])
})

# Seurat object
test_that(".getMatrix.Seurat retrieves and subsets assay matrix", {
  testthat::skip_if_not_installed("Seurat")
  counts <- matrix(1:12,
                   nrow = 3,
                   dimnames = list(c("gene1", "gene2", "gene3"),
                                   c("cell1", "cell2", "cell3", "cell4")))
  suppressWarnings(object <- Seurat::CreateSeuratObject(counts = counts))
  output <- .getMatrix(object = object,
                       use_features = c("gene1", "gene3"),
                       use_cells = c("cell2", "cell4"),
                       verbose = TRUE)
  expect_equal(as.matrix(output),
               counts[c("gene1", "gene3"), c("cell2", "cell4")])
})

# ---------------------------------------------------------------------------
# .requirePackage
# ---------------------------------------------------------------------------

test_that(".requirePackage returns 0 for installed package when load is FALSE", {
  output <- .requirePackage("utils", load = FALSE)
  expect_equal(output, 0)
})

test_that(".requirePackage errors clearly for missing package", {
  expect_error(.requirePackage("definitelyNotARealPackage123", load = FALSE),
               "Required package")
})

test_that(".requirePackage includes CRAN install instructions when source is cran", {
  expect_error(.requirePackage("definitelyNotARealPackage123",
                               load = FALSE,
                               source = "cran"),
               'install.packages\\("definitelyNotARealPackage123"\\)')
})

test_that(".requirePackage includes Bioconductor install instructions when source is bioc", {
  expect_error(.requirePackage("definitelyNotARealPackage123",
                               load = FALSE,
                               source = "bioc"),
               'BiocManager::install\\("definitelyNotARealPackage123"\\)')
})

test_that(".requirePackage errors for unrecognized package source", {
  expect_error(.requirePackage(
    "definitelyNotARealPackage123",
    load = FALSE,
    source = "github"),
    "Unrecognized package source")
})
