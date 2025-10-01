# ---------------------------------------------------------------------------
# Helper Functions
# ---------------------------------------------------------------------------

# Retrieve data from an object
#
# object A 'Seurat' or 'SingleCellExperiment' object
# type A string indicating data type to retrieve
# name A string under which data are stored
# use_cells (Optional) A string vector of cell names to subset
.retrieveData <- function(object, type, name, use_cells = NULL) {
  if (methods::is(object, "Seurat")) {
    # for a Seurat object
    if (type == "cell_metadata") {
      output_data <- object@meta.data[if (is.null(use_cells)) TRUE else use_cells, name]
    }
  } else if (methods::is(object, "SingleCellExperiment")) {
    # for a SingleCellExperiment object
    if (type == "cell_metadata") {
      output_data <- object@colData[if (is.null(use_cells)) TRUE else use_cells, name]
    }
  }
  return(output_data)
}

# Retrieve cell IDs ---------------------------
#
# Extract cell IDs/names from provided object
#
# object -- An object of class Seurat or SingleCellExperiment
# use_assay -- For Seurat objects, character string/vector indicating assay to use
.getCellIDs <- function(object,
                        use_assay = NULL) {
  # By object type
  if (methods::is(object, "Seurat")) {
    if (is.null(use_assay)) {
      use_assay <- Seurat::DefaultAssay(object)
    }
    if (length(use_assay) > 1) {
      # If multiple assays, check that cell IDs are identical
      cell_IDs <- colnames(object[[use_assay[1]]])
      for (i in 2:length(use_assay)) {
        cell_IDs_i <- colnames(object[[use_assay[i]]])
        if (!identical(cell_IDs, cell_IDs_i)) {
          stop("Cell IDs do not match across provided assays indicated by parameter 'use_assay'. Please supply valid input!")
        }
      }
    } else {
      cell_IDs <- colnames(object[[use_assay]])
    }
  } else if (methods::is(object, "SingleCellExperiment")) {
    cell_IDs <- rownames(object@colData)
  }
  return(cell_IDs)
}

# Retrieve matrix ---------------------------
#
# Extract a matrix from provided object
#
# object -- An object of class Seurat or SingleCellExperiment
# use_matrix -- If there is a user-supplied matrix, do not retrieve a matrix from the object
# use_assay -- For Seurat or SingleCellExperiment objects, a character string indicating the assay to use
# use_layer -- For Seurat objects, a character string indicating the layer/slot to use
# use_features -- A vector of feature names to use to subset the matrix
# exclude_features -- A vector of feature names to exclude from the matrix
# use_cells -- A vector of cell IDs to use to subset the matrix
# verbose -- A boolean value indicating whether to use verbose output during the execution of this function
.getMatrix <- function(object = NULL,
                       use_matrix = NULL,
                       use_assay = NULL,
                       use_layer = NULL,
                       use_features = NULL,
                       exclude_features = NULL,
                       use_cells = NULL,
                       verbose) {
  # By object type
  if (is.null(object) | methods::is(object, "Seurat")) {
    use_matrix <- .getMatrix.Seurat(object,
                                    use_matrix,
                                    use_assay,
                                    use_layer,
                                    use_features,
                                    exclude_features,
                                    use_cells,
                                    verbose)
  } else if (methods::is(object, "SingleCellExperiment")) {
    use_matrix <- .getMatrix.SingleCellExperiment(object,
                                                  use_matrix,
                                                  use_assay,
                                                  use_layer,
                                                  use_features,
                                                  exclude_features,
                                                  use_cells,
                                                  verbose)
  }
  # Return matrix
  return(use_matrix)
}

# Extract a matrix from a Seurat object
#
# object -- A Seurat object
# use_matrix -- If there is a user-supplied matrix, do not retrieve a matrix from the object
# use_assay -- a string indicating the assay to use
# use_layer -- a string indicating the layer/slot to use
# use_features -- A vector of feature names to use to subset the matrix
# exclude_features -- A vector of feature names to exclude from the matrix
# use_cells -- A vector of cell IDs to use to subset the matrix
# verbose -- A boolean value indicating whether to use verbose output during the executio
.getMatrix.Seurat <- function(object,
                              use_matrix,
                              use_assay,
                              use_layer,
                              use_features,
                              exclude_features,
                              use_cells,
                              verbose) {
  # If matrix is not provided as input
  if (is.null(use_matrix)) {
    # Get assay
    if (is.null(use_assay)) {
      use_assay <- Seurat::DefaultAssay(object)
    } else {
      # Check that input assay is present in object
      .validInput(use_assay, "use_assay", list(object, FALSE, NULL))
    }
    # Determine which layer to use
    if (is.null(use_layer)) {
      use_layer <- "counts"
    }
    # Check that selected layer is present within selected assay in object
    .validInput(use_layer, "use_layer", list(object, use_assay, FALSE, NULL))
    # Extract matrix
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Fetching feature x cell matrix using '", use_assay, "' assay and '", use_layer, "' layer..")
    if ("Assay5" %in% methods::is(object[[use_assay]])) {
      use_matrix <- object[[use_assay]]@layers[[use_layer]]
      colnames(use_matrix) <- colnames(object[[use_assay]])
      rownames(use_matrix) <- rownames(object[[use_assay]])
    } else {
      use_matrix <- methods::slot(object[[use_assay]], name = use_layer)
    }
  } else {
    # If use_assay is not NULL
    if (!is.null(use_assay)) {
      if (verbose) warning("Input for parameter 'use_assay' is not used when a matrix is provided for parameter 'use_matrix'.")
    }
    # If use_layer is not NULL
    if (!is.null(use_layer)) {
      if (verbose) warning("Input for parameter 'use_layer' is not used when a matrix is provided for parameter 'use_matrix'.")
    }
  }

  # If matrix has no row names
  if (is.null(rownames(use_matrix))) {
    # Stop if trying to subset features
    if (!is.null(use_features) | !is.null(exclude_features)) {
      stop("Provided 'use_matrix' has no row names, therefore, input for parameters 'use_features' and 'exclude_features' cannot be used.")
    }
    rownames(use_matrix) <- seq(1, nrow(use_matrix))
  }
  # Subset matrix by selected features
  if (!is.null(use_features)) {
    unused_features <- use_features[!(use_features %in% rownames(use_matrix))]
    if (verbose & (length(unused_features) > 0)) {
      warning("Could not find the following ",
              length(unused_features),
              " features provided by 'use_features' in 'use_matrix': \n",
              unused_features)
    }
  } else {
    use_features <- rownames(use_matrix)
  }
  use_features <- use_features[!(use_features %in% exclude_features)]
  if (length(use_features) < 1) {
    stop("No remaining features in matrix. Please check input to 'use_features' and/or 'exclude_features'!")
  }

  # If matrix has no column names
  if (is.null(colnames(use_matrix))) {
    # Stop if trying to subset cells
    if (!is.null(use_cells)) {
      stop("Provided 'use_matrix' has no column names, therefore, input for parameter 'use_cells' cannot be used.")
    }
    colnames(use_matrix) <- seq(1, ncol(use_matrix))
  }
  # Subset matrix by selected cells
  if (!is.null(use_cells)) {
    unused_cells <- use_cells[!(use_cells %in% colnames(use_matrix))]
    if (verbose & (length(unused_cells) > 0)) {
      warning("Could not find the following ",
              length(unused_cells),
              " cells provided by 'use_cells' in 'use_matrix': \n",
              unused_cells)
    }
  } else {
    use_cells <- colnames(use_matrix)
  }

  use_matrix <- use_matrix[use_features, use_cells]
  return(use_matrix)
}

.getMatrix.SingleCellExperiment <- function(object,
                                            use_matrix,
                                            use_assay,
                                            use_layer,
                                            use_features,
                                            exclude_features,
                                            use_cells,
                                            verbose) {
  # If matrix is not provided as input
  if (is.null(use_matrix)) {
    # Get assay
    if (is.null(use_assay)) {
      use_assay <- "counts"
      .validInput(use_assay, "use_assay", list(object, FALSE, NULL))
    }
    use_matrix <- object@assays@data[[use_assay]]
    if (verbose) message(format(Sys.time(), "%Y-%m-%d %X"), " : Fetching feature x cell matrix using '", use_assay, "' assay..")
  } else {
    # If assay is not NULL
    if (!is.null(use_assay)) {
      if (verbose) warning("Input for parameter 'use_assay' is not used when a matrix is provided for parameter 'use_matrix'.")
    }
  }
  # If use_layer is not NULL
  if (!is.null(use_layer)) {
    if (verbose) warning("Input for parameter 'use_layer' is not used when input for parameter 'object' is of type 'SingleCellExperiment'.")
  }

  # If matrix has no row names
  if (is.null(rownames(use_matrix))) {
    # Stop if trying to subset features
    if (!is.null(use_features) | !is.null(exclude_features)) {
      stop("Provided 'use_matrix' has no row names, therefore, input for parameters 'use_features' and 'exclude_features' cannot be used.")
    }
    rownames(use_matrix) <- seq(1, nrow(use_matrix))
  }
  # Subset matrix by selected features
  if (!is.null(use_features)) {
    unused_features <- use_features[!(use_features %in% rownames(use_matrix))]
    if (verbose & (length(unused_features) > 0)) {
      warning("Could not find the following ",
              length(unused_features),
              " features provided by 'use_features' in 'use_matrix': \n",
              unused_features)
    }
  } else {
    use_features <- rownames(use_matrix)
  }
  use_features <- use_features[!(use_features %in% exclude_features)]
  if (length(use_features) < 1) {
    stop("No remaining features in matrix. Please check input to 'use_features' and/or 'exclude_features'!")
  }

  # If matrix has no column names
  if (is.null(colnames(use_matrix))) {
    # Stop if trying to subset cells
    if (!is.null(use_cells)) {
      stop("Provided 'use_matrix' has no column names, therefore, input for parameter 'use_cells' cannot be used.")
    }
    colnames(use_matrix) <- seq(1, ncol(use_matrix))
  }
  # Subset matrix by selected cells
  if (!is.null(use_cells)) {
    unused_cells <- use_cells[!(use_cells %in% colnames(use_matrix))]
    if (verbose & (length(unused_cells) > 0)) {
      warning("Could not find the following ",
              length(unused_cells),
              " cells provided by 'use_cells' in 'use_matrix': \n",
              unused_cells)
    }
  } else {
    use_cells <- colnames(use_matrix)
  }

  use_matrix <- use_matrix[use_features, use_cells]
  return(use_matrix)
}

# Require package ---------------------------
#
# Load required package
# Adapted from ArchR code, Jeffrey Granja & Ryan Corces
#
# x -- Name of package
# load -- Whether to load package
# installInfo -- Installation info
# source -- cran/bioc, etc.
.requirePackage <- function(x = NULL, load = TRUE, installInfo = NULL, source = NULL){
  if(x %in% rownames(utils::installed.packages())){
    if(load){
      suppressPackageStartupMessages(require(x, character.only = TRUE))
    }else{
      return(0)
    }
  }else{
    if (!is.null(source) & is.null(installInfo)) {
      if (tolower(source) == "cran") {
        installInfo <- paste0('install.packages("',x,'")')
      } else if (tolower(source) == "bioc"){
        installInfo <- paste0('BiocManager::install("',x,'")')
      } else {
        stop("Unrecognized package source, available are cran/bioc!")
      }
    }
    if (!is.null(installInfo)) {
      stop(paste0("Required package : ", x, " is not installed/found!\n  Package Can Be Installed : ", installInfo))
    } else {
      stop(paste0("Required package : ", x, " is not installed/found!"))
    }
  }
}

# Startup ---------------------------
#
# Adapted from ArchR code, Jeffrey Granja & Ryan Corces

.onAttach <- function(libname, pkgname){
  # ASCII permuteDE logo
  packageStartupMessage("                                                         _       ___   _____
 ========.  .=======!\\   _ _     __ _ __ _ __ ___  _   _| |_ ___|  _ \\| ____|
 --,-,--. \\/ .-,-,--|/  | '  \\ / = \\ '__| '_ ` _ \\| | | | __/ = \\ | | |  _|
   ! !  /\\ \\/  ! !      | |)  |  __/ |  | | | | | | |_| | ||  __/ |_| | |___
 ======' /\\ `=======!\\  | .__/ \\___|_|  |_| |_| |_|\\____|\\__\\___|____/|_____|
 -------'  `--------|/  |_|
 ")
  # package startup
  v <- utils::packageVersion("permuteDE")
  packageStartupMessage("permuteDE : Version ", v,
                        "\nIf you encounter a bug please report : https://github.com/CorcesLab/permuteDE/issues")
}
