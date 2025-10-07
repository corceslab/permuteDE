<br>
<a href ="<LINK>"><img src="man/figures/logo.png" width="200px" align="right" /></a>

<!-- badges: start -->
<!-- badges: end -->

**permuteDE** is an R package intended to help users assess which differential expression comparisons can be trusted and which should be met with skepticism. 

Differential expression analyses are susceptible to false positives, and it can be difficult to prioritize the most robust results for further study. permuteDE uses permutation testing to identify which comparisons have a higher number of significant differentially expressed features than would be expected by chance.

<br>

## Installation

permuteDE is designed to be run on Unix-based operating systems such as macOS and linux.

permuteDE installation currently requires `remotes` and `BiocManager` for installation of GitHub and Bioconductor packages. Run the following commands to install the various dependencies used by permuteDE:

First, install remotes (for installing GitHub packages) if it isn’t already installed:
```
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
```

Then, install BiocManager (for installing bioconductor packages) if it isn’t already installed:
```
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```

Then, install permuteDE:
```
remotes::install_github("corceslab/permuteDE", ref="main", repos = BiocManager::repositories(), upgrade = "never")
```

Notes:

* Installation should complete in under 2 minutes.
* This package is supported for macOS and Linux. 
* Package dependencies can be found in the "DESCRIPTION" file.

## Usage

Please see the [tutorial](<LINK>). It takes less than 10 minutes to run on a standard laptop.

## How permuteDE works

permuteDE takes as input a Seurat object, SingleCellExperiment object, or matrix. permuteDE enables the user to run several different methods for differential expression analysis using the `runDE()` function. We encourage users to run pseudobulk differential expression methods to decrease the number of likely false positives, but permuteDE is also compatible with cell-level tests. 

These results are then passed to the `permuteDE()` function. Across a number of iterations, permuteDE shuffles the group labels for a given set of replicates and then performs pseudobulk differential expression between those permuted groups. The metric used by the permutation test is the number of significantly differentially expressed features for the true group comparison vs. the comparisons that use permuted group labels. 

See the [tutorial](<LINK>) for further details.

<hr>

<p align="left"><a href ="https://www.corceslab.com/"><img src="man/figures/CorcesLab_logo.png" alt="" width="300"></a></p>

permuteDE is developed and maintained by the [Corces Lab](https://www.corceslab.com/) at the Gladstone Institutes.
