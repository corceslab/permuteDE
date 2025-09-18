
# permuteDE

<!-- badges: start -->
<!-- badges: end -->

This package performs permutation testing for differential expression analysis by shuffling the group labels for a given set of replicates and then performing pseudobulk differential expression between those permuted groups. This can help to assess how many false positive significant differentially expressed features can be expected by chance. The metric used by the permutation test is the number of significantly differentially expressed features for the true group comparison vs. the comparisons that use permuted group labels.

## Installation

Because this is a private repository, you will need to create an access token at: github.com/settings/tokens

Then provide your GitHub credentials prior to installing the package, using the following command:

``` r
credentials::set_github_pat()
```

You can then install the package with the following command:

``` r
remotes::install_github("cathrinesant/permuteDE")
```
