mitology R package
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src=./vignettes/figures/puzzle_mitology.png width="200" align="right" />

Mitology is an R package to dissect mitochondrial activity from gene
expression data. We provide in R environment MitoCarta3.0 pathways and
we exploited Reactome pathways and GO hierarchies to derive and
re-organize mitochondrial-specific gene sets.

## Installation

The `mitology` package can be installed from Bioconductor:

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("mitology")
```

Also you can install the latest version of `mitology` from GitHub with:

``` r
library(remotes)
remotes::install_github("CaluraLab/mitology")
```
