# README


This is the landing page for COMP683 Course Project.

## Getting setup

``` r
# unless you are starting an R session
# from within this directory, source
# the `.Rprofile` from the project root.
source(".Rprofile")

# If running for the first time use the code below to setup SLICER by installing dependencies
# ## SLICER and some other dependencies are no longer on CRAN, thus they must be installed manually
box::use(mods/setup_SLICER[...])

# If you have run the above, future work may be called by doing the above again or the below:
library(SLICER) # traditional
box::use(SLICER[...]) # using `box::use()`
```

# Project Proposal

<!-- include Quarto doc in README -->

## DE with SLICER

> **Note**
>
> This is a work in progress

## Group Members

-   Justin Landis
-   Victor Adediwura

## Abstract

SLICER is a method to select features (genes) to build a trajectory of
cells. In Single Cell RNAseq (scRNAseq), this method may be helpful in
the context of cell differentiation analyses. The goal of our project is
to investigate Differential Expression (DE) approaches to the features
selected by SLICER.

## Formal Statement of the Problem

While SLICER automatically selects genes that are important for defining
a trajectory among the data, it does not associate which features are
most important to defined cell types. <!-- Not sure if this is true -->

## Related Work

TBD

## Contributions

## Datasets

We will be using data sets from [Single-cell dattasets for temporal gene
expression integration](https://zenodo.org/records/6587903),
specifically utilizing a few Hematopoiesis differentiation dataset (as
there are 2).

## Intended Experiments

TBD

## Expected Challenges

Immediate challenges will be the disparity between softwares. Data is
stored in a `h5ad` format that can be read into memory via
`scanpy.read_h5ad(...)`. However `SLICER` is implemented in R and will
need to be locally installed as it was [removed from
CRAN](https://cran.r-project.org/web/packages/SLICER/index.html) in
`2022`.

## Implementation

Since `SLICER` is implemented in R, we will be implementing our DE in R
as well. Our code will be posted on
[GitHub](https://github.com/jtlandis/Comp683-Proj)

## Preliminary Results

TBD
