
# `h5ad2SE`= h5ad --> SummarizedExperiment (SE)
#
#  Intended to read a h5ad object downloaded from
#  [here](https://zenodo.org/records/6587903) and
#  convert it into something usable in R.
#
#
#

#' @export
box::use(SummarizedExperiment[...],
         Matrix[...])

box::use(reticulate[import],
         methods[as])

# import scanpy
scanpy <- import("scanpy")

#' @export
h5ad2SE <- function(file) {

  # read the file
  data <- scanpy$read_h5ad(file)

  col_data <- data$obs
  row_data <- data$var
  mat <- data$layers[["matrix"]]

  # construct SE object
  SummarizedExperiment(
    #traditionally, the feature are the rows
    # in SE objects
    list(counts = t(mat)),
    rowData = as(row_data, "DataFrame"),
    colData = as(col_data, "DataFrame")
    )

}
