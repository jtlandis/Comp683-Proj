
box::use(./h5ad2SE[...])
box::use(./setup_SLICER[...],
         lle[lle],
         tibble[tibble, as_tibble],
         ggplot2[...],
         DESeq2[...])


#' @export
read_h5ad <- function(path) {
  stopifnot("path should be a character"=is.character(path),
            "path should exist"=file.exists(path))
  h5ad2SE(path)
}

# check to see if `x` is like an
# integer
is_int <- function(x) {
  is.integer(x) || all((x %% 1)==0)
}

#' @param data Summarized Exeriment object
#' @param gene_subset indicies to subset `data`. If NULL, then
#' SLICER::select_genes(...) will be run
#' @param lle_dim number of dimensions to reduce to after
#' local linear embedding procedure
#' @param knn Number of nearest neighbors for trajectory
#' @export
do_SLICER <- function(data,
                      gene_subset = NULL,
                      select_k = 15,
                      lle_dim = 2,
                      knn = 5) {
  counts <- t(assay(data, "counts"))
  if (is.null(gene_subset)) {
    gene_subset <- SLICER::select_genes(counts)
  } else {
    stopifnot("`gene_subset` should be `NULL` or an integer vector"=is_int(gene_subset))
  }
  traj <- counts[,gene_subset]
  if (is.null(select_k)) {
    select_k <- SLICER::select_k(traj)
  } else {
    stopifnot("`select_k` should be `NULL` or a scalar integer vector"=is_int(select_k)&&(length(select_k)==1))
  }
  stopifnot("`lle_dim` should be a scalar integer vector"=is_int(lle_dim)&&(length(lle_dim)==1))
  stopifnot("`knn` should be a scalar integer vector"=is_int(knn)&&(length(knn)==1))
  traj_lle <- traj |>
    as.matrix() |>
    lle(m = lle_dim, select_k) |>
    _$Y
  traj_graph <- conn_knn_graph(traj_lle, knn)
  ends <- find_extreme_cells(traj_graph, traj_lle)
  start <- ends[1]
  cells_ordered <- cell_order(traj_graph, start)
  branches <- assign_branches(traj_graph, start)

  out <- as_tibble(traj_lle,
            .name_repair = ~paste0("Dim", seq_along(.x))) |>
    dplyr::mutate(
      branch = factor(.env$branches),
      order = order(.env$cells_ordered)
    ) |>
    dplyr::bind_cols(
      as_tibble(colData(data),
                rownames = ".names")
    )
  class(out) <- c("SLICER_traj_df", class(out))
  out
}

#' @export
plot_SLICER <- function(traj_df, x, y) {
  traj_df |>
    dplyr::arrange(order) |>
    ggplot(aes({{x}}, {{y}})) +
    geom_point(aes(color = branch, group = order)) +
    labs(title = "Cell Ordering") +
    theme_bw()
}

#' @export
do_deseq <- function(data, design) {
  assay(data, 1) <- as.matrix(assay(data, 1))
  dds <- DESeqDataSet(data, design = design)
  dds_fit <- DESeq(dds)
  res_names <- resultsNames(dds_fit)[-1]
  m <- regexec("^(.*)_([^_]+)_vs_(.*)$", res_names)
  m <- regmatches(res_names, m)
  col_name <- gsub("_", " ", m[[1]][2])
  against <- m[[1]][4]
  out <- mapply(function(name, group, x) {
    out <- results(x, name = name)
    attr(out, "group") <- factor(group)
    out
  },
  name = res_names,
  group = vapply(m, `[`, 3, FUN.VALUE = character(1)),
  MoreArgs = list(x = dds_fit),
  SIMPLIFY = F)
  attr(out, "against") <- against
  attr(out, "col_name") <- col_name
  out
}

#' @export
plot_deseq <- function(deseq_res, log2FC = 2, padj = 0.05) {
  against <- attr(deseq_res, "against")
  col_name <- attr(deseq_res, "col_name")
  data <- lapply(deseq_res,
                 function(x) {
                   group <- attr(x, "group")
                   data <- tibble::as_tibble(x, rownames = "genes")
                   data$group <- group
                   data
                 }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      padj = {padj[is.na(padj)] <- 1; padj},
      type = dplyr::case_when(
        abs(.data$log2FoldChange)>.env$log2FC & .data$padj < .env$padj ~ "Significant",
        TRUE ~ ""
      )
    )
  ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(data = ~subset(.x, type==""),
               color = "grey",
               alpha = .3) +
    geom_point(data = ~subset(.x, type!=""),
               aes(color = type)) +
    theme_bw() +
    facet_wrap(~group) +
    labs(title = sprintf("DE of '%s' against '%s'", col_name, against),
         caption = sprintf("Significance was determined if a gene had abs(log2FoldChange)>%s and padj<%s",
                           format(log2FC), format(padj))) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = "darkred")
}
