#' @export
box::use(./h5ad2SE[...],
         lle[lle])
box::use(./setup_SLICER[...],
         tibble[tibble, as_tibble],
         ggplot2[...],
         DESeq2[...],
         stats[relevel])


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
                      knn = 5,
                      start = NULL) {

  counts <- t(assay(data, "counts"))
  if (is.null(gene_subset)) {
    cat("selecting genes...\n")
    gene_subset <- SLICER::select_genes(counts)
  } else {
    stopifnot("`gene_subset` should be `NULL` or an integer vector"=is_int(gene_subset))
  }
  traj <- counts[,gene_subset]
  if (is.null(select_k)) {
    cat("finding optimal k...\n")
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
  if (is.null(start)) {
    ends <- find_extreme_cells(traj_graph, traj_lle)
    start <- ends[1]
  } else {
    stopifnot("`start` should be a scalar integer vector"=is_int(knn)&&(length(knn)==1))
  }

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
  attr(out, "traj_graph") <- traj_graph
  attr(out, "start_node") <- start
  attr(out, "gene_subset") <- gene_subset
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
get_paths <- function(traj_df) {
  graph <- attr(traj_df, "traj_graph")
  start <- attr(traj_df, "start_node")
  vpaths <- igraph::shortest_paths(graph,
                         from = start,
                         to = igraph::V(graph))$vpath
  branches <- traj_df$branch
  n_branches <- length(levels(branches))
  paths <- lapply(vpaths,
                  function(i, x) rle(as.integer(x[i]))$values,
                  x = branches)
  ##
  upaths <- unique(paths)
  len <- vapply(upaths, length, 1L)
  upaths <- upaths[order(len)]
  n <- length(upaths)
  is_contained <- logical(n)

  for (i in seq_along(upaths)[-n]) {
    cur_path <- upaths[[i]]
    cur_len <- length(cur_path)
    cur_seq <- 1:cur_len
    for (j in (i+1):n) {
      other_path <- upaths[[j]]
      comp <- other_path[cur_seq] == cur_path
      if (all(comp)) {
        is_contained[i] <- TRUE
        break
      }
    }
  }
  p <- upaths[!is_contained]

  g <- igraph::make_empty_graph(n = n_branches) |>
    igraph::add_edges(
      lapply(p, as_edge_sequence) |>
        unlist(recursive = F) |>
        unique() |>
        unlist()
    )
  attr(p, "branch_graph") <- g
  attr(p, "traj_graph") <- graph
  class(p) <- "branch_paths"
  p
}

#' @export
as_edge_sequence <- function(x) {
  ln <- length(x)

  lapply(seq_along(x)[-ln],
         function(i, v) v[i:(i+1)],
         v = x)
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

get_unique_comp <- function(paths) {
  edges <- lapply(paths, as_edge_sequence) |>
    unlist(recursive = F) |>
    unique()
  ln <- length(edges)
  dup <- logical(ln)
  for (i in seq_along(edges)[-ln]) {
    cur_edge <- edges[[i]]
    for (j in (i+1):ln) {
      oth_edge <- edges[[j]]
      if (all(cur_edge %in% oth_edge)) {
        dup[i] <- TRUE
        break
      }
    }
  }
  edges[!dup]
}

#' @export
do_deseq_along <- function(data, traj_df) {
  data$branch <- traj_df$branch
  paths <- get_paths(traj_df)
  comps <- get_unique_comp(paths)
  ln <- length(comps)
  out <- vector("list", ln)
  for (i in seq_along(out)) {
    to_comp <- comps[[i]]
    data_sub <- data[,data$branch %in% to_comp]
    data_sub$branch <- factor(relevel(data_sub$branch, ref = to_comp[1]))
    cat(srpintf("starting comparison %i of %i\n", i, ln))
    out[[i]] <- do_deseq(data_sub, ~ branch)
  }
  against <- vapply(out, attr, FUN.VALUE = character(1), which = "against")
  out <- unlist(out)
  class(out) <- "deseq_along"
  attr(out, "comps") <- vapply(comps, function(x) sprintf("%s/%s", x[2], x[1]), character(1))
  out
}

#' @export
plot_deseq_along <- function(deseq_res, log2FC = 2, padj = 0.05) {
  comps <- attr(deseq_res, "comps")
  data <- purrr::map2(deseq_res,
                      comps,
                 function(x, a) {
                   group <- attr(x, "group")
                   data <- tibble::as_tibble(x, rownames = "genes")
                   data$comp <- a
                   data
                 }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      padj = {padj[is.na(padj)] <- 1; padj},
      type = dplyr::case_when(
        abs(.data$log2FoldChange)>.env$log2FC & .data$padj < .env$padj ~ "Significant",
        TRUE ~ ""
      ),
      comp = factor(comp, levels = .env$comps)
    )
  ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(data = ~subset(.x, type==""),
               color = "grey",
               alpha = .3) +
    geom_point(data = ~subset(.x, type!=""),
               aes(color = type)) +
    theme_bw() +
    facet_wrap(~comp) +
    labs(title = "DE along trajectories",
         caption = sprintf("Significance was determined if a gene had abs(log2FoldChange)>%s and padj<%s",
                           format(log2FC), format(padj))) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = "darkred")
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
