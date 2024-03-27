
box::use(mods/h5ad2SE[...])
box::use(mods/setup_SLICER[...],
         lle[lle],
         tibble[tibble, as_tibble])

# as_seq <- function(vec) {
#   f <- factor(vec)
#   levels(f) <- seq_along(levels(f))
#   as.integer(f)
# }
# calc_vertex_freqs <- function(path_counts, path_len) {
#   n <- length(path_len)
#   vertex_freqs <- matrix(0, n, n)
#   for (i in seq_len(n)) {
#     vertex_freqs[path_len[i], i] <- path_counts[i]
#   }
#   vertex_freqs
# }
#
# calc_path_counts <- function(vpath, which_ind) {
#   counts <- integer(length(which_ind))
#   for (i in seq_along(vpath)) {
#     index <- which_ind[vpath[[i]]]
#     counts[index] <- counts[index] + 1L
#   }
#   counts
# }

data <- h5ad2SE(box::file("data/adata_olsson.h5ad"))
ptraj <- t(assay(data, "counts"))

genes <- readRDS(box::file("data/olsson_SLICER_select_genes.rds"))#select_genes(ptraj)
k <- 15 #select_k(traj[,genes], kmin=5) --> 15
traj_lle <- lle(as.matrix(ptraj[,genes]), m=2, k)$Y
traj_graph <- conn_knn_graph(traj_lle, 5)
ends <- find_extreme_cells(traj_graph, traj_lle)
start <- ends[1]
cells_ordered <-cell_order(traj_graph, start)
branches <- assign_branches(traj_graph,start)
d <- tibble(x = traj_lle[,1],
            y = traj_lle[,2],
            # z = traj_lle[,3],
            branch = factor(branches),
            order = order(cells_ordered),
            cell_type = colData(data)$cell_type) |>
  dplyr::left_join(y = cell_o, by = c("i"="cell_i"))
aricode::NMI(d$branch, d$cell_type)
table(d$branch, d$cell_type)
box::use(
  ggplot2[...],
  gganimate[...]
)
p <-  d |>
  dplyr::arrange(order) |>
  ggplot(aes(x, y)) +
  geom_point(aes(color = branch, group = order), size = 6) +
  geom_point(aes(color = branch, group = order), size = 2) +
  labs(x = "Manifold 1", y = "Manifold 2", title = "Cell Ordering") +
  theme_bw()

anim <- p + transition_components(order, enter_length = 50L, exit_length = 50L) +
  shadow_mark(exclude_layer = 1) +
  enter_grow() +
  enter_fade() +
  exit_shrink(size = .2) +
  exit_fade(alpha = .3)

anim_save(box::file("prelim_data.gif"),
          animation = anim, duration = 20, fps = 30)


box::use(DESeq2[...])
assay(data, 1) <- as.matrix(assay(data, 1))
dds <- DESeqDataSet(data, design = ~ cell_type)
dds_fit <- DESeq(dds)
res <- results(dds_fit, name = "cell_type_GMP_vs_CMP")
as_tibble(res, rownames = "gene") |>
  ggplot(aes(x = log2FoldChange, y = -log(padj))) + geom_point()
