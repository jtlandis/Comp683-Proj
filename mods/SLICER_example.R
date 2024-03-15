
# This is from [SLICER Github](https://github.com/jw156605/SLICER)
box::use(./setup_SLICER[...],
         lle[lle])
genes <- select_genes(traj)
k <- select_k(traj[,genes], kmin=5)
traj_lle <- lle(traj[,genes], m=2, k)$Y
traj_graph <- conn_knn_graph(traj_lle,5)
ends <- find_extreme_cells(traj_graph, traj_lle)
start <- 1
cells_ordered <- cell_order(traj_graph, start)
branches <- assign_branches(traj_graph,start)

# taken from SLICER::compute_geodesic_entropy
#  this function omits automatic plotting.
geodesic_entropy <- function(traj_graph, start){
  smart_table = function(x, n) {
    tab = rep(0, n)
    tab[x] = 1
    return(tab)
  }
  n = length(igraph::V(traj_graph))
  paths = igraph::shortest_paths(traj_graph, from = start)
  path_table = sapply(paths$vpath, smart_table, n)
  path_counts = rowSums(path_table)
  vertex_freqs = matrix(0, n, n)
  for (i in 1:length(paths$vpath)) {
    for (j in 1:length(paths$vpath)) {
      if (length(paths$vpath[[j]]) >= i) {
        vertex_freqs[i, paths$vpath[[j]][i]] = vertex_freqs[i,
                                                            paths$vpath[[j]][i]] + 1
      }
    }
  }
  vertex_probs = vertex_freqs/rowSums(vertex_freqs)
  entropy_i = function(vertex_probs, i) {
    return(-sum(vertex_probs[i, which(vertex_probs[i, ] >
                                        0)] * log2(vertex_probs[i, which(vertex_probs[i,
                                        ] > 0)])))
  }
  max_path_length = max(sapply(paths$vpath, function(x) length(x)))
  entropies = rep(Inf, max_path_length)
  for (i in 1:max_path_length) {
    entropies[i] = entropy_i(vertex_probs, i)
  }
  entropies
}

compute_geodesic_entropy(traj_graph, start)
compute_geodesic_entropy(traj_graph, 496)

box::use(tibble[tibble],
         ggplot2[...])

p <- tibble(x = traj_lle[,1],
       y = traj_lle[,2],
       branch = branches,
       order = cells_ordered) |>
  dplyr::arrange(order) |>
  ggplot(aes(x, y, color = branches)) +
  geom_point() +
  geom_path()

p + gganimate::transition_reveal(order, keep_last = TRUE)
