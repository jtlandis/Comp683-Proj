

select_genes2 <- function(embedding) {
  browser()
  k = SLICER:::min_conn_k(embedding)
  n = nrow(embedding)
  m = ncol(embedding)
  traj_dist = as.matrix(dist(embedding))
  adj_mat = adaptive_knn_graph2(traj_dist, rep(k, n))
  sel_vals = sapply(1:m, selection_val, embedding, adj_mat)
  genes = which(sel_vals > 1)
  return(genes)
}

adaptive_knn_graph2 <- function (traj_dist, k)
{
  adj_mat = matrix(0, nrow = nrow(traj_dist), ncol = ncol(traj_dist))
  knn = t(apply(traj_dist, 1, order))
  for (i in 1:nrow(traj_dist)) {
    adj_mat[i, knn[i, 2:(k[i] + 1)]] = traj_dist[i, knn[i,
                                                        2:(k[i] + 1)]]
  }
  return(adj_mat)
}

# copy_fn <- function(fun) {
#   rlang::new_function(
#     args = formals(fun),
#     body = body(fun),
#     env = environment(fun)
#   )
# }
#
# geo_entro <- copy_fn(SLICER::assign_branches) |>
#   ggside:::mod_fun_at(quote(browser()), 1)

# This is from [SLICER Github](https://github.com/jw156605/SLICER)
box::use(./setup_SLICER[...],
         lle[lle])
genes <- select_genes(traj)
k <- 15 #select_k(traj[,genes], kmin=5) --> 15
traj_lle <- lle(traj[,genes], m=2, k)$Y
traj_graph <- conn_knn_graph(traj_lle,5)
ends <- find_extreme_cells(traj_graph, traj_lle)
start <- 1
cells_ordered <- cell_order(traj_graph, start)
branches <- assign_branches(traj_graph,start)

# taken from SLICER::compute_geodesic_entropy
#  this function omits automatic plotting.
#  Only here to inspect

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
         ggplot2[...],
         gganimate[...])

p <- tibble(x = traj_lle[,1],
       y = traj_lle[,2],
       branch = factor(branches),
       order = cells_ordered) |>
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

animate(anim, duration = 20, fps = 30)
