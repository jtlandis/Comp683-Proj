
box::use(mods/h5ad2SE[...])
box::use(mods/setup_SLICER[...],
         lle[lle])

as_seq <- function(vec) {
  f <- factor(vec)
  levels(f) <- seq_along(levels(f))
  as.integer(f)
}

data <- h5ad2SE("data/adata_olsson.h5ad")
ptraj <- t(assay(data, "counts"))

genes <- select_genes(ptraj)
k <- 15 #select_k(traj[,genes], kmin=5) --> 15
traj_lle <- lle(as.matrix(ptraj[,genes]), m=6, k)$Y
traj_graph <- conn_knn_graph(traj_lle, 5)
ends <- find_extreme_cells(traj_graph, traj_lle)
start <- 1 # ends[1]
cells_ordered <- cell_order(traj_graph, start)
# branches <- assign_branches(traj_graph,start)

calc_vertex_freqs <- function(path_counts, path_len) {
  n <- length(path_len)
  vertex_freqs <- matrix(0, n, n)
  for (i in seq_len(n)) {
    vertex_freqs[path_len[i], i] <- path_counts[i]
  }
  vertex_freqs
}

calc_path_counts <- function(vpath, which_ind) {
  counts <- integer(length(which_ind))
  for (i in seq_along(vpath)) {
    index <- which_ind[vpath[[i]]]
    counts[index] <- counts[index] + 1L
  }
  counts
}

# This updated function is much faster, as should
# fix some bugs associated with SLICER package.
assign_branches <- function (traj_graph, start, min_branch_len = 10, cells = igraph::V(traj_graph))
{
  # browser()
  smart_table = function(x, y, n) {
    tab = rep(0, n)
    tab[y[x]] = 1
    return(tab)
  }
  n = length(cells)
  which_ind <- integer(length(igraph::V(traj_graph)))
  which_ind[cells] <- 1:n
  paths = igraph::shortest_paths(traj_graph, from = start, to = cells)
  # path_table = sapply(paths$vpath, smart_table, which_ind, n)
  # path_counts = rowSums(path_table)
  path_counts <- calc_path_counts(paths$vpath, which_ind)
  path_len <- vapply(paths$vpath, length, FUN.VALUE = integer(1))
  vertex_freqs <- calc_vertex_freqs(path_counts, path_len)
  vertex_probs = vertex_freqs/rowSums(vertex_freqs)
  max_path_length = max(path_len)
  entropies <- vapply(1:max_path_length,
                      function(i, vertex_probs) {
                        probs <- vertex_probs[i,]
                        val <- probs[probs > 0]
                        -sum( val * log2(val))
                      }, FUN.VALUE = numeric(1),
                      vertex_probs = vertex_probs)

  branch_point = which(entropies > 1)[1] - 1
  num_branches = round(2^(entropies[branch_point + 1]))
  if (is.na(num_branches) | is.nan(num_branches)) {
    return(rep(1, n))
  }
  crit_verts = order(vertex_probs[branch_point + 1, ], decreasing = TRUE)[1:num_branches]
  num_branches <- num_branches - sum(vertex_freqs[branch_point + 1, crit_verts] < min_branch_len)
  # for (i in crit_verts) {
  #   if (vertex_freqs[branch_point + 1, i] < min_branch_len) {
  #     num_branches = num_branches - 1
  #   }
  # }t
  if (num_branches < 2) {
    return(rep(1, n))
  }
  assign_cell = function(x, branch_point, crit_verts) {
    if (length(x) <= branch_point) {
      return(1)
    }
    for (i in 1:length(crit_verts)) {
      if (which_ind[x[branch_point + 1]] == crit_verts[i]) {
        return(i + 1)
      }
    }
    return(1)
  }
  path_long <- path_len >= branch_point + 1
  cell_assignments <- rep(1L, n)
  cell_assignments[path_long] <- cell_assignments[path_long] +
    vapply(paths$vpath[path_long], function(path, table, i) {
      match(path[i], table = table, nomatch = 0L)
    }, FUN.VALUE = integer(1), table = crit_verts, i = branch_point + 1)
  branch_assignments = sapply(paths$vpath, assign_cell, branch_point,crit_verts) |>
    as_seq()
  if (max(branch_assignments)-1 > num_branches) {
    # browser()
    # compute once
    cell_dist <- igraph::distances(traj_graph, v = cells, to = cells)

    while ((max(branch_assignments)- 1) > num_branches) {
      # min_dist <- Inf
      # min_i <- 0
      # min_j <- 0
      branch_splits <- split(1:n, branch_assignments)
      combos <- combn(unique(branch_assignments), 2, simplify = F) |>
        Filter(f = function(x) {! 1 %in% x}, x = _)
      mean_dists <- vapply(
        combos,
        FUN = function(id, dist, split_id){
          mean(dist[split_id[[id[1]]],
                    split_id[[id[2]]]],)
          },
        FUN.VALUE = numeric(1),
        dist = cell_dist,
        split_id = branch_splits)
      selected <- combos[[which.min(mean_dists)]]
      branch_assignments[branch_assignments==max(selected)] <- min(selected)
      branch_assignments <- as_seq(branch_assignments)
    }

  }
  while ((max(branch_assignments) - 1) > num_branches) {
    min_dist = Inf
    min_i = 0
    min_j = 0
    browser()
    for (i in 3:max(branch_assignments)) {
      for (j in 2:(i - 1)) {
        mean_dist = mean(igraph::distances(traj_graph,
                                           v = cells[branch_assignments == i],
                                           to = cells[branch_assignments == j]))
        # if (is.nan(mean_dist)) {
        #   browser()
        # }
        cat(mean_dist, " < ", min_dist, "\n")
        if (!is.nan(mean_dist) && mean_dist < min_dist) {
          min_dist = mean_dist
          min_i = i
          min_j = j
        }
      }
    }
    branch_assignments[branch_assignments == min_i] = min_j
  }
  branch_assignments <- as_seq(branch_assignments)
  geodesic_dists = process_distance(traj_graph, start)
  for (i in 1:num_branches) {
    branch_i = which(branch_assignments == (i + 1))
    recurse_br = assign_branches(traj_graph,
                                 start = cells[branch_i[which.min(geodesic_dists[cells[branch_i]])]],
                                 min_branch_len, cells[branch_i])
    rec_num_branches = max(recurse_br)
    if (rec_num_branches > 1) {
      add_this = max(branch_assignments)
      #
      not_one <- recurse_br > 1
      recurse_br[not_one] <- recurse_br[not_one] + add_this
      # for (j in 2:rec_num_branches) {
      #   recurse_br[recurse_br == j] = add_this + j -
      #     1
      # }
      recurse_br[recurse_br == 1] = max(branch_assignments)
      branch_assignments[branch_i] = recurse_br
    }
  }
  as_seq(branch_assignments)
}
branches <- assign_branches(traj_graph,start)

p <- tibble(x = traj_lle[,1],
            y = traj_lle[,2],
            z = traj_lle[,3],
            branch = factor(branches),
            order = cells_ordered,
            cell_type = colData(data)$cell_type) |>
  dplyr::arrange(order) |>
  ggplot(aes(x, y)) +
  geom_point(aes(color = cell_type, group = order), size = 6) +
  geom_point(aes(color = cell_type, group = order), size = 2) +
  labs(x = "Manifold 1", y = "Manifold 2", title = "Cell Ordering") +
  theme_bw()

anim <- p + transition_components(order, enter_length = 50L, exit_length = 50L) +
  shadow_mark(exclude_layer = 1) +
  enter_grow() +
  enter_fade() +
  exit_shrink(size = .2) +
  exit_fade(alpha = .3)

animated <- animate(anim, duration = 20, fps = 30)
