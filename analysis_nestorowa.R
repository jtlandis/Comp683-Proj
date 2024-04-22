box::use(mods/analysis_tools[...],
         ggplot2[...])
box::use(aricode[NMI])

data2 <- read_h5ad(box::file("data/nestorowa/adata_nestorowa.h5ad"))
.gsub <- readRDS(box::file("data/nestorowa/nestorowa_SLICER_select_genes.rds"))
traj_df1 <- do_SLICER(
  data = data2,
  gene_subset = .gsub,
  select_k = 5,
  start = 812
)

find_extreme_cells(
  attr(traj_df1, "traj_graph"),
  attr(traj_df1, "lle")
)

traj_df2 <- do_SLICER(
  data = data2,
  gene_subset = .gsub,
  select_k = 5,
  start = 1695
)

NMI(traj_df1$branch, traj_df2$branch)
NMI(traj_df1$branch, traj_df1$cell_types_IDs)
NMI(traj_df2$branch, data2$cell_types_IDs)
NMI(traj_df2$branch, data2$cell_types_broad)
NMI(traj_df2$branch, data2$cell_types_fine)
NMI(traj_df2$branch, data2$cell_types_broad_cleaned)

paths <- get_paths(traj_df1)
plot(attr(paths, "branch_graph"))
plot_SLICER(traj_df2, Dim1, Dim2)
plot_SLICER(traj_df2, Dim1, Dim2, color = cell_types_IDs)
table(traj_df$branch, traj_df$cell_types_IDs)

# ## Animation ----
#
animate_SLICER(traj_df1, Dim1, Dim2, branch = branch,
               path = box::file("figs/nestorowa-traj-rev.gif"))

animate_SLICER(traj_df2, Dim1, Dim2, branch = branch,
               path = box::file("figs/nestorowa-traj.gif"))
animate_SLICER(traj_df2, Dim1, Dim2, branch = cell_types_broad,
               path = box::file("figs/nestorowa-traj-celltype.gif"), exclude_legend = FALSE)

## DESeq2 ----
## choosing to not do DE against one of the groups like we did in
## olsson data set. I think doing DE along the trajectory makes more
## sense.

dds_along <- do_deseq_along(data, traj_df2)
saveRDS(dds_along, box::file("data/nestorowa/nestorowa_dds_along.rds"))
plot_along <- plot_deseq_along(dds_along, padj = 1e-10)
olsson <- readRDS(box::file("data/olsson/olsson_dds_along.rds"))
olsson_genes <- plot_deseq_along(olsson, padj = 1e-10) |>
  _$data |>
  subset(type!="") |>
  _$genes |> unique()
unique_genes <- plot_along$data |>
  subset(type!="") |>
  dplyr::pull(genes) |>
  unique()

ggsave(filename = box::file("figs/nestorowa-DE-along.png"),
       plot_along, dpi = 500)


## Save Sig-genes ----

# readr::write_csv(subset(p1$data, type!=""), file = box::file("data/output/nestorowa_cell_type_DE_data.csv"))
# readr::write_csv(subset(p2$data, type!=""), file = box::file("data/output/nestorowa_SLICER_branch_DE_data.csv"))
