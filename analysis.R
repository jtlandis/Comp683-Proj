

box::use(mods/analysis_tools[...])

data <- read_h5ad(box::file("data/olsson/adata_olsson.h5ad"))

traj_df <- do_SLICER(
  data = data,
  gene_subset = readRDS(box::file("data/olsson/olsson_SLICER_select_genes.rds"))
)

plot_SLICER(traj_df, Dim1, Dim2)

aricode::NMI(traj_df$branch, traj_df$cell_type)
table(traj_df$branch, traj_df$cell_type)

p <-  plot_SLICER(traj_df, Dim1, Dim2)

## Animation ----

# box::use(gganimate[...])
# anim <- p + transition_components(order, enter_length = 50L, exit_length = 50L) +
#   shadow_mark(exclude_layer = 1) +
#   enter_grow() +
#   enter_fade() +
#   exit_shrink(size = .2) +
#   exit_fade(alpha = .3)
#
# anim_save(box::file("prelim_data.gif"),
#           animation = anim, duration = 20, fps = 30)

## DESeq2 ----
##

data$cell_type <- relevel(factor(data$cell_type), "LSK")
dds_cell_type <- do_deseq(data, ~ cell_type)
p1 <- plot_deseq(dds_cell_type, padj = 1e-10)
data$branches <- relevel(factor(branches), "12")
dds_branch <- do_deseq(data, ~ branches)
p2 <- plot_deseq(dds_branch, padj = 1e-10)


(p1 + theme(legend.position = "none", plot.caption = element_blank(), axis.title.x = element_blank())) |>
  cowplot::plot_grid(p2, align = "v", rel_heights = c(1,2.3), ncol = 1) -> cowplot1

## Save Sig-genes ----

readr::write_csv(subset(p1$data, type!=""), file = box::file("data/output/cell_type_DE_data.csv"))
readr::write_csv(subset(p2$data, type!=""), file = box::file("data/output/SLICER_branch_DE_data.csv"))
