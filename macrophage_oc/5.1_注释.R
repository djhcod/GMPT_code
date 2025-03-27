# 加载包 -------------------------------------------------------------------------------

cat("\014")
rm(list = ls())
{
  library(qs)
  library(Seurat)
  library(tidyverse)
  library(GPTCelltype)
  library(patchwork)
}

# 载入数据 ------------------------------------------------------------------------------

seurat_clustered <- qread("output/seurat_obj/seurat_clustered.qs")
seurat_clustered
head(seurat_clustered)
levels(seurat_clustered)

# 载入颜色集
load("output/color_palette.rdata")

# 加载分群函数
source("function/cluster_functions.r")


# 自动注释 ------------------------------------------------------------------

# 寻找各cluster的差异基因
cluster_degs <- FindAllMarkers(
  seurat_clustered,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  only.pos = TRUE
)
head(cluster_degs)
unique(cluster_degs$gene) %>% length()
write_csv(cluster_degs, "output/cluster_degs.csv")

# 生成用于命令ChatGPT进行细胞类型注释的提示词
gpt_prompt <- gptcelltype(
  cluster_degs,
  tissuename = "human ovary"
)
gpt_prompt
qsave(gpt_prompt, "output/gpt_prompt.qs")

# gpt-3.5-turbo-0125模型的注释结果
cell_types_gpt3.5_turbo <- tribble(
  ~cluster_id, ~cell_type,
  "0", "T cells",
  "1", "Macrophages",
  "2", "Ovarian surface epithelial cells",
  "3", "Ovarian epithelial cells",
  "4", "T cells",
  "5", "Fibroblasts",
  "6", "Fibroblasts",
  "7", "Fibroblasts",
  "8", "Endothelial cells",
  "9", "Smooth muscle cells",
  "10", "B cells",
  "11", "Fibroblasts",
  "12", "Macrophage",
  "13", "Fibroblasts",
  "14", "B cells",
  "15", "Epithelial cells",
  "16", "Proliferating cells",
  "17", "Natural killer cells"
)

# gpt-4-turbo-2024-04-09模型的注释结果
cell_types_gpt4_turbo <- tribble(
  ~cluster_id, ~cell_type,
  "0", "T cell",
  "1", "Macrophages",
  "2", "Epithelial cells",
  "3", "Epithelial cells",
  "4", "T cells",
  "5", "Fibroblasts",
  "6", "Fibroblasts",
  "7", "Fibroblasts",
  "8", "Endothelial cells",
  "9", "Smooth muscle cells",
  "10", "B cell",
  "11", "Fibroblasts",
  "12", "Macrophages",
  "13", "Fibroblasts",
  "14", "B cells",
  "15", "Epithelial cells",
  "16", "Proliferating cells",
  "17", "Natural killer cells"
)



# gpt-4o-2024-05-13模型的注释结果
cell_types_gpt4o <- tribble(
  ~cluster_id, ~cell_type,
  "0", "T cells (cytotoxic subset)",
  "1", "Macrophages",
  "2", "Epithelial cells",
  "3", "Epithelial cells (ovarian)",
  "4", "T cells",
  "5", "Fibroblasts",
  "6", "Fibroblasts",
  "7", "Fibroblasts",
  "8", "Endothelial cells",
  "9", "Smooth muscle cells",
  "10", "Plasma cells",
  "11", "Fibroblasts",
  "12", "Monocytes",
  "13", "Fibroblasts",
  "14", "B cells",
  "15", "Epithelial cells (stromal)",
  "16", "Proliferating cells",
  "17", "Natural killer (NK) cells"
)
cell_types_gpt4o













# 根据传统marker基因进行手动注释 --------------------------------------------------------

# 载入包含了细胞注释marker gene的表格
marker_genes <- read_csv("data/marker_genes.csv")
marker_genes %>% print(n = Inf)
# 定义需要注释的细胞类型
cell_type_to_check <- unique(na.omit(marker_genes$cell_type))
cell_type_to_check

# 构建marker gene列表
genes_to_check <- map(cell_type_to_check, function(x) {
  marker_genes %>%
    filter(cell_type == x) %>%
    pull(marker_gene)
})
genes_to_check <- set_names(genes_to_check, cell_type_to_check)
genes_to_check


# 绘制气泡图
p1 <- DotPlot(
  seurat_clustered,
  features = genes_to_check,
  cols = mycolors
) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.title = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p1
ggsave(width = 18.5, height = 6.2, "figures/seurat_pipeline/dot_plot.pdf")




# 确定最终的注释

new.cluster.ids <- c(
  "0" = "T/NK cells",
  "1" = "Macrophages",
  "2" = "Cancer/Epithelial cells",
  "3" = "Cancer/Epithelial cells",
  "4" = "T/NK cells",
  "5" = "Fibroblasts",
  "6" = "Fibroblasts",
  "7" = "Fibroblasts",
  "8" = "Endothelial cells",
  "9" = "Smooth muscle cells",
  "10" = "B cells/plasma cells",
  "11" = "Fibroblasts",
  "12" = "Monocytes",
  "13" = "Fibroblasts",
  "14" = "B cells/plasma cells",
  "15" = "Stromal cells (undefined)",
  "16" = "Proliferating cells",
  "17" = "T/NK cells"
)

# 合并自动-手动注释结果
anno_compare <- tibble(
  cluster_id = cell_types_gpt3.5_turbo$cluster_id,
  cell_types_gpt3.5_turbo = cell_types_gpt3.5_turbo$cell_type,
  cell_types_gpt4_turbo = cell_types_gpt4_turbo$cell_type,
  cell_types_gpt4o = cell_types_gpt4o$cell_type,
  cell_types_manual = new.cluster.ids
)
anno_compare
write_csv(anno_compare, "output/annotation_compare.csv")




seurat_annotated <- RenameIdents(
  seurat_clustered,
  new.cluster.ids
)
table(Idents(seurat_annotated))

# 把细胞注释信息保存到meta.data中
seurat_annotated$cell_type <- Idents(seurat_annotated)
head(seurat_annotated, 3)
table(seurat_annotated$RNA_snn_res.0.3, seurat_annotated$cell_type)


# 绘制注释后的UMAP图
plot_grid(
  sc.DimPlot(
    seurat_clustered,
    reduction = "umap.integrated",
    colors = mycolors_many,
    pt.size = 0.00001,
    legend.title = "Cell cluster"
  ),
  ggDimPlot(
    seurat_annotated,
    reduction = "umap.integrated",
    colors = mycolors,
    pt.size = 0.00001,
    legend.title = "Cell type"
  ),
  nrow = 1,
  rel_widths = c(85, 100),
  labels = "AUTO"
)
ggsave(
  "figures/seurat_pipeline/umap_plot_annotated.pdf",
  width = 13.7,
  height = 5.4
  )






# marker基因及cluster差异基因可视化验证 --------------------------------------------------------

# 绘制marker gene表达情况的feature plot
marker.FeaturePlot <- function(
    seurat_obj,
    reduction = "umap.integrated",
    features,
    ncol = 5
    ) {
  map(features, ~ FeaturePlot(
    seurat_obj,
    features = .x,
    reduction = reduction,
    cols = c("gray", mycolors[2]),
    pt.size = 0.01,
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) +
    theme_bw()) %>%
    wrap_plots(ncol = ncol)
}

{
  features <- marker_genes %>%
    filter(cell_type == "Cancer/ Epithelial cells" | cell_type == "Immune cells") %>%
    pull(marker_gene)
  marker.FeaturePlot(seurat_annotated, features = features)
  ggsave(
    width = 21.98,
    height = 3.8 * ceiling(length(features) / 5),
    "figures/seurat_pipeline/marker_gene_featureplot/epi_immu.png"
  )


  features <- marker_genes %>%
    filter(cell_type == "T/NK cells") %>%
    pull(marker_gene)
  marker.FeaturePlot(seurat_annotated, features = features)
  ggsave(
    width = 21.98,
    height = 3.8 * ceiling(length(features) / 5),
    "figures/seurat_pipeline/marker_gene_featureplot/t_nk.png"
  )


  features <- marker_genes %>%
    filter(cell_type == "B cells/plasma cells") %>%
    pull(marker_gene)
  marker.FeaturePlot(seurat_annotated, features = features)
  ggsave(
    width = 21.98,
    height = 3.8 * ceiling(length(features) / 5),
    "figures/seurat_pipeline/marker_gene_featureplot/b_plasma.png"
  )


  features <- marker_genes %>%
    filter(cell_type == "Monocytes/Macrophages") %>%
    pull(marker_gene)
  marker.FeaturePlot(seurat_annotated, features = features)
  ggsave(
    width = 21.98,
    height = 3.8 * ceiling(length(features) / 5),
    "figures/seurat_pipeline/marker_gene_featureplot/mono_macro.png"
  )

  features <- marker_genes %>%
    filter(cell_type == "Fibroblasts") %>%
    pull(marker_gene)
  marker.FeaturePlot(seurat_annotated, features = features)
  ggsave(
    width = 21.98,
    height = 3.8 * ceiling(length(features) / 5),
    "figures/seurat_pipeline/marker_gene_featureplot/fib.png"
  )


  features <- marker_genes %>%
    filter(cell_type == "Smooth muscle cells") %>%
    pull(marker_gene)
  marker.FeaturePlot(seurat_annotated, features = features)
  ggsave(
    width = 21.98,
    height = 3.8 * ceiling(length(features) / 5),
    "figures/seurat_pipeline/marker_gene_featureplot/smc.png"
  )


  features <- marker_genes %>%
    filter(cell_type == "Endothelial cells") %>%
    pull(marker_gene)
  marker.FeaturePlot(seurat_annotated, features = features)
  ggsave(
    width = 21.98,
    height = 3.8 * ceiling(length(features) / 5),
    "figures/seurat_pipeline/marker_gene_featureplot/endo.png"
  )


  features <- marker_genes %>%
    filter(cell_type == "Proliferating cells") %>%
    pull(marker_gene)
  marker.FeaturePlot(seurat_annotated, features = features)
  ggsave(
    width = 21.98,
    height = 3.8 * ceiling(length(features) / 5),
    "figures/seurat_pipeline/marker_gene_featureplot/proliferating.png"
  )
}








# 寻找各亚群的TOP差异基因
all_degs <- FindAllMarkers(
  seurat_annotated,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1
)
head(all_degs)
write.csv(all_degs, "output/annotated_cluster_degs.csv")


# 获取每个cluster的前10个差异基因
top_degs <- all_degs %>%
  group_by(cluster) %>%
  slice_head(n = 10)
top_degs
write.csv(top_degs, "output/top_degs.csv", row.names = FALSE)


# 绘制TOP差异基因的热图
ScaleData(seurat_annotated, features = top_degs$gene) %>%
  subset(downsample = 2000) %>%
  DoHeatmap(
    features = top_degs$gene,
    group.colors = mycolors,
    size = 2,
    raster = FALSE
  ) +
  scale_fill_gradientn(colors = mycolors[1:3]) +
  theme(axis.text = element_text(size = 6))
ggsave(
  width = 8.62,
  height = 8.16,
  "figures/seurat_pipeline/top_degs_heatmap.png"
  )



qsave(seurat_annotated, "output/seurat_obj/seurat_annotated.qs")














# 保存sessionInfo -----------------------------------------------------------

sink(
  "sessionInfo/5.1_Cell_annotation.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
