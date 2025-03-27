# 加载包 ---------------------------------------------------------------------

{
  rm(list = ls())
  cat("\014")

  library(qs)
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(clustree)
  library(scatterplot3d) # 绘制3D UMAP图
  library(gridGraphics) # 将base-R plots转换为ggplot2绘图框架

  # 载入颜色集
  load("output/color_palette.rdata")

  # 加载分群相关函数
  source("function/cluster_functions.r")
}



# 载入数据 ------------------------------------------------------------------------------

seurat_macro <- qread("output/seurat_obj/seurat_macro.qs")
seurat_macro
head(seurat_macro)
levels(seurat_macro)
table(seurat_macro$sample_id)

# 去掉细胞数过少(<100)的样本
samples <- setdiff(
  unique(seurat_macro$sample_id),
  c("BT1306", "GSM5599220", "GSM5599223", "scrSOL007")
)

seurat_macro <- subset(
  seurat_macro,
  sample_id %in% samples
)
table(seurat_macro$sample_id)




# 归一化 ---------------------------------------------------------------------

seurat_phase <- SCTransform(seurat_macro, verbose = FALSE)
seurat_phase



# PCA
seurat_phase <- RunPCA(seurat_phase)



## 分析细胞周期和线粒体基因的影响---------------------------------------------

# 加载细胞周期标记基因
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# 为每个细胞进行细胞周期打分
seurat_phase <- CellCycleScoring(
  seurat_phase,
  g2m.features = g2m.genes,
  s.features = s.genes
)
head(seurat_phase)
table(seurat_phase$Phase)





# 查看线粒体基因的表达情况
quart_mito <- summary(seurat_phase$mitoRatio)
# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase$mitoFr <- cut(
  seurat_phase$mitoRatio,
  breaks = c(-Inf, quart_mito[2], quart_mito[3], quart_mito[5], Inf),
  labels = c("Low", "Medium", "Medium high", "High")
)
table(seurat_phase$mitoFr)


# Plot the PCA colored by cell cycle phase
p1 <- DimPlot(
  seurat_phase,
  reduction = "pca",
  group.by = "Phase",
  cols = mycolors,
  raster = FALSE # 禁止栅格化
) +
  labs(color = "Cell cycle phase") +
  theme_light() +
  theme(plot.title = element_blank())
p2 <- DimPlot(
  seurat_phase,
  reduction = "pca",
  group.by = "Phase",
  split.by = "Phase",
  cols = mycolors,
  raster = FALSE
) +
  theme_bw() +
  theme(plot.title = element_blank())
p_cell_circ <- plot_grid(
  p1 + NoLegend(),
  p2 + NoLegend(),
  get_legend(p1),
  nrow = 1,
  rel_widths = c(3, 7, 1)
)
p_cell_circ




# 评估线粒体基因的影响
# 根据各细胞线粒体基因的比例信息绘制PCA
p1 <- DimPlot(
  seurat_phase,
  reduction = "pca",
  group.by = "mitoFr",
  cols = mycolors,
  raster = FALSE
) +
  labs(color = "Mitochondrial gene ratio") +
  theme_bw() +
  theme(plot.title = element_blank())
p2 <- DimPlot(
  seurat_phase,
  reduction = "pca",
  group.by = "mitoFr",
  split.by = "mitoFr",
  cols = mycolors,
  raster = FALSE
) +
  theme_bw() +
  theme(plot.title = element_blank())
p_mito <- plot_grid(
  p1 + NoLegend(),
  p2 + NoLegend(),
  get_legend(p1),
  nrow = 1,
  rel_widths = c(3, 7, 1)
)
p_mito
plot_grid(p_cell_circ, p_mito, nrow = 2)
ggsave(
  width = 16.5, height = 8,
  filename = "figures/macro/cell_circle_mito_PCA.pdf"
)





## 分割layer，重新归一化-------------------------------------------------------

seurat_split <- seurat_phase
seurat_split[["RNA"]] <- split(
  seurat_split[["RNA"]],
  f = seurat_split$sample_id
)
seurat_split



# 在每个layer中执行归一化
seurat_split <- SCTransform(seurat_split, verbose = FALSE)
seurat_split

qsave(seurat_split, file = "output/macro/seurat_split.qs")





# 整合 ----------------------------------------------------------------------

seurat_split <- qread("output/macro/seurat_split.qs")

seurat_split <- RunPCA(seurat_split)

# Run UMAP
seurat_split <- seurat_split %>%
  RunUMAP(
    dims = 1:30,
    reduction = "pca",
    reduction.name = "umap.unintegrated",
    verbose = FALSE
  ) %>%
  FindNeighbors(dims = 1:30, reduction = "pca") %>%
  FindClusters(cluster.name = "unintegrated_clusters")

head(seurat_split)
table(seurat_split$unintegrated_clusters)


# 检查整合前的分群情况

p1 <- sc.DimPlot(
  seurat_split,
  reduction = "umap.unintegrated",
  group.by = "sample_id",
  colors = mycolors_many,
  legend.title = "Sample ID"
)
p1

p2 <- sc.DimPlot(
  seurat_split,
  reduction = "umap.unintegrated",
  group.by = "study_id",
  colors = mycolors,
  legend.title = "Study ID"
)
p2

p3 <- sc.DimPlot(
  seurat_split,
  reduction = "umap.unintegrated",
  group.by = "sample_type",
  colors = mycolors,
  legend.title = "Sample type"
)
p3

plot_grid(p1, p2, p3, nrow = 1, labels = "AUTO", rel_widths = c(100, 80, 75))
ggsave("figures/macro/UMAP_before_integration.pdf", height = 5, width = 19)





# 整合

seurat_split@reductions


seurat_integrated <- IntegrateLayers(
  seurat_split,
  method = CCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "integrated.reduction",
  verbose = FALSE
)
# 整合后重新合并layers
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])
seurat_integrated


qsave(seurat_integrated, "output/macro/seurat_integrated.qs")




# 分群 ----------------------------------------------------------------------

seurat_integrated <- qread("output/macro/seurat_integrated.qs")

# 重新非线性降维、聚类
seurat_integrated <- FindNeighbors(
  seurat_integrated,
  dims = 1:30,
  reduction = "integrated.reduction"
) %>%
  RunUMAP(
    reduction = "integrated.reduction",
    reduction.name = "umap.integrated",
    dims = 1:30,
    min.dist = 0.05,
    verbose = FALSE
  )
DimPlot(seurat_integrated, reduction = "umap.integrated")



# 确定分辨率
resolutions <- c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.8, 1)
seurat_integrated <- FindClusters(
  seurat_integrated,
  resolution = resolutions,
  verbose = FALSE
)
head(seurat_integrated)
select(
  seurat_integrated@meta.data,
  starts_with("SCT_snn_res.")
) %>%
  lapply(levels) %>%
  print()


# 批量绘制不同分辨率下的UMAP图
plot_umap_res <- sc.DimPlot.res(seurat_integrated, prefix = "SCT_snn_res.")


# 绘制聚类树展示不同分辨率下的细胞分群情况及相互关系
plot_tree <- clustree(
  seurat_integrated,
  prefix = "SCT_snn_res."
) +
  scale_color_manual(values = mycolors) +
  theme(
    legend.box.background = element_rect(color = "black", fill = NA),
    plot.margin = margin(0, 20, 0, 0) # 增加右侧边距
  )
plot_tree


plot_grid(plot_tree, plot_umap_res, nrow = 1, labels = "AUTO")
ggsave("figures/macro/resolution_exploration.pdf", width = 15, height = 10)

qsave(seurat_integrated, "output/macro/seurat_res.qs")





# 设定分辨率、绘制最终的UMAP图
seurat_integrated <- qread("output/macro/seurat_res.qs")

Idents(seurat_integrated) <- "SCT_snn_res.0.2"

# 绘制最终的UMAP图
sc.DimPlot(
  seurat_integrated,
  reduction = "umap.integrated",
  colors = mycolors,
  label = TRUE,
  legend.title = "Cluster ID"
)
ggsave("figures/macro/umap_plot.pdf", width = 8, height = 5)



# 在整合后的UMAP图上映射样本来源信息
p1 <- sc.DimPlot(
  seurat_integrated,
  reduction = "umap.integrated",
  group.by = "sample_id",
  colors = mycolors_many,
  legend.title = "Sample ID"
)
p1

p2 <- sc.DimPlot(
  seurat_integrated,
  reduction = "umap.integrated",
  group.by = "study_id",
  colors = mycolors_many,
  legend.title = "Study ID"
)
p2

p3 <- sc.DimPlot(
  seurat_integrated,
  reduction = "umap.integrated",
  group.by = "sample_type",
  colors = mycolors_many,
  legend.title = "Sample type"
)
p3

plot_grid(p1, p2, p3, nrow = 1, labels = "AUTO", rel_widths = c(100, 80, 75))
ggsave("figures/macro/UMAP_after_integration.pdf", height = 5, width = 19)


qsave(seurat_integrated, "output/macro/seurat_clustered.qs")






# 亚群命名 --------------------------------------------------------------------

seurat_clustered <- qread("output/macro/seurat_clustered.qs")
seurat_clustered
head(seurat_clustered)
levels(seurat_clustered)


# 寻找各cluster的差异基因
seurat_clustered <- PrepSCTFindMarkers(seurat_clustered)
cluster_degs <- FindAllMarkers(
  seurat_clustered,
  min.pct = 0.2,
  logfc.threshold = 0.25,
  only.pos = TRUE
)
head(cluster_degs)
unique(cluster_degs$gene) %>% length()
write_csv(cluster_degs, "output/macro/cluster_degs.csv")


# 获取每个cluster的top差异基因
top_degs <- cluster_degs %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
top_degs
write.csv(top_degs, "output/macro/top_degs.csv", row.names = FALSE)

# 绘制TOP差异基因的热图
ScaleData(seurat_clustered, features = top_degs$gene) %>%
  subset(downsample = 1000) %>%
  DoHeatmap(
    features = top_degs$gene,
    group.colors = mycolors,
    size = 2,
    raster = FALSE,
    disp.min = 0
  ) +
  scale_fill_gradientn(colors = mycolors[1:3]) +
  theme(axis.text = element_text(size = 6))
ggsave(width = 8.62, height = 8.16, "figures/macro/top_degs_heatmap.png")




# 重命名亚群
new.cluster.ids <- str_c("Macro ", levels(seurat_clustered))
names(new.cluster.ids) <- levels(seurat_clustered)
new.cluster.ids

seurat_annotated <- RenameIdents(
  seurat_clustered,
  new.cluster.ids
)
table(Idents(seurat_annotated))

# 把细胞注释信息保存到meta.data中
seurat_annotated$cell_type <- Idents(seurat_annotated)
head(seurat_annotated, 3)
table(seurat_annotated$SCT_snn_res.0.2, seurat_annotated$cell_type)




## 绘制最终的自定义UMAP图-------------------------------------------------------
umap_plot <- ggDimPlot(
  seurat_annotated,
  reduction = "umap.integrated",
  colors = mycolors,
  legend.title = "Cluster ID"
)
umap_plot

# 绘制3D UMAP图
seurat_annotated <- RunUMAP(
  seurat_annotated,
  reduction = "integrated.reduction",
  reduction.name = "umap.3d",
  dims = 1:30,
  min.dist = 0.05,
  n.components = 3,
  verbose = FALSE
)

umap_3d_tab <- seurat_annotated@reductions$umap.3d@cell.embeddings %>%
  as.data.frame() %>%
  mutate(
    cell_type = seurat_annotated$cell_type,
    type_color = mycolors[as.numeric(as.factor(seurat_annotated$cell_type))]
  )
head(umap_3d_tab)


scatterplot3d(
  x = umap_3d_tab[, 1],
  y = umap_3d_tab[, 2],
  z = umap_3d_tab[, 3],
  color = umap_3d_tab$type_color,
  pch = ".",
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  zlab = "UMAP_3"
)
umap_plot_3d <- recordPlot()

plot_grid(umap_plot_3d, umap_plot, rel_widths = c(100, 95), labels = "AUTO")
ggsave("figures/macro/umap_plot_annotated.pdf", width = 11, height = 5)


qsave(seurat_annotated, "output/macro/seurat_annotated.qs")
















# 保存sessionInfo -----------------------------------------------------------

sink(
  "sessionInfo/6.1_Macrophages_clustering.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
