# 加载包 -------------------------------------------------------------------------------

cat("\014")
rm(list = ls())
{
  library(qs)
  library(clustree)
}
# 加载分群相关函数
source("function/cluster_functions.r")


# 载入数据 ------------------------------------------------------------------------------

seurat_integrated <- qread("output/seurat_obj/seurat_integrated.qs")
seurat_integrated
head(seurat_integrated)

# 载入颜色集
load("output/color_palette.rdata")





# 决定后续分析的主成分 ------------------------------------------------------------------------

# Explore heatmap of PCs
png(
  "figures/seurat_pipeline/pc_heatmap.png",
  width = 10,
  height = 10,
  units = "in",
  res = 300
  )
DimHeatmap(
  seurat_integrated,
  dims = 1:15,
  cells = 2000
)
dev.off()


# Plot the elbow plot
plot_elbow <- ElbowPlot(seurat_integrated, ndims = 50) +
  theme_bw()
plot_elbow

# 绘制标注了理想主成分数量的自定义肘图
plot_cum_elbow <- sc.ElbowPlot(seurat_integrated)

# 合并两幅肘图
plot_grid(plot_elbow, plot_cum_elbow, nrow = 1, labels = "AUTO")
ggsave("figures/seurat_pipeline/elbow_plot.pdf", width = 11, height = 5)




# 根据确定的主成分数量重新非线性降维、聚类
seurat_integrated <- FindNeighbors(
  seurat_integrated,
  dims = 1:pc_aim,
  reduction = "integrated.reduction"
) %>%
  RunUMAP(
    dims = 1:pc_aim,
    reduction = "integrated.reduction",
    reduction.name = "umap.integrated",
    verbose = FALSE
  )




# 确定分辨率 -----------------------------------------------------------------------------
resolutions <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.8)
seurat_integrated <- FindClusters(
  seurat_integrated,
  resolution = resolutions,
  verbose = FALSE
)
head(seurat_integrated)
select(
  seurat_integrated@meta.data,
  starts_with("RNA_snn_res.")
) %>%
  lapply(levels) %>% print()


# 批量绘制不同分辨率下的UMAP图
plot_umap_res <- sc.DimPlot.res(seurat_integrated)


# 绘制聚类树展示不同分辨率下的细胞分群情况及相互关系
plot_tree <- clustree(
  seurat_integrated,
  prefix = "RNA_snn_res."
) +
  scale_color_manual(values = mycolors) +
  theme(
    legend.box.background = element_rect(color = "black", fill = NA),
    plot.margin = margin(0, 20, 0, 0) # 增加右侧边距
  )
plot_tree


plot_grid(plot_tree, plot_umap_res, nrow = 1, labels = "AUTO")
ggsave("figures/seurat_pipeline/resolution_exploration.pdf", width = 15, height = 10)




# 设定分辨率、绘制最终的UMAP图 ---------------------------------------------------
Idents(seurat_integrated) <- "RNA_snn_res.0.3"

# 绘制最终的UMAP图
sc.DimPlot(seurat_integrated, reduction = "umap.integrated", colors = mycolors_many)
ggsave("figures/seurat_pipeline/umap_plot.pdf", width = 6.5, height = 5.4)


# 在整合后的UMAP图上映射样本来源信息
p1 <- DimPlot(
  seurat_integrated,
  reduction = "umap.integrated",
  group.by = "sample_id",
  cols = mycolors_many,
  raster = FALSE
) +
  labs(color = "Sample ID") +
  theme_bw()
p2 <- DimPlot(
  seurat_integrated,
  reduction = "umap.integrated",
  group.by = "study_id",
  cols = mycolors,
  raster = FALSE
) +
  labs(color = "Study ID") +
  theme_bw()
p3 <- DimPlot(
  seurat_integrated,
  reduction = "umap.integrated",
  group.by = "sample_type",
  cols = mycolors,
  raster = FALSE
) +
  labs(color = "Sample type") +
  theme_bw()

plot_grid(p1, p2, p3, nrow = 1, labels = "AUTO", rel_widths = c(90, 80, 75))
ggsave("figures/seurat_pipeline/UMAP_after_integration.pdf", height = 5, width = 19)




qsave(seurat_integrated, "output/seurat_obj/seurat_clustered.qs")






# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/4_Clustering.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
