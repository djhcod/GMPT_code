# 加载包 -------------------------------------------------------------------------------

cat("\014")
{
  library(Seurat)
  library(qs)
  library(tidyverse)
  library(beepr)
  library(cowplot)
}


# 数据读取 ------------------------------------------------------------------------------

seurat_filtered <- qread("output/seurat_obj/seurat_combined.qs")
seurat_filtered
head(seurat_filtered)




# 标准化、找高变基因、归一化 -------------------------------------------------------

seurat_phase <- NormalizeData(seurat_filtered)
seurat_phase


# Identify the most variable genes
seurat_phase <- FindVariableFeatures(
  seurat_phase,
  selection.method = "vst",
  nfeatures = 2000
)

# 归一化
seurat_phase <- ScaleData(seurat_phase)
seurat_phase



# PCA
{
  seurat_phase <- RunPCA(seurat_phase)
  beep()
}








# 分析细胞周期和线粒体基因的影响------------------------------------------------------------

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




# 载入颜色集
load("output/color_palette.rdata")


# Plot the PCA colored by cell cycle phase
p1 <- DimPlot(
  seurat_phase,
  reduction = "pca",
  group.by = "Phase",
  cols = mycolors,
  raster = FALSE # 禁止栅格化
) +
  theme_bw()
p2 <- DimPlot(
  seurat_phase,
  reduction = "pca",
  group.by = "Phase",
  split.by = "Phase",
  cols = mycolors,
  raster = FALSE
) +
  theme_bw()
p_cell_circ <- plot_grid(p1, p2, ncol = 2, labels = "AUTO")
p_cell_circ

# Visualize the distribution of cell cycle markers across
Idents(seurat_phase) <- "Phase" # 更改idents为细胞周期便于绘制RidgePlot
levels(seurat_phase)
RidgePlot(
  seurat_phase,
  features = c("TUBB4B", "SMC4", "CKS1B", "HMGB2"),
  cols = mycolors,
  ncol = 2
)
ggsave(width = 9.2, height = 5.8, "figures/seurat_pipeline/cell_cycle_marker_exp.pdf")
Idents(seurat_phase) <- "orig.ident"



# 评估线粒体基因的影响
# 根据各细胞线粒体基因的比例信息绘制PCA
p1 <- DimPlot(
  seurat_phase,
  reduction = "pca",
  group.by = "mitoFr",
  cols = mycolors,
  raster = FALSE
) +
  theme_bw()
p2 <- DimPlot(
  seurat_phase,
  reduction = "pca",
  group.by = "mitoFr",
  split.by = "mitoFr",
  cols = mycolors,
  raster = FALSE
) +
  theme_bw()
p_mito <- plot_grid(p1, p2, ncol = 2, labels = c("C", "D"))
p_mito
plot_grid(p_cell_circ, p_mito, nrow = 2)
ggsave(
  width = 4500,
  height = 3000,
  units = "px",
  filename = "figures/seurat_pipeline/cell_circle_mito_PCA.pdf"
  )





# 分割layer，重新归一化 ------------------------------------------------------------------------------

seurat_split <- seurat_phase
seurat_split[["RNA"]] <- split(
  seurat_split[["RNA"]],
  f = seurat_split$sample_id
)
seurat_split




# 在每个layer中查找高变基因
seurat_split <- FindVariableFeatures(seurat_split, verbose = FALSE)

# 在每个layer中执行归一化
{
  seurat_split <- ScaleData(seurat_split)
  beep()
}
seurat_split

qsave(seurat_split, file = "output/seurat_obj/seurat_split.qs")





# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/2_ScaleData.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
