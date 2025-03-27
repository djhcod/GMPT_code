# 加载包 -------------------------------------------------------------------------------

rm(list = ls())
cat("\014")
{
  library(Seurat)
  library(tidyverse)
  library(qs)
  library(beepr)
}

# 加载数据 ------------------------------------------------------------------------------

seurat_split <- qread("output/seurat_obj/seurat_split.qs")
seurat_split
seurat_split@assays
head(seurat_split)

# 不进行整合时检验细胞分群情况 -------------------------------------------------------------

# PCA
{
  seurat_split <- RunPCA(seurat_split)
  beep()
}

# Run UMAP
{
  seurat_split <- seurat_split %>%
    RunUMAP(
      dims = 1:10,
      reduction = "pca",
      reduction.name = "umap.unintegrated",
      verbose = FALSE
    ) %>%
    FindNeighbors(dims = 1:10, reduction = "pca") %>%
    FindClusters(cluster.name = "unintegrated_clusters")
  beep()
}
head(seurat_split)
table(seurat_split$unintegrated_clusters)



# 载入颜色集
load("output/color_palette.rdata")

# 载入自定义UMAP图绘制函数
source("function/cluster_functions.r")

# 检查整合前的分群情况
p1 <- DimPlot(
  seurat_split,
  reduction = "umap.unintegrated",
  group.by = "sample_id",
  cols = mycolors_many,
  raster = FALSE
) +
  labs(color = "Sample ID") +
  theme_bw()
p2 <- DimPlot(
  seurat_split,
  reduction = "umap.unintegrated",
  group.by = "study_id",
  cols = mycolors,
  raster = FALSE
) +
  labs(color = "Study ID") +
  theme_bw()
p3 <- DimPlot(
  seurat_split,
  reduction = "umap.unintegrated",
  group.by = "sample_type",
  cols = mycolors,
  raster = FALSE
) +
  labs(color = "Sample type") +
  theme_bw()

plot_grid(p1, p2, p3, nrow = 1, labels = "AUTO", rel_widths = c(90, 80, 75))
ggsave("figures/seurat_pipeline/UMAP_before_integration.pdf", height = 5, width = 19)









# 整合 --------------------------------------------------------------------------------

seurat_split@reductions

if (TRUE) {
  library(SeuratWrappers)
  seurat_integrated <- IntegrateLayers(
    seurat_split,
    method = FastMNNIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.reduction"
  )
  # 整合后重新合并layers
  seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])
  beep()
}


if (FALSE) {
  seurat_integrated <- IntegrateLayers(
    seurat_split,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.reduction",
    verbose = FALSE
  )
  # 整合后重新合并layers
  seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])
  beep()
}
seurat_integrated


qsave(seurat_integrated, "output/seurat_obj/seurat_integrated.qs")









# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/3_Integration.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
