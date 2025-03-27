library(Seurat)
library(tidyverse)
library(ggrepel)
library(beepr)
library(patchwork)
library(cowplot)
library(ggunchull) # 绘制UMAP图的边界
library(ggh4x) # 更改UMAP图的坐标轴
library(RColorBrewer)




# 自定义肘图绘制函数 ---------------------------------------------------------------------------


sc.ElbowPlot <- function(seurat_obj) {
  # Determine percent of variation associated with each PC
  pc_sd <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100

  # Calculate cumulative percents for each PC
  cumu_sd <- cumsum(pc_sd)

  # Determine which PC exhibits cumulative percent greater than 90%
  # and % variation associated with the PC as less than 5
  co1 <- which(cumu_sd > 90 & pc_sd < 5)[1]

  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(
    which((pc_sd[1:length(pc_sd) - 1] - pc_sd[2:length(pc_sd)]) > 0.1),
    decreasing = T
  )[1] + 1

  # Minimum of the two calculation
  pc_aim <<- min(co1, co2)
  str_glue(
    "---------------------------------------------------------------------------------\n",
    "理想的主成分数量【pc_aim】为【{pc_aim}】\n",
    "---------------------------------------------------------------------------------"
  ) %>% print()


  # Create a table with values
  plot_df <- tibble(
    pc_sd = pc_sd,
    cumu_sd = cumu_sd,
    n_pc = 1:length(pc_sd)
  )
  plot_df

  # Elbow plot to visualize
  plot_cum_elbow <- ggplot(
    plot_df,
    aes(
      x = cumu_sd,
      y = pc_sd,
      label = n_pc,
      color = if_else(n_pc > pc_aim, "PCs below ideal PC", "PCs above ideal PC")
    )
  ) +
    geom_point() +
    scale_color_manual(values = mycolors) +
    geom_label_repel(
      data = subset(plot_df, n_pc == pc_aim),
      aes(label = n_pc),
      label.size = 0.5, # 标注框的粗细
      box.padding = 2, # 标注框与标注对象的间距
      point.padding = 0.3, # 指示线和目标点之间的间距（默认是0）
      nudge_x = 4,
      nudge_y = 0.05,
      arrow = arrow(
        angle = 10,
        length = unit(0.15, "inches"),
        ends = "last",
        type = "closed"
      ),
      show.legend = FALSE
    ) +
    geom_vline(xintercept = 90, color = "gray50", lty = 2) +
    geom_hline(yintercept = min(pc_sd[pc_sd > 5]), color = "gray50", lty = 2) +
    xlab("Cumulative contribution of each PC to the standard deviation (%)") +
    ylab("Contribution of each PC to the standard deviation (%)") +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      legend.direction = "horizontal",
      legend.position = "inside",
      legend.position.inside = c(0.5, 0.93),
      legend.box.background = element_rect(color = "black", fill = NA)
    )
  print(plot_cum_elbow)
  return(plot_cum_elbow)
}










# UMAP/PCA图绘制函数 --------------------------------------------------------------

sc.DimPlot <- function(
    seurat_obj,
    reduction,
    group.by = "ident",
    split.by = NULL,
    colors = brewer.pal(name = "Paired", n = 12),
    pt.size = 0.1,
    label = FALSE,
    label.size = 6,
    legend.title = "Cell type") {
  # 定义截短坐标轴，便于绘制短坐标轴UMAP图
  axis <- guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(3, "cm")
  )

  # 修改坐标轴标题
  if (str_detect(reduction, "pca")) {
    xlab <- "PCA_1"
    ylab <- "PCA_2"
  } else {
    xlab <- "UMAP_1"
    ylab <- "UMAP_2"
  }

  p <- DimPlot(
    seurat_obj,
    reduction = reduction,
    group.by = group.by,
    split.by = split.by,
    cols = colors,
    pt.size = pt.size,
    label = label,
    label.size = label.size,
    repel = TRUE,
    raster = FALSE
  ) +
    xlab(xlab) +
    ylab(ylab) +
    labs(color = legend.title) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      axis.line = element_line(
        # 定义坐标轴的箭头样式
        arrow = arrow(type = "closed", length = unit(0.1, "inches"), angle = 20)
      ),
      axis.title = element_text(hjust = 0.05, face = "italic", size = 12),
      legend.title = element_text(size = 12),
      plot.title = element_blank()
    ) +
    scale_x_continuous(breaks = NULL) + # 取消坐标轴刻度线
    scale_y_continuous(breaks = NULL) +
    guides(x = axis, y = axis) # 定义短坐标轴
  return(p)
}




# 基于ggplot的自定义UMAP图，圈出了各个亚群的边界
ggDimPlot <- function(
    seurat_obj,
    reduction,
    colors = brewer.pal(name = "Paired", n = 12),
    pt.size = 0.1,
    legend.title = "Cell type") {
  # 提取UMAP坐标信息及分群信息
  df <- seurat_obj@reductions[[reduction]]@cell.embeddings %>%
    bind_cols(Idents(seurat_obj))
  colnames(df) <- c("umap_1", "umap_2", "cluster")

  # 定义截短坐标轴，便于绘制短坐标轴UMAP图
  axis <- guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(3, "cm")
  )

  p <- ggplot(df, aes(umap_1, umap_2, color = cluster)) +
    geom_point(size = pt.size, show.legend = FALSE) +
    # 圈出各个亚群
    stat_unchull(
      aes(color = cluster, fill = cluster),
      alpha = 0.25,
      linewidth = 0.25,
      lty = 2,
      delta = 0.3 # 圈的边距
    ) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    xlab("UMAP_1") +
    ylab("UMAP_2") +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      axis.line = element_line(
        # 定义坐标轴的箭头样式
        arrow = arrow(
          type = "closed",
          length = unit(0.1, "inches"),
          angle = 20
        )
      ),
      axis.title = element_text(hjust = 0.05, face = "italic"),
      legend.title = element_text(size = 12)
    ) +
    scale_x_continuous(breaks = NULL) + # 取消坐标轴刻度线
    scale_y_continuous(breaks = NULL) +
    # 重新绘制坐标轴和图例
    guides(
      fill = guide_legend(
        override.aes = list(alpha = 1), # 重新指定图例的透明度
        title = legend.title
      ),
      color = guide_legend(
        override.aes = list(lty = 1), # 重新指定图例的线型
        title = legend.title
      ),
      x = axis, # 定义短坐标轴
      y = axis
    )
  return(p)
}





# 细胞亚群质量评估的小提琴图 ---------------------------------------------------------------------

sc.cluster.VlnPlot <- function(seurat_obj) {
  plot_violin <- VlnPlot(
    seurat_obj,
    features = c("nCount_RNA", "nFeature_RNA"),
    cols = mycolors,
    ncol = 1,
    raster = FALSE
  ) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.title = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )
  print(plot_violin)
  return(plot_violin)
}





# 不同分辨率下的UMAP图的绘制 -------------------------------------------------------------------

sc.DimPlot.res <- function(seurat_obj, colors = mycolors_many, prefix = "RNA_snn_res") {
  plot_umap_res <- map(
    grep(prefix, colnames(seurat_obj@meta.data), value = TRUE),
    function(res) {
      Idents(seurat_obj) <- res
      DimPlot(
        seurat_obj,
        reduction = "umap.integrated",
        cols = colors,
        label = TRUE,
        label.size = 4,
        repel = TRUE,
        raster = FALSE
      ) +
        ggtitle(res) +
        xlab("UMAP_1") +
        ylab("UMAP_2") +
        theme_bw() +
        theme(legend.position = "none")
    }
  ) |>
    wrap_plots(ncol = 3)
  print(plot_umap_res)
  return(plot_umap_res)
}
