# 加载包及数据 ---------------------------------------------------------------------

{
  rm(list = ls())
  cat("\014")

  library(Seurat)
  library(qs)
  library(tidyverse)
  library(cowplot)
  library(ggalign)

  # 载入颜色集
  load("output/color_palette.rdata")
}

# 载入细胞比例信息
cell_prop <- qread("output/BayesPrism/BayesPrism_cell_propotions.qs") %>%
  as.data.frame()
head(cell_prop)

# 规范化列名，便于引用
colnames(cell_prop) <- str_remove(colnames(cell_prop), "\\+") %>%
  str_replace_all(" |-", "_")

cell_prop <- cell_prop %>%
  rownames_to_column("sample_id") %>%
  tibble()
cell_prop
summary(cell_prop)

# 检查每个样本的细胞比例总和是否都为1
rowSums(cell_prop[, -1]) %>% summary()




# 提取此前的亚群顺序，以保持绘图时颜色分配的统一
seurat_annotated <- qread("output/macro/seurat_bio_annotated.qs")

labels <- unique(
  seurat_annotated$cell_annotation[order(seurat_annotated$cell_type)]
)
labels

label_levels <- labels %>%
  str_remove("\\+") %>%
  str_replace_all(" |-", "_")
label_levels




# 重新排序cell_prop的列
cell_prop <- cell_prop %>%
  relocate(c(sample_id, all_of(label_levels)))


save(
  cell_prop, label_levels, labels,
  file = "output/BayesPrism/cell_prop_and_label_levels.rdata"
)



# 可视化细胞比例 -----------------------------------------------------------------

# 细胞比例密度图
p_density <- cell_prop %>%
  pivot_longer(
    cols = colnames(cell_prop)[-1],
    names_to = "cluster",
    values_to = "prop"
  ) %>%
  mutate(
    cluster = factor(
      cluster,
      levels = label_levels,
      labels = labels
    )
  ) %>%
  ggplot(aes(x = log10(prop), color = cluster, fill = cluster)) +
  geom_density(alpha = 0.2) +
  scale_fill_manual(values = mycolors, name = "Cluster") +
  scale_color_manual(values = mycolors, name = "Cluster") +
  xlab("Log10 (proportion)") +
  ylab("Sample density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(
    legend.title = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12)
  )
p_density




# 细胞占比甜甜圈图
p_donut <- tibble(
  cluster = colnames(cell_prop[, -1]),
  sum_prop = cell_prop %>%
    select(-sample_id) %>%
    colSums()
) %>%
  arrange(-sum_prop) %>%
  mutate(
    cluster = factor(cluster, levels = label_levels, labels = labels),
    fraction = sum_prop / sum(sum_prop),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, -1)),
    label = str_c(cluster, ": ", sprintf("%.2f", fraction * 100), "%")
  ) %>%
  ggplot(aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = cluster)) +
  geom_rect() +
  scale_fill_manual(values = mycolors, name = "Cluster") +
  geom_text(
    aes(x = 0, y = 0, label = str_c(label, collapse = "\n")),
    hjust = 0.5,
    vjust = 0.7,
    color = "black",
    size = 2
  ) +
  coord_polar(theta = "y") +
  xlim(c(-1, 4)) +
  theme_void() +
  theme(legend.title = element_text(size = 12))
p_donut


# 细胞占比箱型图
p_box <- cell_prop %>%
  pivot_longer(
    cols = colnames(cell_prop)[-1],
    names_to = "cluster",
    values_to = "prop"
  ) %>%
  mutate(
    cluster = factor(
      cluster,
      levels =
        label_levels,
      labels = labels
    )
  ) %>%
  mutate(prop = log10(prop)) %>%
  ggplot(aes(x = cluster, y = prop, fill = cluster)) +
  geom_boxplot(outliers = FALSE) +
  scale_fill_manual(values = mycolors) +
  geom_jitter(color = "black", size = 0.1, alpha = 0.5) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 11)
  ) +
  xlab("Cluster") +
  ylab("Log10 (proportion)") +
  scale_y_continuous(expand = c(0, .5)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p_box

plot_grid(
  p_density,
  p_box,
  p_donut,
  labels = "AUTO",
  nrow = 1,
  rel_widths = c(100, 70, 70),
  align = "h"
)
ggsave(
  width = 15,
  height = 4,
  "figures/macro/Visualization_of_deconvolution_results.pdf"
)





# 各样本细胞占比热图

# 整理热图绘制数据
heatmap_dat <- cell_prop %>%
  column_to_rownames("sample_id") %>%
  log10() %>%
  `colnames<-`(labels)
head(heatmap_dat)

# 提取层次聚类后的cluster顺序，便于绘制热图时保证顶部条形图顺序和热图列名一致
cluster_order <- heatmap_dat %>%
  ggheatmap() +
  hmanno(position = "top") +
  align_dendro()
cluster_order <- cluster_order@data %>%
  colnames()
cluster_order


p_heatmap <- heatmap_dat %>%
  ggheatmap() +
  scale_fill_gradientn(colors = mycolors[3:1]) +
  guides(fill = guide_colorbar(title = "Log10 (proportion)")) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y = element_blank(),
    plot.margin = unit(c(0, 0, 0, 1), "cm")
  ) +
  hmanno(
    position = "top", # 激活哪个方向的注释
    size = 0.5 # 热图注释的相对高度
  ) +
  align_dendro(size = 0.3) + # 根据层次聚类或树状图排序行或列（由hmanno的position参数确定）
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  align_gg( # 添加其他ggplot图形，类似ggplot函数的作用
    data = cell_prop[,-1] %>%
      `colnames<-`(labels) %>%
      pivot_longer(
        cols = everything(),
        names_to = "cluster",
        values_to = "prop"
      ) %>%
      group_by(cluster) %>%
      summarise(median = median(prop)) %>%
      mutate(cluster = factor(cluster, levels = cluster_order)) %>%
      arrange(cluster)
  ) +
  geom_bar(aes(y = log10(median), fill = cluster), stat = "identity") +
  scale_fill_manual(values = mycolors) +
  ylab("Log10 (proportion)") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.title = element_text(size = 10))
p_heatmap

ggsave(width = 7.3, height = 8.9, "figures/macro/Heatmap_of_deconvolution_results.pdf")











# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/7.3_BayesPrism_cell_prop_visualization.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
