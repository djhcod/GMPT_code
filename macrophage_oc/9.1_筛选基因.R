# 加载包及数据 ---------------------------------------------------------------------

{
  rm(list = ls())
  cat("\014")

  library(qs)
  library(tidyverse)
  library(Seurat)
  library(ggpubr) # 绘制45度火山图
  library(ggrepel) # 图片标注

  # 载入颜色集
  load("output/color_palette.rdata")

  # 指定输出目录
  out_dir <- "output/identify_candidate_genes/"
  fig_dir <- "figures/identify_candidate_genes/"
}


# 载入完成了巨噬细胞细分注释的Seurat对象
seurat_annotated <- qread("output/macro/seurat_bio_annotated.qs", nthreads = 10)
seurat_annotated
table(Idents(seurat_annotated))
colnames(seurat_annotated@meta.data)
table(seurat_annotated$sample_type)



# 指定需要分析的目标细胞亚群
aim_cluster <- "CCL3L1+ Inflam-TAMs"


# 确定差异基因筛选的阈值
log2FC_cut <- 0.25
adjp_cut <- 0.05







# 寻找目标巨噬细胞亚群和其他巨噬细胞亚群之间的差异基因------------------------------

deg_tab <- FindMarkers(
  seurat_annotated,
  ident.1 = aim_cluster,
  logfc.threshold = 0,
  min.pct = 0.1,
  only.pos = FALSE
)
head(deg_tab)



deg_tab <- deg_tab %>%
  rownames_to_column("gene") %>%
  tibble() %>%
  mutate(
    status = case_when(
      avg_log2FC < -log2FC_cut & p_val_adj < adjp_cut ~ "Down",
      avg_log2FC > log2FC_cut & p_val_adj < adjp_cut ~ "Up",
      .default = "Unchanged"
    ) %>%
      factor()
  )
summary(deg_tab$status)

write_csv(
  deg_tab,
  str_c(out_dir, "DEGs_between_", aim_cluster, "_and_other_macrophage_subtypes.csv")
)


top_degs <- deg_tab %>%
  filter(status != "Unchanged") %>%
  group_by(status) %>%
  slice_max(abs(avg_log2FC), n = 10)
top_degs


# 火山图
p_vol <- ggscatterhist(
  deg_tab,
  x = "pct.1",
  y = "pct.2",
  color = "status",
  palette = c(mycolors[1], "gray", mycolors[2]),
  size = 1,
  alpha = 0.8,
  ggtheme = theme_bw(),
  margin.plot = "histogram",
  margin.params = list(fill = "gray", color = NA),
  xlab = "Pct.1",
  ylab = "Pct.2"
)
p_vol$sp <- p_vol$sp +
  geom_label_repel(
    data = top_degs,
    aes(x = pct.1, y = pct.2, label = gene, color = status),
    size = 2,
    box.padding = 0.5,
    label.padding = 0.1,
    label.size = 0.2,
    max.overlaps = 100,
    show.legend = FALSE
  ) +
  geom_hline(yintercept = 0.01, color = mycolors[3]) +
  geom_vline(xintercept = 0.01, color = mycolors[3])
p_vol

ggsave(
  width = 6,
  height = 6,
  str_c(fig_dir, "DEGs_between_", aim_cluster, "_and_other_macrophage_subtypes.pdf")
)

markers <- deg_tab %>% filter(status == "Up") %>% pull(gene)


write_csv(
  deg_tab %>% filter(status == "Up"),
  str_c(out_dir, "markers.csv"),
  col_names = FALSE
)


qsave(markers, str_c(out_dir, "markers.qs"))












# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/9.1_Identify_candidate_genes.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()

