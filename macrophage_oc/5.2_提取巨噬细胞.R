# 加载包 ---------------------------------------------------------------------

rm(list = ls())
cat("\014")
{
  library(qs)
  library(Seurat)
  library(tidyverse)
  library(cowplot)
}
# 载入颜色集
load("output/color_palette.rdata")

# 载入数据 ------------------------------------------------------------------------------

seurat_annotated <- qread("output/seurat_obj/seurat_annotated.qs")
seurat_annotated
head(seurat_annotated)
levels(seurat_annotated)




# 统计良、恶性样本间细胞类型及比例的差异 -----------------------------------------------------

p1 <- seurat_annotated@meta.data %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = unique(sample_id[order(sample_type, study_id)])
    )
  ) %>%
  ggplot(
    aes(y = sample_id, fill = cell_type)
  ) +
  geom_bar(width = 0.8, position = "fill") +
  scale_fill_manual(values = mycolors) +
  xlab("Ratio") +
  ylab("Sample ID") +
  labs(fill = "Cell type") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  facet_wrap(~sample_type)
p1

p2 <- ggplot(
  data = seurat_annotated@meta.data,
  aes(x = sample_type, fill = cell_type)
) +
  geom_bar(width = 0.8, position = "fill") +
  scale_fill_manual(values = mycolors) +
  xlab("Sample type") +
  ylab("Ratio") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
p2

p3 <- seurat_annotated@meta.data %>%
  filter(is.na(tumor_stage) == FALSE) %>%
  ggplot(aes(x = tumor_stage, fill = cell_type)) +
  geom_bar(width = 0.8, position = "fill") +
  scale_fill_manual(values = mycolors) +
  xlab("FIGO stage") +
  ylab("Ratio") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
p3

plot_grid(
  p1 + NoLegend(),
  p2 + NoLegend(),
  p3 + NoLegend(),
  get_legend(p1),
  nrow = 1,
  rel_widths = c(4, 1, 1.5, 1)
)

ggsave("figures/seurat_pipeline/cell_type_ratio.pdf", width = 14, height = 6)




# 提取巨噬细胞 ------------------------------------------------------------------

seurat_macro <- subset(seurat_annotated, cell_type == "Macrophages")
seurat_macro
colnames(seurat_macro@meta.data)

# 去除不需要的meta.data列
new_meta.data <- select(seurat_macro@meta.data, -c(unintegrated_clusters:cell_type))

# 重新构建Seurat对象
seurat_macro <- CreateSeuratObject(
  counts = seurat_macro[["RNA"]]$counts,
  meta.data = new_meta.data
)
seurat_macro
head(seurat_macro)

qsave(seurat_macro, "output/seurat_obj/seurat_macro.qs")





# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/5.2_Extract_macrophages.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()

