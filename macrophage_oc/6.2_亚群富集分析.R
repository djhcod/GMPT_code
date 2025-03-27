# 加载包 ---------------------------------------------------------------------

{
  rm(list = ls())
  cat("\014")

  library(clusterProfiler)
  library(msigdbr) # 注释基因集
  library(Seurat)
  library(GSVA)
  library(GOplot) # 绘制弦图
  library(tidyverse)
  library(pheatmap)
  library(cowplot)
  library(patchwork)
  library(grid) # 保存pheatmap绘制的热图
  library(scatterplot3d) # 绘制3D UMAP图
  library(qs)
  library(ggpubr) # 给ggplot图形添加统计指标

  # 载入颜色集
  load("output/color_palette.rdata")
}



# 普通ORA富集分析 ----------------------------------------------------------

## 差异基因准备 ---------------------------------------------------------

# 载入差异基因列表
degs <- read_csv("output/macro/cluster_degs.csv") %>% as_tibble()
degs
degs <- rename(degs, SYMBOL = gene)

# 筛选上调的degs
up_degs <- degs %>%
  filter(p_val_adj < 0.05, avg_log2FC >= 0.25)
up_degs
# 统计每个cluster的degs数量
up_degs %>%
  group_by(cluster) %>%
  summarise(n())

# 将差异基因转换为列表，列表中的每个对象为每个cluster
up_degs_list <- split(up_degs$SYMBOL, f = up_degs$cluster)
str(up_degs_list)




## GO富集分析 --------------------------------------------------------------

go <- compareCluster(
  up_degs_list,
  fun = "enrichGO",
  OrgDb = "org.Hs.eg.db",
  keyType = "SYMBOL", # 指定输入的基因ID类型为SYMBOL
  ont = "ALL", # 一次性得到"BP", "MF", 和 "CC"
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

# 提取富集结果
go_result <- go@compareClusterResult %>% as_tibble()
go_result
write_csv(go_result, "output/macro/GO_enrichment_result.csv")

# 提取每个cluster的前20个富集通路
top_go_result <- go_result %>%
  group_by(Cluster, ONTOLOGY) %>%
  slice_max(order_by = Count, n = 20) %>%
  slice_min(
    order_by = p.adjust,
    n = 20,
    with_ties = FALSE # 并列的数据只取前n个
  ) %>%
  ungroup()
top_go_result
# 统计各个cluster的top富集通路数量
table(top_go_result$Cluster, top_go_result$ONTOLOGY)
write_csv(top_go_result, "output/macro/Top_GO_enrichment_result.csv")






## 可视化 -----------------------------------------------------------------

### 气泡图-----------------------------------------------------
# 富集分析气泡图自编函数
enrich.dotplot <- function(dat, title = "Enrichment result") {
  # 让文字过长时换行显示
  dat$Description <- str_wrap(dat$Description, width = 50)

  ggplot(dat, aes(x = Count, y = reorder(Description, Count))) +
    geom_point(aes(size = Count, color = -log10(p.adjust))) +
    scale_colour_gradient2(
      low = mycolors[1],
      mid = mycolors[2],
      high = mycolors[3],
      midpoint = median(-log10(dat$p.adjust))
    ) +
    scale_x_continuous(expand = c(0.2, 0)) +
    scale_y_discrete(expand = c(0.05, 0)) +
    labs(
      x = "Gene Numbers",
      y = NULL,
      size = "Count"
    ) +
    ggtitle(title) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 10),
      axis.text.y = element_text(size = 6)
    )
}

# 批量绘制并合并各cluster的top富集通路气泡图
p_bp <- map(
  unique(top_go_result$Cluster),
  function(cluster) {
    top_go_result %>%
      filter(Cluster == cluster, ONTOLOGY == "BP") %>%
      enrich.dotplot(title = str_c("Cluster ", cluster))
  }
) %>%
  wrap_plots(ncol = 5)
p_bp
ggsave(width = 28, height = 13, "figures/macro/GO_BP_dotplots.pdf")

p_cc <- map(
  unique(top_go_result$Cluster),
  function(cluster) {
    top_go_result %>%
      filter(Cluster == cluster, ONTOLOGY == "CC") %>%
      enrich.dotplot(title = str_c("Cluster ", cluster))
  }
) %>%
  wrap_plots(ncol = 5)
p_cc
ggsave(width = 28, height = 13, "figures/macro/GO_CC_dotplots.pdf")

p_mf <- map(
  unique(top_go_result$Cluster),
  function(cluster) {
    top_go_result %>%
      filter(Cluster == cluster, ONTOLOGY == "MF") %>%
      enrich.dotplot(title = str_c("Cluster ", cluster))
  }
) %>%
  wrap_plots(ncol = 5)
p_mf
ggsave(width = 28, height = 13, "figures/macro/GO_MF_dotplots.pdf")



### 弦图（chord plot） -------------------------------------------------

# 富集结果弦图绘制自编函数
enrich.chordplot <- function(
    dat, # 富集结果
    degs, # 差异基因列表
    group = 0, # 绘制哪个cluster的富集结果
    ontology = "BP", # 基因集
    colors = mycolors_many) {
  # 按要求整理富集结果数据
  enrich_tab <- dat %>%
    filter(Cluster == group, ONTOLOGY == ontology) %>%
    select(
      Category = ONTOLOGY,
      ID,
      Term = Description,
      Genes = geneID,
      adj_pval = p.adjust
    ) %>%
    mutate(
      Genes = str_replace_all(Genes, "/", ", ")
    )

  # 按要求整理差异基因列表
  chord_genelist <- degs %>%
    filter(cluster == group) %>%
    select(ID = SYMBOL, logFC = avg_log2FC)

  # 需要绘制的基因（这里所有基因都需要绘制）
  chord_genes <- chord_genelist %>% as.data.frame()

  # 需要绘制的感兴趣的通路（这里绘制所有通路）
  chord_process <- enrich_tab$Term

  # 生成弦图绘制的数据
  circ <- circle_dat(enrich_tab, chord_genelist)
  chord <- chord_dat(circ, chord_genes)

  p_chord <- GOChord(
    chord,
    space = 0.02, # 圆圈各部分的间距
    gene.order = "logFC",
    gene.space = 0.2, # 基因离圆圈的距离
    gene.size = 3,
    lfc.col = colors[2:3], # logFC的颜色映射
    ribbon.col = colors[1:length(chord_process)], # 各个通路的颜色
    border.size = 0, # 描边粗细
    process.label = 3 # 图例大小
  ) +
    theme(legend.text = element_text(size = 8))
  return(p_chord)
}


# 批量绘制并合并各cluster的top富集通路弦图
p_chord_bp <- map(
  unique(top_go_result$Cluster),
  ~ enrich.chordplot(
    top_go_result,
    degs = up_degs,
    group = .x,
    ontology = "BP"
  )
) %>%
  wrap_plots(ncol = 5)
ggsave(
  width = 8 * 10, height = 17 * 2, "figures/macro/GO_BP_chord_plots.pdf",
  limitsize = FALSE
)







# GSVA富集分析 ----------------------------------------------------------------

## 整理注释背景基因集-----------------------------------------------------
if (FALSE) {
  msigdbr_species() # 查看所有的物种
  msigdbr_collections() %>% print(n = Inf) # 查看所有的注释基因集信息
}

# 提取人类GO注释基因集（C5）
go_geneset <- msigdbr(species = "Homo sapiens", category = "C5")
head(go_geneset)
unique(go_geneset$gs_subcat)

# 按要求整理注释基因列表的自编函数
genelist.prepare <- function(geneset, subcat) {
  geneset <- geneset %>%
    filter(gs_subcat == subcat) %>%
    select(gs_name, gene_symbol) %>%
    mutate(
      gs_name = str_remove(subcat, ":") %>%
        str_c("^", ., "_") %>%
        str_remove(gs_name, pattern = .) %>%
        str_replace_all("_", " ") %>%
        str_to_sentence()
    )
  str_c(
    "提取的", subcat, "基因集",
    "-----------------------------------------------------------------------"
  ) %>% print()
  geneset %>% print()
  genelist <- split(geneset$gene_symbol, geneset$gs_name)
  return(genelist)
}

# 批量整理得到"BP", "MF", 和 "CC"的注释基因列表
go_genelist <- map(
  c("GO:BP", "GO:CC", "GO:MF"),
  ~ genelist.prepare(
    go_geneset,
    subcat = .x
  )
)
names(go_genelist) <- c("GO:BP", "GO:CC", "GO:MF")



## 构建单细胞pseudobulk表达矩阵--------------------------------------------

# 载入Seurat对象
seurat_annotated <- qread("output/macro/seurat_annotated.qs")
seurat_annotated
levels(seurat_annotated)

# 执行pseudobulk
pseudobulk <- AggregateExpression(
  seurat_annotated,
  group.by = "cell_type"
)[["RNA"]]

pseudobulk[1:4, 1:4]
dim(pseudobulk)
summary(pseudobulk[, 1])


# 去除在所有cluster中表达量均为0的基因
pseudobulk <- pseudobulk[rowSums(pseudobulk) > 0, ]
dim(pseudobulk)



# GSVA自编函数
sc.gsva <- function(geneset) {
  gsva_result <- gsva(gsvaParam(pseudobulk, geneset))

  gsva_result <- as.data.frame(gsva_result)
  str_c(
    deparse(substitute(geneset)),
    "的GSVA结果数据部分展示",
    "-----------------------------------------------------------------------"
  ) %>% print()
  gsva_result[1:4, 1:4] %>% print()

  str_c(
    deparse(substitute(geneset)),
    "的GSVA结果数据行列数",
    "-----------------------------------------------------------------------"
  ) %>% print()
  dim(gsva_result) %>% print()

  return(gsva_result)
}

# 批量进行"BP", "MF", 和 "CC"的GSVA分析
gsva_result <- map(
  names(go_genelist),
  ~ sc.gsva(go_genelist[[.x]])
)
names(gsva_result) <- names(go_genelist)






## 可视化-----------------------------------------------------------------

# GSVA TOP富集通路热图展示自编函数
gsva.pheatmap <- function(gsva_result, top_n = 5) {
  top_pathway <- c()
  for (i in 1:ncol(gsva_result)) {
    top_pathway <- gsva_result[order(-gsva_result[, i]), ] %>%
      .[1:top_n, ] %>%
      rownames() %>%
      c(top_pathway, .)
  }
  top_pathway <- unique(top_pathway)

  p <- gsva_result[rownames(gsva_result) %in% top_pathway, ] %>%
    pheatmap(
      scale = "column",
      cluster_cols = FALSE,
      show_colnames = T,
      angle_col = "45",
      border_color = "white"
    )
  return(p)
}





# 保存pheatmap图片的函数
save_pheatmap_pdf <- function(x, filename, width = 7, height = 7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid.newpage()
  grid.draw(x$gtable)
  dev.off()
}


walk(
  names(gsva_result),
  function(x) {
    p <- gsva.pheatmap(gsva_result[[x]], top_n = 5)
    save_pheatmap_pdf(
      p,
      width = 11,
      height = 10,
      filename = str_replace(x, ":", "_") %>%
        str_c("figures/macro/GSVA_", ., "_heatmap.pdf")
    )
  }
)






# 亚群生物学命名 -----------------------------------------------------------------

abbrev <- tribble(
  ~abbreviation, ~full_name,
  "LP-TAMs", "Lipid-associated TAMs",
  "Reg-TAMs", "Immune regulatory TAMs",
  "Inflam-TAMs", "Inflammatory cytokine-enriched TAMs",
  "Prolif-TAMs", "Proliferating TAMs",
  "MT-TAMs", "Metallothionein-associated TAMs",
  "RTM-TAMs", "Resident-tissue macrophages-like TAMs",
  "HSP-TAMs", "Heat shock protein-enriched TAMs",
  "RB-TAMs", "Ribosome biogenesis-associated TAMs",
  "IFN-TAMs", "Interferon-primed TAMs"
)
abbrev
write_csv(abbrev, "output/macro/Macrophage_Subclusters_Abbreviations.csv")

seurat_annotated <- RenameIdents(
  seurat_annotated,
  "Macro 0" = "APOE+ LP-TAMs",
  "Macro 1" = "HLA-DPB1+ Reg-TAMs",
  "Macro 2" = "CCL3L1+ Inflam-TAMs",
  "Macro 3" = "TOP2A+ Prolif-TAMs",
  "Macro 4" = "MT1H+ MT-TAMs",
  "Macro 5" = "LYVE1+ RTM-TAMs",
  "Macro 6" = "TNFAIP6+ Inflam-TAMs",
  "Macro 7" = "HSPA6+ HSP-TAMs",
  "Macro 8" = "LILRA4+ RB-TAMs",
  "Macro 9" = "CXCL10+ IFN-TAMs"
)
levels(seurat_annotated)

# 把细胞注释信息保存到meta.data中
seurat_annotated$cell_annotation <- Idents(seurat_annotated)
head(seurat_annotated, 3)
table(seurat_annotated$cell_type, seurat_annotated$cell_annotation)

qsave(seurat_annotated, "output/macro/seurat_bio_annotated.qs")





## 绘制UMAP图 -----------------------------------------------------------------

# 加载分群相关函数
source("function/cluster_functions.r")

umap_plot <- ggDimPlot(
  seurat_annotated,
  reduction = "umap.integrated",
  colors = mycolors,
  legend.title = "Cluster ID"
)
umap_plot

# 绘制3D UMAP图
umap_3d_tab <- seurat_annotated@reductions$umap.3d@cell.embeddings %>%
  as.data.frame() %>%
  mutate(
    cell_type = seurat_annotated$cell_annotation,
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
ggsave("figures/macro/umap_plot_bio_annotated.pdf", width = 11, height = 5)








# 细胞比例可视化----------------------------------------------------------------

seurat_annotated <- qread("output/macro/seurat_bio_annotated.qs")

# 统计各巨噬细胞亚群的占比情况

p_1 <- seurat_annotated@meta.data %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = unique(sample_id[order(sample_type, study_id)])
    )
  ) %>%
  ggplot(
    aes(y = sample_id, fill = cell_annotation)
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
p_1

p_2 <- seurat_annotated@meta.data %>%
  group_by(study_id, cell_annotation) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(study_id) %>%
  mutate(proportion = count / sum(count)) %>%
  ggplot(aes(x = study_id, y = proportion, fill = cell_annotation)) +
  geom_bar(stat = "identity", width = 0.8, position = "fill") +
  scale_fill_manual(values = mycolors) +
  xlab("Study ID") +
  ylab("Ratio") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
p_2

p_3 <- ggplot(
  data = seurat_annotated@meta.data,
  aes(x = sample_type, fill = cell_annotation)
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
p_3

plot_grid(
  p_1 + NoLegend(),
  p_2 + NoLegend(),
  p_3 + NoLegend(),
  get_legend(p_1),
  nrow = 1,
  rel_widths = c(4, 1.5, 1.5, 1),
  axis = "b",
  align = "h"
)

ggsave("figures/macro/cell_type_ratio.pdf", width = 14, height = 6)




# 统计良、恶性样本间细胞比例的差异

p_prop <- seurat_annotated@meta.data %>%
  as_tibble() %>%
  group_by(sample_id, sample_type, cell_annotation) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(total_count = sum(count)) %>%
  mutate(proportion = count / total_count) %>%
  ggplot(aes(x = cell_annotation, y = proportion * 100, fill = sample_type)) +
  geom_boxplot() +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(0, 3)) +
  xlab("Cell type") +
  ylab("Proportion (%)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  stat_compare_means(
    aes(group = sample_type),
    method = "wilcox.test",
    vjust = -1
  )
p_prop

ggsave(
  width = 17.2,
  height = 5.3,
  "figures/macro/cell_type_ratio_between_normal_and_OC_samples.pdf"
)




# 保存sessionInfo -----------------------------------------------------------

sink(
  "sessionInfo/6.2_Macro_subclusters_DEGs_enrichment.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()

