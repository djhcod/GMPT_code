# 加载包及数据 ---------------------------------------------------------------------

{
  rm(list = ls())
  cat("\014")

  library(BayesPrism)
  library(Seurat)
  library(qs)
  library(tidyverse)
}

# 载入 bulk RNA-seq 原始计数矩阵
bk.dat <- qread("output/tcga/TCGA_raw_count_matrix.qs")
# 载入完成了巨噬细胞细分注释的Seurat对象
seurat_bio_annotated <- qread("output/macro/seurat_bio_annotated.qs")







# 去卷积分析数据准备 ---------------------------------------------------------------

# 提取单细胞表达矩阵
sc.dat <- seurat_bio_annotated[["RNA"]]$counts
sc.dat[1:5, 1:5]
# rownames 是样本ID，colnames 是基因名称
sc.dat <- t(sc.dat)
dim(sc.dat)


# 提取 bulk RNA-seq 原始计数矩阵
bk.dat[1:4, 1:4]
# rownames 是样本ID，colnames 是基因名称
bk.dat <- t(bk.dat)
dim(bk.dat)



# 整理细胞注释标签(长度等于sc.dat的行数)
cell_type_info <- seurat_bio_annotated@meta.data %>%
  select(cell_annotation) %>%
  rownames_to_column("barcode") %>%
  as_tibble()
cell_type_info
# 转换为字符向量
cell.type.labels <- as.character(cell_type_info$cell_annotation)

head(cell.type.labels)
table(cell.type.labels)
Idents(seurat_bio_annotated) %>% table()







# 质量检查及过滤 -----------------------------------------------------------------

# 细胞类型质量检查，plotting the pairwise correlation matrix between cell types
plot.cor.phi(
  sc.dat,
  input.labels = cell.type.labels,
  title = "cell type correlation",
  cexRow = 0.5,
  cexCol = 0.5,
)


# 过滤离群基因
sc.stat <- plot.scRNA.outlier(
  sc.dat,
  cell.type.labels = cell.type.labels,
  species = "hs",
  return.raw = TRUE # return the data used for plotting
)
head(sc.stat)

bk.stat <- plot.bulk.outlier(
  bk.dat,
  sc.input = sc.dat,
  cell.type.labels = cell.type.labels,
  species = "hs",
  return.raw = TRUE
)
head(bk.stat)

# 过滤基因
sc.dat.filtered <- cleanup.genes(
  sc.dat,
  input.type = "count.matrix",
  species = "hs",
  gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
  exp.cells = 5
)
dim(sc.dat.filtered)

# 检查不同类型基因表达的一致性
plot.bulk.vs.sc(sc.dat.filtered, bulk.input = bk.dat)

# 提取蛋白编码基因
sc.dat.filtered.pc <- select.gene.type(sc.dat.filtered, gene.type = "protein_coding")





# 执行BayesPrism -------------------------------------------------------------

# 构建prism对象
myPrism <- new.prism(
  reference = sc.dat.filtered.pc,
  mixture = bk.dat,
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = NULL,
  key = NULL, # cell.type.label 中恶性细胞类型对应的字符向量
  outlier.cut = 0.05,
  outlier.fraction = 0.1,
)

# 执行BayesPrism
bp.res <- run.prism(myPrism, n.cores = 10)


bp.res
slotNames(bp.res)
theta <- get.fraction(
  bp.res,
  which.theta = "final",
  state.or.type = "type"
)
head(theta)


# 保存细胞比例信息
qsave(theta, file = "output/BayesPrism/BayesPrism_cell_propotions.qs")












# 保存sessionInfo---------------------------------------------------------
sink(
  "sessionInfo/7.2_BayesPrism.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
