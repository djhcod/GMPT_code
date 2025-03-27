{
  library(ArrayExpress)
  library(tidyverse)
  library(janitor)
  library(data.table)
  library(Seurat)
  library(qs)
}

# 开全局代理运行
if (FALSE) {
  getAE(
    "E-MTAB-8107",
    path = "data/E-MTAB-8107",
    type = "processed"
  )
}

# SRDF文件（临床信息）是以制表符分隔的表格
srdf <- read_delim(
  "data/E-MTAB-8107/metadata/E-MTAB-8107.sdrf.txt",
  delim = "\t"
) %>%
  clean_names() # 规范列名
srdf %>% glimpse()
# 选择需要提取的信息
srdf <- srdf %>% select(
  source_name, characteristics_individual, characteristics_organism_part,
  characteristics_sampling_site, characteristics_metastatic_site
)
map(srdf, unique)
table(srdf$characteristics_organism_part)
# 筛选其中代表卵巢癌的样本
srdf <- srdf %>%
  filter(characteristics_organism_part == "ovary") %>%
  select(!characteristics_organism_part)
srdf
# 统计每个患者的样本情况
srdf %>%
  group_by(characteristics_individual, characteristics_sampling_site) %>%
  summarise(n())
# 保存SRDF文件
qsave(srdf, "output/seurat_obj/EMTAB8107_srdf.qs")



# 读取IDF文件（数据集基本信息）
idf <- read_delim(
  "data/E-MTAB-8107/metadata/E-MTAB-8107.idf.txt",
  col_names = FALSE,
  delim = "\t"
)
idf %>% print(n = Inf)

# 提取对实验总体设计的描述
idf %>%
  filter(X1 == "Experiment Description") %>%
  select(X2) %>%
  as.character()
# 保存IDF文件
qsave(idf, "output/seurat_obj/EMTAB8107_idf.qs")



# 读取卵巢癌的表达量数据
files <- list.files("data/E-MTAB-8107")
files
# 根据SRDF文件的“source_name”一列可以从这些表达量数据中筛选代表卵巢癌的文件
ov_files <- map(srdf$source_name, ~ str_detect(files, .x) %>% files[.])
ov_files
# 试验性读取其中一个文件
test <- fread(str_c("data/E-MTAB-8107/", ov_files[1]))
test[1:4, 1:4]
# 发现需要把第一列转换成行名
rm(test)

# 批量读取表达矩阵并构建Seurat对象
seurat_list <- map(ov_files, function(ov_files) {
  print(ov_files)
  count <- fread(
    str_c("data/E-MTAB-8107/", ov_files),
    data.table = FALSE # 不要以data.table格式读取，因为需要添加行名（data.table不支持行名）
  )
  rownames(count) <- count[, 1]
  count[, 1] <- NULL
  seurat <- CreateSeuratObject(
    counts = count,
    project = str_split(ov_files, "\\.", simplify = TRUE)[, 1],
    min.cells = 5
  )
  return(seurat)
})
head(seurat_list[[1]])

merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
)
merged_seurat

merged_seurat <- JoinLayers(merged_seurat)
merged_seurat
head(merged_seurat)
table(merged_seurat$orig.ident)



# 根据此前的SRDF文件将额外的样本信息添加到Seurat的meta.data中
head(srdf)
new_meta.data <- left_join(
  merged_seurat@meta.data,
  srdf,
  join_by(orig.ident == source_name)
)
# 检查合并后的数据的行数是否等于细胞数
nrow(new_meta.data) == ncol(merged_seurat)
# 添加行名
row.names(new_meta.data) <- rownames(merged_seurat@meta.data)
head(new_meta.data)
# 重命名列及添加新列
new_meta.data <- new_meta.data %>%
  mutate(
    study_id = "E-MTAB-8107", # 添加研究编号，便于合并后识别
    sample_id = orig.ident, # 添加样本标识
    patient_id = str_replace(characteristics_individual, "^.", "E-MTAB-8107-"),
    sample_type = characteristics_sampling_site,
    metastatic = if_else(
      characteristics_metastatic_site == "not applicable",
      "Local",
      "omentum/peritoneum metastases"
    )
  ) %>%
  select(!c(
    characteristics_individual,
    characteristics_sampling_site,
    characteristics_metastatic_site
  ))
head(new_meta.data)
table(new_meta.data$sample_id)
table(new_meta.data$sample_type)
table(new_meta.data$metastatic)

# 将样本信息添加回meta.data中
merged_seurat@meta.data <- new_meta.data

qsave(merged_seurat, "output/seurat_obj/EMTAB8107_merged_seurat.qs")





# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/0.1_E-MTAB-8107.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
