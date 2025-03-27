{
  library(tidyverse)
  library(Seurat)
  library(janitor)
  library(qs)
}


files <- list.files("data/GSE154600_RAW", pattern = "^GSM")
files
dir_names <- str_split(files, "_", simplify = TRUE)[, 1] %>%
  unique()
dir_names

dir_names <- str_c("data/GSE154600_RAW/", dir_names)
dir_names
walk(dir_names, dir.create)

walk(files, function(files) {
  # 复制各个10X文件到其对应的样本文件夹内
  new_dir <- str_split(files, "_", simplify = TRUE)[, 1] %>%
    str_c("data/GSE154600_RAW/", .)
  file.copy(
    from = str_c("data/GSE154600_RAW/", files),
    to = new_dir
  )
  # 重命名3个10X文件
  new_filename <- case_when(
    str_detect(files, "barcodes") ~ "barcodes.tsv.gz",
    str_detect(files, "genes") ~ "features.tsv.gz",
    .default = "matrix.mtx.gz"
  )
  file.rename(
    from = file.path(new_dir, files),
    to = file.path(new_dir, new_filename)
  )
})



# 构建 Seurat 对象 ----------------------------------------------------------------------

seurat_list <- map(dir_names, function(dir_names) {
  print(dir_names)
  sce <- CreateSeuratObject(
    counts = Read10X(dir_names),
    project = str_split(dir_names, "/", simplify = TRUE)[, 3],
    min.cells = 5
  )
  return(sce)
})
head(seurat_list[[1]])

merged_seurat <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = str_split(dir_names, "/", simplify = T)[, 3]
)
merged_seurat

merged_seurat <- JoinLayers(merged_seurat)
merged_seurat
head(merged_seurat)

# 将样本标识添加到meta.data的新的一列“sample”中:
merged_seurat$sample_id <- merged_seurat$orig.ident
# 添加研究编号，便于合并后识别
merged_seurat$study_id <- "GSE154600"

# 添加SRA Run Selector中提供的样本信息
add_meta.data <- read_csv("data/GSE154600_RAW/SraRunTable.txt") %>%
  clean_names() # 规范列名
add_meta.data %>% glimpse()
add_meta.data <- add_meta.data %>%
  select(
    sample_id = sample_name,
    chemotherapy_response,
    tumor_stage)
add_meta.data

new_meta.data <- left_join(
  merged_seurat@meta.data,
  add_meta.data,
  by = "sample_id"
)
# 添加样本类型列
new_meta.data <- new_meta.data %>%
  mutate(
    sample_type = "Cancer",
    .after = "sample_id"
  )
head(new_meta.data)
# 检查合并后的数据的行数是否等于细胞数
nrow(new_meta.data) == ncol(merged_seurat)
# 添加行名
rownames(new_meta.data) <- rownames(merged_seurat@meta.data)
# 将样本信息添加回meta.data中
merged_seurat@meta.data <- new_meta.data
head(merged_seurat)
table(merged_seurat$chemotherapy_response)
table(merged_seurat$tumor_stage)

qsave(merged_seurat, "output/seurat_obj/GSE154600_merged_seurat.qs")



# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/0.2_GSE154600.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
