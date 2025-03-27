# 加载包 -------------------------------------------------------------------------------

{
  library(Seurat)
  library(tidyverse)
  library(janitor)
  library(qs)
}



# 数据读取 ------------------------------------------------------------------------------

gse <- "GSE173682"
data_path <- str_c("data/", gse, "_RAW")

files <- list.files(data_path, pattern = "^GSM")
files

dir_names <- str_split(files, "_", simplify = TRUE)[, 1] %>%
  unique()
dir_names

dir_names <- file.path(data_path, dir_names)
dir_names
walk(dir_names, dir.create)

walk(files, function(file) {
  # 复制各个10X文件到其对应的样本文件夹内
  new_dir <- str_split(file, "_", simplify = TRUE)[, 1] %>%
    file.path(data_path, .)
  file.copy(
    from = file.path(data_path, file),
    to = new_dir
  )
  # 重命名3个10X文件
  new_filename <- case_when(
    str_detect(file, "barcode") ~ "barcodes.tsv.gz",
    str_detect(file, "gene|feature") ~ "features.tsv.gz",
    .default = "matrix.mtx.gz"
  )
  file.rename(
    from = file.path(new_dir, file),
    to = file.path(new_dir, new_filename)
  )
})


# 批量读取
seurat_list <- map(dir_names, function(dir_name) {
  print(dir_name)
  sce <- CreateSeuratObject(
    counts = Read10X(dir_name),
    project = str_split(dir_name, "/", simplify = TRUE)[, 3],
    min.cells = 5
  )
  return(sce)
})
seurat_list[[1]]
head(seurat_list[[1]])

merged_seurat <- merge(
  seurat_list[[1]],
  seurat_list[-1],
  add.cell.ids = str_split(dir_names, "/", simplify = TRUE)[, 3]
)
merged_seurat
head(merged_seurat)

# 合并layers
merged_seurat <- JoinLayers(merged_seurat)

# 将样本标识添加到meta.data的新的一列“sample_id”中:
merged_seurat$sample_id <- merged_seurat$orig.ident
# 添加研究编号，便于合并后识别
merged_seurat$study_id <- gse








# 添加SRA Run Selector中提供的样本信息（这里只下载了上面两个样本对应的SraRunTable文件）
add_meta.data <- read_csv(file.path(data_path, "SraRunTable.csv")) %>%
  clean_names() # 规范列名
add_meta.data %>% glimpse()
add_meta.data <- add_meta.data %>%
  select(
    sample_id = library_name,
    sample_type = histological_type,
  ) %>%
  distinct()
add_meta.data

# 根据原始文献修改和添加列
add_meta.data <- add_meta.data %>%
  mutate(
    age = if_else(str_detect(sample_id, "3BAE2L"), 61, 59),
    tumor_stage = if_else(str_detect(sample_id, "3BAE2L"), "II", "III"),
    sample_id = if_else(str_detect(sample_id, "^Patient_8"), "GSM5276940", "GSM5276943")
  )
add_meta.data

new_meta.data <- left_join(
  merged_seurat@meta.data,
  add_meta.data,
  by = "sample_id"
)
head(new_meta.data)
tail(new_meta.data)
# 检查合并后的数据的行数是否等于细胞数
nrow(new_meta.data) == ncol(merged_seurat)
# 添加行名
rownames(new_meta.data) <- rownames(merged_seurat@meta.data)

# 将样本信息添加回meta.data中
merged_seurat@meta.data <- new_meta.data
head(merged_seurat)
tail(merged_seurat)
table(merged_seurat$tumor_stage)
table(merged_seurat$sample_id)

qsave(merged_seurat, str_c("output/seurat_obj/", gse, "_merged_seurat.qs"))




# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/0.3_GSE173682.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
