# 加载环境 -------------------------------------------------------------------------------

cat("\014")
# 加载质控自编函数
source("function/qc_functions.r")


# 数据读取 ------------------------------------------------------------------------------

gse <- c("GSE184880", "EMTAB8107", "GSE154600", "GSE173682")
map(gse, function(gse) {
  sce <- str_c("output/seurat_obj/", gse, "_merged_seurat.qs") %>%
    qread()
  assign(gse, sce, envir = .GlobalEnv)
})

# 检查各数据集的临床信息
map(gse, function(gse_i) {
  sce <- get(gse_i)
  map(sce@meta.data[4:ncol(sce@meta.data)], table)
})

# 修改meta.data中不合适的变量
GSE184880@meta.data <- GSE184880@meta.data %>%
  mutate(
    sample_type = if_else(
      sample_type == "High-grade serous ovarian cancer tissue", "HGSOC", "Normal"
    ),
    tumor_stage = str_extract(tumor_stage, "^I+")
  )
EMTAB8107@meta.data <- EMTAB8107@meta.data %>%
  mutate(
    sample_type = if_else(sample_type == "neoplasm", "OC", "Normal"),
    metastatic = if_else(metastatic == "Local", "No", "Yes")
  )
GSE154600@meta.data <- GSE154600@meta.data %>%
  mutate(
    tumor_stage = if_else(tumor_stage == "IIIc", "III", "IV")
  )
# 再次检查修改后的数据集的临床信息
map(c("GSE184880", "EMTAB8107", "GSE154600"), function(gse_i) {
  sce <- get(gse_i)
  map(sce@meta.data[4:ncol(sce@meta.data)], table)
})

# 统计待评价的所有数据集的样本总数
all_sample_id <- c()
for (i in gse) {
  sample_id <- get(i)@meta.data[["sample_id"]] %>% unique()
  all_sample_id <- c(all_sample_id, sample_id)
}
length(all_sample_id)
rm(i, sample_id)


# 合并质控前的数据集，便于后续绘制质控前后的细胞数量变化总图
seurat_list <- map(gse, get)
seurat_combined_before_qc <- merge(
  seurat_list[[1]],
  seurat_list[-1]
)
seurat_combined_before_qc




# 定义调色板
head(palettes_d_names)
palettes_d_names %>%
  filter(length >= length(all_sample_id)) %>%
  arrange(length) %>%
  print(n = 20)
mycolors <- c(
  paletteer_d("PrettyCols::Lively"),
  paletteer_d("MapPalettes::green_machine")[3:5],
  paletteer_d("ltc::mterese")
)
mycolors_many <- c(
  paletteer_d("PrettyCols::Lively"),
  paletteer_d("trekcolors::lcars_series")
)
pdf(width = 11, height = 8, "figures/mycolors.pdf")
show_col(mycolors)
dev.off()

pdf(width = 11, height = 8, "figures/mycolors_many.pdf")
show_col(mycolors_many)
dev.off()

save(mycolors, mycolors_many, file = "output/color_palette.rdata")



# 建立空表，用于统计质控信息
qc_info <- tibble(
  study_id = NA,
  nCount_cut = NA,
  nFeature_cut = NA,
  mitoRatio_cut = NA,
  decont_cut = NA,
  removed_sample_count = NA,
  cells_count_original = NA,
  cells_count_after_filtering = NA
)


# 质量评价 --------------------------------------------------------------------------------

sce <- gse[1]

# 指定质控阈值
nCount_cut <- c(500, 80000)
nFeature_cut <- c(200, 8000)
mitoRatio_cut <- 0.3
decont_cut <- 0.2 # 环境RNA污染去除阈值


# 运用函数，检查质控指标
sc.qc(sce)


# 过滤细胞
sc.filter(sce, sample.filt = "GSM5599228")


# 重新质量评价
sc.qc(str_c(sce, "_filtered"))


# 评估和去除环境游离RNA污染
sc.decont(str_c(sce, "_filtered"), decont.cut = decont_cut)

# 保存过滤后的Seurat对象
sc.save(str_c(sce, "_low_con"))

# 记录过滤信息
qc_info <- tibble(
  study_id = sce,
  nCount_cut = str_flatten(nCount_cut, collapse = "~"),
  nFeature_cut = str_flatten(nFeature_cut, collapse = "~"),
  mitoRatio_cut,
  decont_cut,
  removed_sample_count = length(unique(get(sce)$sample_id)) -
    length(unique(get(str_c(sce, "_low_con"))$sample_id)),
  cells_count_original = ncol(get(sce)),
  cells_count_after_filtering = ncol(get(str_c(sce, "_low_con")))
) %>%
  bind_rows(qc_info, .)


# 过滤其他数据集--------------------------------------------------------------------------

sce <- gse[2]
nCount_cut <- c(300, 20000)
nFeature_cut <- c(200, 6000)
mitoRatio_cut <- 0.15
decont_cut <- 0.3

sc.qc(sce)
sc.filter(sce, sample.filt = "BT1303")
sc.qc(str_c(sce, "_filtered"))
sc.decont(str_c(sce, "_filtered"), decont.cut = decont_cut)
sc.save(str_c(sce, "_low_con"))

qc_info <- tibble(
  study_id = sce,
  nCount_cut = str_flatten(nCount_cut, collapse = "~"),
  nFeature_cut = str_flatten(nFeature_cut, collapse = "~"),
  mitoRatio_cut,
  decont_cut,
  removed_sample_count = length(unique(get(sce)$sample_id)) -
    length(unique(get(str_c(sce, "_low_con"))$sample_id)),
  cells_count_original = ncol(get(sce)),
  cells_count_after_filtering = ncol(get(str_c(sce, "_low_con")))
) %>%
  bind_rows(qc_info, .)





sce <- gse[3]
nCount_cut <- c(500, 20000)
nFeature_cut <- c(300, 6500)
mitoRatio_cut <- 0.15
decont_cut <- 0.2

sc.qc(sce)
sc.filter(sce)
sc.qc(str_c(sce, "_filtered"))
sc.decont(str_c(sce, "_filtered"), decont.cut = decont_cut, sample.filt = "GSM4675274")
sc.save(str_c(sce, "_low_con"))

qc_info <- tibble(
  study_id = sce,
  nCount_cut = str_flatten(nCount_cut, collapse = "~"),
  nFeature_cut = str_flatten(nFeature_cut, collapse = "~"),
  mitoRatio_cut,
  decont_cut,
  removed_sample_count = length(unique(get(sce)$sample_id)) -
    length(unique(get(str_c(sce, "_low_con"))$sample_id)),
  cells_count_original = ncol(get(sce)),
  cells_count_after_filtering = ncol(get(str_c(sce, "_low_con")))
) %>%
  bind_rows(qc_info, .)





sce <- gse[4]
nCount_cut <- c(500, 80000)
nFeature_cut <- c(300, 8000)
mitoRatio_cut <- 0.25
decont_cut <- 0.2

sc.qc(sce)
sc.filter(sce)
sc.qc(str_c(sce, "_filtered"))
sc.decont(str_c(sce, "_filtered"), decont.cut = decont_cut)
sc.save(str_c(sce, "_low_con"))

qc_info <- tibble(
  study_id = sce,
  nCount_cut = str_flatten(nCount_cut, collapse = "~"),
  nFeature_cut = str_flatten(nFeature_cut, collapse = "~"),
  mitoRatio_cut,
  decont_cut,
  removed_sample_count = length(unique(get(sce)$sample_id)) -
    length(unique(get(str_c(sce, "_low_con"))$sample_id)),
  cells_count_original = ncol(get(sce)),
  cells_count_after_filtering = ncol(get(str_c(sce, "_low_con")))
) %>%
  bind_rows(qc_info, .)



# 保存质控信息
qc_info[-1, ] %>% print(width = Inf)
write_csv(qc_info[-1, ], "output/QC_info.csv")


# 合并Seurat --------------------------------------------------------------------------

seurat_list <- map(gse, function(gse_i) {
  str_c(gse_i, "_low_con") %>% get()
})
seurat_list

seurat_combined <- merge(
  seurat_list[[1]],
  seurat_list[-1]
)
seurat_combined
head(seurat_combined)




# 绘制质控前后的细胞数量变化总图
cells_before_qc <- seurat_combined_before_qc@meta.data %>%
  group_by(sample_id, study_id) %>%
  summarise(count = n()) %>%
  mutate(status = "Original")
cells_after_qc <- seurat_combined@meta.data %>%
  group_by(sample_id, study_id) %>%
  summarise(count = n()) %>%
  mutate(status = "After filtering") %>%
  left_join(
    cells_before_qc[, c("sample_id", "study_id")],
    .,
    by = "sample_id",
    keep = FALSE
  ) %>%
  select(!study_id.y) %>%
  rename(study_id = study_id.x) %>%
  mutate(
    count = if_else(is.na(count), 0, count),
    status = "After filtering"
  )
cells_data <- bind_rows(cells_before_qc, cells_after_qc) %>%
  mutate(
    status = factor(status, levels = c("Original", "After filtering")),
    sample_id = fct_relevel(sample_id, unique(sample_id[order(study_id)]))
  )

ggplot(cells_data, aes(x = sample_id, y = count, fill = status)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_fill_manual(values = mycolors) +
  ggtitle("Number of cells before and after QC") +
  xlab("Samples") +
  ylab("Number of cells") +
  scale_y_continuous(
    labels = label_comma(),
    n.breaks = 15,
    expand = expansion(mult = c(0, 0.05)) # 让Y轴从0开始
  ) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(), # 去除主要和次要垂直网格线
    title = element_text(size = 12),
    legend.title = element_blank(),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
  )
ggsave("figures/QC/number_of_cells_before_and_after_QC.pdf", width = 10.5, height = 5.1)



# 检查样本信息的分布
map(
  c(
    "study_id", "sample_type", "pathology", "tumor_stage",
    "chemotherapy_response", "metastatic"
  ),
  ~ table(seurat_combined[[.x]])
)
# 更改列顺序和重编码变量
seurat_combined@meta.data <- seurat_combined@meta.data %>%
  relocate(
    study_id, sample_id, sample_type,
    tumor_stage, metastatic,
    chemotherapy_response, age,
    .after = nFeature_RNA
  ) %>%
  select(!c(patient_id, pathology)) %>%
  mutate(
    sample_type = if_else(
      str_detect(sample_type, "Normal"), "Normal", "OC"
    )
  )
tail(seurat_combined)

# 检查并保存样本信息
sample_info_1 <- seurat_combined@meta.data %>%
  select(study_id:age) %>%
  distinct() %>%
  arrange(study_id) %>%
  tibble()
sample_info_1


# 统计各样本的细胞数、UMI数、基因数
sample_info_2 <- seurat_combined@meta.data %>%
  group_by(sample_id) %>%
  summarise(
    ncells = n(),
    nUMIs = sum(nCount_RNA),
    ngenes = sum(nFeature_RNA)
  ) %>%
  tibble()
sample_info_2

# 合并上面两个表格
sample_info <- left_join(
  sample_info_1,
  sample_info_2,
  by = "sample_id"
)
sample_info
write_csv(sample_info, file = "output/sample_info.csv")

# 检查良、恶性样本的细胞数
table(seurat_combined$sample_type)



# 合并layers
{
  seurat_combined <- JoinLayers(seurat_combined)
  beep()
}
seurat_combined


qsave(seurat_combined, file = "output/seurat_obj/seurat_combined.qs")







# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/1_QC.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
