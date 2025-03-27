# 加载包 ---------------------------------------------------------------------

{
  rm(list = ls())
  cat("\014")

  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(tidyverse)
  library(qs)
  library(gtsummary)
  library(flextable) # 输出.docx格式文档
}


# 下载GDC项目 -----------------------------------------------------------------

if (FALSE) {
  # 查询所有的GDC项目ID
  all_projects <- getGDCprojects()
  all_projects$project_id

  # 查询卵巢癌的项目信息
  all_projects %>%
    filter(project_id %in% c("TCGA-OV", "CPTAC-2"))

  # 查看项目概要
  getProjectSummary("TCGA-OV")

  # 检索转录组数据
  query.exp <- GDCquery(
    project = "TCGA-OV",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )

  # 下载目标数据
  GDCdownload(query = query.exp, directory = "data")


  # 读取下载的数据并转换成summarizedExperiment对象
  exp_obj <- GDCprepare(
    query = query.exp,
    directory = "data" # 下载的表达量文件所在的目录
  )

  qsave(exp_obj, "output/tcga/summarizedExperiment_obj_TCGA_OV.qs")
}

exp_obj <- qread("output/tcga/summarizedExperiment_obj_TCGA_OV.qs")






# 整理自带的临床数据（Clinical indexed data） ----------------------------------

ncol(exp_obj) # 查看包含的样本数量
exp_obj@metadata[["data_release"]] # 查看数据发布日期
# 查看包含的临床信息
exp_obj@colData@listData %>%
  as_tibble() %>%
  glimpse()

# 记录患者数
tcga_sum <- tibble(
  "stage" = "initial",
  "number_of_patients" = ncol(exp_obj)
)
tcga_sum


# 查看目标变量的分布情况
map(
  exp_obj@colData@listData %>%
    as_tibble() %>%
    select(
      sample_type, primary_diagnosis, year_of_diagnosis, race,
      figo_stage, vital_status
    ),
  ~ factor(.x) %>% summary()
)
map(
  exp_obj@colData@listData %>%
    as_tibble() %>%
    select(age_at_index, days_to_last_follow_up, days_to_death),
  ~ summary(.x)
)


## 整理基础临床数据--------------------------------------------------------
coldata <- exp_obj@colData@listData %>%
  as_tibble() %>%
  filter(
    sample_type == "Primary Tumor",
    primary_diagnosis == "Serous cystadenocarcinoma, NOS"
  ) %>%
  select(
    patient_id = patient,
    sample_id = sample,
    year_of_diagnosis,
    age_at_index,
    race,
    figo_stage,
    vital_status,
    days_to_last_follow_up,
    days_to_death
  )
# 记录患者数
tcga_sum <- tibble(
  "stage" = "Restricted to Primary Tumor and Serous cystadenocarcinoma, NOS",
  "number_of_patients" = nrow(coldata)
) %>%
  bind_rows(tcga_sum, .)
tcga_sum

glimpse(coldata)
map(
  coldata %>%
    select(year_of_diagnosis, race, figo_stage, vital_status),
  ~ factor(.x) %>% summary()
)
map(
  coldata %>%
    select(age_at_index, days_to_last_follow_up, days_to_death),
  ~ summary(.x)
)



coldata <- coldata %>%
  mutate(
    year_of_diagnosis = case_when(
      year_of_diagnosis < 2002 ~ "Before 2002",
      year_of_diagnosis < 2006 ~ "2002~2005",
      year_of_diagnosis >= 2006 ~ "After 2005"
    ) %>%
      factor(levels = c("Before 2002", "2002~2005", "After 2005")),
    age_at_diagnosis = age_at_index,
    race = case_when(
      race == "white" ~ "White",
      race == "black or african american" ~ "Black",
      race == "asian" ~ "Asian",
      .default = "Others"
    ) %>%
      factor(levels = c("White", "Black", "Others")),
    figo_stage = case_when(
      str_detect(figo_stage, "^Stage IV") ~ 4,
      str_detect(figo_stage, "^Stage III") ~ 3,
      str_detect(figo_stage, "Stage II[A-C]?") ~ 2,
      str_detect(figo_stage, "Stage I[A-C]?") ~ 1,
    ) %>%
      factor(),
    vital_status = if_else(vital_status == "Dead", 1, 0) %>%
      factor(levels = c(0, 1), labels = c("Alive", "Dead")),
    surv_time = if_else(
      vital_status == "Dead",
      round(days_to_death / 30, digits = 0),
      round(days_to_last_follow_up / 30, digits = 0)
    ),
    .keep = "unused"
  )
summary(coldata)

qsave(coldata, "output/tcga/TCGA_colData.qs")






## 处理治疗信息 ----------------------------------------------------------------

coldata <- qread("output/tcga/TCGA_colData.qs")

treat_dat <- exp_obj@colData@listData[["treatments"]] %>%
  list_rbind() %>%
  as_tibble() %>%
  mutate(
    patient_id = str_remove(submitter_id, "_treatment.*"),
    .before = 1
  ) %>%
  filter(patient_id %in% coldata$patient_id)
glimpse(treat_dat)
# 查看目标变量的分布情况
map(
  treat_dat %>%
    select(treatment_type, treatment_or_therapy),
  ~ factor(.x) %>% summary()
)

treat_dat <- treat_dat %>%
  select(patient_id, submitter_id, treatment_type, treatment_or_therapy)

# 提取包含多个治疗信息的患者ID
ids <- treat_dat %>%
  group_by(patient_id) %>%
  summarise(count = n()) %>%
  filter(count > 2)
ids
# 查看这些患者的治疗信息
treat_dat %>%
  filter(patient_id %in% ids$patient_id) %>%
  arrange(submitter_id) %>%
  print(n = Inf)


# 整理放疗信息
radiotherapy <- treat_dat %>%
  filter(treatment_type == "Radiation Therapy, NOS") %>%
  select(patient_id, radiotherapy = treatment_or_therapy) %>%
  distinct()
radiotherapy
# 检查每个患者是否有重复的放疗信息
duplicated(radiotherapy$patient_id) %>% any()


# 添加到coldata中
coldata <- coldata %>%
  left_join(
    x = coldata,
    y = radiotherapy,
    by = "patient_id"
  ) %>%
  mutate(
    radiotherapy = case_when(
      radiotherapy == "no" ~ 0,
      radiotherapy == "yes" ~ 1,
      radiotherapy == "not reported" ~ NA
    ) %>%
      factor(levels = c(0, 1), labels = c("No", "Yes"))
  )
coldata %>% print(width = Inf)
summary(coldata$radiotherapy)

qsave(coldata, "output/tcga/TCGA_colData.qs")











# 下载和处理额外临床信息 -------------------------------------------------------------

if (FALSE) {
  # 检索额外临床信息
  query.clinical <- GDCquery(
    project = "TCGA-OV",
    data.category = "Clinical",
    data.format = "bcr xml"
  )
  # 下载额外临床数据
  GDCdownload(query.clinical, directory = "data")

  # 解析患者基本信息
  patient_dat <- GDCprepare_clinic(
    query.clinical,
    clinical.info = "patient",
    directory = "data"
  ) %>%
    as_tibble()
  qsave(patient_dat, "output/tcga/summarizedExperiment_patient_dat.qs")

  # 解析药物治疗信息
  drug_dat <- GDCprepare_clinic(
    query.clinical,
    clinical.info = "drug",
    directory = "data"
  ) %>%
    as_tibble()
  qsave(drug_dat, "output/tcga/summarizedExperiment_drug_dat.qs")

  # 解析复发信息
  recurrence_dat <- GDCprepare_clinic(
    query.clinical,
    clinical.info = "new_tumor_event",
    directory = "data"
  ) %>%
    as_tibble()
  qsave(recurrence_dat, "output/tcga/summarizedExperiment_recurrence_dat.qs")
}



## 整理患者基本信息-----------------------------------------------------------

patient_dat <- qread("output/tcga/summarizedExperiment_patient_dat.qs")
coldata <- qread("output/tcga/TCGA_colData.qs")

# 仅保留coldata中存在的患者
patient_dat <- patient_dat %>%
  filter(bcr_patient_barcode %in% coldata$patient_id)
nrow(patient_dat) == nrow(coldata)

# 查看目标变量的分布情况
patient_dat %>% glimpse()
map(
  patient_dat %>%
    select(
      history_of_neoadjuvant_treatment,
      neoplasm_histologic_grade,
      tumor_residual_disease,
      anatomic_neoplasm_subdivision
    ),
  ~ factor(.x) %>% summary()
)

patient_dat <- patient_dat %>%
  mutate(
    patient_id = bcr_patient_barcode,
    grade = case_when(
      str_detect(neoplasm_histologic_grade, "G3|4") ~ "High grade",
      str_detect(neoplasm_histologic_grade, "G1|2") ~ "Low grade",
      .default = NA
    ) %>%
      factor(levels = c("High grade", "Low grade")),
    residual_disease = factor(
      tumor_residual_disease,
      levels = c("No Macroscopic disease", "1-10 mm", "11-20 mm", ">20 mm")
    ),
    anatomic_subdivision = case_when(
      anatomic_neoplasm_subdivision == "Bilateral" ~ "Bilateral",
      str_detect(anatomic_neoplasm_subdivision, "(Left)|(Right)") ~ "Unilateral",
      is_empty(anatomic_neoplasm_subdivision) ~ NA
    ) %>%
      factor(levels = c("Bilateral", "Unilateral")),
    .keep = "none"
  )
summary(patient_dat)



# 根据“patient_id"匹配“coldata”和"patient_dat"
combined_dat <- left_join(
  x = coldata,
  y = patient_dat,
  by = "patient_id"
)
glimpse(combined_dat)

qsave(combined_dat, "output/tcga/TCGA_combined_clinical_dat.qs")






## 整理药物治疗信息--------------------------------------------------------

drug_dat <- qread("output/tcga/summarizedExperiment_drug_dat.qs")
combined_dat <- qread("output/tcga/TCGA_combined_clinical_dat.qs")

# 仅保留coldata中存在的患者
drug_dat <- drug_dat %>%
  filter(bcr_patient_barcode %in% coldata$patient_id)

unique(drug_dat$bcr_patient_barcode) %>% length()
nrow(coldata)
drug_dat %>% glimpse()

# 查看任意一个患者的药物治疗信息
drug_dat %>%
  filter(bcr_patient_barcode == "TCGA-04-1332") %>%
  print(width = Inf)

# 查看目标变量的分布情况
map(
  drug_dat %>%
    select(
      total_dose_units, therapy_types, drug_name, regimen_indication,
      route_of_administrations
    ),
  ~ factor(.x) %>% summary()
)




### 整理辅助化疗信息----------------------------------
adjuvant_chem <- drug_dat %>%
  filter(
    therapy_types == "Chemotherapy",
    regimen_indication == "ADJUVANT"
  ) %>%
  select(
    patient_id = bcr_patient_barcode,
    number_cycles,
    drug_name
  ) %>%
  mutate(drug_name = str_to_title(drug_name)) %>% # 所有药名首字母大写
  distinct()
adjuvant_chem

# 统计药物的种类和数量
adjuvant_chem %>%
  group_by(drug_name) %>%
  summarise(count = n()) %>%
  arrange(drug_name) %>%
  select(drug_name, count) %>%
  print(n = Inf)

# 合并同义药物名称
adjuvant_chem <- adjuvant_chem %>%
  mutate(
    drug_name_convert = case_when(
      drug_name == "" ~ NA,
      drug_name == "Taxol/Carboplatin" ~ "Paclitaxel + Carboplatin",
      drug_name == "Topotecan/Carboplatin" ~ "Topotecan + Carboplatin",
      str_detect(
        drug_name,
        "Carb.*"
      ) ~ "Carboplatin",
      str_detect(
        drug_name,
        "(Cisplatin.*)|(Ciplastin)"
      ) ~ "Cisplatin",
      str_detect(
        drug_name,
        "Doc|xetaxel"
      ) ~ "Docetaxel",
      str_detect(
        drug_name,
        "(Doxil)|(Doxorubicin.*)"
      ) ~ "Doxorubicin",
      str_detect(
        drug_name,
        "(Gemcita|ibine.*)|(Gemicitabine)|(Gemz.*)"
      ) ~ "Gemcitabine",
      str_detect(
        drug_name,
        "(Taxol.*)|(Paclitaxel.*)|(Pacliatxel)|(Paciltaxe|al)|(Taxotere.*)"
      ) ~ "Paclitaxel",
      .default = drug_name
    )
  )
adjuvant_chem

adjuvant_chem %>%
  group_by(drug_name_convert) %>%
  summarise(count = n()) %>%
  arrange(drug_name_convert) %>%
  select(drug_name_convert, count) %>%
  print(n = Inf)




# 统计药物方案的种类和数量
adjuvant_chem %>%
  group_by(patient_id) %>%
  arrange(drug_name) %>%
  summarise(drug_regimen = str_flatten(drug_name_convert, collapse = " + ")) %>%
  group_by(drug_regimen) %>%
  summarise(count = n()) %>%
  arrange(-count) %>%
  print(n = Inf)


# 整理得到每个患者的化疗方案及周期
adjuvant_chem_tab <- map(
  unique(adjuvant_chem$patient_id),
  function(id) {
    case <- adjuvant_chem[adjuvant_chem$patient_id == id, ]
    drug_regimen <- str_flatten(case$drug_name_convert, collapse = " + ")
    cycles <- max(case$number_cycles)
    adjuvant_chem_drug <- case_when(
      length(unique(case$drug_name_convert)) > 2 ~ "More than 2 drugs",
      str_detect(drug_regimen, "Paclitaxel") &
        str_detect(drug_regimen, "Carboplatin") ~ "Paclitaxel + Carboplatin",
      str_detect(drug_regimen, "Paclitaxel") &
        str_detect(drug_regimen, "Cisplatin") ~ "Paclitaxel + Cisplatin",
      .default = "Others"
    ) %>%
      factor(
        levels = c(
          "Paclitaxel + Carboplatin",
          "Paclitaxel + Cisplatin",
          "More than 2 drugs",
          "Others"
        )
      )
    tibble(patient_id = id, adjuvant_chem_drug, cycles)
  }
) %>%
  list_rbind()
adjuvant_chem_tab
summary(adjuvant_chem_tab)

# 合并到临床数据中
combined_dat <- adjuvant_chem_tab %>%
  mutate(adjuvant_chem = 1) %>%
  left_join(
    x = combined_dat,
    y = .,
    by = "patient_id"
  ) %>%
  mutate(
    adjuvant_chem = if_else(is.na(adjuvant_chem), 0, 1) %>%
      factor(levels = c(1, 0), labels = c("Yes", "No"))
  )
summary(combined_dat)



# 重编码化疗周期
combined_dat <- combined_dat %>%
  mutate(
    cycles = case_when(
      cycles < 6 ~ "<6",
      cycles == 6 ~ "6",
      cycles > 6 ~ ">6"
    ) %>%
      factor(levels = c("6", "<6", ">6")),
    .after = cycles
  )
summary(combined_dat$cycles)

combined_dat <- combined_dat %>%
  mutate(
    cycles_combined = if_else(adjuvant_chem == "No", "0", cycles) %>%
      factor(levels = c("0", levels(combined_dat$cycles))),
    .after = cycles
  )
summary(combined_dat$cycles_combined)


# 添加记录化疗综合信息的变量
combined_dat <- combined_dat %>%
  mutate(
    chemo_combined = if_else(
      adjuvant_chem == "No",
      "No chemotherapy",
      adjuvant_chem_drug
    ) %>%
      factor(),
    .after = adjuvant_chem
  )
summary(combined_dat$chemo_combined)

combined_dat <- combined_dat %>%
  mutate(
    chemo_combined = factor(
      chemo_combined,
      levels = c(
        "Paclitaxel + Carboplatin",
        "No chemotherapy",
        "Paclitaxel + Cisplatin",
        "More than 2 drugs",
        "Others"
      )
    )
  )
summary(combined_dat$chemo_combined)









### 整理激素治疗信息----------------------------------------------------
hormone <- drug_dat %>%
  filter(therapy_types == "Hormone Therapy") %>%
  select(
    patient_id = bcr_patient_barcode,
    drug_name
  ) %>%
  mutate(drug_name = str_to_title(drug_name)) %>% # 所有药名首字母大写
  distinct()
hormone


# 统计药物的种类和数量
hormone %>%
  group_by(drug_name) %>%
  summarise(count = n()) %>%
  arrange(drug_name) %>%
  select(drug_name, count) %>%
  print(n = Inf)

# 筛选接受了激素治疗的患者
patients_hormone <- hormone %>%
  filter(drug_name != "")


# 整理得到每个患者是否接受激素治疗
combined_dat <- combined_dat %>%
  mutate(
    hormone_therapy = if_else(
      patient_id %in% patients_hormone$patient_id, "Yes", "No"
    ) %>%
      factor(levels = c("No", "Yes"))
  )
summary(combined_dat$hormone_therapy)

qsave(combined_dat, "output/tcga/TCGA_combined_clinical_dat.qs")




### 整理靶向治疗信息----------------------------------------------------
combined_dat <- qread("output/tcga/TCGA_combined_clinical_dat.qs")

target_therapy <- drug_dat %>%
  filter(therapy_types == "Targeted Molecular therapy") %>%
  select(
    patient_id = bcr_patient_barcode,
    drug_name
  ) %>%
  mutate(drug_name = str_to_title(drug_name)) %>% # 所有药名首字母大写
  distinct()
target_therapy


# 统计药物的种类和数量
target_therapy %>%
  group_by(drug_name) %>%
  summarise(count = n()) %>%
  arrange(drug_name) %>%
  select(drug_name, count) %>%
  print(n = Inf)

# 筛选接受了Bevacizumab治疗的患者
patients_target <- drug_dat %>%
  filter(
    str_detect(drug_name, "(B|b)evacizum.*|(A|a)vastin.*")
  ) %>%
  select(patient_id = bcr_patient_barcode) %>%
  distinct()
patients_target %>% print(n = Inf)


# 整理得到每个患者是否接受Bevacizumab治疗
combined_dat <- combined_dat %>%
  mutate(
    bevacizumab = if_else(
      patient_id %in% patients_target$patient_id, "Yes", "No"
    ) %>%
      factor(levels = c("No", "Yes"))
  )
summary(combined_dat$bevacizumab)

qsave(combined_dat, "output/tcga/TCGA_combined_clinical_dat.qs")






## 整理复发信息 ------------------------------------------------------------------

recurrence_dat <- qread("output/tcga/summarizedExperiment_recurrence_dat.qs")
combined_dat <- qread("output/tcga/TCGA_combined_clinical_dat.qs")

# 仅保留coldata中存在的患者
recurrence_dat <- recurrence_dat %>%
  filter(bcr_patient_barcode %in% coldata$patient_id)
glimpse(recurrence_dat)


# 查看目标变量的分布情况
recurrence_dat %>%
  select(
    days_to_new_tumor_event_after_initial_treatment,
    new_neoplasm_event_type
  ) %>%
  mutate(
    new_neoplasm_event_type = factor(new_neoplasm_event_type)
  ) %>%
  summary()

# 整理得到复发信息
recurrence_dat <- recurrence_dat %>%
  filter(new_neoplasm_event_type == "Recurrence") %>%
  select(
    patient_id = bcr_patient_barcode,
    recurrence_time = days_to_new_tumor_event_after_initial_treatment
  ) %>%
  mutate(recurrence = 1)
recurrence_dat
summary(recurrence_dat$recurrence_time)

# 检查是否有重复的患者id
duplicated(recurrence_dat$patient_id) %>% any()


combined_dat <- combined_dat %>%
  left_join(
    y = recurrence_dat,
    by = "patient_id"
  ) %>%
  mutate(
    recurrence = if_else(is.na(recurrence), 0, recurrence) %>%
      factor(levels = c(0, 1), labels = c("Recurrence_free", "Recurrence")),
    recurrence_time = ceiling(recurrence_time / 30)
  )
summary(combined_dat)

qsave(combined_dat, "output/tcga/TCGA_combined_clinical_dat.qs")





# 整理最终的临床数据 --------------------------------------------------------------

combined_dat <- qread("output/tcga/TCGA_combined_clinical_dat.qs")
colnames(combined_dat)

combined_dat <- combined_dat %>%
  relocate(
    patient_id:year_of_diagnosis,
    age_at_diagnosis,
    race,
    figo_stage,
    grade,
    anatomic_subdivision,
    residual_disease,
    radiotherapy,
    adjuvant_chem,
    adjuvant_chem_drug,
    chemo_combined,
    cycles,
    cycles_combined,
    hormone_therapy,
    bevacizumab,
    vital_status,
    surv_time,
    recurrence,
    recurrence_time
  ) %>%
  mutate(
    figo_stage = if_else(
      figo_stage == 1 | figo_stage == 2, "I/II", "III/IV"
    ) %>%
      factor(levels = c("III/IV", "I/II"))
  ) %>%
  filter(is.na(surv_time) == FALSE) # 删除生存时间有缺失值的个案
summary(combined_dat)

# 记录患者数
tcga_sum <- tibble(
  "stage" = "Exclude unknown survival time",
  "number_of_patients" = nrow(combined_dat)
) %>%
  bind_rows(tcga_sum, .)
tcga_sum
qsave(tcga_sum, "output/tcga/TCGA_clinical_information_screening_process.qs")

qsave(combined_dat, "output/tcga/TCGA_combined_clinical_dat.qs")










# 绘制临床信息三线表------------------------------------------------------------

combined_dat <- qread("output/tcga/TCGA_combined_clinical_dat.qs")

tab_base <- combined_dat %>%
  select(
    -c(
      patient_id, sample_id, adjuvant_chem_drug,
      cycles, recurrence, recurrence_time
    )
  ) %>%
  tbl_summary(
    label = list( # 修改变量名
      year_of_diagnosis ~ "Year of diagnosis",
      age_at_diagnosis ~ "Age at diagnosis (years)",
      race ~ "Race",
      figo_stage ~ "FIGO stage",
      grade ~ "Grade",
      anatomic_subdivision ~ "Anatomic subdivision",
      residual_disease ~ "Residual disease",
      radiotherapy ~ "Adjuvant radiotherapy",
      adjuvant_chem ~ "Adjuvant chemotherapy",
      chemo_combined ~ "Chemotherapy drugs",
      cycles_combined ~ "Total chemotherapy cycles",
      hormone_therapy ~ "Hormone therapy",
      bevacizumab ~ "Bevacizumab therapy",
      vital_status ~ "Vital status",
      surv_time ~ "Survival time (months)"
    ),
    digits = list(
      all_continuous() ~ 1, # 所有连续变量保留1位小数
      all_categorical() ~ c(0, 1) # 分类变量频数保留整数，百分比保留1位小数
    ),
    type = list(all_dichotomous() ~ "categorical") # 让二分类变量的每个水平都展示
  ) %>%
  bold_labels() # 加粗标签文字
tab_base

# 将基线特征表输出成word
as_flex_table(tab_base) %>%
  save_as_docx(
    tab_base,
    align = "left",
    path = "output/tcga/TCGA_baseline_table.docx"
  )













# 整理表达矩阵 ------------------------------------------------------------------

if (FALSE) {
  exp <- exp_obj@assays@data@listData[["tpm_unstrand"]] %>% as.data.frame()
  exp[1:5, 1:5]
  dim(exp)

  # 添加基因注释
  exp_tab <- exp %>%
    mutate(
      gene_name = exp_obj@rowRanges@elementMetadata@listData[["gene_name"]],
      gene_type = exp_obj@rowRanges@elementMetadata@listData[["gene_type"]],
      .before = 1
    ) %>%
    filter(gene_type == "protein_coding") %>%
    select(-gene_type)
  exp_tab[1:5, 1:5]
  dim(exp_tab)


  # 相同基因取表达量的平均值
  exp_tab_filt <- aggregate(. ~ gene_name, mean, data = exp_tab)
  rownames(exp_tab_filt) <- exp_tab_filt$gene_name
  exp_tab_filt$gene_name <- NULL
  dim(exp_tab_filt)


  # 仅保留在一半以上样本里表达的基因
  exp_tab_filt <- exp_tab_filt %>%
    filter(rowSums(exp_tab_filt > 0) > (ncol(exp_tab_filt) / 2))
  dim(exp_tab_filt)




  # 添加样本名
  colnames(exp_tab_filt) <- exp_obj@colData@rownames

  # 简化样本名称
  combined_dat$sample_id %>% head()
  colnames(exp_tab_filt) %>% head()
  colnames(exp_tab_filt) <- str_sub(colnames(exp_tab_filt), 1, 16)
  colnames(exp_tab_filt) %>% head()

  # 只保留临床数据中存在的样本
  exp_tab_filt <- exp_tab_filt[, combined_dat$sample_id]
  dim(exp_tab_filt)
  combined_dat$sample_id %>% length()

  qsave(exp_tab_filt, "output/tcga/TCGA_raw_count_matrix.qs")
}





exp_CorOutliers <- TCGAanalyze_Preprocessing(
  exp_obj,
  filename = "figures/TCGA_PreprocessingOutput.png"
  )

# normalization of genes
exp_Norm <- TCGAanalyze_Normalization(
  tabDF = exp_CorOutliers,
  geneInfo = geneInfoHT
)

# quantile filter of genes
exp_Filt <- TCGAanalyze_Filtering(
  tabDF = exp_Norm,
  method = "quantile",
  qnt.cut = 0.25
)

exp_Filt[1:5, 1:5]

# 制作基因注释文件
anno <- exp_obj@rowRanges@elementMetadata@listData %>%
  as.data.frame() %>%
  select(gene_id, gene_name, gene_type) %>%
  mutate(gene_id = str_remove_all(gene_id, "\\..*$")) %>%
  distinct()
head(anno)

# 注释表达量矩阵
exp_tab <- exp_Filt %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  left_join(
    anno,
    by = "gene_id"
  ) %>%
  filter(gene_type == "protein_coding") %>%
  select(-c(gene_type, gene_id))


# 相同基因取表达量的平均值
exp_tab <- aggregate(. ~ gene_name, mean, data = exp_tab)
rownames(exp_tab) <- exp_tab$gene_name
exp_tab$gene_name <- NULL
dim(exp_tab)


# 简化样本名称
combined_dat$sample_id %>% head()
colnames(exp_tab) %>% head()
colnames(exp_tab) <- str_sub(colnames(exp_tab), 1, 16)
colnames(exp_tab) %>% head()

# 只保留临床数据中存在的样本
exp_tab_filt <- exp_tab[, combined_dat$sample_id]
dim(exp_tab_filt)
combined_dat$sample_id %>% length()

qsave(exp_tab_filt, "output/tcga/TCGA_raw_count_matrix.qs")





















# 保存sessionInfo -----------------------------------------------------------

sink(
  "sessionInfo/7.1_TCGA_download_sessionInfo.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
