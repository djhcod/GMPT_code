# 加载包及数据 ---------------------------------------------------------------------

{
  rm(list = ls())
  cat("\014")

  library(tidyverse)
  library(qs)
  library(survminer)
  library(survival)
  library(cowplot)
  library(patchwork)
  library(randomForestSRC)
  library(interp) # 调参可视化数据准备（contour plot）

  # 载入颜色集
  load("output/color_palette.rdata")

  fig_dir <- "figures/modeling/"
  out_dir <- "output/modeling/"
}

# 载入共同差异基因
markers <- qread("output/identify_candidate_genes/markers.qs")


# 预后分析数据准备 ------------------------------------------------------------

# 载入TCGA表达量矩阵
exp <- qread("output/tcga/TCGA_raw_count_matrix.qs")
exp[1:4, 1:4]

# 载入TCGA临床信息
clinic <- qread("output/tcga/TCGA_clinical_dat_with_cell_propotations.qs")
glimpse(clinic)



# 提取共同差异基因的表达量信息
exp <- exp[rownames(exp) %in% markers, ]
nrow(exp)

# 整理表达量矩阵
exp <- t(exp) %>%
  as.data.frame() %>%
  rownames_to_column("patient_id") %>%
  as_tibble() %>%
  mutate(patient_id = str_remove(patient_id, "-01.$"))
exp
colnames(exp)

# log2转换
exp[, 2:ncol(exp)] <- apply(exp[, 2:ncol(exp)], 2, function(x) {(log2(x + 1))})
summary(exp[, 2:5])

# 将基因名中的“-”替换成“_”，防止在后续构建公式时引用错误
colnames(exp) <- colnames(exp) %>%
  str_replace_all("-", "_")



# 合并TCGA表达量矩阵和临床信息
tcga_dat <- left_join(
  clinic,
  exp,
  by = "patient_id"
)
tcga_dat




# 拆分训练集和验证集

set.seed(123)

# 创建训练集和验证集的索引
train_indices <- sample(seq_len(nrow(tcga_dat)), size = 0.7 * nrow(tcga_dat))

# 根据索引分割数据集
tcga_dat$split <- "test"
tcga_dat$split[train_indices] <- "train"

table(tcga_dat$split)









# 获取所有共同marker基因的截断值
cut_list <- map(
  colnames(exp[, -1]),
  function(gene) {
    cut_off <- surv_cutpoint(
      tcga_dat,
      time = "surv_time",
      event = "vital_status",
      variables = gene,
      minprop = 0.2
    )
    cut_off <- as.numeric(cut_off[["cutpoint"]][1, 1])

    cut_tab <- data.frame(
      gene_name = gene,
      cut_off = cut_off
    )

    # 根据截断值将各基因的表达量分为低和高
    group <<- if_else(pull(tcga_dat, gene) < cut_off, 0, 1) %>%
      factor(levels = c(0, 1), labels = c("Low", "High"))

    # 将高低分组信息返回到TCGA数据中
    tcga_dat <<- tcga_dat %>%
      mutate(
        !!str_c(gene, "_group") := if_else(
          pull(tcga_dat, gene) < cut_off, "Low", "High"
        ) %>%
          factor(levels = c("Low", "High"))
      )
    return(cut_tab)
  }
)
cutoff_table <- list_rbind(cut_list)
head(cutoff_table)
write_csv(cutoff_table, str_c(out_dir, "gene_cutoff_values.csv"))
summary(tcga_dat$ZEB2_group)

vars <- tcga_dat %>%
  select(ends_with("_group")) %>%
  colnames()
vars






# 随机生存森林分析 -------------------------------------------------------------

# 定义随机森林模型的方程
fun_forest <- vars %>%
  str_c(collapse = " + ") %>%
  str_c("Surv(surv_time, vital_status) ~ ", .) %>%
  as.formula()
fun_forest


rsf_dat <- tcga_dat %>% filter(split == "train")
rsf_dat[, vars] <- lapply(
  rsf_dat[, vars],
  function(x) {as.numeric(x) - 1}
)
summary(rsf_dat$ZEB2_group)


# 通过grid search进行随机森林调参，寻找最优nodesize和mtry参数组合
set.seed(123)
tune_forest <- tune(
  fun_forest,
  data = rsf_dat,
  ntreeTry = 500, # 进行调优的树的数量
  mtryStart = 1, # mtry迭代的起始值
  doBest = TRUE,
  trace = FALSE # 是否显示每次迭代过程
)

# the optimized forest
print(tune_forest$rf)




# 调参结果可视化
# visualize the nodesize/mtry OOB surface
plot.tune <- function(tune_forest, linear = TRUE) {
  x <- tune_forest$results[, "nodesize"] # 获取grid search后的nodesize向量
  y <- tune_forest$results[, "mtry"] # 获取grid search后的mtry向量
  z <- tune_forest$results[, "err"] # 获取grid search后的OBB error向量

  best.nodesize <- x[which.min(z)] # 最小OBB error对应的最优nodesize
  best.mtry <- y[which.min(z)] # 最小OBB error对应的y，即最优mtry

  so <- interp(
    x = x,
    y = y,
    z = z,
    linear = linear
  )

  filled.contour(
    x = so$x,
    y = so$y,
    z = so$z,
    xlim = range(so$x, finite = TRUE) + c(-1, 1),
    ylim = range(so$y, finite = TRUE) + c(-0.1, 0.1),
    color.palette = colorRampPalette(mycolors[3:1]),
    plot.axes = {
      axis(1) # 绘制水平坐标轴
      axis(2) # 绘制垂直坐标轴
      # 绘制grid search中的每个参数点
      points(x, y, pch = 16, cex = .3, col = "white")
      # 标注最优参数点
      points(
        x = best.nodesize,
        y = best.mtry,
        pch = "x",
        cex = 1,
        font = 2
      )
    },
    xlab = "Minimum terminal node size",
    ylab = "Mtry (No. of variables tried at each split)",
    main = paste0(
      "Optimal nodesize: ", best.nodesize,
      "; Optimal mtry: ", best.mtry
    ),
    cex.main = 0.8,
    key.title = title(main = "OOS error", cex.main = 0.8)
  )
}

# plot the surface
plot.tune(tune_forest)
p_tune <- recordPlot() %>% ggdraw()


# 提取最优nodesize和mtry参数
grid.search_result <- as_tibble(tune_forest$results) %>%
  arrange(err)
grid.search_result

best.nodesize <- grid.search_result$nodesize[1]
best.nodesize
best.mtry <- grid.search_result$mtry[1]
best.mtry





# 根据最优mtry和nodesize组合构建随机生存森林模型
set.seed(123)
rfsrc <- rfsrc(
  fun_forest,
  data = rsf_dat,
  ntree = 3000,
  nodesize = best.nodesize,
  mtry = best.mtry,
  block.size = 1,
  bootstrap = "by.root",
  samptype = "swor",
  importance = "permute",
  na.action = "na.impute",
  nimpute = 5
)
rfsrc


# 绘制不同ntrees对应的累积OOB error rate的折线图，确定ntree
rfsrc_data <- tibble(
  ntrees = c(1:length(rfsrc$err.rate)),
  err.rate = rfsrc$err.rate
)
rfsrc_data

p_oob <- ggplot(rfsrc_data, aes(x = ntrees, y = err.rate)) +
  geom_line(lwd = 0.5, color = mycolors[2]) +
  labs(x = "Number of trees", y = "Tree cumulative OOS error rate") +
  theme_bw() +
  theme(
    plot.margin = ggplot2::margin(40, 1, 40, 1),
    aspect.ratio = 1
  )
p_oob






# 根据最优参数组合构建最终的随机生存森林模型
ntree <- 1000

set.seed(123)
rfsrc_final <- rfsrc(
  fun_forest,
  data = rsf_dat,
  ntree = ntree,
  nodesize = best.nodesize,
  mtry = best.mtry,
  block.size = 1,
  bootstrap = "by.root",
  samptype = "swor",
  importance = "permute",
  na.action = "na.impute",
  nimpute = 5
)
rfsrc_final





# 绘制变量重要性图
importance <- tibble(
  vars = names(rfsrc_final$importance) %>% str_remove("_group$"),
  vimp = rfsrc_final$importance
)
importance

p_vimp <- ggplot(
  importance,
  aes(x = vimp, y = reorder(vars, vimp))
) +
  geom_bar(
    aes(fill = vimp),
    stat = "identity",
    width = 0.9
  ) +
  scale_fill_gradientn(colors = mycolors[1:3]) +
  xlab("Variable Importance (VIMP)") +
  ylab("Variables") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_blank()
  )
p_vimp

# 合并随机森林调参图和变量重要性图
plot_grid(p_tune, p_oob, p_vimp, labels = "AUTO", nrow = 1)
ggsave(width = 17.3, height = 5.1, str_c(fig_dir, "rfsrc_gene_filt.pdf"))


# 最终的模型基因
aim_gene <- importance %>%
  slice_max(vimp, n = 10) %>%
  pull(vars)
aim_gene

write_csv(
  importance %>%
    mutate(vars = str_remove(vars, "_group$")),
  str_c(out_dir, "vimp_table.csv")
  )



# 只保留TCGA数据中的模型基因
tcga_filt <- tcga_dat[, colnames(tcga_dat) %in% c(colnames(clinic), aim_gene, "split")]

















# 建模 ----------------------------------------------------------------------

fun_model <- aim_gene %>%
  str_c(collapse = " + ") %>%
  str_c("Surv(surv_time, vital_status) ~ ", .) %>%
  as.formula()
fun_model





# 通过grid search进行随机森林调参，寻找最优nodesize和mtry参数组合
set.seed(123)
tune_forest <- tune(
  fun_model,
  data = rsf_dat,
  ntreeTry = 500, # 进行调优的树的数量
  mtryStart = 1, # mtry迭代的起始值
  doBest = TRUE,
  trace = FALSE # 是否显示每次迭代过程
)

# the optimized forest
print(tune_forest$rf)




# 调参结果可视化
# visualize the nodesize/mtry OOB surface
plot.tune <- function(tune_forest, linear = TRUE) {
  x <- tune_forest$results[, "nodesize"] # 获取grid search后的nodesize向量
  y <- tune_forest$results[, "mtry"] # 获取grid search后的mtry向量
  z <- tune_forest$results[, "err"] # 获取grid search后的OBB error向量

  best.nodesize <- x[which.min(z)] # 最小OBB error对应的最优nodesize
  best.mtry <- y[which.min(z)] # 最小OBB error对应的y，即最优mtry

  so <- interp(
    x = x,
    y = y,
    z = z,
    linear = linear
  )

  filled.contour(
    x = so$x,
    y = so$y,
    z = so$z,
    xlim = range(so$x, finite = TRUE) + c(-1, 1),
    ylim = range(so$y, finite = TRUE) + c(-0.1, 0.1),
    color.palette = colorRampPalette(mycolors[3:1]),
    plot.axes = {
      axis(1) # 绘制水平坐标轴
      axis(2) # 绘制垂直坐标轴
      # 绘制grid search中的每个参数点
      points(x, y, pch = 16, cex = .3, col = "white")
      # 标注最优参数点
      points(
        x = best.nodesize,
        y = best.mtry,
        pch = "x",
        cex = 1,
        font = 2
      )
    },
    xlab = "Minimum terminal node size",
    ylab = "Mtry (No. of variables tried at each split)",
    main = paste0(
      "Optimal nodesize: ", best.nodesize,
      "; Optimal mtry: ", best.mtry
    ),
    cex.main = 0.8,
    key.title = title(main = "OOS error", cex.main = 0.8)
  )
}

# plot the surface
plot.tune(tune_forest)
p_tune <- recordPlot() %>% ggdraw()


# 提取最优nodesize和mtry参数
grid.search_result <- as_tibble(tune_forest$results) %>%
  arrange(err)
grid.search_result

best.nodesize <- grid.search_result$nodesize[1]
best.nodesize
best.mtry <- grid.search_result$mtry[1]
best.mtry





# 根据最优mtry和nodesize组合构建随机生存森林模型
set.seed(123)
rfsrc <- rfsrc(
  fun_model,
  data = rsf_dat,
  ntree = 2000,
  nodesize = best.nodesize,
  mtry = best.mtry,
  block.size = 1,
  bootstrap = "by.root",
  samptype = "swor",
  importance = "permute",
  na.action = "na.impute",
  nimpute = 5
)
rfsrc


# 绘制不同ntrees对应的累积OOB error rate的折线图，确定ntree
rfsrc_data <- tibble(
  ntrees = c(1:length(rfsrc$err.rate)),
  err.rate = rfsrc$err.rate
)
rfsrc_data

p_oob <- ggplot(rfsrc_data, aes(x = ntrees, y = err.rate)) +
  geom_line(lwd = 0.5, color = mycolors[2]) +
  labs(x = "Number of trees", y = "Tree cumulative OOS error rate") +
  theme_bw() +
  theme(
    plot.margin = ggplot2::margin(40, 1, 40, 1),
    aspect.ratio = 1
  )
p_oob






# 根据最优参数组合构建最终的随机生存森林模型
ntree <- 1000

set.seed(123)
model <- rfsrc(
  fun_model,
  data = rsf_dat,
  ntree = ntree,
  nodesize = best.nodesize,
  mtry = best.mtry,
  block.size = 1,
  bootstrap = "by.root",
  samptype = "swor",
  importance = "permute",
  na.action = "na.impute",
  nimpute = 5
)
model





# 绘制变量重要性图
importance <- tibble(
  vars = names(model$importance) %>% str_remove("_group$"),
  vimp = model$importance
)
importance

p_vimp <- ggplot(
  importance,
  aes(x = vimp, y = reorder(vars, vimp))
) +
  geom_bar(
    aes(fill = vimp),
    stat = "identity",
    width = 0.9
  ) +
  scale_fill_gradientn(colors = mycolors[1:3]) +
  xlab("Variable Importance (VIMP)") +
  ylab("Variables") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
p_vimp

# 合并随机森林调参图和变量重要性图
plot_grid(p_tune, p_oob, p_vimp, labels = "AUTO", nrow = 1)
ggsave(width = 17.3, height = 5.1, str_c(fig_dir, "final_rfsrc_model.pdf"))

save(model, aim_gene, tcga_filt, file = str_c(out_dir, "gene_model.rdata"))


# 导出模型数据用于网页构建
saveRDS(model, "output/model_dat_for_web/gene_model.rds")








# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/9.2_Gene_model.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
