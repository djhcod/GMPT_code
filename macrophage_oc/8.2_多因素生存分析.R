# 加载包及数据 ---------------------------------------------------------------------

{
  rm(list = ls())
  cat("\014")

  library(qs)
  library(tidyverse)
  library(Seurat)
  library(cowplot)
  library(survival)
  library(MatchThem)
  library(forestplot)
  library(glmnet)
  library(randomForestSRC)
  library(interp) # 调参可视化数据准备（contour plot）
  library(ggpubr) # 给ggplot图形添加统计指标

  # 载入颜色集
  load("output/color_palette.rdata")
}



# 载入细胞标签顺序
load("output/BayesPrism/cell_prop_and_label_levels.rdata")


# 载入包含了细胞比例信息的TCGA临床数据
tcga <- qread("output/cluster_surve/TCGA_clinical_dat_with_cell_prop_group.qs")
glimpse(tcga)


# 载入多重插补后匹配的数据集
mimids_list <- qread("output/cluster_surve/mimids_list.qs")
names(mimids_list)

# 载入完成了巨噬细胞细分注释的Seurat对象
seurat_macro <- qread("output/macro/seurat_bio_annotated.qs")
levels(seurat_macro)








# 多因素COX回归 ----------------------------------------------------------------

# 在匹配后的数据集中进行Cox Regression with a robust variance estimator
pool_list <- list()
{
  pool_list[["APOE_LP_TAMs_group"]] <- with(
    mimids_list[["APOE_LP_TAMs_group"]],
    coxph(
      Surv(surv_time, vital_status == 1) ~ APOE_LP_TAMs_group +
        year_of_diagnosis + age_at_diagnosis + race + figo_stage + grade +
        anatomic_subdivision + residual_disease + radiotherapy + chemo_combined +
        cycles_combined + hormone_therapy + bevacizumab,
      cluster = mimids_list[["APOE_LP_TAMs_group"]][["models"]][[1]][["subclass"]],
      weights = mimids_list[["APOE_LP_TAMs_group"]][["models"]][[1]][["weights"]]
    )
  ) %>%
    pool()

  pool_list[["HLA_DPB1_Reg_TAMs_group"]] <- with(
    mimids_list[["HLA_DPB1_Reg_TAMs_group"]],
    coxph(
      Surv(surv_time, vital_status == 1) ~ HLA_DPB1_Reg_TAMs_group +
        year_of_diagnosis + age_at_diagnosis + race + figo_stage + grade +
        anatomic_subdivision + residual_disease + radiotherapy + chemo_combined +
        cycles_combined + hormone_therapy + bevacizumab,
      cluster = mimids_list[["HLA_DPB1_Reg_TAMs_group"]][["models"]][[1]][["subclass"]],
      weights = mimids_list[["HLA_DPB1_Reg_TAMs_group"]][["models"]][[1]][["weights"]]
    )
  ) %>%
    pool()

  pool_list[["CCL3L1_Inflam_TAMs_group"]] <- with(
    mimids_list[["CCL3L1_Inflam_TAMs_group"]],
    coxph(
      Surv(surv_time, vital_status == 1) ~ CCL3L1_Inflam_TAMs_group +
        year_of_diagnosis + age_at_diagnosis + race + figo_stage + grade +
        anatomic_subdivision + residual_disease + radiotherapy + chemo_combined +
        cycles_combined + hormone_therapy + bevacizumab,
      cluster = mimids_list[["CCL3L1_Inflam_TAMs_group"]][["models"]][[1]][["subclass"]],
      weights = mimids_list[["CCL3L1_Inflam_TAMs_group"]][["models"]][[1]][["weights"]]
    )
  ) %>%
    pool()

  pool_list[["TOP2A_Prolif_TAMs_group"]] <- with(
    mimids_list[["TOP2A_Prolif_TAMs_group"]],
    coxph(
      Surv(surv_time, vital_status == 1) ~ TOP2A_Prolif_TAMs_group +
        year_of_diagnosis + age_at_diagnosis + race + figo_stage + grade +
        anatomic_subdivision + residual_disease + radiotherapy + chemo_combined +
        cycles_combined + hormone_therapy + bevacizumab,
      cluster = mimids_list[["TOP2A_Prolif_TAMs_group"]][["models"]][[1]][["subclass"]],
      weights = mimids_list[["TOP2A_Prolif_TAMs_group"]][["models"]][[1]][["weights"]]
    )
  ) %>%
    pool()

  pool_list[["MT1H_MT_TAMs_group"]] <- with(
    mimids_list[["MT1H_MT_TAMs_group"]],
    coxph(
      Surv(surv_time, vital_status == 1) ~ MT1H_MT_TAMs_group +
        year_of_diagnosis + age_at_diagnosis + race + figo_stage + grade +
        anatomic_subdivision + residual_disease + radiotherapy + chemo_combined +
        cycles_combined + hormone_therapy + bevacizumab,
      cluster = mimids_list[["MT1H_MT_TAMs_group"]][["models"]][[1]][["subclass"]],
      weights = mimids_list[["MT1H_MT_TAMs_group"]][["models"]][[1]][["weights"]]
    )
  ) %>%
    pool()

  pool_list[["LYVE1_RTM_TAMs_group"]] <- with(
    mimids_list[["LYVE1_RTM_TAMs_group"]],
    coxph(
      Surv(surv_time, vital_status == 1) ~ LYVE1_RTM_TAMs_group +
        year_of_diagnosis + age_at_diagnosis + race + figo_stage + grade +
        anatomic_subdivision + residual_disease + radiotherapy + chemo_combined +
        cycles_combined + hormone_therapy + bevacizumab,
      cluster = mimids_list[["LYVE1_RTM_TAMs_group"]][["models"]][[1]][["subclass"]],
      weights = mimids_list[["LYVE1_RTM_TAMs_group"]][["models"]][[1]][["weights"]]
    )
  ) %>%
    pool()

  pool_list[["TNFAIP6_Inflam_TAMs_group"]] <- with(
    mimids_list[["TNFAIP6_Inflam_TAMs_group"]],
    coxph(
      Surv(surv_time, vital_status == 1) ~ TNFAIP6_Inflam_TAMs_group +
        year_of_diagnosis + age_at_diagnosis + race + figo_stage + grade +
        anatomic_subdivision + residual_disease + radiotherapy + chemo_combined +
        cycles_combined + hormone_therapy + bevacizumab,
      cluster = mimids_list[["TNFAIP6_Inflam_TAMs_group"]][["models"]][[1]][["subclass"]],
      weights = mimids_list[["TNFAIP6_Inflam_TAMs_group"]][["models"]][[1]][["weights"]]
    )
  ) %>%
    pool()

  pool_list[["HSPA6_HSP_TAMs_group"]] <- with(
    mimids_list[["HSPA6_HSP_TAMs_group"]],
    coxph(
      Surv(surv_time, vital_status == 1) ~ HSPA6_HSP_TAMs_group +
        year_of_diagnosis + age_at_diagnosis + race + figo_stage + grade +
        anatomic_subdivision + residual_disease + radiotherapy + chemo_combined +
        cycles_combined + hormone_therapy + bevacizumab,
      cluster = mimids_list[["HSPA6_HSP_TAMs_group"]][["models"]][[1]][["subclass"]],
      weights = mimids_list[["HSPA6_HSP_TAMs_group"]][["models"]][[1]][["weights"]]
    )
  ) %>%
    pool()

  pool_list[["LILRA4_RB_TAMs_group"]] <- with(
    mimids_list[["LILRA4_RB_TAMs_group"]],
    coxph(
      Surv(surv_time, vital_status == 1) ~ LILRA4_RB_TAMs_group +
        year_of_diagnosis + age_at_diagnosis + race + figo_stage + grade +
        anatomic_subdivision + residual_disease + radiotherapy + chemo_combined +
        cycles_combined + hormone_therapy + bevacizumab,
      cluster = mimids_list[["LILRA4_RB_TAMs_group"]][["models"]][[1]][["subclass"]],
      weights = mimids_list[["LILRA4_RB_TAMs_group"]][["models"]][[1]][["weights"]]
    )
  ) %>%
    pool()

  pool_list[["CXCL10_IFN_TAMs_group"]] <- with(
    mimids_list[["CXCL10_IFN_TAMs_group"]],
    coxph(
      Surv(surv_time, vital_status == 1) ~ CXCL10_IFN_TAMs_group +
        year_of_diagnosis + age_at_diagnosis + race + figo_stage + grade +
        anatomic_subdivision + residual_disease + radiotherapy + chemo_combined +
        cycles_combined + hormone_therapy + bevacizumab,
      cluster = mimids_list[["CXCL10_IFN_TAMs_group"]][["models"]][[1]][["subclass"]],
      weights = mimids_list[["CXCL10_IFN_TAMs_group"]][["models"]][[1]][["weights"]]
    )
  ) %>%
    pool()
}



cox_tab <- map(
  names(pool_list),
  function(x) {
    sum_pool_cox <- summary(pool_list[[x]], conf.int = TRUE)[1, ]
    sum_pool_hr <- summary(
      pool_list[[x]],
      conf.int = TRUE,
      exponentiate = TRUE # 将回归系数指数化,得到HR值
    )[1, ]
    variables <- x %>%
      str_remove("_group") %>%
      str_replace_all("_", "-")
    coef <- sprintf("%.2f", sum_pool_cox$estimate)
    hr <- sprintf("%.2f", sum_pool_hr$estimate)
    lci <- sprintf("%.2f", sum_pool_hr$`2.5 %`)
    uci <- sprintf("%.2f", sum_pool_hr$`97.5 %`)
    hr_ci <- str_c(hr, " (", lci, ", ", uci, ")")
    p_value <- case_when(
      sum_pool_hr$p.value < 0.001 ~ "< 0.001*",
      sum_pool_hr$p.value < 0.05 ~ sprintf("%.3f", sum_pool_hr$p.value) %>%
        str_c(., "*"),
      .default = sprintf("%.3f", sum_pool_hr$p.value)
    )
    pool_multicox_result <- tibble(
      Clusters = variables,
      Coefficient = coef,
      HR = hr,
      LCI = lci,
      UCI = uci,
      HR_95_CI = hr_ci,
      P_value = p_value
    )
    return(pool_multicox_result)
  }
) %>%
  list_rbind()
cox_tab
write_csv(cox_tab, "output/cluster_surve/Multivariate_COX_results.csv")




# 绘制森林图
forest_data <- cox_tab
forest_data[, c("HR", "LCI", "UCI")] <- lapply(
  forest_data[, c("HR", "LCI", "UCI")],
  as.numeric
)

pdf(
  width = 6.4,
  height = 4.5,
  "figures/cluster_surve/Multivariate_COX_forestplot.pdf"
)
forestplot(
  forest_data,
  labeltext = c(Clusters, Coefficient, HR_95_CI, P_value), # 选择哪些列以文字展示
  align = "llll", # 设置每列文字的对齐方式
  mean = HR, # 定义HR均值列
  lower = LCI, # 定义HR值的95%CI下限列
  upper = UCI, # 定义HR值的5%CI上限列
  graph.pos = 4, # 设置森林图出现的列
  graphwidth = unit(50, "mm"), # 森林图的宽度
  colgap = unit(5, "mm"), # 设置图形中的列间距
  fn.ci_norm = fpDrawDiamondCI, # 中央HR均值点的形状
  boxsize = 0.2, # 中央HR均值点的大小
  lwd.ci = 1.5, # 设置95%CI线的粗细
  ci.vertices = TRUE, # 添加95%CI线两端的小竖线
  ci.vertices.height = 0.1, # 95%CI线两端的小竖线的高度
  zero = 1, # 设置无效线
  lwd.zero = 1, # 设置无效线的粗细
  col = fpColors(
    box = mycolors[2], # OR均值点的颜色
    lines = mycolors[1], # 95%CI线的颜色
    zero = "black"
  ), # 无效线的颜色
  xlog = T, # 转换为对数坐标轴
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.7), # 表格主体文字的大小
    ticks = gpar(cex = 0.7), # 森林图下方的坐标轴的刻度文字大小
    xlab = gpar(cex = 0.7)
  ),
  new_page = FALSE
) %>% # 森林图下方X轴标签文字的大小
  # 添加表头
  fp_add_header(
    Clusters = "Clusters",
    Coefficient = "Coefficient",
    HR_95_CI = "HR (95%CI)",
    P_value = "P value"
  ) %>%
  # 添加横线，制作成三线表
  fp_add_lines(
    h_1 = gpar(lwd = 2, col = "black"),
    h_2 = gpar(lwd = 1, col = "black"),
    h_12 = gpar(lwd = 2, col = "black")
  ) %>%
  fp_set_zebra_style("grey95") # 添加间隔填充
dev.off()












# LASSO回归 -----------------------------------------------------------------

# 数据准备
x <- tcga %>%
  select(APOE_LP_TAMs_group:CXCL10_IFN_TAMs_group) %>%
  colnames() %>%
  str_c(collapse = " + ") %>%
  str_c("~ ", .) %>%
  as.formula() %>%
  model.matrix(tcga) %>%
  .[, -1]
head(x)

y <- Surv(tcga$surv_time, tcga$vital_status == 1) # 定义响应变量



set.seed(123)

# 拟合LASSO模型
lasso <- glmnet(
  x = x,
  y = y,
  family = "cox",
  alpha = 1,
  standardize = T,
  nlambda = 100
)

# LASSO回归系数交叉验证
lasso.cv <- cv.glmnet(
  x = x,
  y = y,
  family = "cox",
  alpha = 1, # 指定为lasso回归
  nfolds = 10, # 交叉验证次数，一般取10
  type.measure = "deviance"
)


## LASSO回归可视化 --------------------------------------------------------------

# 用ggplot手动绘制LASSO交叉验证图

# 提取结果
cv_data <- tibble(
  lambda = log(lasso.cv$lambda), # lambda值取对数
  cvm = lasso.cv$cvm, # 交叉验证误差均值
  cvup = lasso.cv$cvup, # 交叉验证误差的上界
  cvlo = lasso.cv$cvlo, # 交叉验证误差的下界
  nzero = lasso.cv$nzero # 非零变量的数量
)

# 使用ggplot绘制
plot_lasso.cv <- ggplot(cv_data, aes(x = lambda, y = cvm)) +
  geom_errorbar(
    aes(ymin = cvlo, ymax = cvup),
    width = 0.05,
    color = "gray70"
  ) +
  geom_point(color = mycolors[2], size = 2) +
  geom_vline(
    xintercept = log(lasso.cv$lambda.min),
    color = "gray20",
    lty = 2,
    linewidth = 0.5
  ) +
  scale_x_continuous(
    sec.axis = sec_axis(
      ~.,
      breaks = cv_data$lambda[seq(1, length(cv_data$lambda), by = 5)], # 每隔5个lambda值放一个break
      labels = cv_data$nzero[seq(1, length(cv_data$nzero), by = 5)], # 相应的非零系数数量
      name = "Number of Non-Zero Coefficients"
    )
  ) +
  labs(
    x = "Log (Lambda)",
    y = "Partial Likelihood Deviance"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
plot_lasso.cv




# 用ggplot手动绘制LASSO系数图

# 提取结果
lasso_dat <- data.frame("lambda" = log(lasso$lambda), as.matrix(t(lasso$beta))) %>%
  tibble()
lasso_dat

lasso_dat <- lasso_dat %>%
  pivot_longer(
    cols = 2:ncol(lasso_dat),
    names_to = "variable",
    values_to = "coef"
  ) %>%
  mutate(variable = str_remove(variable, "_groupHigh")) %>%
  mutate(variable = factor(variable, levels = label_levels, labels = labels))
lasso_dat

# 提取非零系数数量
nzero_dat <- tibble(lambda = log(lasso$lambda), nzero = lasso$df)
nzero_dat

# 绘图
plot_lasso <- ggplot(
  lasso_dat,
  aes(x = lambda, y = coef, color = variable)
) +
  geom_line(linewidth = 1, alpha = 0.8) +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(
    sec.axis = sec_axis(
      ~.,
      breaks = nzero_dat$lambda[seq(1, length(nzero_dat$lambda), by = 5)], # 每隔5个lambda值放一个break
      labels = nzero_dat$nzero[seq(1, length(nzero_dat$nzero), by = 5)], # 相应的非零系数数量
      name = "Number of Non-Zero Coefficients"
    )
  ) +
  geom_vline(
    xintercept = log(lasso.cv$lambda.min),
    color = "gray20",
    lty = 2,
    linewidth = 0.5
  ) +
  labs(
    x = "Log (Lambda)",
    y = "Coefficient"
  ) +
  guides(color = guide_legend(title = "Variables")) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
plot_lasso







# 以条形图的形式展示各变量的lasso回归系数
coef_tab <- tibble(
  variables = labels,
  value = coef(lasso.cv, s = "lambda.min") %>% as.matrix() %>% c()
)
coef_tab

plot_coef <- ggplot(
  coef_tab,
  aes(
    x = value,
    y = reorder(variables, value)
  )
) +
  geom_bar(
    stat = "identity",
    show.legend = T,
    width = .9,
    aes(fill = value) # 设置根据变量系数的大小来进行渐变填充
  ) +
  scale_fill_gradientn(colours = mycolors[1:2]) +
  guides(fill = "none") +
  geom_text(
    aes(
      label = sprintf("%0.3f", value),
      x = 0,
      hjust = ifelse(value >= 0, 1.3, -0.3)
    ),
    colour = "black",
    size = 3,
    fontface = "italic"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  xlab("Coefficients") +
  ylab(NULL) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )
plot_coef


# 合并LASSO交叉验证图、系数图、系数条形图
plot_grid(
  plot_lasso.cv,
  plot_lasso,
  plot_coef,
  labels = "AUTO",
  rel_widths = c(80, 100, 80),
  nrow = 1
)
ggsave(width = 20, height = 5, "figures/cluster_surve/LASSO-Cox.pdf")







# 随机生存森林 ------------------------------------------------------------------

# 定义随机森林模型的方程
fun_forest <- tcga %>%
  select(ends_with("_TAMs")) %>%
  colnames() %>%
  str_c(collapse = " + ") %>%
  str_c("Surv(surv_time, vital_status) ~ ", .) %>%
  as.formula()
fun_forest


# 通过grid search进行随机森林调参，寻找最优nodesize和mtry参数组合
set.seed(123)
tune_forest <- tune(
  fun_forest,
  data = tcga,
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
    main = paste0("Optimal nodesize: ", best.nodesize,
                  "; Optimal mtry: ", best.mtry),
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





# 根据最优ntree和nodesize组合构建随机生存森林模型
set.seed(123)
rfsrc <- rfsrc(
  fun_forest,
  data = tcga,
  ntree = 4000,
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
set.seed(123)
rfsrc_final <- rfsrc(
  fun_forest,
  data = tcga,
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
rfsrc_final




# 绘制变量重要性图
importance <- gg_vimp(rfsrc_final) %>%
  as_tibble() %>%
  mutate(
    vars = factor(vars, levels = label_levels, labels = labels)
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
  geom_text(
    aes(
      label = sprintf("%0.4f", vimp),
      x = if_else(vimp <= 0, 0.0006, -0.0006)
    ),
    colour = "black",
    size = 3,
    fontface = "italic"
  ) +
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
ggsave(width = 18.4, height = 5.4, "figures/cluster_surve/randomForestSRC.pdf")





# 指定后续分析的目标巨噬细胞亚群
aim_cluster <- "CCL3L1+ Inflam-TAMs"







# 分析单细胞数据中目标细胞亚群在良恶性样本间的比例差异 ----------------------------------------------

p_prop <- seurat_macro@meta.data %>%
  as_tibble() %>%
  group_by(sample_id, sample_type, cell_annotation) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(total_count = sum(count)) %>%
  mutate(proportion = count / total_count) %>%
  filter(cell_annotation %in% aim_cluster) %>%
  ggplot(aes(x = cell_annotation, y = proportion * 100, fill = sample_type)) +
  geom_boxplot() +
  scale_fill_manual(values = mycolors) +
  scale_y_continuous(expand = c(0, 1)) +
  xlab(aim_cluster) +
  ylab("Proportion (%)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  stat_compare_means(
    aes(group = sample_type),
    method = "wilcox.test",
    vjust = -2
  )
p_prop

ggsave(
  width = 5.5,
  height = 6.7,
  str_c(
    "figures/cluster_surve/Comparison of ",
    aim_cluster,
    " proportions between normal and OC samples.pdf"
    )
  )





# 保存sessionInfo -----------------------------------------------------------

sink(
  "sessionInfo/8.2_Multivariate_survial_analyse.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
