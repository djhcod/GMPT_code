# 加载包及数据 ---------------------------------------------------------------------

{
  rm(list = ls())
  cat("\014")

  library(qs)
  library(tidyverse)
  library(survminer)
  library(timeROC)
  library(survival)
  library(cowplot)
  library(patchwork)
  library(ezcox)
  library(glmnet)
  library(rms)
  library(ggDCA)

  # 载入颜色集
  load("output/color_palette.rdata")

  # 指定输出目录
  out_dir <- "output/modeling/"
  fig_dir <- "figures/modeling/"
}

# 载入模型数据
load("output/modeling/gene_model.rdata")

model
aim_gene
nrow(tcga_filt)
glimpse(tcga_filt)

tcga_filt$gene_score <- predict(model, tcga_filt)$predicted
summary(tcga_filt$gene_score)


set.seed(1)

# 创建训练集和验证集的索引
train_indices <- sample(seq_len(nrow(tcga_filt)), size = 0.7 * nrow(tcga_filt))

# 根据索引分割数据集
tcga_filt$split <- "test"
tcga_filt$split[train_indices] <- "train"

table(tcga_filt$split)



# 基因评分预测能力评估 ---------------------------------------------------------

data <- tcga_filt %>%
  filter(split == "train")
timeroc <- timeROC(
  T = data$surv_time, # 指定随访时间列
  delta = data$vital_status, # 指定生存状态列
  marker = data$gene_score, # 指定预测值列
  cause = 1, # 指定感兴趣的结局事件
  weighting = "marginal", # weighting="marginal"为采用Kaplan-Meier估计删失分布
  times = c(1 * 12, 3 * 12, 5 * 12), # 绘制1年、3年和5年的ROC
  ROC = TRUE, # 是否保存敏感度和特异度的预测值
  iid = T
)

timeroc_df <- data.frame(
  TP_1year = timeroc$TP[, 1],
  FP_1year = timeroc$FP[, 1],
  TP_3year = timeroc$TP[, 2],
  FP_3year = timeroc$FP[, 2],
  TP_5year = timeroc$TP[, 3],
  FP_5year = timeroc$FP[, 3]
)




p1 <- ggplot(data = timeroc_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = mycolors[1]) +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = mycolors[2]) +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = mycolors[3]) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
    x = 0.75, y = 0.25, size = 4.5,
    label = paste0(
      "AUC at 1 year = ",
      sprintf("%.3f", timeroc$AUC[[1]])
    ),
    color = mycolors[1]
  ) +
  annotate("text",
    x = 0.75, y = 0.15, size = 4.5,
    label = paste0(
      "AUC at 3 years = ",
      sprintf("%.3f", timeroc$AUC[[2]])
    ),
    color = mycolors[2]
  ) +
  annotate("text",
    x = 0.75, y = 0.05, size = 4.5,
    label = paste0(
      "AUC at 5 years = ",
      sprintf("%.3f", timeroc$AUC[[3]])
    ),
    color = mycolors[3]
  ) +
  labs(x = "1-Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(face = "plain", size = 11, color = "black"),
    axis.title.x = element_text(
      face = "plain", size = 14, color = "black",
      margin = margin(c(15, 0, 0, 0))
    ),
    axis.title.y = element_text(
      face = "plain", size = 14, color = "black",
      margin = margin(c(0, 15, 0, 0))
    )
  )
p1

plotAUCcurve(timeroc, conf.int = T, conf.band = F, col = mycolors[2])
p2 <- recordPlot() %>% ggdraw()
p2




data <- tcga_filt %>%
  filter(split == "test")
timeroc <- timeROC(
  T = data$surv_time, # 指定随访时间列
  delta = data$vital_status, # 指定生存状态列
  marker = data$gene_score, # 指定预测值列
  cause = 1, # 指定感兴趣的结局事件
  weighting = "marginal", # weighting="marginal"为采用Kaplan-Meier估计删失分布
  times = c(1 * 12, 3 * 12, 5 * 12), # 绘制1年、3年和5年的ROC
  ROC = TRUE, # 是否保存敏感度和特异度的预测值
  iid = T
)

timeroc_df <- data.frame(
  TP_1year = timeroc$TP[, 1],
  FP_1year = timeroc$FP[, 1],
  TP_3year = timeroc$TP[, 2],
  FP_3year = timeroc$FP[, 2],
  TP_5year = timeroc$TP[, 3],
  FP_5year = timeroc$FP[, 3]
)



p3 <- ggplot(data = timeroc_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = mycolors[1]) +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = mycolors[2]) +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = mycolors[3]) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
    x = 0.75, y = 0.25, size = 4.5,
    label = paste0(
      "AUC at 1 year = ",
      sprintf("%.3f", timeroc$AUC[[1]])
    ),
    color = mycolors[1]
  ) +
  annotate("text",
    x = 0.75, y = 0.15, size = 4.5,
    label = paste0(
      "AUC at 3 years = ",
      sprintf("%.3f", timeroc$AUC[[2]])
    ),
    color = mycolors[2]
  ) +
  annotate("text",
    x = 0.75, y = 0.05, size = 4.5,
    label = paste0(
      "AUC at 5 years = ",
      sprintf("%.3f", timeroc$AUC[[3]])
    ),
    color = mycolors[3]
  ) +
  labs(x = "1-Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(face = "plain", size = 11, color = "black"),
    axis.title.x = element_text(
      face = "plain", size = 14, color = "black",
      margin = margin(c(15, 0, 0, 0))
    ),
    axis.title.y = element_text(
      face = "plain", size = 14, color = "black",
      margin = margin(c(0, 15, 0, 0))
    )
  )
p3

plotAUCcurve(timeroc, conf.int = T, conf.band = F, col = mycolors[2])
p4 <- recordPlot() %>% ggdraw()


plot_grid(
  p1, p3, p2, p4,
  labels = "AUTO",
  ncol = 2
)

ggsave(width = 10, height = 9, str_c(fig_dir, "ROC_gene_scores.pdf"))




## 高低基因评分患者的生存分析---------------------------------

cut_off <- surv_cutpoint(
  data,
  time = "surv_time",
  event = "vital_status",
  variables = "gene_score",
  minprop = 0.2
)
cut_off <- as.numeric(cut_off[["cutpoint"]][1, 1])

tcga_filt <- tcga_filt %>%
  mutate(
    risk_group = if_else(gene_score <= cut_off, "Low-risk", "High-risk") %>%
      factor(levels = c("Low-risk", "High-risk"))
  )
summary(tcga_filt$risk_group)






surv_data <- tcga_filt %>%
  filter(split == "train")

surv <- survfit(Surv(surv_time, vital_status == 1) ~ risk_group, data = surv_data)

# 手动计算P值
p_value <- surv_pvalue(surv, data = surv_data)$pval
p_value <- if_else(
  p_value < 0.001,
  "<0.001",
  sprintf("%.3f", p_value)
) %>%
  str_c("Log-rank test P value: ", .)
p_value

p1 <- ggsurvplot(
  surv,
  fun = "pct",
  palette = mycolors,
  conf.int = T,
  conf.int.style = "ribbon",
  pval = FALSE,
  ggtheme = theme_bw(),
  legend.title = "Risk group",
  xlab = "Time in months",
  break.time.by = 36,
  size = 0.6
)[["plot"]] +
  annotate(
    "text",
    x = 110,
    y = 100,
    label = p_value,
    size = 4
  )
p1




surv <- survfit(Surv(recurrence_time, recurrence == 1) ~ risk_group, data = surv_data)

# 手动计算P值
p_value <- surv_pvalue(surv, data = surv_data)$pval
p_value <- if_else(
  p_value < 0.001,
  "<0.001",
  sprintf("%.3f", p_value)
) %>%
  str_c("Log-rank test P value: ", .)
p_value

p2 <- ggsurvplot(
  surv,
  fun = "pct",
  palette = mycolors,
  conf.int = T,
  conf.int.style = "ribbon",
  pval = FALSE,
  ggtheme = theme_bw(),
  legend.title = "Risk group",
  xlab = "Time in months",
  break.time.by = 36,
  size = 0.6
)[["plot"]] +
  annotate(
    "text",
    x = 60,
    y = 100,
    label = p_value,
    size = 4
  ) +
  ylab("Recurrence probability (%)")
p2



surv_data <- tcga_filt %>%
  filter(split == "test")

surv <- survfit(Surv(surv_time, vital_status == 1) ~ risk_group, data = surv_data)

# 手动计算P值
p_value <- surv_pvalue(surv, data = surv_data)$pval
p_value <- if_else(
  p_value < 0.001,
  "<0.001",
  sprintf("%.3f", p_value)
) %>%
  str_c("Log-rank test P value: ", .)
p_value

p3 <- ggsurvplot(
  surv,
  fun = "pct",
  palette = mycolors,
  conf.int = T,
  conf.int.style = "ribbon",
  pval = FALSE,
  ggtheme = theme_bw(),
  legend.title = "Risk group",
  xlab = "Time in months",
  break.time.by = 36,
  size = 0.6
)[["plot"]] +
  annotate(
    "text",
    x = 150,
    y = 100,
    label = p_value,
    size = 4
  )
p3



surv <- survfit(Surv(recurrence_time, recurrence == 1) ~ risk_group, data = surv_data)

# 手动计算P值
p_value <- surv_pvalue(surv, data = surv_data)$pval
p_value <- if_else(
  p_value < 0.001,
  "<0.001",
  sprintf("%.3f", p_value)
) %>%
  str_c("Log-rank test P value: ", .)
p_value

p4 <- ggsurvplot(
  surv,
  fun = "pct",
  palette = mycolors,
  conf.int = T,
  conf.int.style = "ribbon",
  pval = FALSE,
  ggtheme = theme_bw(),
  legend.title = "Risk group",
  xlab = "Time in months",
  break.time.by = 36,
  xlim = c(0, 72),
  size = 0.6
)[["plot"]] +
  annotate(
    "text",
    x = 55,
    y = 100,
    label = p_value,
    size = 4
  ) +
  ylab("Recurrence probability (%)")
p4

plot_grid(
  p1, p2, p3, p4,
  nrow = 2,
  labels = "AUTO"
)

ggsave(width = 9, height = 8, str_c(fig_dir, "risk_group_survival_curves.pdf"))







# 单因素cox回归 ----------------------------------------------------------------

summary(tcga_filt)

if (FALSE) {
  age_cut_off <- surv_cutpoint(
    tcga_filt,
    time = "surv_time",
    event = "vital_status",
    variables = "age_at_diagnosis",
    minprop = 0.2
  )
  age_cut_off <- as.numeric(age_cut_off[["cutpoint"]][1, 1])

  tcga_filt <- tcga_filt %>%
    mutate(
      age_at_diagnosis = if_else(
        age_at_diagnosis <= age_cut_off,
        str_c("<=", age_cut_off),
        str_c(">", age_cut_off)
      ) %>%
        factor(levels = c(str_c("<=", age_cut_off), str_c(">", age_cut_off)))
    )
  summary(tcga_filt$age_at_diagnosis)


  tcga_filt <- tcga_filt %>%
    mutate(
      residual_disease = case_when(
        residual_disease == "No Macroscopic disease" ~ "No Macroscopic disease",
        residual_disease == "1-10 mm" ~ "1-10 mm",
        residual_disease == "11-20 mm" | residual_disease == ">20 mm" ~ ">10 mm",
      ) %>%
        factor(levels = c("No Macroscopic disease", "1-10 mm", ">10 mm"))
    )
  summary(tcga_filt$residual_disease)
}

covariates <- c(
  "age_at_diagnosis", "race", "figo_stage", "grade",
  "anatomic_subdivision", "residual_disease", "radiotherapy", "adjuvant_chem",
  "chemo_combined", "cycles_combined", "hormone_therapy", "bevacizumab",
  "gene_score", "risk_group"
)

data <- tcga_filt %>%
  filter(split == "train") %>%
  select(all_of(c(covariates, "vital_status", "surv_time")))


cox_tab <- ezcox(
  data,
  time = "surv_time",
  status = "vital_status",
  covariates = covariates,
  verbose = FALSE
)
cox_tab %>% print(width = Inf, n = Inf)

sig_vars_tab <- cox_tab %>%
  filter(p.value < 0.05)
print(sig_vars_tab, n = Inf)

write_csv(cox_tab, file = str_c(out_dir, "univariate_cox_regression.csv"))


sig_vars <- sig_vars_tab %>%
  pull(Variable) %>%
  unique()
sig_vars







# 多因素Cox回归 ----------------------------------------------------------------

model_vars <- c(
  "age_at_diagnosis", "residual_disease", "chemo_combined",
  "bevacizumab", "gene_score"
)

model_dat <- data %>%
  select(all_of(c(model_vars, "vital_status", "surv_time"))) %>%
  na.omit()

fun <- model_vars %>%
  str_c(collapse = " + ") %>%
  str_c("Surv(surv_time, vital_status) ~ ", .) %>%
  as.formula()
fun


dd <- datadist(model_dat)
options(datadist = "dd")

model <- cph(
  fun,
  x = T,
  y = T,
  data = model_dat,
  surv = T
)
model
summary(model)

cbind("回归系数" = coef(model), confint(model)) # 展示回归方程的系数及其95%CI
cbind("HR" = exp(coef(model)), exp(confint(model))) # 展示OR值（OR=e^β）及其95%CI


model_tab <- data.frame(
  "Coef" = model[["coefficients"]] %>% round(2),
  "SE" = summary(model)[, "S.E."] %>% na.omit() %>% round(2),
  "HR (95%CI)" = str_c(
    round(exp(coef(model)), 2),
    " (",
    round(exp(confint(model))[, 1], 2),
    ", ",
    round(exp(confint(model))[, 2], 2),
    ")"
    )
)
model_tab
write.csv(
  model_tab,
  file = str_c(out_dir, "multivariate_cox_regression.csv")
  )






# 得到模型的预测值
tcga_final <- tcga_filt %>%
  select(all_of(c(model_vars, "vital_status", "surv_time", "split"))) %>%
  na.omit() %>%
  mutate(model_pred = predict(model, as.data.frame(.)))
summary(tcga_final$model_pred)





# 绘制ROC----------------------------------------------------------------------

data <- tcga_final %>% filter(split == "train")

timeroc <- timeROC(
  T = data$surv_time, # 指定随访时间列
  delta = data$vital_status, # 指定生存状态列
  marker = data$model_pred, # 指定预测值列
  cause = 1, # 指定感兴趣的结局事件
  weighting = "marginal", # weighting="marginal"为采用Kaplan-Meier估计删失分布
  times = c(1 * 12, 3 * 12, 5 * 12), # 绘制1年、3年和5年的ROC
  ROC = TRUE, # 是否保存敏感度和特异度的预测值
  iid = T
)

timeroc_df <- data.frame(
  TP_1year = timeroc$TP[, 1],
  FP_1year = timeroc$FP[, 1],
  TP_3year = timeroc$TP[, 2],
  FP_3year = timeroc$FP[, 2],
  TP_5year = timeroc$TP[, 3],
  FP_5year = timeroc$FP[, 3]
)




p1 <- ggplot(data = timeroc_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = mycolors[1]) +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = mycolors[2]) +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = mycolors[3]) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
    x = 0.75, y = 0.25, size = 4.5,
    label = paste0(
      "AUC at 1 year = ",
      sprintf("%.3f", timeroc$AUC[[1]])
    ),
    color = mycolors[1]
  ) +
  annotate("text",
    x = 0.75, y = 0.15, size = 4.5,
    label = paste0(
      "AUC at 3 years = ",
      sprintf("%.3f", timeroc$AUC[[2]])
    ),
    color = mycolors[2]
  ) +
  annotate("text",
    x = 0.75, y = 0.05, size = 4.5,
    label = paste0(
      "AUC at 5 years = ",
      sprintf("%.3f", timeroc$AUC[[3]])
    ),
    color = mycolors[3]
  ) +
  labs(x = "1-Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(face = "plain", size = 11, color = "black"),
    axis.title.x = element_text(
      face = "plain", size = 14, color = "black",
      margin = margin(c(15, 0, 0, 0))
    ),
    axis.title.y = element_text(
      face = "plain", size = 14, color = "black",
      margin = margin(c(0, 15, 0, 0))
    )
  )
p1

plotAUCcurve(timeroc, conf.int = T, conf.band = F, col = mycolors[2])
p2 <- recordPlot() %>% ggdraw()
p2





timeroc <- timeROC(
  T = tcga_final$surv_time, # 指定随访时间列
  delta = tcga_final$vital_status, # 指定生存状态列
  marker = tcga_final$model_pred, # 指定预测值列
  cause = 1, # 指定感兴趣的结局事件
  weighting = "marginal", # weighting="marginal"为采用Kaplan-Meier估计删失分布
  times = c(1 * 12, 3 * 12, 5 * 12), # 绘制1年、3年和5年的ROC
  ROC = TRUE, # 是否保存敏感度和特异度的预测值
  iid = T
)

timeroc_df <- data.frame(
  TP_1year = timeroc$TP[, 1],
  FP_1year = timeroc$FP[, 1],
  TP_3year = timeroc$TP[, 2],
  FP_3year = timeroc$FP[, 2],
  TP_5year = timeroc$TP[, 3],
  FP_5year = timeroc$FP[, 3]
)




p3 <- ggplot(data = timeroc_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = mycolors[1]) +
  geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = mycolors[2]) +
  geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = mycolors[3]) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
    x = 0.75, y = 0.25, size = 4.5,
    label = paste0(
      "AUC at 1 year = ",
      sprintf("%.3f", timeroc$AUC[[1]])
    ),
    color = mycolors[1]
  ) +
  annotate("text",
    x = 0.75, y = 0.15, size = 4.5,
    label = paste0(
      "AUC at 3 years = ",
      sprintf("%.3f", timeroc$AUC[[2]])
    ),
    color = mycolors[2]
  ) +
  annotate("text",
    x = 0.75, y = 0.05, size = 4.5,
    label = paste0(
      "AUC at 5 years = ",
      sprintf("%.3f", timeroc$AUC[[3]])
    ),
    color = mycolors[3]
  ) +
  labs(x = "1-Specificity", y = "Sensitivity") +
  theme(
    axis.text = element_text(face = "plain", size = 11, color = "black"),
    axis.title.x = element_text(
      face = "plain", size = 14, color = "black",
      margin = margin(c(15, 0, 0, 0))
    ),
    axis.title.y = element_text(
      face = "plain", size = 14, color = "black",
      margin = margin(c(0, 15, 0, 0))
    )
  )
p3

plotAUCcurve(timeroc, conf.int = T, conf.band = F, col = mycolors[2])
p4 <- recordPlot() %>% ggdraw()
p4


plot_grid(
  p1, p3, p2, p4,
  labels = "AUTO",
  ncol = 2
)

ggsave(width = 10, height = 9, str_c(fig_dir, "ROC_final_model.pdf"))









# 校准曲线 --------------------------------------------------------------------

calib_curve <- function(data) {
  calibrate_list <- map(
    list(12, 12 * 3, 12 * 5),
    function(x) {
      fit_cal <- cph(
        Surv(surv_time, vital_status == 1) ~ model_pred,
        x = T,
        y = T,
        data = data,
        surv = T,
        time.inc = x
      )
      cal <- calibrate(
        fit_cal,
        u = x, # 要评价的时间节点
        method = "boot",
        cmethod = "KM",
        m = ceiling(nrow(data) / 4), # 设置每多少个对象为一个评估单位，数值越小节点越多
        B = 1000
      )
      plot(
        cal,
        lwd = 2,
        lty = 1,
        errbar.col = mycolors[1],
        xlab = str_c("Predicted probability of ", x / 12, "-year OS (%)"),
        ylab = str_c("Actual ", x / 12, "-year OS (%)"),
        col = mycolors[2],
        mgp = c(2, 1, 0)
      )
      p <- recordPlot() %>% ggdraw()
      return(p)
    }
  )
  calibrate_list %>%
    wrap_plots(nrow = 1) +
    plot_annotation(tag_levels = "A")
}

plot_grid(
  calib_curve(data),
  calib_curve(tcga_final),
  ncol = 1
)
ggsave(width = 13, height = 8, filename = str_c(fig_dir, "calibration_curves.pdf"))










# DCA曲线绘制------------------------------------------------------------------

dca_list <- map(
  list(data, tcga_final),
  function(x) {
    dca <- dca(
      model,
      new.data = as.data.frame(x),
      times = c(12, 12 * 3, 12 * 5)
    )

    p <- ggplot(
      dca,
      color = TRUE,
      linetype = FALSE,
      lwd = 0.8
    ) +
      scale_color_manual(
        values = mycolors,
        label = c(
          "1 year DCA", "3 year DCA", "5 year DCA",
          "ALL-1 year", "ALL-3 year", "ALL-5 year",
          "None"
        )
      ) +
      theme_bw() +
      theme(
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 10)
      )
    return(p)
  }
)

dca_list %>%
  wrap_plots(nrow = 1) +
  plot_annotation(tag_levels = "A")
ggsave(width = 11, height = 6, filename = str_c(fig_dir, "dca.pdf"))




save(list = ls(), file = str_c(out_dir, "final_model.rdata"))



# 导出模型数据用于网页构建
saveRDS(model, "output/model_dat_for_web/final_model.rds")








# 保存sessionInfo---------------------------------------------------------

sink(
  "sessionInfo/9.3_Final_model.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
