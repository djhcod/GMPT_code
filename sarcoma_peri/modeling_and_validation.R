# 加载包及数据 ------------------------------------------------------------------

{
  rm(list = ls())
  cat("\014")

  library(tidyverse)
  library(ggsci)
  library(mice)
  library(rms)
  library(cowplot)
  library(patchwork)
  library(timeROC)
  library(caret)
  library(ggprism)
  library(ggDCA)
  library(survival)

  mycolors <- pal_jco()(9)
}


load("sarcoma_peri/data/data_for_RSF.Rdata")

# 提取插补后的数据集
dat_complete <- complete(imp_surv) %>% as_tibble()
dat_complete
summary(dat_complete)

# 提取原始数据集
dat_original <- dat[, c("age", "grade2", "figo", "chem", "peri", "dead", "time")] %>%
  as_tibble()
dat_original


# 拆分训练集和验证集 --------------------------------------------------------

set.seed(1)

# 创建训练集和验证集的索引
train_indices <- sample(
  seq_len(nrow(dat_complete)),
  size = 0.7 * nrow(dat_complete)
)

# 根据索引分割数据集
dat_complete$split <- "test"
dat_complete$split[train_indices] <- "train"

table(dat_complete$split)






# 建模 ----------------------------------------------------------------------

dat <- filter(dat_complete, split == "train")


dd <- datadist(dat)
options(datadist = "dd")

model <- cph(
  Surv(time, dead == 1) ~ age + grade2 + figo + chem + peri,
  x = T, y = T, data = dat, surv = T
)
model

cbind("回归系数" = coef(model), confint(model)) # 展示回归方程的系数及其95%CI
cbind("HR" = exp(coef(model)), exp(confint(model))) # 展示OR值（OR=e^β）及其95%CI




# 绘制nomogram
surv <- Survival(model)
surv1 <- function(x) surv(1 * 12, lp = x)
surv2 <- function(x) surv(3 * 12, lp = x)
surv3 <- function(x) surv(5 * 12, lp = x)

nom <- nomogram(model,
  fun = list(surv1, surv2, surv3), # 将预测值转换成列线图中的分值
  lp = F, # 是否显示原始预测值
  maxscale = 100, # nomogram记分轴的最大值
  funlabel = c(
    "1-Year Survival probability",
    "3-Year survival probability",
    "5-Year survival probability"
  ),
  fun.at = c("0.9", "0.85", "0.80", "0.70", "0.6", "0.5", "0.4", "0.3", "0.2", "0.1")
)
plot(nom, col.grid = gray(c(0.8, 0.95)))

recordPlot() %>% ggdraw()
ggsave(
  width = 9.6, height = 6.5,
  filename = "sarcoma_peri/figures/nomogram.pdf"
)






# 得到模型的预测值 --------------------------------------------------------------------

# 插补数据中执行模型预测
dat_complete$pred <- predict(model, as.data.frame(dat_complete))

dat_train <- filter(dat_complete, split == "train")
dat_test <- filter(dat_complete, split == "test")

# 在删除所有缺失值后的原始数据中执行模型预测
dat_all <- na.omit(dat_original)
dat_all$pred <- predict(model, as.data.frame(dat_all))





# 绘制ROC----------------------------------------------------------------------

# 批量在各个数据集中进行ROC分析
roc_list <- map(
  list(dat_train, dat_test, dat_all),
  function(x) {
    timeroc <- timeROC(
      T = x$time, # 指定随访时间列
      delta = as.numeric(x$dead) - 1, # 指定生存状态列
      marker = x$pred, # 指定预测值列
      cause = 1, # 指定感兴趣的结局事件
      weighting = "marginal",
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
    p <- ggplot(data = timeroc_df) +
      geom_line(aes(x = FP_1year, y = TP_1year), linewidth = 1, color = mycolors[1]) +
      geom_line(aes(x = FP_3year, y = TP_3year), linewidth = 1, color = mycolors[2]) +
      geom_line(aes(x = FP_5year, y = TP_5year), linewidth = 1, color = mycolors[4]) +
      geom_abline(slope = 1, intercept = 0, color = "grey", linewidth = 1, linetype = 2) +
      theme_bw() +
      annotate("text",
        x = 0.5, y = 0.25, size = 4.5,
        hjust = 0, # 左对齐
        label = paste0(
          "AUC at 1 year = ",
          sprintf("%.3f", timeroc$AUC[[1]])
        ),
        color = mycolors[1]
      ) +
      annotate("text",
        x = 0.5, y = 0.15, size = 4.5,
        hjust = 0,
        label = paste0(
          "AUC at 3 years = ",
          sprintf("%.3f", timeroc$AUC[[2]])
        ),
        color = mycolors[2]
      ) +
      annotate("text",
        x = 0.5, y = 0.05, size = 4.5,
        hjust = 0,
        label = paste0(
          "AUC at 5 years = ",
          sprintf("%.3f", timeroc$AUC[[3]])
        ),
        color = mycolors[4]
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
    return(p)
  }
)

wrap_plots(roc_list, nrow = 1) +
  plot_annotation(tag_levels = "A")
ggsave(
  width = 15, height = 4.9,
  filename = "sarcoma_peri/figures/roc.pdf"
)



timeroc_list <- map(
  list(dat_train, dat_test, dat_all),
  function(x) {
    timeroc <- timeROC(
      T = pull(x, time), # 指定随访时间列
      delta = (pull(x, dead) %>% as.numeric(x$dead)) - 1, # 指定生存状态列
      marker = pull(x, pred), # 指定预测值列
      cause = 1, # 指定感兴趣的结局事件
      weighting = "marginal",
      times = c(1 * 12, 3 * 12, 5 * 12), # 绘制1年、3年和5年的ROC
      ROC = TRUE, # 是否保存敏感度和特异度的预测值
      iid = T
    )
    plotAUCcurve(timeroc, conf.int = T, conf.band = TRUE, col = mycolors[1])
    p <- recordPlot() %>% ggdraw()
    return(p)
  }
)

timeroc_list %>%
  wrap_plots(roc_list, nrow = 1) +
  plot_annotation(tag_levels = "A")
ggsave(
  width = 15, height = 5,
  filename = "sarcoma_peri/figures/time_roc.pdf"
)




## AUC值的重复交叉验证----------------------------------------------------------

set.seed(1996)

# 设定交叉验证参数K和N，生成N*K个新数据集
folds <- createMultiFolds(
  y = dat_complete$dead,
  k = 10, # 设定为10折交叉验证
  times = 200
)


# 计算200次10折交叉验证的1年AUC值

auc_cv <- as.numeric() # 建立空表，放每次的AUC

for (i in 1:2000) { # 需要循环的次数=K*重复次数
  train <- dat_complete[folds[[i]], ] %>% as.data.frame()
  test <- dat_complete[-folds[[i]], ] %>% as.data.frame() # 标记验证集
  cvmodel <- cph(Surv(time, dead == 1) ~ age + grade2 + figo + chem + peri,
    x = T, y = T, data = train, surv = T
  )
  test$pred <- predict(cvmodel, test)
  roc <- timeROC(
    T = test$time, # 指定随访时间列
    delta = as.numeric(test$dead) - 1, # 指定生存状态列
    marker = test$pred, # 指定预测值列
    cause = 1, # 指定感兴趣的结局事件
    weighting = "marginal", # weighting="marginal"为采用Kaplan-Meier估计删失分布
    times = 12, # 绘制1年ROC
    ROC = TRUE, # 是否保存敏感度和特异度的预测值
    iid = F
  )
  auc_cv <- append(auc_cv, as.numeric(roc$AUC[2])) # 提取并汇总每次计算得到的AUC值
}

summary(auc_cv) # 展示1年交叉验证后的AUC的平均值
cv <- data.frame("roc_1year" = auc_cv) # 将所有AUC值以数据框的形式储存在cv表中



# 计算200次10折交叉验证的3年AUC值

auc_cv <- as.numeric()

for (i in 1:2000) {
  train <- dat_complete[folds[[i]], ] %>% as.data.frame()
  test <- dat_complete[-folds[[i]], ] %>% as.data.frame()
  cvmodel <- cph(Surv(time, dead == 1) ~ age + grade2 + figo + chem + peri,
    x = T, y = T, data = train, surv = T
  )
  test$pred <- predict(cvmodel, test)
  roc <- timeROC(
    T = test$time,
    delta = as.numeric(test$dead) - 1,
    marker = test$pred,
    cause = 1,
    weighting = "marginal",
    times = 12 * 3,
    ROC = TRUE,
    iid = F
  )
  auc_cv <- append(auc_cv, as.numeric(roc$AUC[2]))
}

summary(auc_cv)
cv <- data.frame(cv, "roc_3year" = auc_cv)




# 计算200次10折交叉验证的5年AUC值

auc_cv <- as.numeric()

for (i in 1:2000) {
  train <- dat_complete[folds[[i]], ] %>% as.data.frame()
  test <- dat_complete[-folds[[i]], ] %>% as.data.frame()
  cvmodel <- cph(Surv(time, dead == 1) ~ age + grade2 + figo + chem + peri,
    x = T, y = T, data = train, surv = T
  )
  test$pred <- predict(cvmodel, test)
  roc <- timeROC(
    T = test$time,
    delta = as.numeric(test$dead) - 1,
    marker = test$pred,
    cause = 1,
    weighting = "marginal",
    times = 12 * 5,
    ROC = TRUE,
    iid = F
  )
  auc_cv <- append(auc_cv, as.numeric(roc$AUC[2]))
}

summary(auc_cv)
cv <- data.frame(cv, "roc_5year" = auc_cv)




# AUC值交叉验证结果可视化
# 需要用melt()函数将构建的cv表格转换为两列
# 第一列是值标签（auc_1year、auc_3year、auc_5year）
# 第二列是相应的值

cv_plot <- pivot_longer(
  cv,
  cols = starts_with("roc_"), # 需要转换的列名
  names_to = "Groups", # 新数据集的标签列名
  values_to = "AUC" # 新数据集的数值列名
)
cv_plot

# 绘制AUC值200次10折交叉验证的小提琴图
ggplot(cv_plot, aes(Groups, AUC)) +
  geom_violin(aes(fill = Groups)) +
  geom_boxplot(width = 0.1) +
  scale_fill_manual(values = mycolors[c(1, 2, 4)]) +
  theme_prism(base_size = 15, border = T) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.title = element_blank()
  )
ggsave(width = 8.7, height = 5.7, filename = "sarcoma_peri/figures/auc_cv.pdf")




# 绘制校准曲线-----------------------------------------------------------------

calib_curve <- function(data) {
  calibrate_list <- map(
    list(12, 12 * 3, 12 * 5),
    function(x) {
      fit_cal <- cph(
        Surv(time, dead == 1) ~ pred,
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
        xlab = str_c("Nomogram-Predicted Probability of ", x / 12, "-year OS (%)"),
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
  calib_curve(dat_train),
  calib_curve(dat_test),
  calib_curve(dat_all),
  ncol = 1
)
ggsave(width = 15, height = 15, filename = "sarcoma_peri/figures/calibration_curve.pdf")





# DCA曲线绘制------------------------------------------------------------------


dca_list <- map(
  list(dat_train, dat_test, dat_all),
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
ggsave(width = 15, height = 6, filename = "sarcoma_peri/figures/dca.pdf")

# 输出训练集用于后续构建网页计算器
# 网页计算器的构建需要将分类变量转换成数值编码

dat_train %>%
  select(age, grade2, figo, chem, peri, time, dead) %>%
  mutate(
    grade2 = as.numeric(grade2) %>% factor(),
    figo = as.numeric(figo) %>% factor(),
    chem = (as.numeric(chem) - 1) %>% factor(),
    peri = (as.numeric(peri) - 1) %>% factor(),
    dead = (as.numeric(dead) - 1) %>% factor()
  ) %>% 
  saveRDS(file = "web_nomogram/data/sarcoma_peri_train_dataset.rds")
