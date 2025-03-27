# 加载包及数据 ---------------------------------------------------------------------

{
  rm(list = ls())
  cat("\014")

  library(qs)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(survival)
  library(survminer)
  library(DataExplorer) # 探索缺失模式
  library(VIM) # 探索缺失模式
  library(mice)
  library(MatchThem)
  library(cobalt)
  library(RISCA)
  library(forestplot)

  # 载入颜色集
  load("output/color_palette.rdata")
}


# 载入此前整理后的去卷积结果、细胞标签顺序
load("output/BayesPrism/cell_prop_and_label_levels.rdata")


# 载入TCGA临床数据
tcga <- qread("output/tcga/TCGA_combined_clinical_dat.qs")
glimpse(tcga)

# 载入患者数记录
tcga_sum <- qread("output/tcga/TCGA_clinical_information_screening_process.qs")

# 重新编码生存状态变量
tcga <- tcga %>%
  mutate(
    vital_status = if_else(vital_status == "Alive", 0, 1),
    recurrence = if_else(recurrence == "Recurrence_free", 0, 1)
  ) %>%
  filter(surv_time > 0) # 去除生存时间<1个月的患者

# 记录患者数
tcga_sum <- tibble(
  "stage" = "Exclude survival time < 1 month",
  "number_of_patients" = nrow(tcga)
) %>%
  bind_rows(tcga_sum, .)
tcga_sum

# 合并细胞比例信息和TCGA临床信息
tcga <- left_join(
  tcga,
  cell_prop,
  by = "sample_id"
)
tcga %>% glimpse()

qsave(tcga, "output/tcga/TCGA_clinical_dat_with_cell_propotations.qs")





# 初步生存分析 --------------------------------------------------------------------

# 批量进行各个细胞亚群比例的生存分析
map(
  colnames(cell_prop)[-1],
  function(cluster) {
    # 确定各个细胞亚群比例的截断值
    cut_off <- surv_cutpoint(
      tcga,
      time = "surv_time",
      event = "vital_status",
      variables = cluster,
      minprop = 0.2
    )
    cut_off <- as.numeric(cut_off[["cutpoint"]][1, 1])
    print(cut_off)

    # 根据截断值将各亚群细胞比例分为低和高
    group <<- if_else(pull(tcga, cluster) < cut_off, "Low", "High") %>%
      factor(levels = c("Low", "High"))

    # 将分组信息返回到TCGA数据中
    tcga <<- tcga %>%
      mutate(
        !!str_c(cluster, "_group") := if_else(pull(tcga, cluster) < cut_off, "Low", "High") %>%
          factor(levels = c("Low", "High"))
      )

    # 生存分析
    surv <- survfit(
      Surv(surv_time, vital_status == 1) ~ group,
      data = tcga
    )

    p <- ggsurvplot(
      surv,
      fun = "pct",
      palette = mycolors,
      conf.int = T,
      conf.int.style = "ribbon",
      pval = FALSE,
      ggtheme = theme_bw(),
      legend.title = cluster %>% str_to_upper(),
      xlab = "Time in months",
      break.time.by = 36
    )
    # 手动计算P值
    p_value <- surv_pvalue(surv, data = tcga)$pval
    p_value <- if_else(
      p_value < 0.001,
      "<0.001",
      sprintf("%.3f", p_value)
    ) %>%
      str_c("Log-rank test P value: ", .)

    str_glue(
      "-----------------------------------------------------------------\n",
      "{cluster}\n",
      "-----------------------------------------------------------------"
    ) %>%
      print()
    print(surv)
    summary(surv, time = c(12, 12 * 3, 12 * 5)) %>% print()
    print(p_value)

    return(
      p[["plot"]] +
        annotate(
          "text",
          x = 140,
          y = 100,
          label = p_value,
          size = 4
        )
    )
  }
) %>%
  wrap_plots(ncol = 5)

ggsave(width = 22, height = 9, "figures/cluster_surve/OS_curve.pdf")
qsave(tcga, "output/cluster_surve/TCGA_clinical_dat_with_cell_prop_group.qs")


glimpse(tcga)
# 统计各个亚群高低比例的样本数
tcga %>%
  select(ends_with("_group")) %>%
  summary()









# 多重插补---------------------------------------------------------------

## 插补前准备-----------------------------------------------------------
tcga_select <- tcga %>%
  select(
    year_of_diagnosis:radiotherapy,
    chemo_combined, cycles_combined,
    hormone_therapy:surv_time,
    APOE_LP_TAMs_group:CXCL10_IFN_TAMs_group
  )

# 分析数据集的基本情况
plot_intro(tcga_select, ggtheme = theme_bw())

# 统计各变量的缺失比例
profile_missing(tcga_select) %>%
  arrange(., desc(pct_missing)) %>%
  print(n = Inf)

# 可视化各变量的缺失比例
p_miss <- plot_missing(
  tcga_select,
  ggtheme = theme_bw(),
  group_color = list(
    "Good" = mycolors[1],
    "OK" = mycolors[2],
    "Bad" = mycolors[3],
    "Remove" = mycolors[4]
  )
) +
  xlab(NULL) +
  ylab("Total number of missing values") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.1))) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    legend.position = "top",
    plot.margin = unit(c(5, 0, 12, 0), "mm")
  )
p_miss

# 可视化缺失模式
{
  aggr(
    tcga_select,
    prop = FALSE,
    numbers = TRUE,
    sortVars = TRUE,
    cex.axis = 0.5,
    col = mycolors,
    combined = TRUE
  )
  p_aggr <- recordPlot() %>%
    ggdraw()
}

plot_grid(p_aggr, p_miss, rel_widths = c(100, 100))
ggsave(width = 12, height = 7, "figures/cluster_surve/Missing_pattern.pdf")


# 排除缺失值数量>3个的个案
tcga_filt <- tcga_select %>%
  filter(rowSums(is.na(tcga_select)) <= 3)

# 记录患者数
tcga_sum <- tibble(
  "stage" = "Exclusion of cases with >3 missing values",
  "number_of_patients" = nrow(tcga_filt)
) %>%
  bind_rows(tcga_sum, .)
tcga_sum

write_csv(
  tcga_sum,
  "output/cluster_surve/TCGA_clinical_information_screening_process.csv"
)



## 插补-----------------------------------------------------------

# 简单运行插补获取predictorMatrix矩阵和默认插补算法
mice_obj <- mice(tcga_filt, maxit = 0, print = FALSE)
names(mice_obj)

# 获取并更改predictorMatrix
pred <- mice_obj$predictorMatrix
colnames(pred)
pred[, which(colnames(pred) == "APOE_LP_TAMs_group"):ncol(pred)] <- 0
pred[which(rownames(pred) == "APOE_LP_TAMs_group"):nrow(pred), ] <- 0

# 各变量的插补方式
mice_obj$method



# 正式插补
imp_dat <- mice(
  tcga_filt,
  m = 10,
  maxit = 30,
  predictorMatrix = pred,
  seed = 123,
  print = FALSE
)
imp_dat



# 提取有缺失值的变量
vars_miss <- tcga_filt[, is.na(tcga_filt) %>% colSums() > 0] %>% colnames()
vars_miss

# 检查有缺失值的变量的插补情况
summary(tcga_filt[, vars_miss])
summary(complete(imp_dat)[, vars_miss])

# 评估插补数据集的分布情况
# 折线图：评估各变量在各插补数据集及迭代过程中平均值和标注差的变化情况
pdf(
  width = 8,
  height = 6,
  "figures/cluster_surve/Mean_and_SD_distributions_for_imputed_data_sets.pdf"
)
plot(imp_dat, col = mycolors)
dev.off()

# 散点图：评估在插补数据集中各变量的取值
stripplot(
  imp_dat,
  str_c(vars_miss, collapse = "+") %>% str_c(., " ~ .imp") %>% as.formula(),
  pch = 20,
  col = mycolors
)

# 绘制各缺失值变量插补前后的密度图
map(
  vars_miss,
  function(x) {
    str_c("~", x) %>%
      as.formula() %>%
      densityplot(imp_dat, ., col = mycolors) %>%
      ggdraw()
  }
) %>%
  wrap_plots(ncol = 4)




# 多重插补后匹配---------------------------------------------------------------

vars_x <- tcga_filt %>%
  select(-ends_with("_group"), -c(vital_status, surv_time)) %>%
  colnames()
vars_x

vars_y <- tcga_filt %>%
  select(ends_with("_group")) %>%
  colnames()
vars_y


# 批量进行各个因变量的PSM
mimids_list <- map(
  vars_y,
  function(x) {
    # 生成PS计算公式
    fun <- str_c(vars_x, collapse = "+") %>%
      str_c(x, "~", .) %>%
      as.formula()
    print(fun)

    # 匹配
    imp_match <- matchthem(
      fun,
      datasets = imp_dat,
      method = "full",
      distance = "glm",
      approach = "across"
    )
    return(imp_match)
  }
)
names(mimids_list) <- vars_y
complete(mimids_list[[1]]) %>% head()

qsave(mimids_list, "output/cluster_surve/mimids_list.qs")



# 均衡性检验
balance_tab <- map(
  names(mimids_list),
  function(x) {
    tab_bal <- bal.tab(
      mimids_list[[x]],
      un = TRUE, # 是否计算匹配前的均衡性
      binary = "std", # 让二元协变量也输出SMD
      stats = c("mean.diffs", "ks.statistics"),
      thresholds = c(m = .2), # 设定各检验的阈值
      imp.fun = "max",
      imbalanced.only = FALSE
    )
    print(tab_bal)
    tab_bal[["Balance.Across.Imputations"]] %>%
      select(
        Max_SMD_before_PSM = Max.Diff.Un,
        Max_SMD_after_PSM = Max.Diff.Adj
      ) %>%
      rownames_to_column("Features") %>%
      {
        mutate(., Cluster = rep(x, nrow(.)), .before = 1)
      } %>%
      return()
  }
) %>%
  list_rbind()
head(balance_tab)

# 提取未能达到匹配标准的变量
unbalance_var <- balance_tab %>%
  filter(Max_SMD_after_PSM >= 0.2)
unbalance_var

write_csv(balance_tab, "output/cluster_surve/Balance_table.csv")




# 可视化展示各变量在匹配前后的均衡性

# 创建一个映射变量名的数据框，便于在图上显示正式的变量名
var_names <- data.frame(
  old = vars_x,
  new = c(
    "Year of diagnosis", "Age at diagnosis", "Race", "FIGO stage", "Grade",
    "Anatomic subdivision", "Residual disease", "Radiotherapy", "Chemotherapy",
    "Chemotherapy cycles", "Hormone therapy", "Bevacizumab"
  )
)
var_names

# 展示匹配前后SMD变化的气泡图
p_loveplot <- map(
  names(mimids_list),
  ~ love.plot(
    mimids_list[[.x]],
    stats = c("mean.diffs"), # 均衡性计算方法
    abs = TRUE, # 是否对均衡性数值取平均
    agg.fun = "range", # 均衡性汇总方式
    drop.distance = F, # 是否在顶部显示distance即PS或weight在匹配前后的情况
    thresholds = c(mean.diffs = 0.2),
    line = F, # 是否绘制各变量的连线
    var.order = "unadjusted",
    var.names = var_names,
    sample.names = c("Unmatched", "Matched"),
    binary = "std",
    stars = "raw",
    shapes = c("circle filled", "circle"),
    colors = mycolors[1:2],
    position = "top", # 图例的位置
    grid = F
  ) +
    scale_x_continuous(limits = c(0, 0.6)) +
    xlab("Absolute SMD") +
    ggtitle(
      .x %>%
        str_remove("_group") %>%
        str_to_upper() %>%
        str_replace_all("_", "-")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12),
      plot.subtitle = element_blank(),
      legend.position = "none",
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8)
    )
) %>%
  wrap_plots(ncol = 5)
p_loveplot


# 展示PS在匹配前后的密度图
p_balplot <- map(
  names(mimids_list),
  ~ bal.plot(
    mimids_list[[.x]],
    var.name = "distance",
    which = "both",
    which.imp = .none,
    colors = mycolors[1:2],
    disp.means = T,
    sample.names = c("Unmatched", "Matched")
  ) +
    ggtitle(
      .x %>%
        str_remove("_group") %>%
        str_to_upper() %>%
        str_replace_all("_", "-")
    ) +
    xlab("Distance") +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12),
      axis.title = element_text(size = 10)
    )
) %>%
  wrap_plots(ncol = 5)
p_balplot

plot_grid(
  p_loveplot,
  p_balplot,
  ncol = 1,
  rel_heights = c(100, 60),
  labels = "AUTO"
)

ggsave(width = 24, height = 12, "figures/cluster_surve/PSM_plot.pdf")








# 匹配后生存分析 -----------------------------------------------------------------

map(
  vars_y,
  function(x) {
    group <<- pull(tcga_filt, x)
    surv <- survfit(
      Surv(surv_time, vital_status == 1) ~ group,
      data = tcga_filt,
      weights = mimids_list[[x]][["models"]][[1]][["weights"]],
      robust = TRUE, # To request cluster-robust standard errors
      cluster = mimids_list[[x]][["models"]][[1]][["subclass"]]
    )

    p <- ggsurvplot(
      surv,
      fun = "pct",
      palette = mycolors,
      conf.int = T,
      conf.int.style = "ribbon",
      pval = FALSE,
      ggtheme = theme_bw(),
      legend.title = x %>% str_to_upper(),
      xlab = "Time in months",
      break.time.by = 36
    )
    print(surv)
    summary(surv, time = c(12, 12 * 3, 12 * 5)) %>% print()

    # 对于匹配后仍然没有达到均衡的变量通过cox回归得到校正这些变量后的P值
    if (x == "CCL3L1_Inflam_TAMs_group") {
      adj_p <- with(
        imp_dat,
        coxph(
          Surv(surv_time, vital_status == 1) ~ CCL3L1_Inflam_TAMs_group + residual_disease,
          robust = TRUE,
          cluster = mimids_list[[x]][["models"]][[1]][["subclass"]],
          weights = mimids_list[[x]][["models"]][[1]][["weights"]]
        )
      ) %>%
        pool() %>%
        summary() %>%
        .$p.value %>%
        .[1]
      adj_p <- ifelse(adj_p < 0.001, "< 0.001", sprintf("%.3f", adj_p)) %>%
        str_c("*Adjusted P value: ", .)
    } else if (x == "MT1H_MT_TAMs_group") {
      adj_p <- with(
        imp_dat,
        coxph(
          Surv(surv_time, vital_status == 1) ~ MT1H_MT_TAMs_group + cycles_combined,
          robust = TRUE,
          cluster = mimids_list[[x]][["models"]][[1]][["subclass"]],
          weights = mimids_list[[x]][["models"]][[1]][["weights"]]
        )
      ) %>%
        pool() %>%
        summary() %>%
        .$p.value %>%
        .[1]
      adj_p <- ifelse(adj_p < 0.001, "< 0.001", sprintf("%.3f", adj_p)) %>%
        str_c("*Adjusted P value: ", .)
    } else {
      adj_p <- ipw.log.rank(
        times = tcga_filt$surv_time,
        failures = tcga_filt$vital_status,
        variable = pull(tcga_filt, x) %>% as.numeric() - 1,
        weights = mimids_list[[x]][["models"]][[1]][["weights"]]
      )$p.value
      adj_p <- ifelse(adj_p < 0.001, "< 0.001", sprintf("%.3f", adj_p)) %>%
        str_c("Adjusted P value: ", .)
    }
    print(adj_p)

    return(
      p[["plot"]] +
        annotate(
          "text",
          x = 150,
          y = 100,
          label = adj_p,
          size = 4
        )
    )
  }
) %>%
  wrap_plots(ncol = 5)

ggsave(width = 22, height = 9, "figures/cluster_surve/OS_curve_PSM.pdf")














# 保存sessionInfo -----------------------------------------------------------

sink(
  "sessionInfo/8.1_Survival_analysis.txt",
  append = FALSE,
  split = FALSE
)
sessionInfo()
sink()
