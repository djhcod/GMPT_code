# 环境配置-----------------------------------------------------------------------
packagecheck<-function(package.names) {
  CRANpackages <- available.packages()
  result = c()
  for(i in 1:length(package.names)) {
    if(!package.names[i] %in% CRANpackages) {
      result[i] = paste0("【!】", package.names[i],"包无法在CRAN中找到")
    } else if(!require(package.names[i], character.only=T)) {
      install.packages(package.names[i])
      result[i] = paste0(package.names[i],"包已下载并加载")
    } else {
      result[i] = paste0(package.names[i],"包已加载")
    }
  }
  result
}

packages <-c("dplyr", "DataExplorer", "ggsci", "scales", "ggplot2", "VIM","qqplotr",
             "cowplot", "nortest", "gtsummary", "flextable", "mice", "openxlsx2", "rpart", 
             "rpart.plot", "partykit", "survival", "survminer", "cmprsk", "cobalt", 
             "MatchThem", "RISCA", "randomForestSRC", "ggRandomForests", "meta", "forestplot")
packagecheck(packages)
rm(list = ls())



# 加载数据集 -----------------------------------------------------------------------------
load("sarcoma_peri/data/lms_ess_2010_2016.Rdata")
dat <- mydata
rm(mydata)


# 🔴数据检查 --------------------------------------------------------------------------------------
summary(dat)
library(dplyr)

#生成数据探索性分析综合报告
if(F){
  library(DataExplorer)
  configure_report(global_ggtheme = quote(theme_bw())) %>%
    create_report(dat, 
                  y = "peri", 
                  config = ., 
                  output_file = "Data Profiling Report.html", 
                  output_dir = "sarcoma_peri/output")
}


#定义颜色集
library(ggsci)
#pal_jco pal_npg pal_futurama pal_frontiers pal_lancet pal_simpsons pal_flatui
mycolors <- pal_jco()(9)
mycolors
library(scales)
show_col(mycolors)


library(ggplot2)
# 绘制数据基本信息图
plot_intro(dat) + 
  theme_bw() + 
  theme(panel.grid.major.y = element_blank()) #移除垂直网格线
# 绘制缺失值统计图
plot_missing(dat) + 
  scale_x_discrete(expand = expansion(mult = c(0.06, 0.06))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.09))) +
  theme_bw() + 
  theme(panel.grid.major.y = element_blank())
# 绘制缺失模式图
profile_missing(dat) %>% arrange(., desc(pct_missing))
library(VIM)
aggr(dat,
     prop = T,
     numbers = T,
     sortVars = TRUE,
     gap = 0,
     cex.axis = 0.6,
     col = mycolors,
     ylab = c("Histogram of missing data","Pattern"))
aggr(dat,
     prop = F,
     numbers = T,
     sortVars = TRUE,
     gap = 0,
     cex.axis = 0.6,
     col = mycolors,
     ylab = c("Histogram of missing data","Pattern"))


# 🔴绘制基线特征表 -----------------------------------------------------------------------------------
#多维度正态性检验自编函数(直方图、Q-Q图、正态性参数检验)
library(ggplot2)
library(qqplotr) #Q-Q图
library(cowplot) #合并图形
library(nortest) #Kolmogorov-Smirnov正态性检验
library(dplyr) #select()函数
#定义自编函数distr.test()
distr.test <- function(data, vector) {
  for (i in 1:length(vector)) {
    #首先去除待检验变量的缺失值
    dat <- select(data, vector[i]) %>% na.omit()
    #直方图
    binwidth <- (max(dat) - min(dat)) / 15
    p1 <- ggplot(dat, 
                 aes(!!sym(vector[i]))) + #把变量名转换为symbol
      geom_histogram(binwidth = binwidth, #组距
                     colour = "white", 
                     fill = mycolors[2]) +
      geom_density(eval(bquote(aes(y = after_stat(count) * binwidth))), #绘制核密度曲线
                   colour = mycolors[1], 
                   fill = mycolors[1], 
                   alpha = 0.3) + 
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      theme_bw() + 
      theme(panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank())
    labs(x = vector[i], 
         y = "Frequency")
    #Q-Q图
    p2 <- ggplot(data = dat, 
                 mapping = aes(sample = !!sym(vector[i]))) + 
      stat_qq_line() +
      stat_qq_point() +
      stat_qq_band(fill = mycolors[1], alpha = 0.3) +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles")+
      theme_bw()
    #合并两幅图
    p <- plot_grid(p1, p2, 
                   ncol = 2, 
                   nrow = 1,
                   rel_widths = c(1,1))
    print(p)
    #正态性参数检验
    if(nrow(dat) < 5000) {
      test <- shapiro.test(dat[,vector[i]]) 
      p_value <- test[["p.value"]]
      ifelse(p_value < 0.001, print(paste0("【", vector[i], "】的Shapiro–Wilk正态性检验：P < 0.001")), 
             print(paste0("【", vector[i], "】的Shapiro–Wilk正态性检验：P = ", round(p_value, 3))))
    } else {
      test <- lillie.test(dat[,vector[i]]) 
      p_value <- test[["p.value"]]
      ifelse(p_value < 0.001, 
             print(paste0("【", vector[i], "】的Kolmogorov-Smirnov正态性检验：P < 0.001")), 
             print(paste0("【", vector[i], "】的Kolmogorov-Smirnov正态性检验：P = ", round(p_value, 3))))
    }
  }
}

select(dat, where(is.numeric)) %>% colnames() #选择连续型变量
c("age","tumor_size") %>% #定义需要进行正态性检验的连续型变量
  distr.test(dat, .) #应用自编函数进行多维度正态性检验

#指定所有需要统计的变量
vars <- c("year2", "age", "race", "marriage", "income", "grade2", "tumor_size", 
          "his", "T_stage", "N_stage", "M_stage", "figo", "surg", "plnd", "alnd", 
          "lnd", "rad", "chem")

#绘制基线特征表
library(gtsummary)
# 修改分组标签
tab_base <- tbl_summary(dat,
                        by = peri, #定义分组变量
                        include = all_of(vars), #指定需要比较的变量
                        #给变量添加标签
                        label = list(year2 ~ "Year of Diagnosis", 
                                     age ~ "Age (years)",
                                     race ~ "Race",
                                     marriage ~ "Marital status",
                                     income ~ "Family income",
                                     grade2 ~ "Grade",
                                     tumor_size ~ "Tumor size (mm)",
                                     his ~ "Histology",
                                     T_stage ~ "T stage",
                                     N_stage ~ "N stage",
                                     M_stage ~ "M stage",
                                     figo ~ "FIGO stage",
                                     surg ~ "Surgery",
                                     plnd ~ "Pelvic lymph node dissection",
                                     alnd ~ "Para-aortic lymph node dissection",
                                     lnd ~ "Lymphadenectomy",
                                     rad ~ "Radiotherapy",
                                     chem ~ "Chemotherapy"),
                        #定义各变量的表示形式
                        statistic = list(age ~ "{mean} ({sd})", #定义年龄以均数+标准差展示
                                         #定义所有分类变量以频数+百分比的形式展示
                                         all_categorical() ~ c("{n} ({p}%)")),
                        #定义小数位数
                        digits = list(all_continuous() ~ 1, #所有连续变量保留1位小数
                                      all_categorical() ~ c(0, 1)), #分类变量频数保留整数，百分比保留1位小数
                        #定义所有二分类变量为分类变量，从而让每个水平都展示
                        type = list(all_dichotomous() ~ "categorical"), 
                        missing = "ifany", #定义缺失值的展示：“ifany”有缺失值则展示；“no”不统计缺失值
                        #指定百分比是按行还是按列统计，默认为按列统计
                        percent = "column") %>%
  add_overall() %>%
  #添加P值。连续变量默认采用Wilcoxon秩和检验。
  #分类变量若预期单元格计数>=5，默认使用“chisq.test.no.right”，即Pearson卡方检验；
  #若预期单元格计数<5，则默认使用“Fisher.test”，即Fisher精确概率法检验。
  add_p(test = list(age ~ "t.test"), #指定t检验的变量
        pvalue_fun = ~style_pvalue(.x, digits = 3)) %>% #P值保留3位小数
  separate_p_footnotes() %>%
  add_significance_stars(thresholds = c(0.05)) %>%
  # 修改表头文字
  modify_header(p.value = "**P value**") %>%
  # 添加跨列表头
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Peritoneal Cytology**") %>%
  # 加粗标签文字
  bold_labels()
tab_base

#将基线特征表输出成word
library(flextable)
as_flex_table(tab_base) %>% 
  save_as_docx(tab_base, align = "left", path = "sarcoma_peri/output/baseline_table.docx")



# 绘制条形图展示不同分期下的阳性腹膜细胞学比例
p_T <- ggplot(dat, 
              aes(x = T_stage, fill = peri)) +
  geom_bar(position = "fill",
           stat = "count", 
           width = .9) + 
  scale_fill_manual(values = mycolors[1:2]) + 
  ylab("Relative proportions") +  
  xlab("T stage") +
  labs(fill = "Peritoneal cytology")+
  theme_bw() +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        plot.margin = margin(0, 0, 0, 10))
p_N <- ggplot(na.omit(dat[,c("N_stage", "peri")]), 
              aes(x = N_stage, fill = peri)) +
  geom_bar(position = "fill",
           stat = "count", 
           width = .9) + 
  scale_fill_manual(values = mycolors[1:2]) + 
  xlab("N stage")+
  ylab(NULL) +
  theme_bw() +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        plot.margin = margin(0, 0, 0, 20))
p_M <- ggplot(dat, 
              aes(x = M_stage, fill = peri)) +
  geom_bar(position = "fill",
           stat = "count", 
           width = .9) + 
  scale_fill_manual(values = mycolors[1:2]) + 
  xlab("M stage")+
  ylab(NULL) +
  theme_bw() +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        plot.margin = margin(0, 0, 0, 20))
p_figo <- ggplot(dat, 
                 aes(x = figo, fill = peri)) +
  geom_bar(position = "fill",
           stat = "count", 
           width = .9) + 
  scale_fill_manual(values = mycolors[1:2]) + 
  xlab("FIGO stage")+
  ylab(NULL) +
  theme_bw() +
  theme(panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        plot.margin = margin(0, 0, 0, 20))

get_legend(p_T + theme(legend.box.margin = margin(0, 0, 0, -60))) %>%
  plot_grid(p_T + theme(legend.position="none"), 
            p_N + theme(legend.position="none"), 
            p_M + theme(legend.position="none"), 
            p_figo + theme(legend.position="none"),
            .,
            nrow  = 1,
            axis = "tb",
            labels = c("A", "B", "C", "D"))



# 🔴分析影响腹膜细胞学的因素 ------------------------------------------------------------------------------
vars_logi <- c("year2", "age", "race", "marriage", "income", "grade2", "tumor_size", "his", 
               "T_stage", "N_stage", "M_stage", "figo", "peri")
dat_logi <- dat[ ,vars_logi]
dat_logi2 <- dat_logi[rowSums(is.na(dat_logi)) <= 3 & is.na(dat_logi$N_stage) == F, ]
paste0("【排除缺失值数量>3个和N分期未知的个案后】排除：", nrow(dat_logi) - nrow(dat_logi2), 
       "个；剩余：", nrow(dat_logi2), "个") %>% print()

library(DataExplorer)
plot_intro(dat_logi2, ggtheme = theme_bw())
plot_missing(dat_logi2, ggtheme = theme_bw())
profile_missing(dat_logi2) %>% arrange(., desc(pct_missing))
aggr(dat_logi2,
     prop=T,
     numbers=T,
     sortVars=TRUE,
     gap=0,
     cex.axis = 0.6,
     col = mycolors,
     ylab=c("Histogram of missing data","Pattern"))
# library(naniar)
# mcar_test(dat_logi2)
#提取包含缺失值的变量名
vars_miss <- colnames(dat_logi2[,colSums(is.na(dat_logi2)) > 0])


# 🔴多重插补
library(mice)
imp <- mice(dat_logi2, m = 10, maxit = 5, seed = 1996, print = FALSE)
imp
#该图显示了各缺失变量在各插补数据中的平均值(左)和标准差(右)。
#通常，我们希望这些线条混合在一起，并且在以后的迭代中不会出现任何趋势。
plot(imp, col = mycolors)
#绘制条带图
stripplot(imp, marriage + race + tumor_size + N_stage ~ .imp, 
          pch = 20, 
          col = mycolors)
#绘制各缺失值变量插补前后的密度图
imp.densityplot <- function(imp, vars_miss, col) {
  for (i in 1:length(vars_miss)) {
    paste0("~", vars_miss[i]) %>% as.formula() %>%
      densityplot(imp, ., col=col) %>% print()
  }
}
imp.densityplot(imp, vars_miss, mycolors)

# 基于多重插补数据的批量单因素logistic回归分析
pool.glm <- function(imp, x, y) {
  result <- c()
  for (i in 1:length(x)) {
    fit <- with(imp, glm(as.formula(paste0(y,"~",x[i])), family = binomial()))
    pool.fit <- pool(fit)
    sum_pool <- summary(pool.fit, conf.int = TRUE)
    sum_pool_or <- summary(pool.fit, conf.int = TRUE, exponentiate = TRUE)
    Variables <- as.character(sum_pool$term)
    Coefficient <- round(sum_pool$estimate, 2)
    OR <- round(sum_pool_or$estimate, 2)
    LCI <- round(sum_pool_or$`2.5 %`, 2)
    UCI <- round(sum_pool_or$`97.5 %`, 2)
    P_value <- ifelse(sum_pool_or$p.value < 0.001, "<0.001", 
                      round(sum_pool_or$p.value, 3))
    single_result <- data.frame(Variables, Coefficient, OR, LCI, UCI, P_value)[-1,]
    result <- rbind(result, single_result)
  }
  result
}
pool.glm_result <- c("year2", "age", "race", "marriage", "income", "grade2", 
                     "tumor_size", "his", "T_stage", "N_stage", "M_stage", "figo") %>%
  pool.glm(imp, ., "peri")
pool.glm_result
pool.glm_result$OR <- paste0(format(pool.glm_result$OR, nsmall = 2), " (", 
                             format(pool.glm_result$LCI, nsmall = 2), ", ", 
                             format(pool.glm_result$UCI, nsmall = 2), ")")
pool.glm_result <- pool.glm_result[,-c(4,5)]
pool.glm_result

library(openxlsx2)
write_xlsx(pool.glm_result, "sarcoma_peri/output/glm_result.xlsx")


#决策树分析
#筛选P<0.1的变量
pool.glm_result[pool.glm_result$P_value < 0.1, "Variables"] %>% unique()
#提取决策树分析的数据集并移除缺失值
dat_tree <- na.omit(dat[,c("grade2", "tumor_size", "his",
                           "T_stage", "N_stage", "M_stage", 
                           "peri")])
library(rpart)
set.seed(1996)
tree <- rpart(peri ~ grade2 + tumor_size + his + T_stage + N_stage + M_stage, 
              data = dat_tree, 
              method = "class")
tree$cptable
#决策树绘制
library(rpart.plot)
prune <- prune(tree,0.01)
rpart.plot(prune,
           type = 2,
           extra = 1+100,
           box.palette = mycolors[1:4])


library(partykit)
plot(as.party(prune))#展示各组阳性事件占比

#可视化每个变量的重要性
prune$variable.importance %>%
  data.frame(var = names(.), impor = .) %>%
  {
    ggplot(., aes(x = reorder(var,-impor), y = impor))+
      geom_bar(stat = "identity", fill = mycolors[1:nrow(.)])+
      labs(x = "Variables", y = "Importance")+
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + #让Y轴从0开始
      theme_bw() +
      theme(panel.grid.major.x=element_blank()) #移除垂直网格线
  }






# 🔴初步预后分析 -----------------------------------------------------------------------------------
# 计算中位随访时间
summary(dat$time)

# 总人群生存分析
library(survival)
surv <- survfit(Surv(time, dead == 1) ~ 1, data = dat)
#绘制生存曲线
library(survminer)
ggsurvplot(surv,
           fun="pct", palette = mycolors,
           conf.int=T, conf.int.style = "ribbon", conf.int.alpha = 1, 
           surv.median.line = "hv",
           risk.table="abs_pct",risk.table.col="black",risk.table.y.text=F,
           risk.table.height=0.2,
           ggtheme=theme_bw(),
           xlab="Time in months",
           break.time.by=40)
surv
summary(surv, time = c(12, 12*3, 12*5))
# 导出设置500*500


# 竞争风险CIF曲线
library(cmprsk)
cuminc_fit <- cuminc(ftime = dat$time, fstatus = dat$status, cencode = 0)
cuminc_fit 
plot(cuminc_fit, 
     xlab = 'Months', 
     ylab = 'CIF', 
     lwd=3, 
     lty = c(1, 3),
     col = mycolors)
grid(col=c('grey85', "grey95"), lty = 1)
box(lwd=1)
# 导出设置500*500


# 腹膜细胞学对OS的影响
surv <- survfit(Surv(time, dead == 1) ~ peri, data = dat)
ggsurvplot(surv,
           fun="pct", palette = mycolors,
           conf.int=T, conf.int.style = "ribbon", conf.int.alpha = 1,
           pval=TRUE,
           risk.table="abs_pct",risk.table.col="black",risk.table.y.text=F,
           risk.table.height=0.25,
           ncensor.plot=F,
           ggtheme=theme_bw(),
           legend.title="Peritoneal cytology",xlab="Time in months",
           break.time.by=40)
surv
summary(surv, time = c(12, 12*3, 12*5))


# 竞争风险CIF曲线绘制自编函数
compete.risk <- function(ftime, fstatus, cencode, group) {
  cuminc_fit <- cuminc(ftime = ftime, fstatus = fstatus, cencode = cencode, group = group)
  #Fine-Gray检验结果。$Test第一行的"pv"即Fine-Gray检验的P值。表示在控制了竞争风险事件
  #（即第二行计算的统计量和P值）后，两组间CSS的统计学差异
  #$est表示估计的各时间点各组的累计癌症发生率与累计竞争风险事件发生率（分别用1和2来区分，
  #与第一行第二行一致）。
  #$var表示估计的各时间点各组的累计癌症发生率与累计竞争风险事件发生率的方差
  print(cuminc_fit)
  P_value <- round(cuminc_fit$Tests[1,"pv"], 3) #提取P值
  P_value <- ifelse(P_value < 0.001, "<0.001", P_value)
  plot(cuminc_fit, xlab = 'Months', ylab = 'CIF', lwd=3, lty = c(1, 1, 3, 3), 
       col = mycolors)
  grid(col=c('grey85', "grey95"), lty = 1)
  text(60, 0.9, pos = 4, paste0("P value of Gray's test: ", P_value), cex = 0.9)
  box(lwd=1)
}

# 腹膜细胞学CIF曲线
compete.risk(dat$time, dat$status, 0, group = dat$peri)

save(list = ls(), file = "sarcoma_peri/data/data_before_match.Rdata")







# 🔴不同FIGO分期下腹膜细胞学对预后的影响 --------------------------------------------------------------
# FIGO I/II
dat$peri_figo <- ifelse(dat$figo == "I" & dat$peri == "Negtive", 
                        "FIGO I + Negtive peritoneal cytology", 
                        ifelse(dat$figo == "I" & dat$peri == "Malignant", 
                               "FIGO I + Malignant peritoneal cytology", 
                               ifelse(dat$figo == "II" & dat$peri == "Negtive", 
                                      "FIGO II + Negtive peritoneal cytology", 
                                      ifelse(dat$figo == "II" & dat$peri == "Malignant", 
                                             "FIGO II + Malignant peritoneal cytology", NA))))
dat$peri_figo <- factor(dat$peri_figo, levels = c("FIGO I + Negtive peritoneal cytology", 
                                                  "FIGO I + Malignant peritoneal cytology", 
                                                  "FIGO II + Negtive peritoneal cytology",
                                                  "FIGO II + Malignant peritoneal cytology"))
summary(dat$peri_figo)


surv <- survfit(Surv(time, dead == 1) ~ peri_figo, data = dat)
ggsurvplot(surv,
           fun="pct", palette = mycolors,
           conf.int=F, 
           pval=TRUE,
           risk.table="abs_pct",risk.table.col="black",risk.table.y.text=F,
           risk.table.height=0.25,
           ncensor.plot=F,
           ggtheme=theme_bw(),
           legend.title="Peritoneal cytology",xlab="Time in months",
           break.time.by=40)
surv
summary(surv, time = c(12, 12*3, 12*5))

# FIGO I/II vs. FIGO III
dat$peri_figo2 <- ifelse((dat$figo == "I" | dat$figo == "II")  & dat$peri == "Negtive", 
                         "FIGO I/II + Negtive peritoneal cytology", 
                         ifelse((dat$figo == "I" | dat$figo == "II") & dat$peri == "Malignant", 
                                "FIGO I/II + Malignant peritoneal cytology", 
                                ifelse(dat$figo == "III" & dat$peri == "Negtive", 
                                       "FIGO III + Negtive peritoneal cytology", 
                                       ifelse(dat$figo == "III" & dat$peri == "Malignant", 
                                              "FIGO III + Malignant peritoneal cytology", NA))))
dat$peri_figo2 <- factor(dat$peri_figo2, levels = c("FIGO I/II + Negtive peritoneal cytology", 
                                                    "FIGO I/II + Malignant peritoneal cytology", 
                                                    "FIGO III + Negtive peritoneal cytology",
                                                    "FIGO III + Malignant peritoneal cytology"))
summary(dat$peri_figo2)



surv <- survfit(Surv(time, dead == 1) ~ peri_figo2, data = dat)
ggsurvplot(surv,
           fun="pct", palette = mycolors,
           conf.int=F, 
           pval=TRUE,
           risk.table="abs_pct",risk.table.col="black",risk.table.y.text=F,
           risk.table.height=0.25,
           ncensor.plot=F,
           ggtheme=theme_bw(),
           legend.title="Peritoneal cytology",xlab="Time in months",
           break.time.by=40)
surv
summary(surv, time = c(12, 12*3, 12*5))










#🔴匹配（多重插补后匹配）-----------------------------------------------------------------
load("sarcoma_peri/data/data_before_match.Rdata")
# 提取用于生存分析的数据集
dat_surv <- dat[,c(vars, c("peri", "dead", "time"))]

# 分析数据集的基本情况
plot_intro(dat_surv, ggtheme = theme_bw())
plot_missing(dat_surv, ggtheme = theme_bw())
profile_missing(dat_surv) %>% arrange(., desc(pct_missing))
aggr(dat_surv,
     prop=F,
     numbers=T,
     sortVars=TRUE,
     gap=0,
     cex.axis = 0.6,
     col = mycolors)

# 排除缺失值数量>3个的个案和N分期未知的个案
dat_surv2 <- dat_surv[rowSums(is.na(dat_surv)) <= 3 & is.na(dat_surv$N_stage) == F, ]
paste0("【排除缺失值数量>3个和N分期未知的个案后】排除：", 
       nrow(dat_surv) - nrow(dat_surv2), 
       "个；剩余：", 
       nrow(dat_surv2), "个") %>% print()

plot_intro(dat_surv2, ggtheme = theme_bw())
plot_missing(dat_surv2, ggtheme = theme_bw())
aggr(dat_surv2,
     prop = T,
     numbers = T,
     sortVars = TRUE,
     gap = 0,
     cex.axis = 0.6,
     col = mycolors)



# 🔴进行缺失值的多重插补

imp_surv <- mice(dat_surv2, m = 10, maxit = 5, seed = 1996, print = FALSE)
imp_surv

#评估插补数据集的分布情况
plot(imp_surv, col = mycolors)
vars_miss_surv <- colnames(dat_surv2[,colSums(is.na(dat_surv2)) > 0]) #提取有缺失值的变量
vars_miss_surv
stripplot(imp_surv, 
          marriage + race + tumor_size + plnd + alnd + lnd + rad + chem ~ .imp, 
          pch = 20,
          col = mycolors)
imp.densityplot(imp_surv, vars_miss_surv, mycolors)


# 🔴多重插补后匹配
# 分析匹配前各变量的均衡性
fun <- paste(vars, collapse  = "+") %>% paste0("peri~", .) %>% as.formula()
library(cobalt)
bal.tab(fun,
        imp_surv, 
        binary = "std", 
        stats = "ks",
        thresholds = c(m = .1),
        imp.fun = 'max') 

#指定需要匹配的变量
vars_psm <- c("year2", "age", "race", "marriage", "income", "grade2", 
              "tumor_size", "his", "figo", "surg", "lnd", "rad", "chem")
#生成PS计算公式
fun <- paste(vars_psm, collapse  = "+") %>% paste0("peri~", .) %>% as.formula()
#匹配
library(MatchThem)
imp_match <- matchthem(fun, 
                       datasets = imp_surv,
                       approach = 'across',
                       distance = "glm",
                       link =   "logit",
                       method = 'full')

#基于cobalt包的均衡性检验
#输出的表格显示了匹配前后的变量均衡性。
#后缀为“.Un”的列表示匹配前的结果；后缀为“.Adj”的列表示匹配后的结果。结果解读：
#Diff.Adj：standardized mean differences (SMD)
#SMD=变量在治疗组和对照组间的平均值差异除以匹配前数据中（默认）的标准差
#V.Ratio.Adj：variance ratios，只有连续型变量会计算此值。越接近1越好，推荐>0.5~0.2
#KS.Adj：Kolmogorov-Smirnov (KS) statistic。
#表示eCDF（经验累积密度函数，empirical cumulative density functions)
#在治疗组和对照组间的最大差异，衡量组间协变量总体分布的差异。范围在0到1之间，值越接近0表示平衡性越好。
#Kolmogorov-Smirnov statistic可作为SMD的补充，用于均衡性评估
tab_bal <- bal.tab(imp_match, 
                   un = TRUE, #是否计算匹配前的均衡性
                   #二元协变量默认输出比例的原始差值。这里为了让二元协变量也输出SMD，需要指定binary = "std"。
                   binary = "std", 
                   addl = vars, #可以添加未参与匹配的变量，以一并检验其均衡性
                   stats = c("m", #SMD
                             "v", #variance ratios
                             "ks"),#Kolmogorov-Smirnov statistics
                   thresholds = c(m = .1),#可以设定各检验的阈值，最后会统计满足和超出该阈值的变量数量
                   #以及差异最大的变量
                   imbalanced.only = F, #只展示不平衡的变量。imbalanced.only = TRUE
                   imp.fun = 'max') 
tab_bal
write_xlsx(tab_bal[["Balance.Across.Imputations"]], 
           "sarcoma_peri/output/baseline_tab_after_match.xlsx", 
           row.names = T)

#Love plot of the SMD
data.frame(old = vars_psm,
           new = c("Year of diagnosis", "Age", "Race", "Marital status", "Family income", 
                   "Grade", "Tumor size", "Histologic type", "FIGO Stage", "Surgery", 
                   "Lymphadenectomy", "Radiotherapy", "Chemotherapy")) %>%
  love.plot(imp_match, 
            stats = c("m", "ks"), 
            abs = TRUE, #是否采用SMD的绝对值绘图
            drop.distance = F, #是否在顶部显示distance即PS或weight在匹配前后的情况
            thresholds = c(m = .1),
            line = F, #是否绘制各变量的连线
            var.order = "unadjusted", 
            var.names = .,
            sample.names = c("Unmatched", "Matched"),
            binary = "std",
            stars = "raw", #"raw": X轴表示"Standardized Mean Differences"，
            #标星号的变量代表其X轴代表的实际是"Mean Differences"
            shapes = c("circle filled", "circle"), 
            colors = mycolors[1:2],
            position = "top", #图例的位置
            grid = F)
# 导出设置1000*630

#展示PS在匹配前后的密度图
bal.plot(imp_match, 
         var.name = "distance", 
         which = "both", 
         which.imp = .none, 
         colors = mycolors[1:2], 
         disp.means=T, 
         sample.names = c("Unmatched", "Matched"))
#单个变量的均衡性检查。连续型变量默认展示kernel density plots；因子变量默认展示条形图
for (i in 1:length(vars_psm)) {
  bal.plot(imp_match, 
           var.name = vars_psm[i],
           which = "both", 
           which.imp = .none, 
           colors = mycolors[1:2], 
           disp.means=T, #是否展示连续型变量的均数
           sample.names = c("Unmatched", "Matched")) %>%
    print()
}
# 导出尺寸530*280




#🔴匹配后生存分析-------------------------------------------------------------------------
surv <- survfit(Surv(time, dead == 1) ~ peri, 
                data = dat_surv2, 
                weights = imp_match[["models"]][[1]][["weights"]], 
                robust = TRUE, #To request cluster-robust standard errors 
                cluster = imp_match[["models"]][[1]][["subclass"]])
ggsurvplot(surv,
           fun="pct", 
           palette = mycolors,
           conf.int=T, 
           conf.int.style = "ribbon",
           conf.int.alpha = 1,
           pval=F,
           risk.table="abs_pct",
           risk.table.col="black",
           risk.table.y.text=F,
           risk.table.height=0.25,
           ncensor.plot=F,
           ggtheme=theme_bw(),
           legend.title="Histology",
           xlab="Time in months",
           break.time.by=40)
surv
summary(surv, time = c(12, 12*3, 12*5))

library(RISCA)
ipw.log.rank(times = dat_surv2$time,
             failures = as.numeric(dat_surv2$dead) - 1,
             variable = as.numeric(dat_surv2$peri) - 1, 
             weights = imp_match[["models"]][[1]][["weights"]])


# 筛选匹配后SMD<0.1的变量
tab_bal[["Balance.Across.Imputations"]] %>% .[.$Max.Diff.Adj > 0.1, ] %>% rownames()
# 计算校正不平衡变量后的P值
adj.p <- with(imp_match, coxph(Surv(time, dead == 1) ~ peri + year2 + T_stage + alnd,
                               robust = TRUE, #To request cluster-robust SEs
                               cluster = imp_match[["models"]][[1]][["subclass"]], 
                               weights = imp_match[["models"]][[1]][["weights"]])) %>%
  pool() %>% summary() %>% .$p.value %>% .[1]
adj.p <- ifelse(adj.p < 0.001, "< 0.001", round(adj.p, 3))
adj.p




# 🔴基于匹配后多重插补数据的批量单因素Cox回归分析
if(F){
  vars_cox <- c("year2", "age", "race", "marriage", "income", "grade2", "tumor_size", "his", 
                "T_stage","N_stage", "M_stage", "figo", "peri", "surg", "plnd", "alnd", "lnd", 
                "rad", "chem")
  length(vars_cox)
  batch.poolcox_result <- c()
  pool_cox <- with(imp_match, coxph(formula(paste0("Surv(time, dead == 1) ~", vars_cox[19])), 
                                    robust = TRUE, #To request cluster-robust SEs
                                    cluster = imp_match[["models"]][[1]][["subclass"]], 
                                    weights = imp_match[["models"]][[1]][["weights"]])) %>% 
    pool()
  {
    sum_pool_cox <- summary(pool_cox, conf.int = TRUE)
    sum_pool_hr <- summary(pool_cox, conf.int = TRUE, exponentiate = TRUE)
    Variables <- as.character(sum_pool_cox$term)
    Coefficient <- round(sum_pool_cox$estimate, 2)
    HR <- round(sum_pool_hr$estimate, 2)
    LCI <- round(sum_pool_hr$`2.5 %`, 2)
    UCI <- round(sum_pool_hr$`97.5 %`, 2)
    HR_95CI <- paste0(HR, " (", LCI, ", ", UCI, ")")
    P_value <- round(sum_pool_hr$p.value, 3)
    single_result <- data.frame(Variables, Coefficient, HR, LCI, UCI, HR_95CI, P_value)
    batch.poolcox_result <- rbind(batch.poolcox_result, single_result)
    batch.poolcox_result
    }
  # 标注P<0.05的变量
  batch.poolcox_result$sig <- ifelse(batch.poolcox_result$P_value < 0.05, "*", "")
  # 在前面添加一列varnames表示变量名称，对应原来表格的多个哑变量，便于后面根据P值筛选变量
  varnames <- c()
  for (i in vars_cox) {
    if (class(dat_surv2[,i]) == "factor") { 
      varnames <- rep(i, length(levels(dat_surv2[,i])) - 1) %>%
        c(varnames, .)
    } else {
      varnames <- c(varnames, i)
    }
  }
  batch.poolcox_result <- data.frame(varnames, batch.poolcox_result)
  
  save(batch.poolcox_result, file = "sarcoma_peri/output/batch.poolcox_result.Rdata")
}



# 🔴匹配后多因素cox回归分析
load("sarcoma_peri/output/batch.poolcox_result.Rdata")
# 筛选单因素cox回归P<0.05的变量
vars_multicox <- batch.poolcox_result[batch.poolcox_result$P_value < 0.05, "varnames"] %>%
  unique()
vars_multicox
# 去除存在共线性的变量
vars_multicox <- c("age", "grade2", "tumor_size", "his", "figo", "peri", "lnd", "chem")

# Cox Regression with a robust variance estimator
{
  pool_multicox <- with(imp_match, 
                        coxph(as.formula(paste0("Surv(time, dead == 1)", "~", 
                                                paste(vars_multicox, collapse = "+"))),
                              robust = TRUE, #To request cluster-robust SEs
                              cluster = imp_match[["models"]][[1]][["subclass"]], 
                              weights = imp_match[["models"]][[1]][["weights"]])) %>%
    pool()
  sum_pool_cox <- summary(pool_multicox, conf.int = TRUE)
  sum_pool_hr <- summary(pool_multicox, conf.int = TRUE, exponentiate = TRUE)
  Variables <- as.character(sum_pool_cox$term)
  Coefficient <- round(sum_pool_cox$estimate, 2)
  HR <- round(sum_pool_hr$estimate, 2)
  LCI <- round(sum_pool_hr$`2.5 %`, 2)
  UCI <- round(sum_pool_hr$`97.5 %`, 2)
  HR_95CI <- paste0(HR, " (", LCI, ", ", UCI, ")")
  P_value <- ifelse(sum_pool_hr$p.value < 0.001, "< 0.001", 
                    round(sum_pool_hr$p.value, 3))
  pool_multicox_result <- data.frame(Variables, Coefficient, HR, LCI, UCI, HR_95CI, P_value)
  pool_multicox_result
}
# 标注P<0.05的变量
pool_multicox_result$sig <- ifelse(pool_multicox_result$P_value < 0.05, "*", "")

# 合并单因素和多因素cox回归结果
# 添加一列序号，便于merge合并后恢复原来的顺序
batch.poolcox_result$n <- c(1:nrow(batch.poolcox_result)) 
cox_result <- merge(x = batch.poolcox_result,
                    y = pool_multicox_result,#x、y为要合并的数据框或者对象
                    by ="Variables", 
                    all = T)
# 按照原来的顺序重新排序
cox_result <-  arrange(cox_result, n)
# 去掉添加的“varnames”列、“n”列以及单独列出的HR列和95%置信区间列
cox_result <- -c(which(colnames(cox_result) == "varnames"), 
                 which(colnames(cox_result) == "n"),
                 which(colnames(cox_result) == "HR.x"),
                 which(colnames(cox_result) == "LCI.x"),
                 which(colnames(cox_result) == "UCI.x"),
                 which(colnames(cox_result) == "HR.y"),
                 which(colnames(cox_result) == "LCI.y"),
                 which(colnames(cox_result) == "UCI.y")) %>%
  cox_result[,.]
rownames(cox_result) <- cox_result[,1]
cox_result <- cox_result[,-1]
colnames(cox_result) <- c("Univariate analysis", rep("", 3), 
                          "Multivariate analysis", rep("", 3))
cox_result <- rep(c("Coef", "HR (95%CI)", "P value", "Sig"), 2) %>%
  rbind(., cox_result)

write_xlsx(cox_result, "sarcoma_peri/output/cox_result.xlsx", row.names = T)
save(list = ls(), file = "sarcoma_peri/data/data_for_RSF.Rdata")







# 🔴随机生存森林分析 -------------------------------------------------------------------
load("sarcoma_peri/data/data_for_RSF.Rdata")

# 定义随机森林模型的方程
vars_cox <- c("year2", "age", "race", "marriage", "income", "grade2", "tumor_size", "his", 
              "T_stage", "N_stage", "M_stage", "figo", "peri", "surg", "plnd", "alnd", "lnd", 
              "rad", "chem")
fun_forest <- paste(vars_cox, collapse = " + ") %>% 
  paste0("Surv(time, dead) ~ ", .) %>% as.formula()

# 通过grid search进行随机森林调参，寻找最优nodesize和mtry参数组合
{
  library(randomForestSRC)
  set.seed(110)
  tune_forest <- tune(fun_forest, 
                      data = dat_surv2, 
                      ntreeTry = 500, #进行调优的树的数量 500
                      mtryStart = 1, #mtry迭代的起始值
                      doBest = T,
                      trace = T) # 是否显示每次迭代过程
  library(beepr)
  beep(sound = "ping") #脚本运行结束后播放提示声
  }

## the optimized forest 
print(tune_forest$rf)

## visualize the nodesize/mtry OOB surface
if (library("interp", logical.return = TRUE)) {
  ## nice little wrapper for plotting results
  plot.tune <- function(tune_forest, linear = TRUE) {
    x <- tune_forest$results[,1] # 获取grid search后的nodesize向量
    y <- tune_forest$results[,2] # 获取grid search后的mtry向量
    z <- tune_forest$results[,3] # 获取grid search后的OBB error向量
    so <- interp(x = x, 
                 y = y, 
                 z = z, 
                 linear = linear)
    best.nodesize <- x[which.min(z)] # 最小OBB error对应的最优nodesize
    best.mtry <- y[which.min(z)] # 最小OBB error对应的y，即最优mtry
    filled.contour(x = so$x,
                   y = so$y,
                   z = so$z,
                   xlim = range(so$x, finite = TRUE) + c(-2, 2),
                   ylim = range(so$y, finite = TRUE) + c(-2, 2),
                   color.palette = colorRampPalette(c("lightgreen", "darkblue")),
                   xlab = "Minimum terminal node size",
                   ylab = "Mtry (No. of variables tried at each split)",
                   main = paste0("Optimal nodesize: ", best.nodesize, 
                                 "; Optimal mtry: ", best.mtry, 
                                 "; Minimal OOS error: ", round(min(z), 2)),
                   cex.main = 0.8,
                   key.title = title(main = "OOS error", cex.main = 1),
                   plot.axes = {
                     axis(1); # 绘制水平坐标轴
                     axis(2); # 绘制垂直坐标轴
                     # 绘制grid search中的每个参数点
                     points(x, y, 
                            pch = 16, cex = .3, col = "white")
                     # 标注最优参数点
                     points(x = best.nodesize, 
                            y = best.mtry, 
                            pch = "x", cex = 1, font = 2)
                   })
  }
  ## plot the surface
  plot.tune(tune_forest)
}
# 导出大小：540*480


# 提取最优nodesize和mtry参数
grid.search_result <- data.frame(tune_forest$results)
grid.search_result <- grid.search_result[order(grid.search_result$err),]
head(grid.search_result)
best.nodesize <- grid.search_result$nodesize[1]
best.nodesize
best.mtry <- grid.search_result$mtry[1]
best.mtry



# 根据最优参数组合构建随机生存森林模型
{
set.seed(110)
rfsrc <- rfsrc(fun_forest, 
               data = dat_surv2, 
               ntree = 3000,
               nodesize = best.nodesize,
               mtry = best.mtry,
               block.size = 1,
               bootstrap = "by.root",
               samptype = "swor", 
               importance = "permute",
               na.action = "na.impute", 
               nimpute = 5)
beep(sound = "ping")
rfsrc
}


# 绘制不同ntrees对应的累积OOB error rate的折线图
rfsrc_data <- rfsrc$err.rate[,"event.2"]
rfsrc_data <- data.frame(ntrees = c(1:length(rfsrc_data)), err.rate = rfsrc_data)
head(rfsrc_data)

ggplot(data = rfsrc_data,
       aes(x = ntrees, y = err.rate)) +
  geom_line(lwd = 0.5)+
  labs(x="Number of trees", y = "Tree cumulative OOS error rate") +
  theme_bw()
# 导出设置：500*480




# 绘制变量重要性图
library(ggRandomForests)
importance <- gg_vimp(rfsrc)
importance <- subset(importance, set == "event.2")
ggplot(importance,
       aes(x = vimp, 
           y = reorder(vars, vimp))) + 
  geom_bar(stat = "identity",   
           show.legend = T,   
           width = .9,   
           aes(fill = vimp)) +  #设置根据“coef”的大小来进行渐变填充
  scale_fill_gradient2(low = mycolors[1],  #设置条形图的渐变颜色
                       mid = mycolors[6],
                       high = mycolors[9])+
  geom_text(aes(label = sprintf("%0.3f", vimp),
                x= -0.001), 
            colour = "black",
            size = 3)+
  xlab("Variable Importance (VIMP)") +  
  ylab("Variables")+
  theme_bw()+ 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())
# 导出大小：800*440

save(list = ls(), file = "sarcoma_peri/data/data_before_subgroup_analysis.Rdata")






# 🔴亚组cox回归分析----------------------------------------------------------------------------------
load("sarcoma_peri/data/data_before_subgroup_analysis.Rdata")

#转换age为分类变量
set.seed(1996)
tree_age <- rpart(Surv(time, dead == 1) ~ age, data = dat)
tree_age$cptable
# find the complexity parameter correlated with the smallest cross-validation error.
# and use this complexity parameter to prune tree
prune_age <- tree_age$cptable[which.min(tree_age$cptable[,"xerror"]),"CP"] %>%
  prune(tree_age, .)
#决策树绘制
rpart.plot(prune_age,
           type = 2,
           extra = 1+100,
           box.palette = mycolors[1:4])
plot(as.party(prune_age))#展示各组生存曲线
# 导出大小：500*400

dat$age2 <- ifelse(dat$age < 53, "18-52", "≥53")
dat$age2 <- factor(dat$age2, levels = c("18-52", "≥53"))
summary(dat$age2)

dat_surv2$age2 <- ifelse(dat_surv2$age < 53, "18-52", "≥53")
dat_surv2$age2 <- factor(dat_surv2$age2, levels = c("18-52", "≥53"))
summary(dat_surv2$age2)

surv <- survfit(Surv(time, dead == 1) ~ age2, data = dat)
summary(surv, time = 12*5)





#转换tumor_size为分类变量
tree_tumor_size <- rpart(Surv(time, dead == 1) ~ tumor_size, data = dat)
tree_tumor_size$cptable
prune_tumor_size <- tree_tumor_size$cptable[which.min(tree_tumor_size$cptable[,"xerror"]),
                                            "CP"] %>%
  prune(tree_tumor_size, .)
rpart.plot(prune_tumor_size,
           type = 2,
           extra = 1+100,
           box.palette = mycolors[1:4])
plot(as.party(prune_tumor_size))
# 导出大小：500*400

dat$tumor_size2 <- ifelse(dat$tumor_size < 81, "<81",
                          ifelse(dat$tumor_size < 158, "81-157", "≥158"))
dat$tumor_size2 <- factor(dat$tumor_size2, levels = c("<81", "81-157", "≥158"))
summary(dat$tumor_size2)

dat_surv2$tumor_size2 <- ifelse(dat_surv2$tumor_size < 81, "<81", 
                                ifelse(dat_surv2$tumor_size < 158, "81-157", "≥158"))
dat_surv2$tumor_size2 <- factor(dat_surv2$tumor_size2, levels = c("<81", "81-157", "≥158"))
summary(dat_surv2$tumor_size2)

surv <- survfit(Surv(time, dead == 1) ~ tumor_size2, data = dat)
summary(surv, time = 12*5)


# 转换T分期为2分类变量
summary(dat$T_stage)
dat$T_stage2 <- ifelse(dat$T_stage == "T1" | dat$T_stage == "T2", "T1/2", "T3/4")
dat$T_stage2 <- factor(dat$T_stage2)
summary(dat$T_stage2)

dat_surv2$T_stage2 <- ifelse(dat_surv2$T_stage == "T1" | dat_surv2$T_stage == "T2", 
                             "T1/2", "T3/4")
dat_surv2$T_stage2 <- factor(dat_surv2$T_stage2)
summary(dat_surv2$T_stage2)

# 转换FIGO分期为二分类变量
summary(dat$figo)
dat$figo2 <- ifelse(dat$figo == "I" | dat$figo == "II", "I/II", "III/IV")
dat$figo2 <- factor(dat$figo2)
summary(dat$figo2)

dat_surv2$figo2 <- ifelse(dat_surv2$figo == "I" | dat_surv2$figo == "II", "I/II", "III/IV")
dat_surv2$figo2 <- factor(dat_surv2$figo2)
summary(dat_surv2$figo2)

# 只保留手术方式中的TH+BSO和RH/EH两个水平
summary(dat$surg)
dat$surg2 <- ifelse(dat$surg == "TH+BSO", "TH+BSO", 
                    ifelse(dat$surg == "RH/EH", "RH/EH", NA))
dat$surg2 <- factor(dat$surg2, levels = c("TH+BSO", "RH/EH"))
summary(dat$surg2)

dat_surv2$surg2 <- ifelse(dat_surv2$surg == "TH+BSO", "TH+BSO", 
                          ifelse(dat_surv2$surg == "RH/EH", "RH/EH", NA))
dat_surv2$surg2 <- factor(dat_surv2$surg2, levels = c("TH+BSO", "RH/EH"))
summary(dat_surv2$surg2)





# 定义亚组
groups <- c("year2", "age2", "race", "marriage", "income", "grade2", "tumor_size2", "his", 
            "T_stage2","N_stage", "M_stage", "figo2", "surg2", "plnd", "alnd", "lnd", 
            "rad", "chem")

# 批量亚组cox回归分析
subcox_result <- c()
for (i in 1:length(groups)) {
  factor.levels <- levels(dat[,groups[i]])
  sub_result <- c()
  for (j in 1:length(factor.levels)) {
    fit <- coxph(Surv(time, dead == 1) ~ peri, 
                 data = dat, 
                 subset = dat[,groups[i]] == factor.levels[j])
    #提取cox回归结果参数
    cox_list <- summary(fit)
    Subgroups <- factor.levels[j]
    Coefficient <- round(cox_list$coefficients[,"coef"], 2)
    HR <- round(cox_list$conf.int[,"exp(coef)"], 2)
    LCI <- round(cox_list$conf.int[,"lower .95"], 2)
    UCI <- round(cox_list$conf.int[,"upper .95"], 2)
    HR_95_CI <- paste0(HR, " (", LCI, ", ", UCI, ")")
    P_value <- ifelse(cox_list$coefficients[,"Pr(>|z|)"] < 0.001, "<0.001", 
                      round(cox_list$coefficients[,"Pr(>|z|)"], 3))
    #汇总单个亚组的cox回归结果
    single_result <- data.frame(Subgroups, Coefficient, HR, LCI, UCI, HR_95_CI, P_value)
    #合并单个变量所有亚组的cox回归结果
    sub_result <- rbind(sub_result, single_result)
  }
  groups_name <- data.frame(Subgroups = "", 
                            Coefficient = "", 
                            HR = "",
                            LCI = "",
                            UCI = "", 
                            HR_95_CI = "", 
                            P_value = "")
  groups_name[1] <- groups[i]
  subcox_result <- rbind(subcox_result, groups_name, sub_result)
}
head(subcox_result)
#展示P<0.05的亚组
subcox_result[subcox_result$P_value < 0.05, ]

write_xlsx(subcox_result, "sarcoma_peri/output/subcox_result.xlsx", row.names = FALSE)



# 🔴亚组分析森林图
forest_data <- subcox_result
forest_data[,c("HR", "LCI", "UCI")] <- lapply(forest_data[,c("HR", "LCI", "UCI")], 
                                              as.numeric)

library(meta)
meta <- metagen(log(HR), lower = log(LCI), upper = log(UCI), pval = P_value,
                data = forest_data, 
                sm = "HR", 
                studlab = Subgroups, 
                comb.fixed = T,
                comb.random = T)
meta
forest(meta, 
       layout = "JAMA",
       common = F,
       random = T,
       col.study = mycolors[1],
       col.square = mycolors[2],
       col.square.lines = NA,
       col.random = "red",
       header.line = T)

library(forestplot)
forestplot(forest_data,
           labeltext = c(Subgroups, Coefficient, HR_95_CI, P_value), 
           align ="llll",#设置每列文字的对齐方式
           mean = HR, #定义HR均值列
           lower = LCI, #定义HR值的95%CI下限列
           upper = UCI, #定义HR值的5%CI上限列
           graph.pos = 4, #设置森林图出现的列
           graphwidth = unit(30,'mm'), #森林图的宽度
           colgap = unit(3,'mm'), #设置图形中的列间距
           fn.ci_norm = fpDrawNormalCI, #中央HR均值点的形状
           boxsize = 0.4, #中央HR均值点的大小
           lwd.ci = 1.5, #设置95%CI线的粗细
           ci.vertices = TRUE, #添加95%CI线两端的小竖线
           zero = 1, #设置无效线
           lwd.zero = 1, #设置无效线的粗细
           col=fpColors(box = mycolors[1], #OR均值点的颜色
                        lines = mycolors[2], #95%CI线的颜色
                        summary = mycolors[9], #合并HR值的颜色
                        zero = "black"), #无效线的颜色
           xlog = T, #转换为对数坐标轴
           is.summary = forest_data$HR_95_CI == "", #设置为TRUE的行文字以粗体出现
           txt_gp = fpTxtGp(label = gpar(cex = 0.7), #表格主体文字的大小
                            ticks = gpar(cex = 0.7), #森林图下方的坐标轴的刻度文字大小
                            xlab = gpar(cex = 0.7))) %>% #森林图下方X轴标签文字的大小
  #添加表头
  fp_add_header(Subgroups = "Subgroups", 
                Coefficient = "Coef", 
                HR_95_CI = "HR (95%CI)", 
                P_value = "P value")  %>%
  #添加一行，展示合并效应量
  fp_append_row(Subgroups = "Overall", 
                Coefficient = "", 
                HR_95_CI = "3.48 (3.15, 3.83)", 
                mean  = 3.48,
                lower = 3.15,
                upper = 3.83,
                P_value = "<0.001",
                is.summary = T) %>%
  # 添加横线，制作成三线表（最后横线添加到nrow(forest_data) + 3行，这里即h_63）
  fp_add_lines(h_1 = gpar(lwd = 2, col = "black"),
               h_2 = gpar(lwd = 1, col = "black"),
               h_63 = gpar(lwd = 2, col = "black")) %>% 
  fp_set_zebra_style("grey95") # 添加间隔填充
# 导出大小: 500*1000






# 🔴基于匹配后数据的批量亚组cox回归分析-----------------------------------------------------
dat_surv2$cluster <- imp_match[["models"]][[1]][["subclass"]]
dat_surv2$weights <- imp_match[["models"]][[1]][["weights"]]

subcox_result_psm <- c()
for (i in 1:length(groups)) {
  factor.levels <- levels(dat_surv2[,groups[i]])
  sub_result <- c()
  for (j in 1:length(factor.levels)) {
    subdata <- subset(dat_surv2, dat_surv2[,groups[i]] == factor.levels[j])
    fit <- coxph(Surv(time, dead == 1) ~ peri, 
                 data = subdata, 
                 robust = TRUE, #To request cluster-robust SEs
                 cluster = subdata$cluster,
                 weights = subdata$weights)
    #提取cox回归结果参数
    cox_list <- summary(fit)
    Subgroups <- factor.levels[j]
    Coefficient <- round(cox_list$coefficients[,"coef"], 2)
    HR <- round(cox_list$conf.int[,"exp(coef)"], 2)
    LCI <- round(cox_list$conf.int[,"lower .95"], 2)
    UCI <- round(cox_list$conf.int[,"upper .95"], 2)
    HR_95_CI <- paste0(HR, " (", LCI, ", ", UCI, ")")
    P_value <- ifelse(cox_list$coefficients[,"Pr(>|z|)"] < 0.001, "<0.001", 
                      round(cox_list$coefficients[,"Pr(>|z|)"], 3))
    #汇总单个亚组的cox回归结果
    single_result <- data.frame(Subgroups, Coefficient, HR, LCI, UCI, HR_95_CI, P_value)
    #合并单个变量所有亚组的cox回归结果
    sub_result <- rbind(sub_result, single_result)
  }
  groups_name <- data.frame(Subgroups = "", 
                            Coefficient = "", 
                            HR = "",
                            LCI = "",
                            UCI = "", 
                            HR_95_CI = "", 
                            P_value = "")
  groups_name[1] <- groups[i]
  subcox_result_psm <- rbind(subcox_result_psm, groups_name, sub_result)
}
head(subcox_result_psm)
#展示P<0.05的亚组
subcox_result_psm[subcox_result_psm$P_value < 0.05, ]

write_xlsx(subcox_result_psm, "sarcoma_peri/output/subcox_result_psm.xlsx", row.names = FALSE)



# 亚组分析森林图
forest_data_psm <- subcox_result_psm
forest_data_psm[,c("HR", "LCI", "UCI")] <- lapply(forest_data_psm[,c("HR", "LCI", "UCI")], 
                                                  as.numeric)

meta <- metagen(log(HR), lower = log(LCI), upper = log(UCI), pval = P_value,
                data = forest_data_psm, 
                sm = "HR", 
                studlab = Subgroups, 
                comb.fixed = T,
                comb.random = T)
meta
forest(meta, 
       layout = "JAMA",
       common = F,
       random = T,
       col.study = mycolors[1],
       col.square = mycolors[2],
       col.square.lines = NA,
       col.random = "red",
       header.line = T)


forestplot(forest_data_psm,
           labeltext = c(Subgroups, Coefficient, HR_95_CI, P_value), 
           align ="llll",#设置每列文字的对齐方式
           mean = HR, #定义HR均值列
           lower = LCI, #定义HR值的95%CI下限列
           upper = UCI, #定义HR值的5%CI上限列
           graph.pos = 4, #设置森林图出现的列
           graphwidth = unit(30,'mm'), #森林图的宽度
           colgap = unit(3,'mm'), #设置图形中的列间距
           fn.ci_norm = fpDrawNormalCI, #中央HR均值点的形状
           boxsize = 0.4, #中央HR均值点的大小
           lwd.ci = 1.5, #设置95%CI线的粗细
           ci.vertices = TRUE, #添加95%CI线两端的小竖线
           zero = 1, #设置无效线
           lwd.zero = 1, #设置无效线的粗细
           col=fpColors(box = mycolors[1], #OR均值点的颜色
                        lines = mycolors[2], #95%CI线的颜色
                        summary = mycolors[9], #合并HR值的颜色
                        zero = "black"), #无效线的颜色
           xlog = T, #转换为对数坐标轴
           is.summary = forest_data_psm$HR_95_CI == "", #设置为TRUE的行文字以粗体出现
           txt_gp = fpTxtGp(label = gpar(cex = 0.7), #表格主体文字的大小
                            ticks = gpar(cex = 0.7), #森林图下方的坐标轴的刻度文字大小
                            xlab = gpar(cex = 0.7))) %>% #森林图下方X轴标签文字的大小
  #添加表头
  fp_add_header(Subgroups = "Subgroups", 
                Coefficient = "Coef", 
                HR_95_CI = "HR (95%CI)", 
                P_value = "P value")  %>%
  #添加一行，展示合并效应量
  fp_append_row(Subgroups = "Overall", 
                Coefficient = "", 
                HR_95_CI = "2.05 (1.88, 2.24)", 
                mean  = 2.05,
                lower = 1.88,
                upper = 2.24,
                P_value = "<0.001",
                is.summary = T) %>%
  # 添加横线，制作成三线表（最后横线添加到nrow(forest_data_psm) + 3行，这里即h_63）
  fp_add_lines(h_1 = gpar(lwd = 2, col = "black"),
               h_2 = gpar(lwd = 1, col = "black"),
               h_63 = gpar(lwd = 2, col = "black")) %>% 
  fp_set_zebra_style("grey95") # 添加间隔填充
# 导出大小: 500*1000









sessionInfo()









