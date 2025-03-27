# 🔴首先分析删除缺失值后对临床特征的影响-----------------------------------------
complete_data <- read.csv("sarcoma_rad/sarcoma.csv", header = TRUE, stringsAsFactors = F)
# complete_data中的include列记录了患者是否被排除。“excluded”指患者因标准标准中除
# 非缺失值以外的原因被排除；而“0”代表因缺失值而被排除的患者；“1”代表最终纳入的患者
complete_data <- subset(complete_data, include == 0 | include == 1)
complete_data <- complete_data[c(1:2, 4, 6, 8, 10, 12, 15, 18, 20, 22, 24, 26, 28, 30, 34:36)]
dput(names(complete_data))

library(VIM)
# 展示缺失值的比例
aggr(complete_data,
  prop = T, numbers = T,
  sortVars = TRUE,
  gap = 2, ylab = c("Histogram of missing data", "Pattern")
)
# 展示缺失值的数量
aggr(complete_data,
  prop = F, numbers = T,
  sortVars = TRUE,
  gap = 2, ylab = c("Histogram of missing data", "Pattern")
)


# 转换因子变量
complete_data$year <- ifelse(complete_data$year <= 2007, 1,
  ifelse(complete_data$year <= 2011, 2, 3)
)
complete_data$year <- factor(complete_data$year)
complete_data$grade <- ifelse(complete_data$grade <= 2, 1, 2)
complete_data$rad <- ifelse(complete_data$radiotherapy == 0, 0, 1)
complete_data$rad <- factor(complete_data$rad)

complete_data[4:6] <- lapply(complete_data[4:6], factor)
complete_data[8:15] <- lapply(complete_data[8:15], factor)
complete_data$tumor_size <- as.numeric(complete_data$tumor_size)



# 🔴比较缺失值和完整数据
myvars1 <- c(
  "age", "year", "race", "single", "grade", "tumor_size", "his", "T_stage",
  "N_stage", "surgery", "lymphadenectomy", "chemotherapy", "rad"
)
nonvar <- c("age", "tumor_size") # 指定偏态分布的变量
library(tableone)
misscompare_data <- CreateTableOne(
  data = complete_data,
  vars = myvars1, # 指定需要汇总的变量
  strata = "include", # 分组依据
  includeNA = F,
  test = T,
  addOverall = T
) # 是否添加汇总项
misscompare_tab <- print(misscompare_data,
  showAllLevels = T, # 完整显示二分类变量所有水平的数据
  nonnormal = nonvar, # 指定非正态分布的变量
  smd = F, # 展示标准化均数差SMD，SMD大于10%，这通常被认为是变量不平衡
  quote = F, noSpaces = T
)
# write.csv(misscompare_tab, file = "sarcoma_rad/misscompare.csv")



# 🔴正式分析数据-----------------------------------------------------------------
mydata <- subset(complete_data, include == 1)
dput(names(mydata))

aggr(mydata, prop = FALSE, numbers = TRUE)
aggr(mydata, prop = TRUE, numbers = TRUE)

mydata$grade <- as.numeric(mydata$grade)
grade <- as.data.frame(mydata[, 6])
colnames(grade)[1] <- "grade"
grade[is.na(grade) == TRUE] <- 3
mydata$grade <- NULL
mydata <- cbind(mydata, grade)

# 转换因子变量
mydata$age <- ifelse(mydata$age < 49, 1, ifelse(mydata$age < 58, 2, 3))
mydata$age <- factor(mydata$age)
mydata$tumor_size <- ifelse(mydata$tumor_size < 70, 1, ifelse(mydata$tumor_size < 165, 2, 3))
mydata[is.na(mydata) == TRUE] <- 4
mydata$tumor_size <- factor(mydata$tumor_size)
mydata$grade <- as.factor(mydata$grade)



# 🔴绘制基线特征表---------------------------------------------------------------
myvars2 <- c(
  "age", "year", "race", "single", "grade", "tumor_size", "his", "T_stage",
  "N_stage", "surgery", "lymphadenectomy", "chemotherapy"
)

basedata <- CreateTableOne(
  data = mydata,
  vars = myvars2,
  strata = "rad",
  test = T,
  addOverall = T
)
basetab <- print(basedata,
  showAllLevels = T,
  smd = T,
  quote = F, noSpaces = T
)
# write.csv(basetab, file = "sarcoma_rad/baseline_table.csv")



# 🔴sTPTW匹配--------------------------------------------------------------------
psmodel <- glm(
  rad ~ age + year + race + single + grade + tumor_size + his + T_stage + N_stage + surgery +
    lymphadenectomy + chemotherapy,
  data = mydata,
  family = binomial(link = "logit")
)
mydata$ps <- predict(psmodel, type = "response") # 得到倾向性评分PS

# 计算Pt，Pt=治疗组人数/总人数
table(mydata$rad) # rad=1为治疗组
pt <- 586 / (586 + 2285)
mydata$w <- ifelse(mydata$rad == 1, pt / mydata$ps, (1 - pt) / (1 - mydata$ps))

library(survey)
iptw <- svydesign(ids = ~0, data = mydata, weights = ~w)
unmatcheddata <- CreateTableOne(
  data = mydata,
  vars = myvars2,
  strata = "rad",
  addOverall = F, test = T
)
unmatchedtab <- print(unmatcheddata,
  showAllLevels = TRUE,
  smd = TRUE,
  quote = F, noSpaces = T
)

matcheddata <- svyCreateTableOne(
  data = iptw,
  vars = myvars2,
  strata = "rad",
  addOverall = F, test = T
)
matchedtab <- print(matcheddata,
  showAllLevels = TRUE,
  smd = TRUE,
  quote = F, noSpaces = T
)


# 输出匹配前后的基线特征表-------------------------------------------------------
table_psm <- cbind(unmatchedtab, matchedtab)

table_psm <- rbind(
  Group = rep(c(
    "Level", "No radiotherapy", "Radiotherapy", "P value",
    "test method", "SMD"
  ), 2),
  table_psm
) # 插入一行分组

colnames(table_psm) <- c(
  "Level", "Unmatched", NA, NA, NA, NA, "Level", "sIPTW",
  NA, NA, NA, NA
) # 更改列名

print(table_psm, quote = FALSE)
# write.csv(table_psm, file = "sarcoma_rad/baseline_table_sIPTW.csv")


# 绘制匹配前后的SMD图------------------------------------------------------------
dataPlot <- data.frame(
  Variables = rownames(ExtractSmd(unmatcheddata)),
  Unmatched = as.numeric(ExtractSmd(unmatcheddata)),
  sIPTW = as.numeric(ExtractSmd(matcheddata))
)
library(reshape2)
dataPlotMelt <- melt(
  data = dataPlot,
  id.vars = c("Variables"),
  variable.name = "Method",
  value.name = "SMD"
)
varNames <- as.character(dataPlot$Variables)[order(dataPlot$Unmatched)]
dataPlotMelt$Variables <- factor(dataPlotMelt$Variables,
  levels = varNames
)
library(ggplot2)
ggplot(
  data = dataPlotMelt,
  mapping = aes(
    x = Variables, y = SMD,
    group = Method,
    color = Method,
    shape = Method
  )
) +
  geom_line() +
  geom_point(size = 3) +
  geom_hline(
    yintercept = 0.1,
    color = "red",
    lty = 2,
    size = 0.8
  ) +
  coord_flip() +
  theme_bw(base_size = 18)





# 🔴初步生存分析-----------------------------------------------------------------
# 总人群生存分析
library(survival)
surv.all <- survfit(Surv(time, dead == 1) ~ 1, data = mydata)

# 绘制生存曲线
library(survminer) # 包含ggsurvplot()函数
ggsurvplot(surv.all,
  fun = "pct",
  conf.int = F,
  risk.table = "abs_pct", risk.table.col = "black", risk.table.y.text = F,
  risk.table.height = 0.2,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  xlab = "Time in months",
  break.time.by = 40
)
surv.all # 展示样本量、结局时间发生数、中位生存时间
summary(surv.all, time = c(12, 12 * 3, 12 * 5)) # 查看3年、5年总生存率




# 不同组织学类型的生存分析
surv.byhis <- survfit(Surv(time, dead == 1) ~ his, data = mydata)
ggsurvplot(surv.byhis,
  fun = "pct",
  conf.int = F,
  pval = TRUE,
  risk.table = "abs_pct", risk.table.col = "black", risk.table.y.text = F,
  risk.table.height = 0.25,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Histology", xlab = "Time in months",
  break.time.by = 40
)
surv.byhis
summary(surv.byhis, time = c(12, 12 * 3, 12 * 5))






# 放疗对总人群预后的影响
surv.byrad <- survfit(Surv(time, dead == 1) ~ rad, data = mydata)
ggsurvplot(surv.byrad,
  fun = "pct", # fun="pct"，绘制OS曲线；"event"绘制累积事件曲线
  size = 1, # 生存曲线的粗细
  censor.size = 4, # 删失形状的大小
  conf.int = T, # 添加置信区间
  conf.int.style = "ribbon", # 设置置信区间的风格，"ribbon"或"step"
  pval = TRUE, # 在图上添加log-rank检验的p值
  # surv.median.line = "hv", #标注出中位生存时间
  risk.table = "abs_pct", # 在图下方添加风险表
  risk.table.col = "black", # "strata"根据不同的分组为风险表添加不同颜色
  risk.table.y.text = F, # 以图示（F）或文字（T）展示风险表中的组标签
  risk.table.height = 0.25, # 风险表的相对高度
  ncensor.plot = F, # 是否以条形图的形式展示随访过程中不同时间点死亡和删失的情况
  ggtheme = theme_bw(), # 改变图形风格
  palette = "lancet", # 配色风格
  legend.title = "Radiotherapy", # 图例标题
  xlab = "Time in months", # 设置x轴标签
  break.time.by = 40
) # 时间间隔
surv.byrad
summary(surv.byrad, time = c(12, 12 * 3, 12 * 5)) # 查看3年、5年总生存率


# 匹配后生存曲线
surv.byrad_iptw <- survfit(Surv(time, dead == 1) ~ rad, data = mydata, weights = w)

ggsurvplot(surv.byrad_iptw,
  fun = "pct",
  size = 1, censor.size = 3,
  conf.int = F,
  pval = F,
  risk.table = "abs_pct", risk.table.col = "black", risk.table.y.text = T,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiotherapy", xlab = "Time in months",
  break.time.by = 40
)
# 匹配后P值
library(RISCA)
ipw.log.rank(
  times = mydata$time,
  failures = mydata$dead == 1,
  variable = mydata$rad,
  weights = mydata$w
)
surv.byrad_iptw
summary(surv.byrad_iptw, time = c(12, 12 * 3, 12 * 5))





# 不同组织学类型内的匹配---------------------------------------------------------
lms_data <- subset(mydata, his == 1)
ess_data <- subset(mydata, his == 2)
ac_data <- subset(mydata, his == 3)

myvars3 <- c(
  "age", "year", "race", "single", "grade", "tumor_size", "T_stage", "N_stage",
  "surgery", "lymphadenectomy", "chemotherapy"
)

# 子宫平滑肌肉瘤内匹配
unmatched_lms_data <- CreateTableOne(
  data = lms_data, vars = myvars3,
  strata = "rad",
  addOverall = F, test = T
)
unmatched_lms_tab <- print(unmatched_lms_data,
  showAllLevels = TRUE,
  smd = F, quote = F, noSpaces = T
)


psmodel_lms <- glm(
  rad ~ age + year + race + single + grade + tumor_size + T_stage +
    N_stage + surgery + lymphadenectomy + chemotherapy,
  data = lms_data,
  family = binomial(link = "logit")
)
lms_data$ps <- predict(psmodel_lms, type = "response")
table(lms_data$rad)
pt_lms <- 317 / (317 + 1192)
lms_data$w <- ifelse(lms_data$rad == 1, pt_lms / lms_data$ps, (1 - pt_lms) / (1 - lms_data$ps))

iptw_lms <- svydesign(ids = ~0, data = lms_data, weights = ~w)
matched_lms_data <- svyCreateTableOne(
  data = iptw_lms,
  vars = myvars3, strata = "rad",
  addOverall = F, test = T
)
matched_lms_tab <- print(matched_lms_data,
  showAllLevels = TRUE,
  smd = F, quote = F, noSpaces = T
)

table_psm_lms <- cbind(unmatched_lms_tab, matched_lms_tab)
table_psm_lms <- rbind(
  Group = rep(c("Level", "No-RT", "RT", "P", "test method"), 2),
  table_psm_lms
)
colnames(table_psm_lms) <- c(
  "Level", "Unmatched", NA, NA, NA, "Level", "sIPTW",
  NA, NA, NA
)
print(table_psm_lms, quote = FALSE)



# 子宫间质肉瘤内匹配
unmatched_ess_data <- CreateTableOne(
  data = ess_data, vars = myvars3,
  strata = "rad",
  addOverall = F, test = T
)
unmatched_ess_tab <- print(unmatched_ess_data,
  showAllLevels = TRUE,
  smd = F,
  quote = F, noSpaces = T
)

psmodel_ess <- glm(
  rad ~ age + year + race + single + grade + tumor_size + T_stage +
    N_stage + surgery + lymphadenectomy + chemotherapy,
  data = ess_data,
  family = binomial(link = "logit")
)
ess_data$ps <- predict(psmodel_ess, type = "response")
table(ess_data$rad)
pt_ess <- 192 / (192 + 744)
ess_data$w <- ifelse(ess_data$rad == 1, pt_ess / ess_data$ps, (1 - pt_ess) / (1 - ess_data$ps))

iptw_ess <- svydesign(ids = ~0, data = ess_data, weights = ~w)
matched_ess_data <- svyCreateTableOne(
  data = iptw_ess,
  vars = myvars3,
  strata = "rad",
  addOverall = F, test = T
)
matched_ess_tab <- print(matched_ess_data,
  showAllLevels = TRUE,
  smd = F, quote = F, noSpaces = T
)

table_psm_ess <- cbind(unmatched_ess_tab, matched_ess_tab)
table_psm_ess <- rbind(
  Group = rep(c("Level", "No-RT", "RT", "P", "test method"), 2),
  table_psm_ess
)
colnames(table_psm_ess) <- c(
  "Level", "Unmatched", NA, NA, NA, "Level", "sIPTW",
  NA, NA, NA
)
print(table_psm_ess, quote = FALSE)




# 子宫腺肉瘤内匹配
unmatched_ac_data <- CreateTableOne(
  data = ac_data,
  vars = myvars3,
  strata = "rad",
  addOverall = F, test = T
)
unmatched_ac_tab <- print(unmatched_ac_data,
  showAllLevels = TRUE,
  smd = F, quote = F, noSpaces = T
)

psmodel_ac <- glm(
  rad ~ age + year + race + single + grade + tumor_size + T_stage +
    N_stage + surgery + lymphadenectomy + chemotherapy,
  data = ac_data,
  family = binomial(link = "logit")
)
ac_data$ps <- predict(psmodel_ac, type = "response")
table(ac_data$rad)
pt_ac <- 77 / (77 + 349)
ac_data$w <- ifelse(ac_data$rad == 1, pt_ac / ac_data$ps, (1 - pt_ac) / (1 - ac_data$ps))

iptw_ac <- svydesign(ids = ~0, data = ac_data, weights = ~w)
matched_ac_data <- svyCreateTableOne(
  data = iptw_ac,
  vars = myvars3,
  strata = "rad",
  addOverall = F, test = T
)
matched_ac_tab <- print(matched_ac_data,
  showAllLevels = TRUE,
  smd = F, quote = F, noSpaces = T
)

table_psm_ac <- cbind(unmatched_ac_tab, matched_ac_tab)
table_psm_ac <- rbind(
  Group = rep(c("Level", "No-RT", "RT", "P", "test method"), 2),
  table_psm_ac
)
colnames(table_psm_ac) <- c(
  "Level", "Unmatched", NA, NA, NA, "Level", "sIPTW",
  NA, NA, NA
)
print(table_psm_ac, quote = FALSE)

# 放疗在不同组织学类型中对预后的影响---------------------------------------------
# 放疗对子宫平滑肌肉瘤预后的影响
surv_lms <- survfit(Surv(time, dead == 1) ~ rad, data = lms_data)
ggsurvplot(surv_lms,
  fun = "pct",
  conf.int = F,
  pval = TRUE,
  risk.table = "abs_pct", risk.table.col = "black", risk.table.y.text = T,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiotherapy", xlab = "Time in months",
  break.time.by = 40
)
surv_lms
summary(surv_lms, time = c(12, 12 * 3, 12 * 5))


surv_lms_matched <- survfit(Surv(time, dead == 1) ~ rad, data = lms_data, weights = w)
ggsurvplot(surv_lms_matched,
  fun = "pct",
  conf.int = F,
  pval = F,
  risk.table = "abs_pct", risk.table.col = "black", risk.table.y.text = T,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiotherapy", xlab = "Time in months",
  break.time.by = 40
)
ipw.log.rank(
  times = lms_data$time,
  failures = lms_data$dead == 1,
  variable = lms_data$rad,
  weights = lms_data$w
)
surv_lms_matched
summary(surv_lms_matched, time = c(12, 12 * 3, 12 * 5))




# 放疗对子宫间质肉瘤预后的影响
surv_ess <- survfit(Surv(time, dead == 1) ~ rad, data = ess_data)
ggsurvplot(surv_ess,
  fun = "pct",
  conf.int = F,
  pval = TRUE,
  risk.table = "abs_pct", risk.table.col = "black", risk.table.y.text = T,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiotherapy", xlab = "Time in months",
  break.time.by = 40
)
surv_ess
summary(surv_ess, time = c(12, 12 * 3, 12 * 5))


surv_ess_matched <- survfit(Surv(time, dead == 1) ~ rad, data = ess_data, weights = w)
ggsurvplot(surv_ess_matched,
  fun = "pct",
  conf.int = F,
  pval = F,
  risk.table = "abs_pct", risk.table.col = "black", risk.table.y.text = T,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiotherapy", xlab = "Time in months",
  break.time.by = 40
)
ipw.log.rank(
  times = ess_data$time,
  failures = ess_data$dead == 1,
  variable = ess_data$rad,
  weights = ess_data$w
)
surv_ess_matched
summary(surv_ess_matched, time = c(12, 12 * 3, 12 * 5))



# 放疗对子宫腺肉瘤预后的影响
surv_ac <- survfit(Surv(time, dead == 1) ~ rad, data = ac_data)
ggsurvplot(surv_ac,
  fun = "pct",
  conf.int = F,
  pval = TRUE,
  risk.table = "abs_pct", risk.table.col = "black", risk.table.y.text = T,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiotherapy", xlab = "Time in months",
  break.time.by = 40
)
surv_ac
summary(surv_ac, time = c(12, 12 * 3, 12 * 5))


surv_ac_matched <- survfit(Surv(time, dead == 1) ~ rad, data = ac_data, weights = w)
ggsurvplot(surv_ac_matched,
  fun = "pct",
  conf.int = F,
  pval = F,
  risk.table = "abs_pct", risk.table.col = "black", risk.table.y.text = T,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiotherapy", xlab = "Time in months",
  break.time.by = 40
)
ipw.log.rank(
  times = ac_data$time,
  failures = ac_data$dead == 1,
  variable = ac_data$rad,
  weights = ac_data$w
)
surv_ac_matched
summary(surv_ac_matched, time = c(12, 12 * 3, 12 * 5))









# 🔴模型构建---------------------------------------------------------------------
# 🔴变量筛选---------------------------------------------------------------------
data_no_rad <- subset(mydata, radiotherapy == 0)
data_no_rad$rad <- NULL
data_no_rad$radiotherapy <- NULL

# 对于多分类变量建立一个设计矩阵以哑变量（0和1）的形式表示，并去掉第一列（截距项）
x <- model.matrix(
  ~ year + age + race + single + grade + tumor_size + his + T_stage + N_stage +
    surgery + lymphadenectomy + chemotherapy,
  data_no_rad
)[, -1]


y <- Surv(data_no_rad$time, data_no_rad$dead == 1) # 定义响应变量
library(glmnet)
set.seed(1996) # 设定随机种子数
lasso.cv <- cv.glmnet(
  x = x, y = y,
  family = "cox",
  alpha = 1, # 指定为lasso回归
  nfolds = 10, # 交叉验证次数，一般取10
  type.measure = "deviance",
  weights = data_no_rad$w
)
plot(lasso.cv)
abline(v = log(lasso.cv$lambda.min), lty = 1, lwd = 1, col = "black")

lasso.cv$lambda.min # 展示binominal deviance最小时的λ
log(lasso.cv$lambda.min) # 展示binominal deviance最小时的log(λ)
lasso.cv$lambda.1se # 展示binominal deviance 1个标准差时的λ
log(lasso.cv$lambda.1se) # 展示binominal deviance 1个标准差时的log(λ)

lasso <- glmnet(
  x = x, y = y,
  family = "cox",
  alpha = 1,
  standardize = T,
  nlambda = 100,
  weights = data_no_rad$w
)
plot(lasso, xvar = "lambda", label = T, lwd = 1.5)
abline(v = log(lasso.cv$lambda.min), lty = 1, lwd = 1, col = "black")
abline(v = log(lasso.cv$lambda.1se), lty = 3, lwd = 2, col = "black")

coef(lasso.cv, s = "lambda.min")
coef(lasso.cv, s = "lambda.1se")

# 以条形图的形式展示各变量的lasso回归系数----------------------------------------
lasso_coef <- data.frame(
  variables = c(
    "Age >58", "High-grade", "Gx",
    "Tumor size 70-165 mm", "Tumor size >165 mm",
    "ESS", "Adenosarcoma",
    "T2", "T3", "T4",
    "N stage", "Chemotherapy"
  ),
  coef = c(
    0.45524036, 0.89386973, 0.24944302,
    0.10882493, 0.38170697,
    -0.28306627, -0.02247745,
    0.59660811, 0.56347931, 0.93573491,
    0.48481647, 0.42979470
  )
)

ggplot(
  lasso_coef,
  aes(
    x = coef,
    y = reorder(variables, coef)
  )
) +
  geom_bar(
    stat = "identity",
    show.legend = T,
    width = .9,
    aes(fill = coef)
  ) + # 设置根据“coef”的大小来进行渐变填充
  scale_fill_gradient2(
    low = "#004688", # 设置条形图的渐变颜色
    mid = "blue",
    high = "#EA0100"
  ) +
  geom_text(
    aes(
      label = sprintf("%0.3f", coef),
      x = -0.07
    ),
    colour = "black",
    size = 3
  ) +
  xlab("Coefficients") +
  ylab("Variables") +
  theme_bw()






# 🔴cox回归建模------------------------------------------------------------------
library(rms)
data_no_rad$dead <- as.factor(data_no_rad$dead)
dd <- datadist(data_no_rad)
options(datadist = "dd")


model <- cph(
  Surv(time, dead == 1) ~ age + grade + tumor_size + his + T_stage + N_stage +
    chemotherapy,
  x = T, y = T, data = data_no_rad, surv = T, weights = w
)
model
cbind("回归系数" = coef(model), confint(model)) # 展示回归方程的系数及其95%CI
cbind("HR" = exp(coef(model)), exp(confint(model))) # 展示OR值（OR=e^β）及其95%CI




# 绘制cox回归森林图--------------------------------------------------------------
library(gdata)
forestplot <- read.xls("sarcoma_rad/forestplot.xlsx",
  sheet = 1,
  header = FALSE
) # 原始数据中第一行不是变量名称，而是标签

library(forestplot)
forestplot(
  labeltext = as.matrix(forestplot[, 1:7]), # 设置原表格中左边7列作为文本在图中展示
  mean = forestplot$V8, # 定义OR均值列
  lower = forestplot$V9, # 定义OR值的95%CI下限列
  upper = forestplot$V10, # 定义OR值的5%CI上限列
  graph.pos = 6, # 设置森林图出现的列，此处设置为3，则出现在第三列
  is.summary = c(
    T, T, F, F, F, T, F, F, F,
    T, F, F, F, F, T, F, F, F, T, F, F, F,
    F, T, F, F, T, F, F
  ), # 设置为TRUE的行文字以粗体出现,从而突出表头文字
  # 设置每列文字的对齐方式：l=左对齐；r=右对齐；c=居中对齐
  align = c("l", "c", "c", "c", "c", "l", "l"),
  # 设置图形中的行距，可设置为"auto"或者固定值，如unit(5,'mm')
  lineheight = unit(6, "mm"),
  colgap = unit(4, "mm"), # 设置图形中的列间距
  zero = 1, # 设置无效线，此处展示的是HR值，故无效线是1，而不是0
  lwd.zero = 1, # 设置无效线的粗细
  boxsize = 0.3, # 设置中央OR均值点的大小
  # 中央OR均值点的形状：fpDrawNormalCI=方形；fpDrawCircleCI=圆形
  fn.ci_norm = fpDrawCircleCI,
  lwd.ci = 2, # 设置95%CI线的粗细
  lty.ci = 1, # 95%CI线的线型：1=实线；>1=虚线
  clip = c(0.5, 4), # 设置森林图展示的置信区间范围，超过的部分用箭头展示
  xticks = seq(0.5, 4.5, by = 1), # X轴刻度范围（和置信区间范围一致）和精度
  xlog = F, # 转换为对数坐标轴
  ci.vertices = TRUE, # 是否显示95%CI线两端的小竖线
  ci.vertices.height = 0.15, # 设置95%CI线箭头的大小
  lwd.xaxis = 2, # 设置X轴线的粗细
  xlab = "Hazard Ratio", # 设置x轴标签
  col = fpColors( # 使用fpColors()函数定义图形元素的颜色
    box = "#EA0100", # OR均值点的颜色
    lines = "#004688", # 95%CI线的颜色
    zero = "darkgray"
  ), # 无效线的颜色
  txt_gp = fpTxtGp( # 设置所有文本元素的格式
    label = gpar(cex = 0.9), # 表格主体文字的大小
    ticks = gpar(cex = 0.9), # 森林图下方的坐标轴的刻度文字大小
    xlab = gpar(cex = 1.0)
  ), # 森林图下方X轴标签文字的大小
  graphwidth = unit(50, "mm")
) # 森林图的宽度





# 🔴绘制nomogram-----------------------------------------------------------------
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





# 导出模型便于后续构建网页计算器---------------------------------
saveRDS(model, file = "web_nomogram/data/sarcoma_rad_model.rds")
saveRDS(nom, file = "web_nomogram/data/sarcoma_rad_nomogram.rds")











# 🔴绘制ROC----------------------------------------------------------------------
mydata$pred <- predict(model, mydata)

library(timeROC)
timeroc <- timeROC(
  T = mydata$time, # 指定随访时间列
  delta = mydata$dead, # 指定生存状态列
  marker = mydata$pred, # 指定预测值列
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
ggplot(data = timeroc_df) +
  geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#BC3C29FF") +
  geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#0072B5FF") +
  geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
  theme_bw() +
  annotate("text",
    x = 0.75, y = 0.25, size = 4.5,
    label = paste0(
      "AUC at 1 year = ",
      sprintf("%.3f", timeroc$AUC[[1]])
    ),
    color = "#BC3C29FF"
  ) +
  annotate("text",
    x = 0.75, y = 0.15, size = 4.5,
    label = paste0(
      "AUC at 3 years = ",
      sprintf("%.3f", timeroc$AUC[[2]])
    ),
    color = "#0072B5FF"
  ) +
  annotate("text",
    x = 0.75, y = 0.05, size = 4.5,
    label = paste0(
      "AUC at 5 years = ",
      sprintf("%.3f", timeroc$AUC[[3]])
    ),
    color = "#E18727FF"
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

plotAUCcurve(timeroc, conf.int = T, conf.band = F, col = "#0072B5FF")




# 🔴AUC值的重复交叉验证----------------------------------------------------------
set.seed(1996)
# 设定交叉验证参数设定K和N，生成N*K个新数据集
library(caret)
folds <- createMultiFolds(
  y = mydata$dead,
  k = 10, # 设定为10折交叉验证
  times = 200
) # 设定重复次数


# 计算200次10折交叉验证的1年AUC值------------------------------------------------
auc_1year <- as.numeric() # 建立空表，放每次的AUC
for (i in 1:2000) { # 需要循环的次数=K*重复次数
  train <- mydata[folds[[i]], ]
  test <- mydata[-folds[[i]], ] # 标记验证集
  cvmodel <- cph(Surv(time, dead == 1) ~ age + grade + tumor_size + his + T_stage + N_stage + chemotherapy,
    x = T, y = T, data = train, surv = T, weights = w
  )
  test$pred <- predict(cvmodel, test)
  roc_1year <- timeROC(
    T = test$time, # 指定随访时间列
    delta = test$dead, # 指定生存状态列
    marker = test$pred, # 指定预测值列
    cause = 1, # 指定感兴趣的结局事件
    weighting = "marginal", # weighting="marginal"为采用Kaplan-Meier估计删失分布
    times = 12, # 绘制1年ROC
    ROC = TRUE, # 是否保存敏感度和特异度的预测值
    iid = F
  )
  auc_1year <- append(auc_1year, as.numeric(roc_1year$AUC[2])) # 提取并汇总每次计算得到的AUC值
}

summary(auc_1year) # 展示1年交叉验证后的AUC的平均值
cv <- data.frame(auc_1year) # 将所有AUC值以数据框的形式储存在cv表中



# 计算200次10折交叉验证的3年AUC值------------------------------------------------
auc_3year <- as.numeric()
for (i in 1:2000) {
  train <- mydata[folds[[i]], ]
  test <- mydata[-folds[[i]], ]
  cvmodel <- cph(Surv(time, dead == 1) ~ age + grade + tumor_size + his + T_stage + N_stage + chemotherapy,
    x = T, y = T, data = train, surv = T, weights = w
  )
  test$pred <- predict(cvmodel, test)
  roc_3year <- timeROC(
    T = test$time,
    delta = test$dead,
    marker = test$pred,
    cause = 1,
    weighting = "marginal",
    times = 12 * 3,
    ROC = TRUE,
    iid = F
  )
  auc_3year <- append(auc_3year, as.numeric(roc_3year$AUC[2]))
}

summary(auc_3year)
cv$auc_3year <- with(cv, auc_3year) # with()函数将3年AUC添加到之前的cv表中



# 计算200次10折交叉验证的5年AUC值------------------------------------------------
auc_5year <- as.numeric()
for (i in 1:2000) {
  test <- mydata[-folds[[i]], ]
  train <- mydata[folds[[i]], ]
  cvmodel <- cph(Surv(time, dead == 1) ~ age + grade + tumor_size + his + T_stage + N_stage + chemotherapy,
    x = T, y = T, data = train, surv = T, weights = w
  )
  test$pred <- predict(cvmodel, test)
  roc_5year <- timeROC(
    T = test$time,
    delta = test$dead,
    marker = test$pred,
    cause = 1,
    weighting = "marginal",
    times = 12 * 5,
    ROC = TRUE,
    iid = F
  )
  auc_5year <- append(auc_5year, as.numeric(roc_5year$AUC[2]))
}

summary(auc_5year)
cv$auc_5year <- with(cv, auc_5year)




# AUC值交叉验证结果可视化--------------------------------------------------------
# 需要用melt()函数将构建的cv表格转换为两列
# 第一列是值标签（auc_1year、auc_3year、auc_5year）
# 第二列是相应的值
library(reshape2)
cv <- data.frame(cv)
cv_plot <- melt(cv,
  measure.vars = c("auc_1year", "auc_3year", "auc_5year"), # 需要转换的列名
  variable.name = "Groups", # 新数据集的标签列名
  value.name = "AUC"
) # 新数据集的数值列名

# 绘制AUC值200次10折交叉验证的小提琴图
library(ggplot2)
library(ggprism)
ggplot(cv_plot, aes(Groups, AUC)) +
  geom_violin(aes(fill = Groups)) +
  geom_boxplot(width = 0.1) +
  theme_prism(base_size = 15, border = T) +
  theme_bw() +
  theme(legend.position = "none")





# 🔴绘制校准曲线-----------------------------------------------------------------
dd <- datadist(mydata)
options(datadist = "dd")

# 1年OS校准曲线------------------------------------------------------------------
fit_cal <- cph(Surv(time, dead == 1) ~ pred,
  x = T, y = T,
  data = mydata,
  surv = T,
  time.inc = 12,
  weights = w
)
cal_1 <- calibrate(fit_cal,
  u = 12, # 要评价的时间节点
  method = "boot",
  cmethod = "KM",
  m = 717, # 设置每多少个对象为一个评估单位，数值越小节点越多
  B = 1000
) # 迭代次数，通常设置为200或300
plot(cal_1,
  lwd = 2, lty = 1, ## 设置线条形状和尺寸
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
  xlab = "Nomogram-Predicted Probability of 1-year OS (%)",
  ylab = "Actual 1-year OS (%)",
  col = c(rgb(192, 98, 83, maxColorValue = 255)),
  xlim = c(0.6, 1), ylim = c(0.6, 1), ## x轴和y轴范围
  mgp = c(2, 1, 0)
) # 控制坐标轴的位置



# 3年OS校准曲线------------------------------------------------------------------
fit_cal2 <- cph(Surv(time, dead == 1) ~ pred,
  x = T, y = T,
  data = mydata,
  surv = T,
  time.inc = 12 * 3,
  weights = w
)
cal_2 <- calibrate(fit_cal2,
  u = 12 * 3,
  cmethod = "KM",
  method = "boot",
  m = 717,
  B = 1000
)
plot(cal_2,
  lwd = 2, lty = 1,
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
  xlab = "Nomogram-Predicted Probability of 3-year OS (%)",
  ylab = "Actual 3-year OS (%)",
  col = c(rgb(192, 98, 83, maxColorValue = 255)),
  xlim = c(0.4, 1), ylim = c(0.4, 1),
  mgp = c(2, 1, 0)
)



# 5年OS校准曲线------------------------------------------------------------------
fit_cal3 <- cph(Surv(time, dead == 1) ~ pred,
  x = T, y = T,
  data = mydata,
  surv = T,
  time.inc = 12 * 5,
  weights = w
)
cal_3 <- calibrate(fit_cal3,
  u = 12 * 5,
  cmethod = "KM", method = "boot",
  m = 717,
  B = 1000
)
plot(cal_3,
  lwd = 2, lty = 1,
  errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
  xlab = "Nomogram-Predicted Probability of 5-year OS (%)",
  ylab = "Actual 5-year OS (%)",
  col = c(rgb(192, 98, 83, maxColorValue = 255)),
  xlim = c(0.2, 1), ylim = c(0.2, 1),
  mgp = c(2, 1, 0)
)






# 🔴DCA曲线绘制------------------------------------------------------------------
library(ggDCA)
dca1 <- dca(model, new.data = mydata, times = c(12, 12 * 3, 12 * 5))

library(ggprism)
ggplot(dca1,
  linetype = F, # 线型
  lwd = 1
) + # 线宽
  theme_classic() + # 使用直线坐标系
  theme(legend.position = "top") + # 图例放在上方
  scale_x_continuous(
    limits = c(0, 1), # x轴范围并加入小刻度
    guide = "prism_minor"
  ) +
  scale_y_continuous(
    limits = c(-0.1, 0.4), # y轴范围并加入小刻度
    guide = "prism_minor"
  ) +
  scale_colour_prism(
    palette = "prism_dark", # 颜色,lengths(ggprism_data$colour_palettes)查看所有颜色主题
    name = "Cylinders",
    label = c(
      "1 year DCA", "3 year DCA", "5 year DCA",
      "ALL-1 year", "ALL-3 year", "ALL-5 year",
      "None"
    )
  ) +
  labs(title = "DCA based on ggDCA package") # 图形标题

AUDC(dca1) # 报告DCA曲线下面积（Area under Decision Curve, AUDC)

write.csv(dca1, file = "sarcoma_rad/net_benefit.csv") # 导出净收益，确定净获益的阈概率范围






# 🔴预后分层---------------------------------------------------------------------
# 得到模型计算的每个患者的1年、3年、5年生存概率
library(pec)
mydata$survprob <- predictSurvProb(model, newd = mydata, times = c(1 * 12, 3 * 12, 5 * 12))

# 计算每个患者的得分
library(nomogramFormula)
results <- formula_rd(nomogram = nom)
data_no_rad$points <- points_cal(formula = results$formula, rd = data_no_rad)
mydata$points <- points_cal(formula = results$formula, rd = mydata)
head(data_no_rad$points)

write.csv(mydata, "sarcoma_rad/mydata.csv")

# 决策树预后分层
library(rpart)
set.seed(1996)
tree <- rpart(Surv(time, dead == 1) ~ points, data = data_no_rad, weights = w)
tree$cptable
cp <- tree$cptable[which.min(tree$cptable[, "xerror"]), "CP"]
prune <- prune(tree, 0.02)
# 决策树绘制
library(rpart.plot)
rpart.plot(prune, # 决策树对象
  type = 2, # 决策树样式
  extra = 1 + 100, # 展示各组结局事件发生数(extra=1)和比例(+100)
  under = F, box.palette = "auto", shadow.col = "gray"
) # 是否将信息在树叶下方展示
library(partykit)
plot(as.party(prune)) # 展示各组生存曲线


mydata$risk_group <- ifelse(mydata$points < 135, 1,
  ifelse(mydata$points < 221, 2, 3)
)




# 各风险组生存曲线比较
surv.by_risk_group <- survfit(Surv(time, dead == 1) ~ risk_group, data = mydata)
ggsurvplot(surv.by_risk_group,
  fun = "pct",
  conf.int = T, conf.int.style = "step",
  pval = TRUE,
  risk.table = "abs_pct", risk.table.col = "strata", risk.table.y.text = T,
  ncensor.plot = F,
  ggtheme = theme_bw(),
  palette = c("#42B43F", "#00468A", "#EA0100"),
  legend.title = "Risk Group", xlab = "Time in months",
  break.time.by = 40
)
surv.by_risk_group
summary(surv.by_risk_group, time = c(12, 12 * 3, 12 * 5))




# 🔴各预后亚组sIPTW--------------------------------------------------------------
data_group_1 <- subset(mydata, mydata$risk_group == 1)
data_group_2 <- subset(mydata, mydata$risk_group == 2)
data_group_3 <- subset(mydata, mydata$risk_group == 3)

# 低危组sIPTW匹配----------------------------------------------------------------
unmatcheddata_1 <- CreateTableOne(
  data = data_group_1,
  vars = myvars2,
  strata = "rad",
  addOverall = F, test = T
)
unmatchedtab_1 <- print(unmatcheddata_1,
  showAllLevels = TRUE,
  smd = F, quote = F, noSpaces = T
)


psmodel_1 <- glm(
  rad ~ age + year + race + single + grade + his + tumor_size + T_stage + N_stage +
    surgery + lymphadenectomy + chemotherapy,
  data = data_group_1,
  family = binomial(link = "logit")
)
data_group_1$ps <- predict(psmodel_1, type = "response")
table(data_group_1$rad)
pt_1 <- 176 / (176 + 1043)
data_group_1$w <- ifelse(data_group_1$rad == 1, pt_1 / data_group_1$ps,
  (1 - pt_1) / (1 - data_group_1$ps)
)

iptw_1 <- svydesign(ids = ~0, data = data_group_1, weights = ~w)
matcheddata_1 <- svyCreateTableOne(
  data = iptw_1,
  vars = myvars2,
  strata = "rad",
  addOverall = F, test = T
)
matchedtab_1 <- print(matcheddata_1,
  showAllLevels = TRUE,
  smd = F, quote = F, noSpaces = T
)



table_psm_1 <- cbind(unmatchedtab_1, matchedtab_1)
table_psm_1 <- rbind(
  Group = rep(c("Level", "No-RT", "RT", "P", "test method"), 2),
  table_psm_1
)
colnames(table_psm_1) <- c(
  "Level", "Unmatched", NA, NA, NA, "Level", "sIPTW",
  NA, NA, NA
)
print(table_psm_1, quote = FALSE)
write.csv(table_psm_1, file = "sarcoma_rad/baseline_table_sIPTW_group_1.csv")




# 中危组sIPTW匹配----------------------------------------------------------------
unmatcheddata_2 <- CreateTableOne(
  data = data_group_2,
  vars = myvars2,
  strata = "rad",
  addOverall = F, test = T
)
unmatchedtab_2 <- print(unmatcheddata_2,
  showAllLevels = TRUE,
  smd = F, quote = F, noSpaces = T
)


psmodel_2 <- glm(
  rad ~ age + year + race + single + grade + his + tumor_size + T_stage + N_stage +
    surgery + lymphadenectomy + chemotherapy,
  data = data_group_2,
  family = binomial(link = "logit")
)
data_group_2$ps <- predict(psmodel_2, type = "response")
table(data_group_2$rad)
pt_2 <- 288 / (288 + 917)
data_group_2$w <- ifelse(data_group_2$rad == 1, pt_2 / data_group_2$ps,
  (1 - pt_2) / (1 - data_group_2$ps)
)

iptw_2 <- svydesign(ids = ~0, data = data_group_2, weights = ~w)
matcheddata_2 <- svyCreateTableOne(
  data = iptw_2,
  vars = myvars2,
  strata = "rad",
  addOverall = F, test = T
)
matchedtab_2 <- print(matcheddata_2,
  showAllLevels = TRUE,
  smd = F, quote = F, noSpaces = T
)

table_psm_2 <- cbind(unmatchedtab_2, matchedtab_2)
table_psm_2 <- rbind(
  Group = rep(c("Level", "No-RT", "RT", "P", "test method"), 2),
  table_psm_2
)
colnames(table_psm_2) <- c(
  "Level", "Unmatched", NA, NA, NA, "Level", "sIPTW",
  NA, NA, NA
)
print(table_psm_2, quote = FALSE)
write.csv(table_psm_2, file = "sarcoma_rad/baseline_table_sIPTW_group_2.csv")




# 高危组sIPTW匹配----------------------------------------------------------------
unmatcheddata_3 <- CreateTableOne(
  data = data_group_3,
  vars = myvars2,
  strata = "rad",
  addOverall = F, test = T
)
unmatchedtab_3 <- print(unmatcheddata_3,
  showAllLevels = TRUE,
  smd = F, quote = F, noSpaces = T
)


psmodel_3 <- glm(
  rad ~ age + year + race + single + grade + his + tumor_size + T_stage + N_stage +
    surgery + lymphadenectomy + chemotherapy,
  data = data_group_3,
  family = binomial(link = "logit")
)
data_group_3$ps <- predict(psmodel_3, type = "response")
table(data_group_3$rad)
pt_3 <- 122 / (122 + 325)
data_group_3$w <- ifelse(data_group_3$rad == 1, pt_3 / data_group_3$ps,
  (1 - pt_3) / (1 - data_group_3$ps)
)

iptw_3 <- svydesign(ids = ~0, data = data_group_3, weights = ~w)
matcheddata_3 <- svyCreateTableOne(
  data = iptw_3,
  vars = myvars2,
  strata = "rad",
  addOverall = F, test = T
)
matchedtab_3 <- print(matcheddata_3,
  showAllLevels = TRUE,
  smd = F, quote = F, noSpaces = T
)

table_psm_3 <- cbind(unmatchedtab_3, matchedtab_3)
table_psm_3 <- rbind(
  Group = rep(c("Level", "No-RT", "RT", "P", "test method"), 2),
  table_psm_3
)
colnames(table_psm_3) <- c(
  "Level", "Unmatched", NA, NA, NA, "Level", "sIPTW",
  NA, NA, NA
)
print(table_psm_3, quote = FALSE)
write.csv(table_psm_3, file = "sarcoma_rad/baseline_table_sIPTW_group_3.csv")




# 🔴在各预后组中评价放疗效果-----------------------------------------------------
# 低危组放疗效果评价-------------------------------------------------------------
surv.by_rad_risk_group_1 <- survfit(Surv(time, dead == 1) ~ rad, data = data_group_1)
ggsurvplot(surv.by_rad_risk_group_1,
  fun = "pct",
  conf.int = T, conf.int.style = "step",
  pval = TRUE,
  risk.table = "abs_pct", risk.table.col = "strata", risk.table.y.text = F,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiatherapy", xlab = "Time in months",
  break.time.by = 40
)
surv.by_rad_risk_group_1
summary(surv.by_rad_risk_group_1, time = c(12, 12 * 3, 12 * 5))


surv.by_rad_risk_group_1_psm <- survfit(Surv(time, dead == 1) ~ rad,
  data = data_group_1,
  weights = data_group_1$w
)
ggsurvplot(surv.by_rad_risk_group_1_psm,
  fun = "pct",
  conf.int = T, conf.int.style = "step",
  pval = F,
  risk.table = "abs_pct", risk.table.col = "strata", risk.table.y.text = F,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiatherapy", xlab = "Time in months",
  break.time.by = 40
)
ipw.log.rank(
  times = data_group_1$time,
  failures = data_group_1$dead == 1,
  variable = data_group_1$rad,
  weights = data_group_1$w
)
surv.by_rad_risk_group_1_psm
summary(surv.by_rad_risk_group_1_psm, time = c(12, 12 * 3, 12 * 5))


# 中危组放疗效果评价-------------------------------------------------------------
surv.by_rad_risk_group_2 <- survfit(Surv(time, dead == 1) ~ rad, data = data_group_2)
ggsurvplot(surv.by_rad_risk_group_2,
  fun = "pct",
  conf.int = T, conf.int.style = "step",
  pval = TRUE,
  risk.table = "abs_pct", risk.table.col = "strata", risk.table.y.text = F,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiatherapy", xlab = "Time in months",
  break.time.by = 40
)
surv.by_rad_risk_group_2
summary(surv.by_rad_risk_group_2, time = c(12, 12 * 3, 12 * 5))

surv.by_rad_risk_group_2_psm <- survfit(Surv(time, dead == 1) ~ rad,
  data = data_group_2,
  weights = data_group_2$w
)
ggsurvplot(surv.by_rad_risk_group_2_psm,
  fun = "pct",
  conf.int = T, conf.int.style = "step",
  pval = F,
  risk.table = "abs_pct", risk.table.col = "strata", risk.table.y.text = F,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiatherapy", xlab = "Time in months",
  break.time.by = 40
)
ipw.log.rank(
  times = data_group_2$time,
  failures = data_group_2$dead == 1,
  variable = data_group_2$rad,
  weights = data_group_2$w
)
surv.by_rad_risk_group_2_psm
summary(surv.by_rad_risk_group_2_psm, time = c(12, 12 * 3, 12 * 5))



# 高危组放疗效果评价-------------------------------------------------------------
surv.by_rad_risk_group_3 <- survfit(Surv(time, dead == 1) ~ rad, data = data_group_3)
ggsurvplot(surv.by_rad_risk_group_3,
  fun = "pct",
  conf.int = T, conf.int.style = "step",
  pval = TRUE,
  risk.table = "abs_pct", risk.table.col = "strata", risk.table.y.text = F,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiatherapy", xlab = "Time in months",
  break.time.by = 40
)
surv.by_rad_risk_group_3
summary(surv.by_rad_risk_group_3, time = c(12, 12 * 3, 12 * 5))

surv.by_rad_risk_group_3_psm <- survfit(Surv(time, dead == 1) ~ rad,
  data = data_group_3,
  weights = data_group_3$w
)
surv.by_rad_risk_group_3_psm
summary(surv.by_rad_risk_group_3_psm, time = c(12, 12 * 3, 12 * 5))
ggsurvplot(surv.by_rad_risk_group_3_psm,
  fun = "pct",
  conf.int = T, conf.int.style = "step",
  pval = F,
  risk.table = "abs_pct", risk.table.col = "strata", risk.table.y.text = F,
  ncensor.plot = F,
  ggtheme = theme_bw(), palette = "lancet",
  legend.title = "Radiatherapy", xlab = "Time in months",
  break.time.by = 40
)
ipw.log.rank(
  times = data_group_3$time,
  failures = data_group_3$dead == 1,
  variable = data_group_3$rad,
  weights = data_group_3$w
)
surv.by_rad_risk_group_3_psm
summary(surv.by_rad_risk_group_3_psm, time = c(12, 12 * 3, 12 * 5))








group_compare <- CreateTableOne(
  data = mydata,
  vars = myvars1,
  strata = "risk_group",
  test = T,
  addOverall = T
)
group_compare_tab <- print(group_compare,
  showAllLevels = T,
  smd = F,
  quote = F, noSpaces = T
)
write.csv(group_compare_tab, file = "sarcoma_rad/risk_group_compare_table.csv")



# 🟢在不同组织学类型中评价模型---------------------------------------------------
# 🟢ROC
lms_data_2 <- subset(mydata, his == 1)
lms_data_2 <- lms_data_2[, -23]
ess_data_2 <- subset(mydata, his == 2)
ess_data_2 <- ess_data_2[, -23]
as_data_2 <- subset(mydata, his == 3)
as_data_2 <- as_data_2[, -23]

mytimeroc <- function(dat) {
  timeroc <- timeROC(
    T = dat$time, # 指定随访时间列
    delta = dat$dead, # 指定生存状态列
    marker = dat$pred, # 指定预测值列
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
  ggplot(data = timeroc_df) +
    geom_line(aes(x = FP_1year, y = TP_1year), size = 1, color = "#BC3C29FF") +
    geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#0072B5FF") +
    geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
    geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
    theme_bw() +
    annotate("text",
      x = 0.75, y = 0.25, size = 4.5,
      label = paste0(
        "AUC at 1 year = ",
        sprintf("%.3f", timeroc$AUC[[1]])
      ),
      color = "#BC3C29FF"
    ) +
    annotate("text",
      x = 0.75, y = 0.15, size = 4.5,
      label = paste0(
        "AUC at 3 years = ",
        sprintf("%.3f", timeroc$AUC[[2]])
      ),
      color = "#0072B5FF"
    ) +
    annotate("text",
      x = 0.75, y = 0.05, size = 4.5,
      label = paste0(
        "AUC at 5 years = ",
        sprintf("%.3f", timeroc$AUC[[3]])
      ),
      color = "#E18727FF"
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
}
mytimeroc(lms_data_2)
mytimeroc(ess_data_2)
mytimeroc(as_data_2)

# 🟢校准曲线
calcurve <- function(dat, ngroup) {
  dd <- datadist(dat)
  options(datadist = "dd")

  fit_cal <- cph(Surv(time, dead == 1) ~ pred,
    x = T, y = T,
    data = dat,
    surv = T,
    time.inc = 12,
    weights = w
  )
  cal_1 <- calibrate(fit_cal,
    u = 12, # 要评价的时间节点
    method = "boot",
    cmethod = "KM",
    m = nrow(dat) / ngroup, # 设置每多少个对象为一个评估单位，数值越小节点越多
    B = 1000
  ) # 迭代次数，通常设置为200或300
  plot(cal_1,
    lwd = 2, lty = 1, ## 设置线条形状和尺寸
    errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
    xlab = "Nomogram-Predicted Probability of 1-year OS (%)",
    ylab = "Actual 1-year OS (%)",
    col = c(rgb(192, 98, 83, maxColorValue = 255)),
    xlim = c(0.5, 1), ylim = c(0.5, 1), ## x轴和y轴范围
    mgp = c(2, 1, 0)
  ) # 控制坐标轴的位置


  fit_cal2 <- cph(Surv(time, dead == 1) ~ pred,
    x = T, y = T,
    data = dat,
    surv = T,
    time.inc = 12 * 3,
    weights = w
  )
  cal_2 <- calibrate(fit_cal2,
    u = 12 * 3,
    cmethod = "KM",
    method = "boot",
    m = nrow(dat) / ngroup,
    B = 1000
  )
  plot(cal_2,
    lwd = 2, lty = 1,
    errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
    xlab = "Nomogram-Predicted Probability of 3-year OS (%)",
    ylab = "Actual 3-year OS (%)",
    col = c(rgb(192, 98, 83, maxColorValue = 255)),
    xlim = c(0.2, 1), ylim = c(0.2, 1),
    mgp = c(2, 1, 0)
  )



  # 5年OS校准曲线
  fit_cal3 <- cph(Surv(time, dead == 1) ~ pred,
    x = T, y = T,
    data = dat,
    surv = T,
    time.inc = 12 * 5,
    weights = w
  )
  cal_3 <- calibrate(fit_cal3,
    u = 12 * 5,
    cmethod = "KM", method = "boot",
    m = nrow(dat) / ngroup,
    B = 1000
  )
  plot(cal_3,
    lwd = 2, lty = 1,
    errbar.col = c(rgb(0, 118, 192, maxColorValue = 255)),
    xlab = "Nomogram-Predicted Probability of 5-year OS (%)",
    ylab = "Actual 5-year OS (%)",
    col = c(rgb(192, 98, 83, maxColorValue = 255)),
    xlim = c(0.2, 1), ylim = c(0.2, 1),
    mgp = c(2, 1, 0)
  )
}
calcurve(lms_data_2, 4)
calcurve(ess_data_2, 3)
calcurve(as_data_2, 3)




sessionInfo()
