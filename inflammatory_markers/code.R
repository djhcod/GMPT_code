mydata <- read.csv("inflammatory_markers/ovary_data.csv")
library(caret)
set.seed(2021)
training.samples <- createDataPartition(mydata$class,
  p = 0.7,
  list = F
)
lassodata <- mydata[training.samples, ]
train <- mydata[training.samples, ]
test <- mydata[-training.samples, ]


lassodata[3:5] <- NULL
lassodata$menstrual_status <- NULL
lassodata$hgb <- NULL
lassodata$hct <- NULL
lassodata$mcv <- NULL
lassodata$lym <- NULL
lassodata$mono <- NULL
lassodata$neut <- NULL
lassodata$baso <- NULL
lassodata$plt <- NULL
lassodata$mpv <- NULL
lassodata$fib <- NULL
lassodata$alb <- NULL
lassodata$ca125 <- NULL
lassodata$he4 <- NULL
lassodata$nlr <- NULL
lassodata$plr <- NULL
lassodata$mlr <- NULL
lassodata$far <- NULL
lassodata$sii <- NULL

x <- as.matrix(lassodata[, 3:27])
y <- as.matrix(lassodata[, 2])
print(paste0("The training cohort contains ", dim(x)[1], " cases; ", dim(x)[2], " variables"))
colnames(x)


library(glmnet)
set.seed(1)
lasso.cv <- cv.glmnet(
  x = x, y = y, family = "binomial",
  alpha = 1,
  nfolds = 10,
  type.measure = "deviance"
)
plot(lasso.cv)
abline(v = log(lasso.cv$lambda.min), lty = 1, lwd = 2, col = "black")
lasso.cv$lambda.min
log(lasso.cv$lambda.min)
lasso.cv$lambda.1se
log(lasso.cv$lambda.1se)
lasso <- glmnet(
  x = x, y = y, family = "binomial", alpha = 1, upper.limits = 6,
  lower.limits = -3, standardize = T, nlambda = 100
)
plot(lasso, xvar = "lambda", label = T)
abline(v = log(lasso.cv$lambda.min), lty = 1, lwd = 2, col = "black")
abline(v = log(lasso.cv$lambda.1se), lty = 3, lwd = 2, col = "black")
lasso.lamda.min <- glmnet(x = x, y = y, family = "binomial", alpha = 1, lambda = lasso.cv$lambda.min)
lasso.lamda.1se <- glmnet(x = x, y = y, family = "binomial", alpha = 1, lambda = lasso.cv$lambda.1se)
coef(lasso.lamda.min, s = "lambda.min")
coef(lasso.lamda.1se, s = "lambda.1se")



fullfit <- glm(
  class ~ age + ca125_cutoff + he4_cutoff
    + fib_cutoff + mlr_cutoff
    + irregular + bloodflow + solid_areas,
  data = train, family = binomial(), x = TRUE
)
summary(fullfit)
cbind("Coefficients" = coef(fullfit), confint(fullfit))
cbind("OR" = exp(coef(fullfit)), exp(confint(fullfit)))
summary(fullfit)$coefficients[, 3]^2
library(car)
vif(fullfit)


fit1 <- glm(class ~ age + ca125_cutoff + he4_cutoff,
  data = train, family = binomial(), x = TRUE
)
summary(fit1)
cbind("Coefficients" = coef(fit1), confint(fit1))
cbind("OR" = exp(coef(fit1)), exp(confint(fit1)))
summary(fit1)$coefficients[, 3]^2
table(train$ca125_cutoff)
table(train$he4_cutoff)


fit2 <- glm(class ~ age + ca125_cutoff + he4_cutoff + bloodflow + solid_areas,
  data = train, family = binomial(), x = TRUE
)
summary(fit2)
cbind("Coefficients" = coef(fit2), confint(fit2))
cbind("OR" = exp(coef(fit2)), exp(confint(fit2)))
summary(fit2)$coefficients[, 3]^2



fit3 <- glm(
  class ~ age + ca125_cutoff + he4_cutoff
    + fib_cutoff + mlr_cutoff + bloodflow + solid_areas,
  data = train, family = binomial(), x = TRUE
)
summary(fit3)
cbind("Coefficients" = coef(fit3), confint(fit3))
cbind("OR" = exp(coef(fit3)), exp(confint(fit3)))
summary(fit3)$coefficients[, 3]^2
table(train$fib_cutoff)
table(train$mlr_cutoff)





library(forestplot)
forestplot1 <- read.csv("inflammatory_markers/forestplot_data_model_1.csv",
  header = FALSE
)
forestplot(
  labeltext = as.matrix(forestplot1[, 1:7]),
  mean = forestplot1$V8,
  lower = forestplot1$V9,
  upper = forestplot1$V10,
  graph.pos = 6,
  is.summary = c(T, F, T, F, F, T, F, F, T),
  align = c("l", "c", "c", "c", "c", "l", "l"),
  lineheight = unit(8, "mm"),
  colgap = unit(4, "mm"),
  zero = 1,
  lwd.zero = 2,
  boxsize = 0.4,
  fn.ci_norm = fpDrawNormalCI,
  lwd.ci = 2.5,
  lty.ci = 1,
  clip = c(0, 25),
  ci.vertices = TRUE,
  ci.vertices.height = 0.3,
  lwd.xaxis = 2,
  xticks = seq(1, 25, by = 3),
  xlab = "Odds Ratio",
  col = fpColors(
    box = rgb(176, 35, 24, maxColorValue = 255),
    lines = rgb(78, 92, 104, maxColorValue = 255),
    zero = rgb(78, 92, 104, maxColorValue = 255)
  ),
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.9),
    ticks = gpar(cex = 0.8),
    xlab = gpar(cex = 1.0)
  ),
  grid = (
    gp <- gpar(
      col = rgb(78, 92, 104, maxColorValue = 255),
      lty = 3,
      lwd = 1
    )),
  graphwidth = unit(60, "mm")
)

forestplot2 <- read.csv("inflammatory_markers/forestplot_data_model_2.csv",
  header = FALSE
)
forestplot(
  labeltext = as.matrix(forestplot2[, 1:7]),
  mean = forestplot2$V8,
  lower = forestplot2$V9,
  upper = forestplot2$V10,
  graph.pos = 6,
  is.summary = c(T, F, T, F, F, T, F, F, T, F, F, T, F, F, T),
  align = c("l", "c", "c", "c", "c", "l", "l"),
  lineheight = unit(6, "mm"),
  colgap = unit(4, "mm"),
  zero = 1,
  lwd.zero = 2,
  boxsize = 0.4,
  fn.ci_norm = fpDrawNormalCI,
  lwd.ci = 2.5,
  lty.ci = 1,
  clip = c(0, 30),
  ci.vertices = TRUE,
  ci.vertices.height = 0.3,
  lwd.xaxis = 2,
  xticks = seq(1, 31, by = 5),
  xlab = "Odds Ratio",
  col = fpColors(
    box = rgb(176, 35, 24, maxColorValue = 255),
    lines = rgb(78, 92, 104, maxColorValue = 255),
    zero = rgb(78, 92, 104, maxColorValue = 255)
  ),
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.9),
    ticks = gpar(cex = 0.8),
    xlab = gpar(cex = 1.0)
  ),
  grid = (
    gp <- gpar(
      col = rgb(78, 92, 104, maxColorValue = 255),
      lty = 3,
      lwd = 1
    )),
  graphwidth = unit(60, "mm")
)

forestplot3 <- read.csv("inflammatory_markers/forestplot_data_model_3.csv",
  header = FALSE
)
forestplot(
  labeltext = as.matrix(forestplot3[, 1:7]),
  mean = forestplot3$V8,
  lower = forestplot3$V9,
  upper = forestplot3$V10,
  graph.pos = 6,
  is.summary = c(T, F, T, F, F, T, F, F, T, F, F, T, F, F, T, F, F, T, F, F, T),
  align = c("l", "c", "c", "c", "c", "l", "l"),
  lineheight = unit(6, "mm"),
  colgap = unit(4, "mm"),
  zero = 1,
  lwd.zero = 2,
  boxsize = 0.4,
  fn.ci_norm = fpDrawNormalCI,
  lwd.ci = 2.5,
  lty.ci = 1,
  clip = c(0, 19),
  ci.vertices = TRUE,
  ci.vertices.height = 0.3,
  lwd.xaxis = 2,
  xticks = seq(1, 21, by = 5),
  xlab = "Odds Ratio",
  col = fpColors(
    box = rgb(176, 35, 24, maxColorValue = 255),
    lines = rgb(78, 92, 104, maxColorValue = 255),
    zero = rgb(78, 92, 104, maxColorValue = 255)
  ),
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.9),
    ticks = gpar(cex = 0.8),
    xlab = gpar(cex = 1.0)
  ),
  grid = (
    gp <- gpar(
      col = rgb(78, 92, 104, maxColorValue = 255),
      lty = 3,
      lwd = 1
    )),
  graphwidth = unit(60, "mm")
)




train$pred1 <- predict(fit1, newdata = train, type = "response")
train$pred2 <- predict(fit2, newdata = train, type = "response")
train$pred3 <- predict(fit3, newdata = train, type = "response")
library(pROC)
roc1 <- roc(train$class, train$pred1)
roc2 <- roc(train$class, train$pred2)
roc3 <- roc(train$class, train$pred3)
roc.roma <- roc(train$class, train$roma_pred)
roc.cph <- roc(train$class, train$cphi_pred)
roc.rmi <- roc(train$class, train$rmi)
roc.lr2 <- roc(train$class, train$lr2_pred)
roc.rops <- roc(train$class, train$rops)
auc(roc1)
auc(roc2)
auc(roc3)
ci.auc(roc1)
ci.auc(roc2)
ci.auc(roc3)
roc.test(roc1, roc2, method = "delong")
roc.test(roc1, roc3, method = "delong")
roc.test(roc2, roc3, method = "delong")
plot.roc(roc1,
  col = rgb(61, 64, 91, maxColorValue = 255),
  lwd = 1.7, lty = 1,
  main = "ROC",
  legacy.axes = T, xlim = c(1, 0), ylim = c(0, 1), asp = 1
)
plot.roc(roc2,
  add = T,
  col = rgb(58, 134, 255, maxColorValue = 255),
  lwd = 1.7, lty = 1
)
plot.roc(roc3,
  add = T,
  col = rgb(239, 71, 111, maxColorValue = 255),
  lwd = 2, lty = 1
)
legend(0.5, 0.3,
  c(
    paste("AUC of Model 1: ", sprintf("%.3f", auc(roc1))),
    paste("AUC of Model 2: ", sprintf("%.3f", auc(roc2))),
    paste("AUC of Model 3: ", sprintf("%.3f", auc(roc3)))
  ),
  col = c(
    rgb(61, 64, 91, maxColorValue = 255),
    rgb(58, 134, 255, maxColorValue = 255),
    rgb(239, 71, 111, maxColorValue = 255)
  ),
  lty = c(1, 1, 1),
  lwd = c(1.7, 1.7, 2),
  box.lty = 2, cex = 0.5, horiz = F, seg.len = 5
)




library(nricens)
nribin(
  mdl.std = fit1,
  mdl.new = fit2,
  updown = "category",
  cut = 0.393,
  niter = 1000
)
library(PredictABEL)
reclassification(
  data = train,
  cOutcome = 5,
  predrisk1 = train$pred1,
  predrisk2 = train$pred2,
  cutoff = c(0, 0.393, 1)
)
nribin(
  mdl.std = fit1,
  mdl.new = fit3,
  updown = "category",
  cut = 0.18,
  niter = 1000
)
reclassification(
  data = train,
  cOutcome = 5,
  predrisk1 = train$pred1,
  predrisk2 = train$pred3,
  cutoff = c(0, 0.18, 1)
)
nribin(
  mdl.std = fit2,
  mdl.new = fit3,
  updown = "category",
  cut = 0.18,
  niter = 1000
)
reclassification(
  data = train,
  cOutcome = 5,
  predrisk1 = train$pred2,
  predrisk2 = train$pred3,
  cutoff = c(0, 0.18, 1)
)
write.csv(train, file = "inflammatory_markers/train.csv")





auc(roc.rops)
auc(roc.roma)
auc(roc.cph)
auc(roc.rmi)
auc(roc.lr2)
ci.auc(roc.rops)
ci.auc(roc.roma)
ci.auc(roc.cph)
ci.auc(roc.rmi)
ci.auc(roc.lr2)
roc.test(roc3, roc.rops, method = "delong")
roc.test(roc3, roc.roma, method = "delong")
roc.test(roc3, roc.cph, method = "delong")
roc.test(roc3, roc.rmi, method = "delong")
roc.test(roc3, roc.lr2, method = "delong")
roc.test(roc.rmi, roc.lr2, method = "delong")
plot.roc(roc3,
  col = rgb(239, 71, 111, maxColorValue = 255), lwd = 2, lty = 1,
  main = "ROC",
  legacy.axes = T, xlim = c(1, 0), ylim = c(0, 1), asp = 1
)
plot.roc(roc.rops,
  add = T,
  col = rgb(58, 134, 255, maxColorValue = 255), lwd = 1.7, lty = 1
)
plot.roc(roc.roma,
  add = T,
  col = rgb(61, 64, 91, maxColorValue = 255), lwd = 1.7, lty = 2
)
plot.roc(roc.cph, add = T, col = "orange", lwd = 1.7, lty = 1)
plot.roc(roc.rmi,
  add = T,
  col = rgb(129, 178, 154, maxColorValue = 255), lwd = 1.7, lty = 1
)
plot.roc(roc.lr2,
  add = T,
  col = rgb(61, 64, 91, maxColorValue = 255), lwd = 1.7, lty = 1
)
legend(0.5, 0.3,
  c(
    paste("AUC of Model 3: ", sprintf("%.3f", auc(roc3))),
    paste("AUC of R-OPS: ", sprintf("%.3f", auc(roc.rops))),
    paste("AUC of ROMA: ", sprintf("%.3f", auc(roc.roma))),
    paste("AUC of CPH-I: ", sprintf("%.3f", auc(roc.cph))),
    paste("AUC of RMI 4: ", sprintf("%.3f", auc(roc.rmi))),
    paste("AUC of LR2: ", sprintf("%.3f", auc(roc.lr2)))
  ),
  col = c(rgb(239, 71, 111, maxColorValue = 255),
    rgb(58, 134, 255, maxColorValue = 255),
    rgb(61, 64, 91, maxColorValue = 255),
    col = "orange",
    rgb(129, 178, 154, maxColorValue = 255),
    rgb(61, 64, 91, maxColorValue = 255)
  ),
  lty = c(1, 1, 2, 1, 1, 1),
  lwd = c(2, 1.7, 1.7, 1.7, 1.7, 1.7),
  box.lty = 2, cex = 0.5, horiz = F, seg.len = 5
)



# True positive（a）
TP.fit3 <- dim(train[as.numeric(train$class) == 1 & train$pred3 > 0.18, ])[1]
# False positive（b）
FP.fit3 <- dim(train[as.numeric(train$class) == 0 & train$pred3 > 0.18, ])[1]
# False negative（c）
FN.fit3 <- dim(train[as.numeric(train$class) == 1 & train$pred3 <= 0.18, ])[1]
# True negative（d）
TN.fit3 <- dim(train[as.numeric(train$class) == 0 & train$pred3 <= 0.18, ])[1]
SN.fit3 <- TP.fit3 / (TP.fit3 + FN.fit3)
SN.fit3 # SE
SP.fit3 <- TN.fit3 / (TN.fit3 + FP.fit3)
SP.fit3 # SP
TP.fit3 / (TP.fit3 + FP.fit3) # PPV
TN.fit3 / (TN.fit3 + FN.fit3) # NPV
SN.fit3 + SP.fit3 - 1 # YI


TP.rops <- dim(train[as.numeric(train$class) == 1 & train$rops > 330, ])[1]
FP.rops <- dim(train[as.numeric(train$class) == 0 & train$rops > 330, ])[1]
FN.rops <- dim(train[as.numeric(train$class) == 1 & train$rops <= 330, ])[1]
TN.rops <- dim(train[as.numeric(train$class) == 0 & train$rops <= 330, ])[1]
SN.rops <- TP.rops / (TP.rops + FN.rops)
SN.rops # SE
SP.rops <- TN.rops / (TN.rops + FP.rops)
SP.rops # SP
TP.rops / (TP.rops + FP.rops) # PPV
TN.rops / (TN.rops + FN.rops) # NPV
SN.rops + SP.rops - 1 # YI


TP.cph <- dim(train[as.numeric(train$class) == 1 & train$cphi_pred > 0.07, ])[1]
FP.cph <- dim(train[as.numeric(train$class) == 0 & train$cphi_pred > 0.07, ])[1]
FN.cph <- dim(train[as.numeric(train$class) == 1 & train$cphi_pred <= 0.07, ])[1]
TN.cph <- dim(train[as.numeric(train$class) == 0 & train$cphi_pred <= 0.07, ])[1]
SN.cph <- TP.cph / (TP.cph + FN.cph)
SN.cph # SE
SP.cph <- TN.cph / (TN.cph + FP.cph)
SP.cph # SP
TP.cph / (TP.cph + FP.cph) # PPV
TN.cph / (TN.cph + FN.cph) # NPV
SN.cph + SP.cph - 1 # YI



TP.lr2 <- dim(train[as.numeric(train$class) == 1 & train$lr2_pred > 0.1, ])[1]
FP.lr2 <- dim(train[as.numeric(train$class) == 0 & train$lr2_pred > 0.1, ])[1]
FN.lr2 <- dim(train[as.numeric(train$class) == 1 & train$lr2_pred <= 0.1, ])[1]
TN.lr2 <- dim(train[as.numeric(train$class) == 0 & train$lr2_pred <= 0.1, ])[1]
SN.lr2 <- TP.lr2 / (TP.lr2 + FN.lr2)
SN.lr2 # SE
SP.lr2 <- TN.lr2 / (TN.lr2 + FP.lr2)
SP.lr2 # SP
TP.lr2 / (TP.lr2 + FP.lr2) # PPV
TN.lr2 / (TN.lr2 + FN.lr2) # NPV
SN.lr2 + SP.lr2 - 1 # YI




library(rms)
nomdata <- train
nomdata$CA125 <- factor(nomdata$ca125_cutoff, levels = c(0, 1), labels = c("≤ 100", "> 100"))
nomdata$HE4 <- factor(nomdata$he4_cutoff, levels = c(0, 1), labels = c("≤ 64.5", "> 64.5"))
nomdata$FIB <- factor(nomdata$fib_cutoff, levels = c(0, 1), labels = c("≤ 3.28", "> 3.28"))
nomdata$MLR <- factor(nomdata$mlr_cutoff, levels = c(0, 1), labels = c("≤ 0.25", "> 0.25"))
nomdata$Bloodflow <- factor(nomdata$bloodflow, levels = c(0, 1), labels = c("No", "Yes"))
nomdata$Solid_areas <- factor(nomdata$solid_areas, levels = c(0, 1), labels = c("No", "Yes"))
dd <- datadist(nomdata)
options(datadist = "dd")
fit3.nom <- lrm(class ~ age + CA125 + HE4 + FIB + MLR
  + Bloodflow + Solid_areas, data = nomdata, x = T, y = T)
nom <- nomogram(fit3.nom,
  fun = plogis,
  maxscale = 100,
  fun.at = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99),
  lp = F,
  funlabel = "Malignant probability",
  conf.int = F
)
plot(nom,
  xfrac = .3,
  cex.var = 1.5,
  cex.axis = 1.2,
  tcl = -0.5,
  lmgp = 0.6,
  label.every = 1,
  naxes = 10,
  col.conf = c("red", "seagreen"),
  points.label = "Points", total.points.label = "Total Points",
  col.grid = gray(c(0.8, 0.95))
)



# 构建数值版本的列线图用于构建网页列线图 ------------------------------------------------------

web_data <- train
web_data$CA125 <- web_data$ca125_cutoff
web_data$HE4 <- web_data$he4_cutoff
web_data$FIB <- web_data$fib_cutoff
web_data$MLR <- web_data$mlr_cutoff
web_data$Bloodflow <- web_data$bloodflow
web_data$Solid_areas <- web_data$solid_areas
dd <- datadist(web_data)
options(datadist = "dd")

lrm_model <- lrm(
  class ~ age + CA125 + HE4 + FIB + MLR
    + Bloodflow + Solid_areas,
  data = web_data,
  x = T, y = T
)

nom_web <- nomogram(
  lrm_model,
  fun = plogis,
  maxscale = 100
)

model_web <- glm(
  class ~ age + CA125 + HE4 + FIB + MLR
    + Bloodflow + Solid_areas,
  data = web_data, family = binomial(), x = TRUE
)


saveRDS(model_web, "web_nomogram/data/ov_sir_model.rds")
saveRDS(nom_web, file = "web_nomogram/data/ov_nomogram.rds")




library(regplot)
regplot(fit3,
  plots = c("density", "boxes"),
  center = TRUE,
  observation = mydata[126, ],
  points = FALSE,
  odds = FALSE,
  nsamp = 1000,
  showP = T,
  rank = "sd",
  droplines = F,
  interval = "confidence",
  clickable = F
)



library(ResourceSelection)
hoslem.test(train$pred1, fitted(fit1), g = 10)
hoslem.test(train$pred2, fitted(fit2), g = 10)
hoslem.test(train$pred3, fitted(fit3), g = 10)

val.prob(train$pred1, train$class)
val.prob(train$pred2, train$class)
val.prob(train$pred3, train$class)



library(rmda)
dca1 <- decision_curve(class ~ pred1,
  data = train,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = 0.95,
  study.design = "case-control",
  population.prevalence = 0.35
)
dca2 <- decision_curve(class ~ pred2,
  data = train,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = 0.95,
  study.design = "case-control",
  population.prevalence = 0.35
)
dca3 <- decision_curve(class ~ pred3,
  data = train,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = 0.95,
  study.design = "case-control",
  population.prevalence = 0.35
)
List.test <- list(dca1, dca2, dca3)
plot_decision_curve(List.test,
  curve.names = c("Model 1", "Model 2", "Model 3"),
  cost.benefit.axis = FALSE,
  col = c(
    rgb(61, 64, 91, maxColorValue = 255),
    rgb(58, 134, 255, maxColorValue = 255),
    rgb(239, 71, 111, maxColorValue = 255)
  ),
  lty = c(1, 1, 1),
  confidence.intervals = F,
  standardize = T
)


# Analyses in the Validation cohort
test$pred1 <- predict(fit1, newdata = test, type = "response")
test$pred2 <- predict(fit2, newdata = test, type = "response")
test$pred3 <- predict(fit3, newdata = test, type = "response")
roc1.test <- roc(test$class, test$pred1)
roc2.test <- roc(test$class, test$pred2)
roc3.test <- roc(test$class, test$pred3)
auc(roc1.test)
auc(roc2.test)
auc(roc3.test)
ci.auc(roc1.test)
ci.auc(roc2.test)
ci.auc(roc3.test)
roc.test(roc1.test, roc2.test, method = "delong")
roc.test(roc1.test, roc3.test, method = "delong")
roc.test(roc2.test, roc3.test, method = "delong")
plot.roc(roc1.test,
  col = rgb(61, 64, 91, maxColorValue = 255),
  lwd = 1.7, lty = 1,
  main = "ROC in TEST",
  legacy.axes = T, xlim = c(1, 0), ylim = c(0, 1), asp = 1
)
plot.roc(roc2.test,
  add = T,
  col = rgb(58, 134, 255, maxColorValue = 255),
  lwd = 1.7, lty = 1
)
plot.roc(roc3.test,
  add = T,
  col = rgb(239, 71, 111, maxColorValue = 255),
  lwd = 2, lty = 1
)
legend(0.5, 0.3,
  c(
    paste("AUC of Model 1: ", sprintf("%.3f", auc(roc1.test))),
    paste("AUC of Model 2: ", sprintf("%.3f", auc(roc2.test))),
    paste("AUC of Model 3: ", sprintf("%.3f", auc(roc3.test)))
  ),
  col = c(
    rgb(61, 64, 91, maxColorValue = 255),
    rgb(58, 134, 255, maxColorValue = 255),
    rgb(239, 71, 111, maxColorValue = 255)
  ),
  lty = c(1, 1, 1),
  lwd = c(1.7, 1.7, 2),
  box.lty = 2, cex = 0.5, horiz = F, seg.len = 5
)




nribin(
  event = test$class,
  p.std = test$pred1,
  p.new = test$pred2,
  updown = "category",
  cut = 0.435,
  niter = 1000
)
reclassification(
  data = test,
  cOutcome = 5,
  predrisk1 = test$pred1,
  predrisk2 = test$pred2,
  cutoff = c(0, 0.435, 1)
)

nribin(
  event = test$class,
  p.std = test$pred1,
  p.new = test$pred3,
  updown = "category",
  cut = 0.18,
  niter = 1000
)
reclassification(
  data = test,
  cOutcome = 5,
  predrisk1 = test$pred1,
  predrisk2 = test$pred3,
  cutoff = c(0, 0.18, 1)
)

nribin(
  event = test$class,
  p.std = test$pred2,
  p.new = test$pred3,
  updown = "category",
  cut = 0.18,
  niter = 1000
)
reclassification(
  data = test,
  cOutcome = 5,
  predrisk1 = test$pred2,
  predrisk2 = test$pred3,
  cutoff = c(0, 0.18, 1)
)
write.csv(test, file = "test.csv")




roc.roma.test <- roc(test$class, test$roma_pred)
roc.cph.test <- roc(test$class, test$cphi_pred)
roc.rmi.test <- roc(test$class, test$rmi)
roc.lr2.test <- roc(test$class, test$lr2_pred)
roc.rops.test <- roc(test$class, test$rops)
auc(roc.rops.test)
auc(roc.roma.test)
auc(roc.cph.test)
auc(roc.rmi.test)
auc(roc.lr2.test)
ci.auc(roc.rops.test)
ci.auc(roc.roma.test)
ci.auc(roc.cph.test)
ci.auc(roc.rmi.test)
ci.auc(roc.lr2.test)
roc.test(roc3.test, roc.rops.test, method = "delong")
roc.test(roc3.test, roc.roma.test, method = "delong")
roc.test(roc3.test, roc.cph.test, method = "delong")
roc.test(roc3.test, roc.rmi.test, method = "delong")
roc.test(roc3.test, roc.lr2.test, method = "delong")
roc.test(roc.rmi.test, roc.lr2.test, method = "delong")
plot.roc(roc3.test,
  col = rgb(239, 71, 111, maxColorValue = 255), lwd = 2, lty = 1,
  main = "ROC in TEST",
  legacy.axes = T, xlim = c(1, 0), ylim = c(0, 1), asp = 1
)
plot.roc(roc.rops.test,
  add = T,
  col = rgb(58, 134, 255, maxColorValue = 255), lwd = 1.7, lty = 1
)
plot.roc(roc.roma.test,
  add = T,
  col = rgb(61, 64, 91, maxColorValue = 255), lwd = 1.7, lty = 2
)
plot.roc(roc.cph.test, add = T, col = "orange", lwd = 1.7, lty = 1)
plot.roc(roc.rmi.test,
  add = T,
  col = rgb(129, 178, 154, maxColorValue = 255), lwd = 1.7, lty = 1
)
plot.roc(roc.lr2.test,
  add = T,
  col = rgb(61, 64, 91, maxColorValue = 255), lwd = 1.7, lty = 1
)
legend(0.5, 0.3,
  c(
    paste("AUC of Model 3: ", sprintf("%.3f", auc(roc3.test))),
    paste("AUC of R-OPS: ", sprintf("%.3f", auc(roc.rops.test))),
    paste("AUC of ROMA: ", sprintf("%.3f", auc(roc.roma.test))),
    paste("AUC of CPH-I: ", sprintf("%.3f", auc(roc.cph.test))),
    paste("AUC of RMI 4: ", sprintf("%.3f", auc(roc.rmi.test))),
    paste("AUC of LR2: ", sprintf("%.3f", auc(roc.lr2.test)))
  ),
  col = c(rgb(239, 71, 111, maxColorValue = 255),
    rgb(58, 134, 255, maxColorValue = 255),
    rgb(61, 64, 91, maxColorValue = 255),
    col = "orange",
    rgb(129, 178, 154, maxColorValue = 255),
    rgb(61, 64, 91, maxColorValue = 255)
  ),
  lty = c(1, 1, 2, 1, 1, 1),
  lwd = c(2, 1.7, 1.7, 1.7, 1.7, 1.7),
  box.lty = 2, cex = 0.5, horiz = F, seg.len = 5
)




TP.fit3.test <- dim(test[as.numeric(test$class) == 1 & test$pred3 > 0.18, ])[1]
FP.fit3.test <- dim(test[as.numeric(test$class) == 0 & test$pred3 > 0.18, ])[1]
FN.fit3.test <- dim(test[as.numeric(test$class) == 1 & test$pred3 <= 0.18, ])[1]
TN.fit3.test <- dim(test[as.numeric(test$class) == 0 & test$pred3 <= 0.18, ])[1]
SN.fit3.test <- TP.fit3.test / (TP.fit3.test + FN.fit3.test)
SN.fit3.test # SE
SP.fit3.test <- TN.fit3.test / (TN.fit3.test + FP.fit3.test)
SP.fit3.test # SP
TP.fit3.test / (TP.fit3.test + FP.fit3.test) # PPV
TN.fit3.test / (TN.fit3.test + FN.fit3.test) # NPV
SN.fit3.test + SP.fit3.test - 1 # YI


TP.rops.test <- dim(test[as.numeric(test$class) == 1 & test$rops > 330, ])[1]
FP.rops.test <- dim(test[as.numeric(test$class) == 0 & test$rops > 330, ])[1]
FN.rops.test <- dim(test[as.numeric(test$class) == 1 & test$rops <= 330, ])[1]
TN.rops.test <- dim(test[as.numeric(test$class) == 0 & test$rops <= 330, ])[1]
SN.rops.test <- TP.rops.test / (TP.rops.test + FN.rops.test)
SN.rops.test # SE
SP.rops.test <- TN.rops.test / (TN.rops.test + FP.rops.test)
SP.rops.test # SP
TP.rops.test / (TP.rops.test + FP.rops.test) # PPV
TN.rops.test / (TN.rops.test + FN.rops.test) # NPV
SN.rops.test + SP.rops.test - 1 # YI

TP.cph.test <- dim(test[as.numeric(test$class) == 1 & test$cphi_pred > 0.07, ])[1]
FP.cph.test <- dim(test[as.numeric(test$class) == 0 & test$cphi_pred > 0.07, ])[1]
FN.cph.test <- dim(test[as.numeric(test$class) == 1 & test$cphi_pred <= 0.07, ])[1]
TN.cph.test <- dim(test[as.numeric(test$class) == 0 & test$cphi_pred <= 0.07, ])[1]
SN.cph.test <- TP.cph.test / (TP.cph.test + FN.cph.test)
SN.cph.test # SE
SP.cph.test <- TN.cph.test / (TN.cph.test + FP.cph.test)
SP.cph.test # SP
TP.cph.test / (TP.cph.test + FP.cph.test) # PPV
TN.cph.test / (TN.cph.test + FN.cph.test) # NPV
(TP.cph.test + TN.cph.test) / (TP.cph.test + TN.cph.test + FP.cph.test + FN.cph.test) # 诊断准确度
SN.cph.test + SP.cph.test - 1 # YI

TP.rmi.test <- dim(test[as.numeric(test$class) == 1 & test$rmi > 200, ])[1]
FP.rmi.test <- dim(test[as.numeric(test$class) == 0 & test$rmi > 200, ])[1]
FN.rmi.test <- dim(test[as.numeric(test$class) == 1 & test$rmi <= 200, ])[1]
TN.rmi.test <- dim(test[as.numeric(test$class) == 0 & test$rmi <= 200, ])[1]
SN.rmi.test <- TP.rmi.test / (TP.rmi.test + FN.rmi.test)
SN.rmi.test # SE
SP.rmi.test <- TN.rmi.test / (TN.rmi.test + FP.rmi.test)
SP.rmi.test # SP
TP.rmi.test / (TP.rmi.test + FP.rmi.test) # PPV
TN.rmi.test / (TN.rmi.test + FN.rmi.test) # NPV
(TP.rmi.test + TN.rmi.test) / (TP.rmi.test + TN.rmi.test + FP.rmi.test + FN.rmi.test) # 诊断准确度
SN.rmi.test + SP.rmi.test - 1 # YI

TP.lr2.test <- dim(test[as.numeric(test$class) == 1 & test$lr2_pred > 0.1, ])[1]
FP.lr2.test <- dim(test[as.numeric(test$class) == 0 & test$lr2_pred > 0.1, ])[1]
FN.lr2.test <- dim(test[as.numeric(test$class) == 1 & test$lr2_pred <= 0.1, ])[1]
TN.lr2.test <- dim(test[as.numeric(test$class) == 0 & test$lr2_pred <= 0.1, ])[1]
SN.lr2.test <- TP.lr2.test / (TP.lr2.test + FN.lr2.test)
SN.lr2.test # SE
SP.lr2.test <- TN.lr2.test / (TN.lr2.test + FP.lr2.test)
SP.lr2.test # SP
TP.lr2.test / (TP.lr2.test + FP.lr2.test) # PPV
TN.lr2.test / (TN.lr2.test + FN.lr2.test) # NPV
(TP.lr2.test + TN.lr2.test) / (TP.lr2.test + TN.lr2.test + FP.lr2.test + FN.lr2.test) # 诊断准确度
SN.lr2.test + SP.lr2.test - 1 # YI



fit1.test <- glm(class ~ pred1, data = test, family = binomial(), x = TRUE)
fit2.test <- glm(class ~ pred2, data = test, family = binomial(), x = TRUE)
fit3.test <- glm(class ~ pred3, data = test, family = binomial(), x = TRUE)
hoslem.test(test$pred1, fitted(fit1.test), g = 10)
hoslem.test(test$pred2, fitted(fit2.test), g = 10)
hoslem.test(test$pred3, fitted(fit3.test), g = 10)

val.prob(test$pred1, test$class)
val.prob(test$pred2, test$class)
val.prob(test$pred3, test$class)



dca1.test <- decision_curve(class ~ pred1, # 定义结局变量和自变量
  data = test,
  family = binomial(link = "logit"), # 定义分布族
  thresholds = seq(0, 1, by = 0.01), # 设置横坐标阈概率的范围和刻度，—般是0~1；但如果有某种具体情况，大家一致认为阈概率达到某个值以上，比如 40%，则必须采取干预措施，那么0.4以后的研究就没什么意义了，可以设为0~0.4
  confidence.intervals = 0.95,
  study.design = "case-control", # 指定研究类型，设置为'case-control'时需要指定患病率，"cohort"无需指定患病率
  population.prevalence = 0.35
) # 指定患病率
dca2.test <- decision_curve(class ~ pred2,
  data = test,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = 0.95,
  study.design = "case-control",
  population.prevalence = 0.35
)
dca3.test <- decision_curve(class ~ pred3,
  data = test,
  family = binomial(link = "logit"),
  thresholds = seq(0, 1, by = 0.01),
  confidence.intervals = 0.95,
  study.design = "case-control",
  population.prevalence = 0.35
)
List.test <- list(dca1.test, dca2.test, dca3.test)
plot_decision_curve(List.test,
  curve.names = c("Model 1", "Model 2", "Model 3"),
  cost.benefit.axis = FALSE,
  col = c(
    rgb(61, 64, 91, maxColorValue = 255),
    rgb(58, 134, 255, maxColorValue = 255),
    rgb(239, 71, 111, maxColorValue = 255)
  ),
  lty = c(1, 1, 1),
  confidence.intervals = F,
  standardize = T
)



# Subgroup analyses for early-stage ovarian cancer
early_stage_data <- subset(mydata, figo_stage == 1 | class == 0)
table(early_stage_data$class)

early_stage_data$pred1 <- predict(fit1, newdata = early_stage_data, type = "response")
early_stage_data$pred2 <- predict(fit2, newdata = early_stage_data, type = "response")
early_stage_data$pred3 <- predict(fit3, newdata = early_stage_data, type = "response")

roc1.early <- roc(early_stage_data$class, early_stage_data$pred1)
roc2.early <- roc(early_stage_data$class, early_stage_data$pred2)
roc3.early <- roc(early_stage_data$class, early_stage_data$pred3)
roc.roma.early <- roc(early_stage_data$class, early_stage_data$roma_pred)
roc.cph.early <- roc(early_stage_data$class, early_stage_data$cphi_pred)
roc.rmi.early <- roc(early_stage_data$class, early_stage_data$rmi)
roc.lr2.early <- roc(early_stage_data$class, early_stage_data$lr2_pred)
roc.rops.early <- roc(early_stage_data$class, early_stage_data$rops)
auc(roc1.early)
auc(roc2.early)
auc(roc3.early)
ci.auc(roc1.early)
ci.auc(roc2.early)
ci.auc(roc3.early)
roc.test(roc1.early, roc2.early, method = "delong")
roc.test(roc1.early, roc3.early, method = "delong")
roc.test(roc2.early, roc3.early, method = "delong")
plot.roc(roc1.early,
  col = rgb(61, 64, 91, maxColorValue = 255),
  lwd = 1.7, lty = 1,
  main = "ROC of Early Stage",
  legacy.axes = T, xlim = c(1, 0), ylim = c(0, 1), asp = 1
)
plot.roc(roc2.early,
  add = T,
  col = rgb(58, 134, 255, maxColorValue = 255),
  lwd = 1.7, lty = 1
)
plot.roc(roc3.early,
  add = T,
  col = rgb(239, 71, 111, maxColorValue = 255),
  lwd = 2, lty = 1
)
legend(0.5, 0.3,
  c(
    paste("AUC of Nomogram: ", sprintf("%.3f", auc(roc1.early))),
    paste("AUC of Model 2: ", sprintf("%.3f", auc(roc2.early))),
    paste("AUC of Model 3: ", sprintf("%.3f", auc(roc3.early)))
  ),
  col = c(
    rgb(61, 64, 91, maxColorValue = 255),
    rgb(58, 134, 255, maxColorValue = 255),
    rgb(239, 71, 111, maxColorValue = 255)
  ),
  lty = c(1, 1, 1),
  lwd = c(1.7, 1.7, 2),
  box.lty = 2, cex = 0.5, horiz = F, seg.len = 5
)



TP.fit3.early <- dim(early_stage_data[as.numeric(early_stage_data$class) == 1 & early_stage_data$pred3 > 0.18, ])[1]
FP.fit3.early <- dim(early_stage_data[as.numeric(early_stage_data$class) == 0 & early_stage_data$pred3 > 0.18, ])[1]
FN.fit3.early <- dim(early_stage_data[as.numeric(early_stage_data$class) == 1 & early_stage_data$pred3 <= 0.18, ])[1]
TN.fit3.early <- dim(early_stage_data[as.numeric(early_stage_data$class) == 0 & early_stage_data$pred3 <= 0.18, ])[1]
SN.fit3.early <- TP.fit3.early / (TP.fit3.early + FN.fit3.early)
SN.fit3.early # SE
SP.fit3.early <- TN.fit3.early / (TN.fit3.early + FP.fit3.early)
SP.fit3.early # SP


reclassification(
  data = early_stage_data,
  cOutcome = 5,
  predrisk1 = early_stage_data$pred1,
  predrisk2 = early_stage_data$pred2,
  cutoff = c(0, 0.393, 1)
)

reclassification(
  data = early_stage_data,
  cOutcome = 5,
  predrisk1 = early_stage_data$pred1,
  predrisk2 = early_stage_data$pred3,
  cutoff = c(0, 0.18, 1)
)

reclassification(
  data = early_stage_data,
  cOutcome = 5,
  predrisk1 = early_stage_data$pred2,
  predrisk2 = early_stage_data$pred3,
  cutoff = c(0, 0.18, 1)
)




auc(roc.rops.early)
auc(roc.roma.early)
auc(roc.cph.early)
auc(roc.rmi.early)
auc(roc.lr2.early)
ci.auc(roc.rops.early)
ci.auc(roc.roma.early)
ci.auc(roc.cph.early)
ci.auc(roc.rmi.early)
ci.auc(roc.lr2.early)
roc.test(roc3.early, roc.rops.early, method = "delong")
roc.test(roc3.early, roc.roma.early, method = "delong")
roc.test(roc3.early, roc.cph.early, method = "delong")
roc.test(roc3.early, roc.rmi.early, method = "delong")
roc.test(roc3.early, roc.lr2.early, method = "delong")
roc.test(roc.rmi.early, roc.lr2.early, method = "delong")
plot.roc(roc3.early,
  col = rgb(239, 71, 111, maxColorValue = 255), lwd = 2, lty = 1,
  main = "ROC of Early Stage",
  legacy.axes = T, xlim = c(1, 0), ylim = c(0, 1), asp = 1
)
plot.roc(roc.rops.early,
  add = T,
  col = rgb(58, 134, 255, maxColorValue = 255), lwd = 1.7, lty = 1
)
plot.roc(roc.roma.early,
  add = T,
  col = rgb(61, 64, 91, maxColorValue = 255), lwd = 1.7, lty = 2
)
plot.roc(roc.cph.early, add = T, col = "orange", lwd = 1.7, lty = 1)
plot.roc(roc.rmi.early,
  add = T,
  col = rgb(129, 178, 154, maxColorValue = 255), lwd = 1.7, lty = 1
)
plot.roc(roc.lr2.early,
  add = T,
  col = rgb(61, 64, 91, maxColorValue = 255), lwd = 1.7, lty = 1
)
legend(0.5, 0.3,
  c(
    paste("AUC of Model 3: ", sprintf("%.3f", auc(roc3.early))),
    paste("AUC of R-OPS: ", sprintf("%.3f", auc(roc.rops.early))),
    paste("AUC of ROMA: ", sprintf("%.3f", auc(roc.roma.early))),
    paste("AUC of CPH-I: ", sprintf("%.3f", auc(roc.cph.early))),
    paste("AUC of RMI 4: ", sprintf("%.3f", auc(roc.rmi.early))),
    paste("AUC of LR2: ", sprintf("%.3f", auc(roc.lr2.early)))
  ),
  col = c(rgb(239, 71, 111, maxColorValue = 255),
    rgb(58, 134, 255, maxColorValue = 255),
    rgb(61, 64, 91, maxColorValue = 255),
    col = "orange",
    rgb(129, 178, 154, maxColorValue = 255),
    rgb(61, 64, 91, maxColorValue = 255)
  ),
  lty = c(1, 1, 2, 1, 1, 1),
  lwd = c(2, 1.7, 1.7, 1.7, 1.7, 1.7),
  box.lty = 2, cex = 0.5, horiz = F, seg.len = 5
)


# R packages citations
citation()
citation("glmnet")
citation("car")
citation("forestplot")
citation("pROC")
citation("nricens")
citation("PredictABEL")
citation("ResourceSelection")
citation("rmda")
citation("rms")
citation("regplot")
citation("DynNom")
