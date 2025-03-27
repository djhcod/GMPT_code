#🔴导入原始数据
original_data <- read.csv("sarcoma_peri/data/Uterine_cancer_2000_2019.csv")
chr0 <- paste0("原始数据有：", nrow(original_data), "行")
colnames(original_data)

#🔴转换诊断年份
range(original_data$Year.of.diagnosis)
mydata1 <- subset(original_data, Year.of.diagnosis >= 2010 & Year.of.diagnosis <= 2016)
chr1 <- paste0("限定诊断年限后：排除", nrow(original_data)-nrow(mydata1),
               "行，剩余", nrow(mydata1), "行")
mydata1$year <- as.factor(mydata1$Year.of.diagnosis)
summary(mydata1$year)

mydata1$year2 <- ifelse(mydata1$year == "2010" | mydata1$year == "2011", "2010-2011", 
                       ifelse(mydata1$year == "2012" | mydata1$year == "2013", "2012-2013", "2014-2016"))
mydata1$year2 <- factor(mydata1$year2)
summary(mydata1$year2)

#🔴转换年龄
unique(mydata1$Age.recode.with.single.ages.and.100.)
mydata1$age <- gsub(" years", "", mydata1$Age.recode.with.single.ages.and.100.)
mydata1$age <- as.numeric(mydata1$age)
mydata2 <- subset(mydata1, age >= 18)
chr2 <- paste0("限定年龄后：排除", nrow(mydata1)-nrow(mydata2),
               "行，剩余", nrow(mydata2), "行")
summary(mydata2$age)

#🔴转换种族
unique(mydata2$Race.recode..White..Black..Other.)
mydata2$race <- ifelse(mydata2$Race.recode..White..Black..Other. == "Unknown", NA, 
                       ifelse(mydata2$Race.recode..White..Black..Other. == 
                              "Other (American Indian/AK Native, Asian/Pacific Islander)",
                              "Others", mydata2$Race.recode..White..Black..Other.))
mydata2$race <- as.factor(mydata2$race)
mydata2$race <- relevel(mydata2$race, ref = "White")
summary(mydata2$race)

#🔴转换婚姻状态
unique(mydata2$Marital.status.at.diagnosis)
mydata2$marriage <- ifelse(mydata2$Marital.status.at.diagnosis == "Unknown", NA,
                           ifelse(mydata2$Marital.status.at.diagnosis == 
                                  "Married (including common law)",
                                  "Married", 
                                  ifelse(mydata2$Marital.status.at.diagnosis == 
                                         "Single (never married)" |
                                         mydata2$Marital.status.at.diagnosis == 
                                         "Unmarried or Domestic Partner",
                                         "Single/Unmarried", "Others")))
mydata2$marriage <- factor(mydata2$marriage, levels = c("Married", "Single/Unmarried", "Others"))
summary(mydata2$marriage)

#🔴转换肿瘤分级
unique(mydata2$Grade..thru.2017.)
mydata2. <- subset(mydata2, Grade..thru.2017. != "B-cell; pre-B; B-precursor")
mydata2.$grade <- ifelse(mydata2.$Grade..thru.2017. == "Unknown", NA, mydata2.$Grade..thru.2017.)
mydata2.$grade <- as.factor(mydata2.$grade)
mydata2.$grade <- relevel(mydata2.$grade, ref = "Well differentiated; Grade I")
summary(mydata2.$grade)

mydata2.$grade2 <- ifelse(is.na(mydata2.$grade), "Gx",
                         ifelse(mydata2.$grade == "Well differentiated; Grade I" | 
                                 mydata2.$grade == "Moderately differentiated; Grade II", 
                               "Low-grade", "High-grade"))
mydata2.$grade2 <- factor(mydata2.$grade2, levels = c("Low-grade", "High-grade", "Gx"))
summary(mydata2.$grade2)

#🔴转换肿瘤大小
unique(mydata2.$CS.tumor.size..2004.2015.)
unique(mydata2.$Tumor.Size.Summary..2016..)
mydata2.$tumor_size <- ifelse(mydata2.$Year.of.diagnosis <= 2015, 
                             mydata2.$CS.tumor.size..2004.2015.,
                             mydata2.$Tumor.Size.Summary..2016..)
mydata2.$tumor_size <- as.numeric(mydata2.$tumor_size)
range(mydata2.$tumor_size)
mydata2.$tumor_size <- ifelse(mydata2.$tumor_size >=1 & mydata2.$tumor_size <= 988, 
                             mydata2.$tumor_size, NA)
summary(mydata2.$tumor_size)

#🔴限定阳性组织学诊断
unique(mydata2.$Diagnostic.Confirmation)
mydata3 <- subset(mydata2., Diagnostic.Confirmation == "Positive histology" |
                           Diagnostic.Confirmation == 
                           "Pos hist AND immunophenotyping AND/OR pos genetic studies")
chr3 <- paste0("限定阳性组织学诊断后：排除", nrow(mydata2)-nrow(mydata3),
               "行，剩余", nrow(mydata3), "行")

#🔴限定组织学类型
unique(mydata3$ICD.O.3.Hist.behav)
#LMS, 8890/3, 8891/3, 8896/3; 
#ESS, 8930/3, 8931/3, 8935/3
mydata4 <- mydata3[grep("^889[016]|^893[015]", mydata3$ICD.O.3.Hist.behav), ]
chr4 <- paste0("限定组织学类型后：排除", nrow(mydata3)-nrow(mydata4),
               "行，剩余", nrow(mydata4), "行")
unique(mydata4$ICD.O.3.Hist.behav)
mydata4$his <- ifelse(mydata4$ICD.O.3.Hist.behav == 
                        "8930/3: Endometrial stromal sarcoma, NOS" | 
                        mydata4$ICD.O.3.Hist.behav == 
                        "8931/3: Endometrial stromal sarcoma, low grade" |
                        mydata4$ICD.O.3.Hist.behav ==
                        "8935/3: Stromal sarcoma, NOS", "ESS", "LMS")
mydata4$his <- factor(mydata4$his, levels = c("LMS", "ESS"))
summary(mydata4$his)

#🔴限定无恶性肿瘤病史
unique(mydata4$First.malignant.primary.indicator)
mydata5 <- subset(mydata4, First.malignant.primary.indicator == "Yes")
chr5 <- paste0("限定无恶性肿瘤病史后：排除", nrow(mydata4)-nrow(mydata5),
               "行，剩余", nrow(mydata5), "行")

#🔴转换T分期
unique(mydata5$Derived.AJCC.T..7th.ed..2010.2015.)
unique(mydata5$Derived.SEER.Combined.T..2016.2017.)
mydata5$T_stage <- ifelse(mydata5$Year.of.diagnosis <= 2015, 
                          mydata5$Derived.AJCC.T..7th.ed..2010.2015.,
                          mydata5$Derived.SEER.Combined.T..2016.2017.)
unique(mydata5$T_stage)
mydata5$T_stage <- ifelse(mydata5$T_stage == "Blank(s)" , NA,
                          substr(mydata5$T_stage, 1, 2))
mydata5$T_stage <- gsub("p|c", "T", mydata5$T_stage)
mydata5$T_stage <- ifelse(mydata5$T_stage == "TX", NA, mydata5$T_stage)
mydata5$T_stage <- as.factor(mydata5$T_stage)
summary(mydata5$T_stage)

#🔴转换详细T分期
mydata5$T_stage_plus <- ifelse(mydata5$Year.of.diagnosis <= 2015, 
                               mydata5$Derived.AJCC.T..7th.ed..2010.2015.,
                               mydata5$Derived.SEER.Combined.T..2016.2017.)
unique(mydata5$T_stage_plus)
mydata5$T_stage_plus <- ifelse(mydata5$T_stage_plus == "Blank(s)", NA, 
                               mydata5$T_stage_plus)
mydata5$T_stage_plus <- gsub("NOS", "", mydata5$T_stage_plus)
mydata5$T_stage_plus <- tolower(mydata5$T_stage_plus)
mydata5$T_stage_plus <- gsub("^(t|p|c)", "T", mydata5$T_stage_plus)
mydata5$T_stage_plus <- ifelse(mydata5$T_stage_plus == "Tx", NA, mydata5$T_stage_plus)
mydata5$T_stage_plus <- as.factor(mydata5$T_stage_plus)
summary(mydata5$T_stage_plus)

#🔴转换N分期
unique(mydata5$Derived.AJCC.N..7th.ed..2010.2015.)
unique(mydata5$Derived.SEER.Combined.N..2016.2017.)
mydata5$N_stage <- ifelse(mydata5$Year.of.diagnosis <= 2015, 
                          mydata5$Derived.AJCC.N..7th.ed..2010.2015., 
                          mydata5$Derived.SEER.Combined.N..2016.2017.)
unique(mydata5$N_stage)
mydata5$N_stage <- ifelse(mydata5$N_stage == "Blank(s)", NA, 
                          gsub("^(p|c)", "N", mydata5$N_stage))
mydata5$N_stage <- ifelse(mydata5$N_stage == "NX", NA, mydata5$N_stage)
mydata5$N_stage <- as.factor(mydata5$N_stage)
summary(mydata5$N_stage)

#🔴转换M分期
unique(mydata5$Derived.AJCC.M..7th.ed..2010.2015.)
unique(mydata5$Derived.SEER.Combined.M..2016.2017.)
mydata5$M_stage <- ifelse(mydata5$Year.of.diagnosis <= 2015, 
                          mydata5$Derived.AJCC.M..7th.ed..2010.2015.,
                          mydata5$Derived.SEER.Combined.M..2016.2017.)
unique(mydata5$M_stage)
mydata5$M_stage <- ifelse(mydata5$M_stage == "Blank(s)", NA, 
                          gsub("^(p|c)", "M", mydata5$M_stage))
mydata5$M_stage <- as.factor(mydata5$M_stage)
summary(mydata5$M_stage)

#🔴TNM转FIGO
mydata5$figo <- ifelse(mydata5$T_stage == "T1" & 
                         mydata5$N_stage == "N0" & 
                         mydata5$M_stage == "M0", "I",
                       ifelse(mydata5$T_stage == "T2" & 
                                mydata5$N_stage == "N0" & 
                                mydata5$M_stage == "M0", "II",
                              ifelse(mydata5$T_stage == "T3" & 
                                       mydata5$M_stage == "M0", "III",
                                     ifelse(mydata5$T_stage != "T4" & 
                                              mydata5$N_stage == "N1" & 
                                              mydata5$M_stage == "M0", "III", 
                                            ifelse(mydata5$T_stage == "T4", "IV", 
                                                   ifelse(mydata5$M_stage == "M1", 
                                                          "IV", NA))))))
mydata5$figo <- as.factor(mydata5$figo)
summary(mydata5$figo)

#🔴TNM转详细FIGO
if(T){
mydata5$figo_plus <- ifelse(mydata5$T_stage_plus == "T1a" & 
                              mydata5$N_stage == "N0" & 
                              mydata5$M_stage == "M0", "IA",
                            ifelse(mydata5$T_stage_plus == "T1b" & 
                                     mydata5$N_stage == "N0" & 
                                     mydata5$M_stage == "M0", "IB",
                                   ifelse(mydata5$T_stage_plus == "T1c" & 
                                            mydata5$N_stage == "N0" & 
                                            mydata5$M_stage == "M0", "IC (for AS)",
                                          ifelse(mydata5$T_stage_plus == "T3a" & 
                                                   mydata5$N_stage == "N0" & 
                                                   mydata5$M_stage == "M0", "IIIA",
                                                 ifelse(mydata5$T_stage_plus == "T3b" &
                                                          mydata5$N_stage == "N0" & 
                                                          mydata5$M_stage == "M0", "IIIB", 
                                                        ifelse(mydata5$T_stage_plus != "T4" &
                                                                 mydata5$N_stage == "N1" & 
                                                                 mydata5$M_stage == "M0", "IIIC", 
                                                               ifelse(mydata5$T_stage_plus == "T4" &
                                                                        mydata5$M_stage == "M0", "IVA", 
                                                                      ifelse(mydata5$M_stage == "M1", "IVB", NA))))))))
unique(mydata5$figo_plus)
mydata5$figo_plus <- as.factor(mydata5$figo_plus)
summary(mydata5$figo_plus)
}


# 🔴排除分期未知 ---------------------------------------------------------------
mydata6 <- mydata5[is.na(mydata5$figo) == F,]
chr6 <- paste0("排除分期未知后：排除", nrow(mydata5)-nrow(mydata6),
               "行，剩余", nrow(mydata6), "行")

#🔴转换腹膜细胞学记录
unique(mydata6$Peritoneal.Cytology.Recode..2010..)
mydata6$peri <- ifelse(mydata6$Peritoneal.Cytology.Recode..2010.. == 
                         "Peritoneal cytology/washing negative for malignancy", 0,
                       ifelse(mydata6$Peritoneal.Cytology.Recode..2010.. == 
                                "Peritoneal cytology/washing malignant (positive for malignancy)", 
                              1, NA))
mydata6$peri <- factor(mydata6$peri, levels = c(0, 1), labels = c("Negtive", "Malignant"))
summary(mydata6$peri)
mydata7 <- subset(mydata6, is.na(peri) == F)
chr7 <- paste0("排除腹膜细胞学未知后：排除", nrow(mydata6)-nrow(mydata7),
               "行，剩余", nrow(mydata7), "行")

#🔴转换手术方式
table(mydata7$RX.Summ..Surg.Prim.Site..1998..)
mydata8 <- subset(mydata7, RX.Summ..Surg.Prim.Site..1998.. >= 40 & 
                    RX.Summ..Surg.Prim.Site..1998.. < 65)
chr8 <- paste0("限定手术方式后：排除", nrow(mydata7)-nrow(mydata8),
               "行，剩余", nrow(mydata8), "行")
mydata8$surg <- ifelse(mydata8$RX.Summ..Surg.Prim.Site..1998.. == 40, "TH", 
                       ifelse(mydata8$RX.Summ..Surg.Prim.Site..1998.. == 50, "TH+BSO", "RH/EH"))
mydata8$surg <- factor(mydata8$surg, levels = c("TH+BSO", "TH", "RH/EH"))
summary(mydata8$surg)


#转换淋巴结切除术，from "Regional.nodes.examined..1988.."
if(F) {
unique(mydata8$Regional.nodes.examined..1988..)
mydata8$lnd <- ifelse(mydata8$Regional.nodes.examined..1988.. <= 3 |
                        mydata8$Regional.nodes.examined..1988.. == 95 |
                        mydata8$Regional.nodes.examined..1988.. == 96, 0, 
                      ifelse(mydata8$Regional.nodes.examined..1988.. > 3 &
                               mydata8$Regional.nodes.examined..1988.. <= 90, 1,
                             ifelse(mydata8$Regional.nodes.examined..1988.. == 97, 
                                    1, NA)))
mydata8$lnd <- factor(mydata8$lnd, levels = c(0, 1), labels = c("No", "Yes"))
summary(mydata8$lnd)
}

#🔴转换腹主动脉旁淋巴结切除，https://seer.cancer.gov/seerstat/variables/seer/ssdi/
unique(mydata8$Number.of.Examined.Para.Aortic.Nodes.Recode..2010..)
alnd <- vector()
for (i in 1:nrow(mydata8)) {
  if(mydata8$Number.of.Examined.Para.Aortic.Nodes.Recode..2010..[i] == 
     "No para-aortic nodes examined" | 
     mydata8$Number.of.Examined.Para.Aortic.Nodes.Recode..2010..[i] == 
     "No para-aortic lymph nodes removed, but aspiration or core biopsy only") {
    alnd <- rbind(alnd, 0) 
    } else if(mydata8$Number.of.Examined.Para.Aortic.Nodes.Recode..2010..[i] == 
              "Para-aortic nodes examined, number unknown" | 
              mydata8$Number.of.Examined.Para.Aortic.Nodes.Recode..2010..[i] == 
              "Not documented; Indeterminate; Not assessed or unknown if assessed") {
      alnd <- rbind(alnd, NA)
      } else if(as.numeric(mydata8$Number.of.Examined.Para.Aortic.Nodes.Recode..2010..[i]) <= 3) {
        alnd <- rbind(alnd, 0)
        } else if(as.numeric(mydata8$Number.of.Examined.Para.Aortic.Nodes.Recode..2010..[i]) > 3) {
          alnd <- rbind(alnd, 1)
          } else {
            alnd <- rbind(alnd, NA)
          }
}
mydata8 <- cbind(mydata8, alnd)
mydata8$alnd <- factor(mydata8$alnd, levels = c(0, 1), labels = c("No", "Yes"))
summary(mydata8$alnd)

#🔴转换盆腔淋巴结切除
unique(mydata8$Number.of.Examined.Pelvic.Nodes.Recode..2010..)
plnd <- vector()
for (i in 1:nrow(mydata8)) {
  if(mydata8$Number.of.Examined.Pelvic.Nodes.Recode..2010..[i] == 
     "No pelvic lymph nodes examined" | 
     mydata8$Number.of.Examined.Pelvic.Nodes.Recode..2010..[i] == 
     "No pelvic lymph nodes removed, but aspiration or core biopsy only") {
    plnd <- rbind(plnd, 0) 
  } else if(mydata8$Number.of.Examined.Pelvic.Nodes.Recode..2010..[i] == 
            "Pelvic nodes examined, number unknown" | 
            mydata8$Number.of.Examined.Pelvic.Nodes.Recode..2010..[i] == 
            "Not documented; Indeterminate; Not assessed or unknown if assessed") {
    plnd <- rbind(plnd, NA)
  } else if(as.numeric(mydata8$Number.of.Examined.Pelvic.Nodes.Recode..2010..[i]) <= 3) {
    plnd <- rbind(plnd, 0)
  } else if(as.numeric(mydata8$Number.of.Examined.Pelvic.Nodes.Recode..2010..[i]) > 3) {
    plnd <- rbind(plnd, 1)
  } else {
    plnd <- rbind(plnd, NA)
  }
}
mydata8 <- cbind(mydata8, plnd)
mydata8$plnd <- factor(mydata8$plnd, levels = c(0, 1), labels = c("No", "Yes"))
summary(mydata8$plnd)

#🔴合并淋巴结切除
mydata8$lnd <- ifelse(mydata8$alnd == "No" & mydata8$plnd == "No", 0, 
                      ifelse(mydata8$alnd == "Yes" | mydata8$plnd == "Yes", 1, NA))
mydata8$lnd <- factor(mydata8$lnd, levels = c(0, 1), labels = c("No", "Yes"))
summary(mydata8$lnd)

#🔴转换放疗
table(mydata8$Radiation.recode)
table(mydata8$RX.Summ..Surg.Rad.Seq)
mydata8$rad <- ifelse(mydata8$Radiation.recode == "None/Unknown" & 
                        mydata8$RX.Summ..Surg.Rad.Seq == 
                        "No radiation and/or cancer-directed surgery", 0,
                      ifelse(mydata8$Radiation.recode == "Refused (1988+)", 0,
                             ifelse(mydata8$Radiation.recode == 
                                      "Recommended, unknown if administered", NA, 1)))
unique(mydata8[mydata8$rad == 1, ]$Radiation.recode)
unique(mydata8[mydata8$rad == 1, ]$RX.Summ..Surg.Rad.Seq)
unique(mydata8[mydata8$rad == 0, ]$Radiation.recode)
mydata8$rad <- factor(mydata8$rad, levels = c(0, 1), labels = c("No", "Yes"))
summary(mydata8$rad)

#🔴转换化疗
unique(mydata8$Chemotherapy.recode..yes..no.unk.)
unique(mydata8$RX.Summ..Systemic.Sur.Seq)
mydata8$chem <- ifelse(mydata8$Chemotherapy.recode..yes..no.unk. == "No/Unknown" &
                         mydata8$RX.Summ..Systemic.Sur.Seq == 
                         "No systemic therapy and/or surgical procedures", 0,
                       ifelse(mydata8$Chemotherapy.recode..yes..no.unk. == "Yes", 1, NA))
unique(mydata8[is.na(mydata8$chem), ]$Chemotherapy.recode..yes..no.unk.)
table(mydata8[is.na(mydata8$chem), ]$RX.Summ..Systemic.Sur.Seq)
table(mydata8[mydata8$chem == 0, ]$RX.Summ..Systemic.Sur.Seq)
mydata8$chem <- factor(mydata8$chem, levels = c(0, 1), labels = c("No", "Yes"))
summary(mydata8$chem)

#🔴转换生存状态
unique(mydata8$COD.to.site.recode)
unique(mydata8$SEER.cause.specific.death.classification)
mydata8$dead <- ifelse(mydata8$COD.to.site.recode == "Alive", 0, 1)
mydata8$dead <- factor(mydata8$dead)
summary(mydata8$dead)

mydata8$status <- ifelse(mydata8$COD.to.site.recode == "Alive", 0,
                         ifelse(mydata8$SEER.cause.specific.death.classification == 
                                  "Dead (attributable to this cancer dx)", 1, 
                                ifelse(mydata8$SEER.cause.specific.death.classification ==
                                         "Alive or dead of other cause" &
                                         mydata8$COD.to.site.recode != "Alive", 2, NA)))
mydata8$status <- factor(mydata8$status)
summary(mydata8$status)

#🔴转换生存时间
range(mydata8$Survival.months)
mydata8$time <- as.numeric(mydata8$Survival.months)
mydata9 <- subset(mydata8, Survival.months > 1)
chr9 <- paste0("限定生存时间后：排除", nrow(mydata8)-nrow(mydata9),
               "行，剩余", nrow(mydata9), "行")

#🔴转换家庭收入中位数
table(mydata9$Median.household.income.inflation.adj.to.2019)
mydata9$income <- ifelse(mydata9$Median.household.income.inflation.adj.to.2019 == "< $35,000" |
                           mydata9$Median.household.income.inflation.adj.to.2019 == "$35,000 - $39,999" |
                           mydata9$Median.household.income.inflation.adj.to.2019 == "$40,000 - $44,999" |
                           mydata9$Median.household.income.inflation.adj.to.2019 == "$45,000 - $49,999" |
                           mydata9$Median.household.income.inflation.adj.to.2019 == "$50,000 - $54,999" |
                           mydata9$Median.household.income.inflation.adj.to.2019 == "$55,000 - $59,999", 
                         "<$60,000",
                         ifelse(mydata9$Median.household.income.inflation.adj.to.2019 == "$60,000 - $64,999" |
                                mydata9$Median.household.income.inflation.adj.to.2019 == "$65,000 - $69,999" |
                                mydata9$Median.household.income.inflation.adj.to.2019 == "$70,000 - $74,999",
                                "$60,000-$74,999", 
                                ifelse(mydata9$Median.household.income.inflation.adj.to.2019 
                                       == "$75,000+", ">$75,000", NA)))
mydata9$income <- factor(mydata9$income, levels = c("<$60,000", "$60,000-$74,999", ">$75,000"))
summary(mydata9$income)

summary(mydata9$peri)

#🔴提取需要的变量
mydata <- mydata9[, c(1, which(colnames(mydata9) == "year") : ncol(mydata9))]
rownames(mydata) <- mydata[,1]
mydata <- mydata[,-1]
#检查数据（检查变量类型、分布）
summary(mydata)
mydata$T_stage <- factor(mydata$T_stage, levels = c("T1", "T2", "T3", "T4"))

#导出数据
save(mydata, file = "sarcoma_peri/data/lms_ess_2010_2016.Rdata")



