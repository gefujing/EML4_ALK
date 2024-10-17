# 准备环境
rm(list = ls())
options(stringsAsFactors = F)
library(timeROC)
library(survival)
library(survminer)
library(ROCR)
library(hash)
library(qvalue)
library(VennDiagram)
library(survMisc)
library(pROC)
library(ggplot2)
library(plyr)
library(rms)
library(rmda)
library(ggDCA)

load(file = "./Rdata/tcga.alk.signature.Rdata")

# 数据为分类变量
df$age = ifelse(df$age >= median(df$age), "Older", "Young")
df$age = factor(df$age, levels = c("Young", "Older"))

df$race = ifelse(df$race == "white", "White", "Others")
df$race = factor(df$race, levels = c("White", "Others"))

df$N = ifelse(df$N == "N0", "N0", "N1/2")
df$N = factor(df$N, levels = c("N0", "N1/2"))

df$T = ifelse(df$T %in% c("T1","T2"), "T1/2", "T3/4")
df$T = factor(df$T, levels = c("T1/2", "T3/4"))

df$stage = ifelse(df$stage == "stage I", "stage I", "stage II/III/IV")
df$stage = factor(df$stage, levels = c("stage I", "stage II/III/IV"))

df$smoking_history = ifelse(df$smoking_history == 1, 0, 1)

df$diagnosis = factor(df$diagnosis, levels = unique(df$diagnosis))
df$site = factor(df$site, levels = unique(df$site))

# 数据预处理
rt = df[,c(1:14, ncol(df))]
rt = rt[,c(1:6, 8:12, 15)]
rt[1:4, 1:4]
rt = rt %>% select(OS.time, OS, everything())

bc = rt

bc$signature_by2 = ifelse(bc$signiture < median(rt$signiture), "Low-Score", "High-Score")

rm(list = c("df", "rt"))

# 封装数据
dd = datadist(bc)

GetFactors1 = c("age", "gender", "smoking_history", "stage")
GetFactors2 = c("age", "gender", "smoking_history", "stage", "signiture")

fml1 = as.formula(paste0("Surv(OS.time, OS)~", paste0(GetFactors1, collapse = "+")))
fml2 = as.formula(paste0("Surv(OS.time, OS)~", paste0(GetFactors2, collapse = "+")))

CP = cph(fml1, bc)
CP_singniture = cph(fml2, bc)

data  <- dca(CP, CP_singniture)
ggplot(data = data,
       linetype = F)
ggsave(filename = "./Rplot/alk.DCA_curve.pdf", width = 6, height = 3.5)
















