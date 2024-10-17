
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(WGCNA)
library(ggsci)

load(file = "./Rdata/tcga.alk.signature.Rdata")

# 数据为分类变量
df$age = ifelse(df$age >= median(df$age), "Older", "Young")
df$race = ifelse(df$race == "white", "White", "Others")
df$N = ifelse(df$N == "N0", "N0", "N1/2")
df$T = ifelse(df$T %in% c("T1","T2"), "T1/2", "T3/4")
df$stage = ifelse(df$stage == "stage I", "stage I", "stage II/III/IV")
df$smoking_history = ifelse(df$smoking_history == 1, 0, 1)