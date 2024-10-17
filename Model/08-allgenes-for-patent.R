
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(plyr)
library(survminer)
load(file = "./Rdata/alk.rawdata.Rdata")

## 写一个函数，方便后续批量计算
BaSurv = Surv(time = df$OS.time,
              event = df$OS)

UniCox = function(x){
  FML = as.formula(paste0("BaSurv~", x))
  GCox = coxph(FML, data = df)
  GSum = summary(GCox)
  HR = round(GSum$coefficients[,2],2)
  PValue = round(GSum$coefficients[,5],3)
  CI = paste0(round(GSum$conf.int[,3:4],2), collapse = "-")
  Unicox = data.frame("characteristics" = x,
                      "Hazard Ratio" = HR,
                      "CI95" = CI,
                      "P Value" = PValue)
  return(Unicox)
}

## 挑一个基因试一下
UniCox(colnames(df)[15])

# 设置基因
finalgenes = c("IL1A","FGF7","OSMR","VCAN","CD14",
               "IRAK3","TREM1",'CSF3','CXCL6',
               'IL6','CXCL2','TNF','CXCR4','APLN',
               'NAMPT','PTGS2','CXCL1','SAA1','IL1B',
               'SPP1','CSF1','CD274','SOCS3','CXCL3','CXCL5','AIM2')

# 单因素分析
VarNames = finalgenes
UniVar = lapply(VarNames, UniCox)
UniVar = ldply(UniVar, data.frame)
GetFactors = UniVar$characteristics[which(UniVar$P.Value < 0.3)] %>% as.character()
# GetFactors = VarNames













