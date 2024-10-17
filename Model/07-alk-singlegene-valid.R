
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggsignif)
library(survival)
library(survminer)
library(pheatmap)
library(tinyarray)
load(file = "./Rdata/alk.rawdata.Rdata")

# IL6单基因分析
fp = df$IL6
ri = ifelse(fp<median(fp),"IL6-Low","IL6-High")
ri = factor(ri,levels = c("IL6-Low","IL6-High"))

sfit <- survfit(Surv(OS.time, OS)~ri, data = df)

ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =F,xlab ="Time/days", 
           ggtheme =theme_light(), 
           ncensor.plot = F)

# CXCL1单基因分析
fp = df$CXCL1
ri = ifelse(fp<median(fp),"CXCL1-Low","CXCL1-High")
ri = factor(ri,levels = c("CXCL1-Low","CXCL1-High"))

sfit <- survfit(Surv(OS.time, OS)~ri, data = df)

ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =F,xlab ="Time/days", 
           ggtheme =theme_light(), 
           ncensor.plot = F)

# CXCL5单基因分析
fp = df$CXCL5
ri = ifelse(fp<median(fp),"CXCL5-Low","CXCL5-High")
ri = factor(ri,levels = c("CXCL5-Low","CXCL5-High"))

sfit <- survfit(Surv(OS.time, OS)~ri, data = df)

ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =F,xlab ="Time/days", 
           ggtheme =theme_light(), 
           ncensor.plot = F)
