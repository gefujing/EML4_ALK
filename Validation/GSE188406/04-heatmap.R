
# 导入数据
rm(list = ls())
library(tidyverse)
library(pheatmap)

load(file = "./1-data-prepare/GSE188406.matrix.ph.Rdata")

# 数据预处理
exp = exp[c("IL6", "CXCL1", "CXCL5"),]
ph$Sample_characteristics_ch1.1 = ifelse(ph$Sample_characteristics_ch1.1 == "treatment: DMSO", "Parent",
                                         ifelse(ph$Sample_characteristics_ch1.1 == "treatment: Alectinib 3 umol/L", "AlectinibRes", "LorlatinibRes"))

colanno = data.frame(row.names = ph$Sample_geo_accession,
                     group = factor(ph$Sample_characteristics_ch1.1, levels = unique(ph$Sample_characteristics_ch1.1)))
anncolor = list(group = c(Parent = "#293462", AlectinibRes = "#a64942", LorlatinibRes = "#009900")) 

pheatmap(mat = exp,
         scale = "row",
         show_colnames = F,
         cluster_cols = F,
         annotation_col = colanno,
         annotation_colors = anncolor,
         gaps_col = c(3, 6),
         breaks = seq(-2,2,length.out = 100),
         cellwidth = 20, cellheight = 20)
