
# 准备环境
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(pheatmap)

# 读取数据
exp = read.csv(file = "./0-RawData/GSE186767.matrix.tpm.csv")
exp = exp[!duplicated(exp$Genename),]
rownames(exp) = exp$Genename
exp = exp[,-1]

grouplist = c(rep("Parent", 3), rep("Crizoleres", 3))
grouplist = factor(grouplist, levels = c("Parent", "Crizoleres"))

save(exp, grouplist, file = "./1-data-prepare/GSE186767.matrix.ph.Rdata")

# 热图
exp = exp[c("IL6", "CXCL1", "CXCL5"),]
colanno = data.frame(row.names = colnames(exp),
                     group = factor(grouplist, levels = c("Parent", "Crizoleres")))
anncolor = list(group = c(Parent = "#293462", Crizoleres = "#a64942")) 

pheatmap(mat = exp,
         scale = "row",
         show_colnames = F,
         cluster_cols = F,
         annotation_col = colanno,
         annotation_colors = anncolor,
         gaps_col = c(3),
         breaks = seq(-1,1,length.out = 100),
         cellwidth = 20, cellheight = 20)
