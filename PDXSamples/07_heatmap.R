
# 导入数据
rm(list = ls())
load("./Rdata/Brigatinib_DEG.Rdata")

table(Group)
library(pheatmap)
library(tinyarray)
dat = log2(exp + 1)

# 基因选择
cg = read.csv(file = "./String/mcode.csv")
cg = cg$name


# 热图
h1 = draw_heatmap(dat[cg,], Group, n_cutoff = 2,
                  legend = T, annotation_legend= T,
                  show_rownames = T,
                  cluster_cols = F)