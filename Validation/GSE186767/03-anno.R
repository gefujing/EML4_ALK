
# 导入数据
rm(list = ls())
library(DESeq2)
library(limma)
library(edgeR)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(enrichplot)
library(msigdbr)
library(GSEABase)
library(GseaVis)
load(file = "./2-degs/1-parent.crizol.lorla.Rdata")

# ID转换
CrizolDEG$SYMBOL = rownames(CrizolDEG)
s2e = bitr(CrizolDEG$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
CrizolDEG = inner_join(CrizolDEG, s2e, by = "SYMBOL")

# 准备差异基因列表1
data(geneList)
geneList = CrizolDEG$logFC
geneList = as.numeric(geneList)
names(geneList) = CrizolDEG$ENTREZID
geneList = sort(geneList, decreasing = T)

# gsea
egmt = gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, ont = "ALL", minGSSize = 100, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
egmtd = egmt@result
rownames(egmtd) = 1:109
write.csv(egmtd, file = "./3-anno/1-CrizolDEG.GSEA.csv")
setID = egmtd$ID[c(15,25,62,105,91)]
gseaNb(object = egmt, geneSetID = setID, curveCol = brewer.pal(5, "Set2"), subPlot = 2, addPval = T)



