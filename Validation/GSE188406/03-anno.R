
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
load(file = "./2-degs-pca/1-parent.alc.lorla.Rdata")

# ID转换
AleDEG$SYMBOL = rownames(AleDEG)
s2e = bitr(AleDEG$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Aledeg = inner_join(AleDEG, s2e, by = "SYMBOL")

LorDEG$SYMBOL = rownames(LorDEG)
s2e = bitr(LorDEG$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
Lordeg = inner_join(LorDEG, s2e, by = "SYMBOL")

# GO富集
Alecg = Aledeg[Aledeg$change != "NOT", "ENTREZID"]
Lorcg = Lordeg[Lordeg$change != "NOT", "ENTREZID"]
cg = intersect(Alecg, Lorcg)

ego = enrichGO(gene = cg, OrgDb = org.Hs.eg.db, ont = "ALL",pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
ego.kk = ego@result
rownames(ego.kk) = 1:2069

GO = ego[c(47,90,114,157,206,222,228,355),c(2,3,1,6,9)]
GO$geneID = str_replace_all(GO$geneID,"/",",") ### 修改geneID这一列的分隔符号
names(GO) = c("ID","Term","category","adj_pval","genes")

genedata = data.frame(ID = Aledeg$SYMBOL, logFC = Aledeg$logFC)
circ = circle_dat(GO, genedata)

# 弦图
gene = circ %>% group_by(term) %>% top_n(5, logFC)
chord = chord_dat(circ, genes = gene$genes)
GOChord(chord, gene.order = "logFC", ribbon.col = brewer.pal(8, "Set2"))

ggsave(filename = "./3-anno/1-DEGs-go-circ2.pdf", width = 3.7, height = 4.7)

# 准备差异基因列表1
data(geneList)
geneList = Aledeg$logFC
geneList = as.numeric(geneList)
names(geneList) = Aledeg$ENTREZID
geneList = sort(geneList, decreasing = T)

# gsea
egmt = gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, ont = "ALL", minGSSize = 100, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
egmtd = egmt@result
rownames(egmtd) = 1:545
write.csv(egmtd, file = "./3-anno/1-Aledeg.GSEA.csv")
setID = egmtd$ID[c(176, 199, 378, 412, 89)]
gseaNb(object = egmt, geneSetID = setID, curveCol = brewer.pal(5, "Set2"), subPlot = 2, addPval = T)

# 准备差异基因列表2
data(geneList)
geneList = Lordeg$logFC
geneList = as.numeric(geneList)
names(geneList) = Lordeg$ENTREZID
geneList = sort(geneList, decreasing = T)

# gsea
egmt = gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, ont = "ALL", minGSSize = 100, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
egmtd = egmt@result
rownames(egmtd) = 1:448
setID = egmtd$ID[c(171, 155, 397, 271, 37)]
gseaNb(object = egmt, geneSetID = setID, curveCol = brewer.pal(5, "Set2"), subPlot = 2, addPval = T)
write.csv(egmtd, file = "./3-anno/1-Lordeg.GSEA.csv")
