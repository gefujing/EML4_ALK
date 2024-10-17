
# 导入数据
rm(list = ls())
library(DESeq2)
library(limma)
library(edgeR)
library(FactoMineR)
library(factoextra)
library(ggplot2)
load(file = "./1-data-prepare/GSE186767.matrix.ph.Rdata")
table(grouplist)

# 差异基因CrizolDEG
design = model.matrix(~0+grouplist)
colnames(design) = levels(grouplist)
rownames(design) = colnames(exp)

dge = DGEList(counts = exp)
dge = calcNormFactors(dge)

v = voom(dge,design, normalize = "quantile")
fit = lmFit(v, design)

constrasts = "Crizoleres-Parent"
cont.matrix = makeContrasts(contrasts = constrasts, levels = design) 
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

CrizolDEG = topTable(fit2, coef = constrasts, n = Inf)
CrizolDEG = na.omit(CrizolDEG)

# 差异基因统计分析
pvalue_t = 0.05
logFC_t = log2(2)

k1 = (CrizolDEG$P.Value < pvalue_t)&(CrizolDEG$logFC < -logFC_t);table(k1)
k2 = (CrizolDEG$P.Value < pvalue_t)&(CrizolDEG$logFC > logFC_t);table(k2)
CrizolDEG$change = ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))

table(CrizolDEG$change)
write.csv(x = CrizolDEG, file = "./2-degs/1-crizoltinib.degs.csv")
save(CrizolDEG, file = "./2-degs/1-parent.crizol.lorla.Rdata")

CrizolDEG$sym = rownames(CrizolDEG)
