
# 导入数据
rm(list = ls())
library(DESeq2)
library(limma)
library(edgeR)
library(FactoMineR)
library(factoextra)
library(ggplot2)

load(file = "./1-data-prepare/GSE188406.matrix.ph.Rdata")
grouplist = ifelse(str_detect(string = ph$Sample_characteristics_ch1.1, pattern = "DMSO"), "Parent",
                   ifelse(str_detect(ph$Sample_characteristics_ch1.1, pattern = "Alectinib"), "AlectinibRes", "LorlatinibRes"))
table(grouplist)
grouplist = factor(grouplist, levels = c("Parent", "AlectinibRes", "LorlatinibRes"))

# PCA
dat = as.data.frame(t(exp))
dat.pca = PCA(dat, graph = FALSE)
pca_plot = fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = grouplist, # color by groups
                         palette = c("#293462", "#a64942", "#009900"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)

ggsave(pca_plot, filename = "./2-degs-pca/1-PARENT_RESIS_PCA.pdf", width = 4, height = 2.5)


# 差异基因AleDEG
design = model.matrix(~0+grouplist)
colnames(design) = levels(grouplist)
rownames(design) = colnames(exp)

dge = DGEList(counts = exp)
dge = calcNormFactors(dge)

v = voom(dge,design, normalize = "quantile")
fit = lmFit(v, design)

constrasts = "AlectinibRes-Parent"
cont.matrix = makeContrasts(contrasts = constrasts, levels = design) 
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

AleDEG = topTable(fit2, coef = constrasts, n = Inf)
AleDEG = na.omit(AleDEG)

# 差异基因统计分析
pvalue_t = 0.05
logFC_t = log2(1.5)

k1 = (AleDEG$P.Value < pvalue_t)&(AleDEG$logFC < -logFC_t);table(k1)
k2 = (AleDEG$P.Value < pvalue_t)&(AleDEG$logFC > logFC_t);table(k2)
AleDEG$change = ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))

table(AleDEG$change)
write.csv(x = AleDEG, file = "./2-degs-pca/1-alectinib.degs.csv")

# 差异基因LorDEG
design = model.matrix(~0+grouplist)
colnames(design) = levels(grouplist)
rownames(design) = colnames(exp)

dge = DGEList(counts = exp)
dge = calcNormFactors(dge)

v = voom(dge,design, normalize = "quantile")
fit = lmFit(v, design)

constrasts = "LorlatinibRes-Parent"
cont.matrix = makeContrasts(contrasts = constrasts, levels = design) 
fit2 = contrasts.fit(fit, cont.matrix)
fit2 = eBayes(fit2)

LorDEG = topTable(fit2, coef = constrasts, n = Inf)
LorDEG = na.omit(LorDEG)

# 差异基因统计分析
pvalue_t = 0.05
logFC_t = log2(1.5)

k1 = (LorDEG$P.Value < pvalue_t)&(LorDEG$logFC < -logFC_t);table(k1)
k2 = (LorDEG$P.Value < pvalue_t)&(LorDEG$logFC > logFC_t);table(k2)
LorDEG$change = ifelse(k1, "DOWN", ifelse(k2, "UP", "NOT"))

table(LorDEG$change)
write.csv(x = LorlaDEG, file = "./2-degs-pca/2-lorlatinib.degs.csv")

save(AleDEG, LorDEG, file = "./2-degs-pca/1-parent.alc.lorla.Rdata")


