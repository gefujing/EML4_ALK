
# 导入数据
rm(list = ls())
load("./Rdata/exp_Group.RData")
table(Group)

# PCA
dat=as.data.frame(t(exp))
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = Group, # color by groups
                         palette = c("#293462", "#a64942", "#fe5f55"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)

ggsave(pca_plot, filename = "./Rplot/Brigatinib_PCA.pdf",
       width = 4.5,
       height = 3.5)


exp = log2(exp+1)

# limma-voom
library(limma)
library(edgeR)

design <- model.matrix(~0+Group)
colnames(design)=levels(Group)
rownames(design)=colnames(exp)

dge <- DGEList(counts = exp)
# dge <- calcNormFactors(dge)
v <- voom(dge, design, plot=TRUE)

contr.matrix <- makeContrasts(
  Brigatinib_15vsNC = Brigatinib_15 - Control, 
  Brigatinib_30vsNC = Brigatinib_30 - Control,
  levels = colnames(design))
contr.matrix

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)
colnames(efit)
summary(decideTests(efit))

B15vsNC <- topTreat(efit, coef=1, n=Inf)
B30vsNC <- topTreat(efit, coef=2, n=Inf)

head(B15vsNC)
head(B15vsNC)  

# 标记上下调
pvalue_t = 0.05
logFC_t = 2

k1 = (B15vsNC$P.Value < pvalue_t)&(B15vsNC$logFC < -logFC_t);table(k1)
k2 = (B15vsNC$P.Value < pvalue_t)&(B15vsNC$logFC > logFC_t);table(k2)
B15vsNC$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(B15vsNC$change)
head(B15vsNC)

k1 = (B30vsNC$P.Value < pvalue_t)&(B30vsNC$logFC < -logFC_t);table(k1)
k2 = (B30vsNC$P.Value < pvalue_t)&(B30vsNC$logFC > logFC_t);table(k2)
B30vsNC$change = ifelse(k1,"DOWN",ifelse(k2,"UP","NOT"))
table(B30vsNC$change)
head(B30vsNC)

# merge
tj = data.frame(b15vsNC = as.integer(table(B15vsNC$change)),
                b30vsNC = as.integer(table(B30vsNC$change)),
                row.names = c("DOWN","STABLE","UP")
);tj


save(exp, Group, B15vsNC, B30vsNC, tj, file = "./Rdata/Brigatinib_DEG.Rdata")
