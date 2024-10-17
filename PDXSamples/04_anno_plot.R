
# 导入数据
rm(list = ls())
load("./Rdata/Brigatinib_DEG.Rdata")
table(Group)

library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(stringr)
library(enrichplot)
library(patchwork)
library(org.Hs.eg.db)

# ID转换
## 15
B15vsNC$SYMBOL = rownames(B15vsNC)
s2e <- bitr(B15vsNC$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
dim(B15vsNC)

deg15 = inner_join(B15vsNC, s2e, by="SYMBOL")
dim(deg15)
length(unique(deg15$SYMBOL))

## 30
B30vsNC$SYMBOL = rownames(B30vsNC)
s2e <- bitr(B30vsNC$SYMBOL, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
dim(B30vsNC)

deg30 = inner_join(B30vsNC, s2e, by="SYMBOL")
dim(deg30)
length(unique(deg30$SYMBOL))

# 两组取交集
UP=function(df){
  df$ENTREZID[df$change=="UP"]
}
DOWN=function(df){
  df$ENTREZID[df$change=="DOWN"]
}

gene_up = intersect(UP(deg15),UP(deg30))
gene_down = intersect(DOWN(deg15),DOWN(deg30))

save(Group, deg15, deg30, gene_up, gene_down,
     file = "./Rdata/Brigatinib_anno.Rdata")


# GO富集分析
rm(list = ls())  
load(file = "./Rdata/Brigatinib_anno.Rdata")

## 上调基因
ego_up <- enrichGO(gene = gene_up,
                   OrgDb= org.Hs.eg.db,
                   ont = "ALL", #ont参数：One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                   readable = TRUE)
## 下调基因
ego_down <- enrichGO(gene = gene_down,
                OrgDb= org.Hs.eg.db,
                ont = "ALL", #ont参数：One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
                readable = TRUE)

save(ego_up, ego_down, file = "./Rdata/Brigatinib_GO.Rdata")

# 可视化
rm(list = ls())  
load("./Rdata/Brigatinib_GO.Rdata")
go_up = dotplot(ego_up)
go_down = dotplot(ego_down)
go_up + go_down

go_up2 = dotplot(ego_up, split = "ONTOLOGY", font.size = 10, 
         showCategory = 5) + 
         facet_grid(ONTOLOGY ~ ., scale = "free") + 
         scale_y_discrete(labels = function(x) str_wrap(x, width = 45))
go_down2 = dotplot(ego_down, split = "ONTOLOGY", font.size = 10, 
                 showCategory = 5) +
           facet_grid(ONTOLOGY ~ ., scale = "free") + 
           scale_y_discrete(labels = function(x) str_wrap(x, width = 45))

ggsave(go_down, filename = "./Rplot/Brigatinib_GO.pdf",
       width = 9.5,
       height = 5)

egodown = ego_down@result
write.csv(egodown, file = "./RawData/ALK_Brigtinib_GO_down.csv")

# KEGG pathway analysis
# 输入数据
rm(list = ls())  
load(file = "./Rdata/Brigatinib_anno.Rdata")

# 对上调/下调/所有差异基因进行富集分析
kk.up <- enrichKEGG(gene = gene_up,
                    organism = 'hsa')
kk.down <- enrichKEGG(gene = gene_down,
                      organism = 'hsa')
save(kk.down,kk.up,file = "./Rdata/Brigatinib_KEGG.Rdata")

# 可视化
rm(list = ls()) 
load("./Rdata/Brigatinib_KEGG.Rdata")

table(kk.up@result$p.adjust<0.05)
table(kk.down@result$p.adjust<0.05)

head(kk.up)[,1:6]
head(kk.down)[,1:6]
browseKEGG(kk.down,"hsa04657")

kegg_down = barplot(kk.down)
kegg_down

ggsave(kegg_down, filename = "./Rplot/Brigatinib_KEGG.pdf",
       width = 9.5,
       height = 5)

write.csv(deg15, "./RawData/Brigatinib_15_deg.csv")
write.csv(deg30, "./RawData/Brigatinib_30_deg.csv")