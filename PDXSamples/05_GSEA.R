
# 清空环境变量
rm(list = ls()) 
options(stringsAsFactors = F)

# 加载包
library(GSEABase)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(ggplot2)
library(stringr)
library(msigdbr)
library(enrichplot)
library(ggsci)

#准备差异基因列表
DEG = read.csv("./RawData/Brigatinib_15_deg.csv")
data(geneList)
geneList = DEG$logFC
geneList = as.numeric(geneList)
names(geneList) = DEG$ENTREZID
geneList = sort(geneList, decreasing = T)

#准备参考基因集
geneset = msigdbr(species = "Homo sapiens")
geneset = geneset[ , c(3,5)]

geneset$gs_name = str_remove(geneset$gs_name, "HALLMARK_")
geneset$gs_name = str_replace_all(geneset$gs_name, "[_]", " ")
geneset$gs_name = str_to_title(geneset$gs_name)

# gsea
egmt <- GSEA(geneList, TERM2GENE = geneset, verbose = F)
#> Warning in fgsea(pathways = geneSets, stats = geneList, nperm = nPerm, minSize = minGSSize, : There are ties in the preranked stats (0.05% of the list).
#> The order of those tied genes will be arbitrary, which may produce unexpected results.
egmt2 <- setReadable(egmt, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
egmtd <- data.frame(egmt)
save(egmt, egmt2, egmtd, file = "./Rdata/Brgatinib_15_GSEA.Rdata")


# 可视化
rm(list = ls()) 
options(stringsAsFactors = F)
load(file = "./Rdata/Brgatinib_15_GSEA.Rdata")

#气泡图，展示geneset被激活还是抑制
dotplot(egmt2, split=".sign") + facet_grid(~.sign)
#山峦图，展示每个geneset的基因logFC分布
ridgeplot(egmt, showCategory = 30)


#选择单个gene set作图
table(str_detect(string = egmtd$Description, pattern = "ANTIMICROBIAL"))
k = str_detect(string = egmtd$Description, pattern = "ANTIMICROBIAL")
egmtd$ID[k]

p = gseaplot2(egmt, geneSetID = "GOBP_AXONEME_ASSEMBLY", 
              title = "GOBP_AXONEME_ASSEMBLY",
              color = "#121b74",
              pvalue_table = F,
              ES_geom = "line")
ggsave(p, filename = "./Rplot/GOBP_ACUTE_INFLAMMATORY_RESPONSE.pdf", width = 5.7, height = 4.4)


#多个gnenset合并展示
gseaplot2(egmt, geneSetID = 1:3, pvalue_table = T)


# 气泡图和山峦图
#山峦图，展示每个geneset的基因logFC分布
k = c("GOBP_LIPOPOLYSACCHARIDE_MEDIATED_SIGNALING_PATHWAY",
      "GOBP_RESPONSE_TO_MOLECULE_OF_BACTERIAL_ORIGIN",
      "GOBP_CELL_CHEMOTAXIS",
      "GOBP_LEUKOCYTE_CHEMOTAXIS",
      "GOBP_GRANULOCYTE_CHEMOTAXIS",
      "GOBP_ANTIMICROBIAL_HUMORAL_RESPONSE",
      "GOBP_GRANULOCYTE_MIGRATION",
      "GOBP_ACUTE_INFLAMMATORY_RESPONSE",
      "GOBP_ANTIMICROBIAL_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_ANTIMICROBIAL_PEPTIDE")
egmtd2 = egmtd[egmtd$ID %in% k, ]
colnames(egmtd2)
egmtd2$LP = -log10(egmtd2$pvalue)
egmtd2$LQ = -log10(egmtd2$qvalues)

p = ggplot(data = egmtd2)+
    geom_point(mapping = aes(x = NES, 
                             y = ID, 
                             size = LP,
                             color = LQ))+
    scale_color_gradient(low = "#6b48ff",high = "#ff502f")+
    theme_bw()
ggsave(p, filename = "./Rplot/Brigatinib_15_GESA_NES.pdf", width = 9.34, height = 4.44)

# NES和FDR提取
write.csv(egmtd2, file = "./RawData/Brigatinib_15_GSEA_NES.csv")
