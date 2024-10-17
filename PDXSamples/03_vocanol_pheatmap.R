
# 导入数据
rm(list = ls())
load("./Rdata/Brigatinib_DEG.Rdata")

table(Group)
library(pheatmap)
library(tinyarray)
dat = log2(exp + 1)

cg1 = rownames(B15vsNC)[B15vsNC$change !="NOT"]
cg2 = rownames(B30vsNC)[B30vsNC$change !="NOT"]
cg1 = cg1[cg1 %in% cg2]
cg2 = cg2[cg2 %in% cg1]
cg = cg1

# 热图
h1 = draw_heatmap(dat[cg1,], Group, n_cutoff = 2,
                  legend = T, annotation_legend= T,
                  cluster_cols = F)
ggsave(h1, filename = "./Rplot/Brigatinib_NC_all_heatmap.pdf",
       width = 5.9)

# 火山图
library(EnhancedVolcano)
colnames(B15vsNC)

vo1 = EnhancedVolcano(B15vsNC,
                lab = rownames(B15vsNC),
                selectLab = c("CXCL1", "CXCL2", "CXCL8", 
                              "IL6", "IL1B", "IL1A",
                              "CSF1", "CSF3", "TNF"),
                x = 'logFC',
                y = 'P.Value',
                ylim = c(0, 6),
                title = 'Brigatinib_15 vs. NC',
                pCutoff = 0.05,
                FCcutoff = 2,
                legendLabels = NULL,
                legendIconSize = 0,
                subtitle = NULL,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'grey30')

vo2 = EnhancedVolcano(B30vsNC,
                      lab = rownames(B30vsNC),
                      selectLab = c("CXCL1", "CXCL2", "CXCL8", 
                                    "IL6", "IL1B", "IL1A",
                                    "CSF1", "CSF3", "TNF"),
                      x = 'logFC',
                      y = 'P.Value',
                      ylim = c(0, 6),
                      title = 'Brigatinib_30 vs. NC',
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      legendLabels = NULL,
                      legendIconSize = 0,
                      subtitle = NULL,
                      drawConnectors = TRUE,
                      widthConnectors = 0.2,
                      colConnectors = 'grey30')

library(patchwork)
Brigatinib_vo = vo1 + vo2
ggsave(Brigatinib_vo, filename = "./Rplot/Brigatinib_NC_all_vocanol.pdf",
       width = 11)

# 两组差异基因比较
library(dplyr)
UP=function(df){
  rownames(df)[df$change=="UP"]
}
DOWN=function(df){
  rownames(df)[df$change=="DOWN"]
}

up = intersect(UP(B15vsNC),UP(B30vsNC))
down = intersect(DOWN(B15vsNC),DOWN(B30vsNC))

down_re = c("CXCL1", "CXCL2", "CXCL8", 
            "IL6", "IL1B", "IL1A",
            "CSF1", "CSF3", "TNF")

hp = draw_heatmap(dat[c(down[1:50]),],
                  Group,
                  n_cutoff = 2, 
                  legend = T, 
                  annotation_legend = T,
                  show_rownames = T)

library(patchwork)
Brigatinib_50down = hp
ggsave(Brigatinib_50down, filename = "./Rplot/Brigatinib_50down_heatmap.pdf",
       width = 8)

hp = draw_heatmap(dat[down_re,],
                  Group,
                  n_cutoff = 2, 
                  legend = T, 
                  annotation_legend = T,
                  show_rownames = T)
library(patchwork)
Brigatinib_9down = hp
ggsave(Brigatinib_9down, filename = "./Rplot/Brigatinib_9down_heatmap.pdf",
      width = 5.1,
      height = 3.5)



#上调、下调基因分别画维恩图
up_genes = list(Brigatinib_15 = UP(B15vsNC),
                Brigatinib_30 = UP(B30vsNC))

down_genes = list(Brigatinib_15 = DOWN(B15vsNC),
                  Brigatinib_30 = DOWN(B30vsNC))

up.plot <- draw_venn(up_genes,"UPgene")
down.plot <- draw_venn(down_genes,"DOWNgene")


# 韦恩图



# 韦恩图
library(UpSetR)

## 上调基因
genevene_up = list("Brigatinib_15" = rownames(B15vsNC)[B15vsNC$change == "UP"],
                   "Brigatinib_30" = rownames(B30vsNC)[B30vsNC$change == "UP"])
lapply(genevene_up,head,3)
genevene_up_p = upset(fromList(genevene_up), order.by = "freq",
                      line.size = 1,
                      point.size = 2,
                      text.scale = 2,
                      matrix.dot.alpha = 1,
                      queries = list(list(query = intersects,
                                          params = list("Brigatinib_15","Brigatinib_30"),
                                          color = "darkred",
                                          active = T)))
genevene_up_p

## 下调基因
genevene_down = list("Brigatinib_15" = rownames(B15vsNC)[B15vsNC$change == "DOWN"],
                     "Brigatinib_30" = rownames(B30vsNC)[B30vsNC$change == "DOWN"])
lapply(genevene_down,head,3)
genevene_down_p = upset(fromList(genevene_down), order.by = "freq",
                        line.size = 1,
                        point.size = 2,
                        text.scale = 2,
                        matrix.dot.alpha = 1,
                        queries = list(list(query = intersects,
                                            params = list("Brigatinib_15","Brigatinib_30"),
                                            color = "darkblue",
                                            active = T)))
genevene_down_p 

# 拼图
library(ggplotify)
genevene_up_p <- as.ggplot(genevene_up_p)
genevene_down_p <- as.ggplot(genevene_down_p)
genevene = genevene_up_p + genevene_down_p
ggsave(genevene, filename = "./Rplot/Brigatinib_genevenep.pdf",
       width = 10,
       height = 6)
