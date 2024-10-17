
# 导入数据
rm(list = ls())
library(tidyverse)
library(pheatmap)
library(ggsci)
library(ggsignif)
library(ggrepel)
load(file = "./2-degs/1-parent.crizol.lorla.Rdata")
load(file = "./1-data-prepare/GSE186767.matrix.ph.Rdata")

# 火山图
pdata = CrizolDEG
pdata$symbol = rownames(pdata)
for_label = pdata[c("IL6", "CXCL1", "CXCL5", "CSF2", "IL8", "IL1B", "PTGS2", "CCL5", "LIF"),]

ggplot(data = pdata) + 
  geom_point(aes(x = logFC, y = -log10(P.Value), color = logFC, size = -log10(P.Value))) + 
  geom_text_repel(data =  for_label, aes(x = logFC, y = -log10(P.Value), label = symbol),
                  nudge_x = 0.5, nudge_y = 0.2, segment.curvature = -0.1,
                  segment.ncp = 3, direction = "y", hjust = "left", max.overlaps = 200 ) + 
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  scale_fill_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                       values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(log2(1)), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) + 
  theme_bw() + 
  ggtitle(label = "Volcano Plot")
ggsave(filename = "./2-degs/1-vocanal.pdf", width = 4.6, height = 3.5)

# 热图
pdata = pdata[order(pdata$logFC),]
cg = pdata$symbol[pdata$change != "NOT"]
cg = cg[!str_detect(string = cg, pattern = "-")]
cg = cg[!str_detect(string = cg, pattern = "[.]")]
cg = cg[c(1:10, 987:996)]
cg = c(cg, "IL6", "CXCL1", "CXCL5")
cg = unique(cg)

exp = exp[cg,]

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
         breaks = seq(-1.5,1.5,length.out = 100),
         cellwidth = 20, cellheight = 20)

