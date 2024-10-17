

# 导入数据
rm(list = ls())
library(tidyverse)
library(pheatmap)
library(ggsci)
library(ggsignif)
library(ggrepel)
load(file = "./2-degs-pca/1-parent.alc.lorla.Rdata")
load(file = "./1-data-prepare/GSE188406.matrix.ph.Rdata")

# 散点图
kk = intersect(rownames(AleDEG), rownames(LorDEG))

AleDEG = AleDEG[kk,]
LorDEG = LorDEG[kk,]

cg = ifelse(AleDEG$change == "UP" & LorDEG$change == "UP", "UP",
            ifelse(AleDEG$change == "DOWN" & LorDEG$change == "DOWN", "DOWN", "NOT"))

pdata = data.frame(Alec = AleDEG$logFC,
                   LorLa = LorDEG$logFC,
                   group = cg)
pdata$group = factor(pdata$group, levels = c("UP", "NOT", "DOWN"))
rownames(pdata) = rownames(AleDEG)

lable = pdata[c("TAGLN", "FILIP1L", "TSPAN2", "SFTPB", "NOX5", "ACPP"),]
lable$ID = rownames(lable)

ggplot(data = pdata, mapping = aes(x = Alec, y = LorLa, color = group))+
  geom_point()+
  geom_vline(xintercept = 0, lty = 4)+
  geom_hline(yintercept = 0, lty = 4)+
  geom_label_repel(aes(label = ID),data = lable, color="black")+
  theme_bw()+
  scale_color_manual(values = c("#F15A24", "#C3CFE2", "#20B2AA"))
  
ggsave(filename = "./2-degs-pca/2-degs.pdf", width = 3.9, height = 2.6)

# 热图
pdata = pdata[order(pdata$Alec),]
cg = rownames(pdata)[pdata$group != "NOT"]
cg = cg[c(1:10, 80:89)]

exp = exp[cg,]

colanno = data.frame(group = c(rep("Parent", 3), rep("AlectinibRes",3), rep("LorlatinibRes", 3)))
rownames(colanno) = colnames(exp)

pheatmap(exp,
         scale = "column",
         annotation_col = colanno,
         show_colnames = F,
         breaks = seq(-1.5,1.5,length.out = 100),
         cellwidth = 15, cellheight = 15)

