
# 导入数据
rm(list = ls())  
load(file = "./Rdata/Brigatinib_DEG.Rdata")


# 两组取交集
UP=function(df){
  rownames(df)[df$change=="UP"]
}
DOWN=function(df){
  rownames(df)[df$change=="DOWN"]
}

gene_up = intersect(UP(B15vsNC),UP(B30vsNC))
gene_down = intersect(DOWN(B15vsNC),DOWN(B30vsNC))


# 1.制作string的输入数据
write.table(gene_down,
            file="./String/Brigatinib_gene_down.txt",
            row.names = F,
            col.names = F,
            quote = F)
# 从string网页获得string_interactions.tsv

# 2.准备cytoscape的输入文件
tsv = read.table("./String/string_interactions.tsv",
                 comment.char = "!",
                 header = T)
tsv2 = tsv[,c(1,2)]
head(tsv2)

p = B30vsNC[B30vsNC$change == "DOWN",
            c("logFC","P.Value")]
p = p[rownames(p) %in% gene_down, ]
p$symbol = rownames(p)
p = p[,c("symbol", "logFC","P.Value")]
head(p)
write.table(p,
            file = "./String/deg.txt",
            sep = "\t",
            quote = F,
            row.names = F)
