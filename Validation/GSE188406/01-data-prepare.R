
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(tinyarray)
library(AnnoProbe)
library(GEOquery)

# 读取数据
exp = read.csv(file = "./0-RawData/GSE188406.matrix.csv")
rownames(exp) = exp$ID_REF
exp = exp[,-1]

ph = read.csv(file = "./0-RawData/GSE188406.ph.csv")
rownames(ph) = ph$Sample_geo_accession

# ID转换
GPL = getGEO("GPL23159")
anno = GPL@dataTable@table
ids = anno[,c(1,10)]

genename = str_split(string = ids$SPOT_ID, pattern = "//", simplify = T)[,c(1,3)]
genename = as.data.frame(genename)
genename$ID = anno$ID

table(str_detect(string = genename$V2, pattern = "[(]"))
genename = genename[str_detect(string = genename$V2, pattern = "[(]"),]

genename2 = str_split(string = genename$V2, pattern = "[(]", simplify = T)
genename2 = as.data.frame(genename2)

genename2$V5 = ifelse(genename2$V6 == "", genename2$V5, genename2$V6)
genename2$V4 = ifelse(genename2$V5 == "", genename2$V4, genename2$V5)
genename2$V3 = ifelse(genename2$V4 == "", genename2$V3, genename2$V4)
genename2$V2 = ifelse(genename2$V3 == "", genename2$V2, genename2$V3)

genename2 = str_split(string = genename2$V2, pattern = "[)]", simplify = T)

ids = data.frame(ID = genename$ID, GeneSymbol = genename2[,1])
ids = ids[!str_detect(string = ids$GeneSymbol, pattern = " "),]

exp = exp[rownames(exp) %in% ids$ID, ]
exp = trans_array(exp = exp, ids = ids, from = "ID", to = "GeneSymbol")

# 取A925L
ph = ph[ph$Sample_characteristics_ch1 == "cell line: A925L", ]
exp = exp[,rownames(ph)]

# 保存数据
save(exp, ph, file = "./1-data-prepare/GSE188406.matrix.ph.Rdata")
