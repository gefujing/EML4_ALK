

# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)

# 表型整理
## 读取表型信息
luad.phe = read.csv(file = "./Rawdata/TCGA-LUAD.phenotype.csv")
luad.sur = read.csv(file = "./Rawdata/TCGA-LUAD.survival.csv")
lusc.phe = read.csv(file = "./Rawdata/TCGA-LUSC.phenotype.csv")
lusc.sur = read.csv(file = "./Rawdata/TCGA-LUSC.survival.csv")

## 合并
phe = rbind(luad.phe, lusc.phe)
sur = rbind(luad.sur, lusc.sur)
rm(list = c("luad.phe", "luad.sur", "lusc.phe", "lusc.sur"))

## 储存数据
write.csv(phe, file = "./Table/TCGA-LUAD-LUSC-phe.csv")
write.csv(sur, file = "./Table/TCGA-LUAD-LUSC-sur.csv")

## 筛选过程在excel完成


# 数据读取
## 读取表达量
rm(list = ls())
load("./Rawdata/TCGA-LUAD-LUSC-exp.Rdata")
phe = read.csv(file = "./Rawdata/TCGA-LUAD-LUSC-phe-alk.csv")
sur = read.csv(file = "./Table/TCGA-LUAD-LUSC-sur.csv")

# 整理样本名

## exp样本名
colnames(exp) = str_replace_all(string = colnames(exp), pattern = "[.]", replacement = "-")

## phe样本名
temp = as.data.frame(str_split(string = phe$id, pattern = "-", simplify = T))
temp[,4] = str_remove_all(string = temp[,4], pattern = "A")

id = rnorm(nrow(temp))
for(i in 1:nrow(temp)) {
  id[i] = paste(temp[i,], sep = "-", collapse = "-")
}

phe$id = id

## sur样本名
temp = as.data.frame(str_split(string = sur$sample, pattern = "-", simplify = T))
temp[,4] = str_remove_all(string = temp[,4], pattern = "A")
temp[,4] = str_remove_all(string = temp[,4], pattern = "B")
temp[,4] = str_remove_all(string = temp[,4], pattern = "C")
temp[,4] = str_remove_all(string = temp[,4], pattern = "D")
temp[,4] = str_remove_all(string = temp[,4], pattern = "E")
temp[,4] = str_remove_all(string = temp[,4], pattern = "F")
temp[,4] = str_remove_all(string = temp[,4], pattern = "G")
temp[,4] = str_remove_all(string = temp[,4], pattern = "H")
table(temp[,4])

id = rnorm(nrow(temp))
for(i in 1:nrow(temp)) {
  id[i] = paste(temp[i,], sep = "-", collapse = "-")
}

sur$sample = id

# 取交集和取子集

## 取交集
id = intersect(colnames(exp), intersect(phe$id, sur$sample))

## 重命名
rownames(phe) = phe$id
sur = sur[!duplicated(sur$sample), ]
rownames(sur) = sur$sample

## 取子集
exp = exp[,id]
phe = phe[id,]
sur = sur[id,]

# 表型合并
identical(rownames(phe), rownames(sur))
phe = cbind(phe, sur[,3:4])

# 测序数据过滤
## ID转换
load(file = "./Rdata/gtf_gene.Rdata")
an = gtf_gene[,c("gene_name","gene_id","gene_type")]
rownames(exp) = str_split(rownames(exp), pattern = "[.]", simplify = T) [,1]
an$gene_id = str_split(an$gene_id, pattern = "[.]", simplify = T) [,1]

exp = exp[rownames(exp) %in% an$gene_id,]
an = an[match(rownames(exp),an$gene_id),]
identical(an$gene_id,rownames(exp))

k = !duplicated(an$gene_name);table(k)
an = an[k,]
exp = exp[k,]
rownames(exp) = an$gene_name

## 数据过滤
dim(exp)
exp = exp[apply(exp, 1, function(x) sum(x>0) > 0.5*ncol(exp)), ]
dim(exp)
exp[1:4,1:4]

# 数据保存
save(exp, phe, file = "./Rdata/TCGA-LUAD-LUSC-alk.Rdata")

