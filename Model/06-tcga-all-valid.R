
# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)

# 数据读取
## 读取表达量
load("./Rawdata/TCGA-LUAD-LUSC-exp.Rdata")
phe = read.csv(file = "./Rawdata/TCGA-LUAD-LUSC-phe-some.csv")
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
save(exp, phe, file = "./Rdata/TCGA-LUAD-LUSC-all.Rdata")


# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggsignif)
library(survival)
library(survminer)
library(pheatmap)
library(tinyarray)
load(file = "./Rdata/TCGA-LUAD-LUSC-all.Rdata")


# 数据预处理
## 去除不需要的表型信息
phe = phe[, -1]

## 转置表达量矩阵
exp = as.data.frame(t(exp))
identical(rownames(exp), rownames(phe))

## 联合表型
df = cbind(phe, exp)
df = df[,!duplicated(colnames(df))]

## 保存数据
save(df, file = "./Rdata/all.rawdata.Rdata")


# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggsignif)
library(survival)
library(survminer)
library(pheatmap)
library(tinyarray)

load(file = "./Rdata/all.rawdata.Rdata")
load(file = "./Rdata/tcga.model2.Rdata")
Final_Getfactors = c("CXCL1", "CXCL5", "IL6")


## 模型预测
fp <- predict(model, df)
fivenum(fp)
df$signiture = fp
save(df, file = "./Rdata/tcga.all.signature.Rdata")


## 划分高低风险
names(fp) = df$sample
ri = ifelse(fp<median(fp),"lowrisk","highrisk")
ri = factor(ri,levels = c("lowrisk","highrisk"))

sfit <- survfit(Surv(OS.time, OS)~ri, data = df)

ggsurvplot(sfit,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =F,xlab ="Time/days", 
           ggtheme =theme_light(), 
           ncensor.plot = F)


## 风险因子三图
fp_dat = data.frame(patientid=1:length(fp),
                    fp=as.numeric(sort(fp)),
                    ri= ri[order(fp)])

sur_dat = data.frame(patientid=1:length(fp),
                     time=df[order(fp),"OS.time"] ,
                     event=df[order(fp),"OS"]) 

sur_dat$event = ifelse(sur_dat$event==0,'Alive','Death')
sur_dat$event = factor(sur_dat$event,levels = c("Death","Alive"))
exp_dat = t(df[order(fp),Final_Getfactors])

### 第一个图-风险因子
p1 = ggplot(fp_dat,aes(x=patientid,y=fp))+
  geom_point(aes(color = ri))+
  scale_color_manual(values = c("blue","red"))+
  geom_vline(xintercept = 0.5*nrow(sur_dat),lty = 2)+
  scale_x_continuous(expand=c(0.01,0))+
  labs(x = "",y = "risk score")+
  theme_bw()

### 第二个图-结局事件图
p2 = ggplot(sur_dat,aes(x=patientid,y=time))+
  geom_point(aes(col=event))+
  scale_color_manual(values = c("red","blue"))+
  geom_vline(xintercept = 0.5*nrow(sur_dat),lty = 2)+
  scale_x_continuous(expand=c(0.01,0))+
  labs(x = "")+
  theme_bw()

## 第三个图-基因表达图
n = scale(t(exp_dat))
n[n>3] = 3
n[n<(-3)] = -3
p3 = ggheat(n,fp_dat$ri,show_rownames = F,legend_color= c("blue","red"),color = c("blue","white","red"))+
  theme(axis.text = element_text(size = 8))
p3

### 拼图
p1 /p2 /p3 + plot_layout(design = 'A
                         B
                         C
                         C
                         C
                         C')
ggsave(filename = "./Rplot/alL.modle.risk.pdf", width = 5, height = 6)


# TimeROC
df$sample = rownames(df)
dat = cbind(df[ ,c("sample", "OS", "OS.time")],
            fp = fp)

library(survminer)
library(survival)
library(timeROC)
result <-with(dat, timeROC(T=OS.time,
                           delta=OS,
                           marker=fp,
                           cause=1,
                           times=c(12*30,24*30,36*30),
                           iid = TRUE))
plot_dat = data.frame(fpr = as.numeric(result$FP),
                      tpr = as.numeric(result$TP),
                      time = rep(as.factor(c(12*30,24*30,36*30)),each = nrow(result$TP)))

library(ggplot2)
ggplot() + 
  geom_line(data = plot_dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,2,3),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
ggsave(filename = "./Rplot/alL.tcga.modle.roc.pdf", width = 9, height = 4)



