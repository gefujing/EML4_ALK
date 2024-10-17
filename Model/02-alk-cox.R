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
load(file = "./Rdata/TCGA-LUAD-LUSC-alk.Rdata")


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
save(df, file = "./Rdata/alk.rawdata.Rdata")


# 设置环境
rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(tidyverse)
library(plyr)
library(survminer)
load(file = "./Rdata/alk.rawdata.Rdata")

## 写一个函数，方便后续批量计算
BaSurv = Surv(time = df$OS.time,
              event = df$OS)

UniCox = function(x){
  FML = as.formula(paste0("BaSurv~", x))
  GCox = coxph(FML, data = df)
  GSum = summary(GCox)
  HR = round(GSum$coefficients[,2],2)
  PValue = round(GSum$coefficients[,5],3)
  CI = paste0(round(GSum$conf.int[,3:4],2), collapse = "-")
  Unicox = data.frame("characteristics" = x,
                      "Hazard Ratio" = HR,
                      "CI95" = CI,
                      "P Value" = PValue)
  return(Unicox)
}

## 挑一个基因试一下
UniCox(colnames(df)[15])

# 设置基因
finalgenes = read.csv("./Rawdata/mcode.csv")
finalgenes = finalgenes$name

# 单因素分析
VarNames = finalgenes
UniVar = lapply(VarNames, UniCox)
UniVar = ldply(UniVar, data.frame)
GetFactors = UniVar$characteristics[which(UniVar$P.Value < 0.3)] %>% as.character()
# GetFactors = VarNames

# 多因素回归
fml = as.formula(paste0("BaSurv~", paste0(GetFactors, collapse = "+")))
Multicox = coxph(fml, data = df)
MultiSum = summary(Multicox)


MultiName = as.character(GetFactors)
MHR = round(MultiSum$coefficients[,2],2)
MPValue = round(MultiSum$coefficients[,5],3)
MCIL = round(MultiSum$conf.int[,3],2)
MCIU = round(MultiSum$conf.int[,4],2)
MCI = paste0(MCIL, "-", MCIU)
Mulcox = data.frame("characteristics" = MultiName,
                    "Hazard Ratio" = MHR,
                    "CI95" = MCI,
                    "P Value" = MPValue)

## 根据多因素COX回归中p值小于0.2选择
Final = merge.data.frame(UniVar, Mulcox, by = "characteristics", all = T, sodf = T)
Final_Getfactors = Final$characteristics[which(Final$P.Value.y <= 0.2)] %>% as.character()
# Final_Getfactors = VarNames
# Final_Getfactors = c("CXCL5", "IL6")

## 构建模型
fml = as.formula(paste0("BaSurv~", paste0(Final_Getfactors, collapse = "+")))
model = coxph(formula =fml, data = df)
ggforest(model, data = df)
ggsave(filename = "./Rplot/tcga.alk.genes.modle.pdf", width = 6, height = 4)
save(model, file = "./Rdata/tcga.model2.Rdata")

## 模型预测
fp <- predict(model, df)
fivenum(fp)
df$signiture = fp
save(df, file = "./Rdata/tcga.alk.signature.Rdata")

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
ggsave(filename = "./Rplot/alk.modle.risk.sur.pdf", width = 4.5, height = 7.5)


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
ggsave(filename = "./Rplot/alk.modle.risk.pdf", width = 5, height = 6)


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
ggsave(filename = "./Rplot/alk.tcga.modle.roc.pdf", width = 9, height = 4)



