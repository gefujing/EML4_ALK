
# 准备环境
rm(list = ls())
options(stringsAsFactors = F)
library(forestplot)
library(survival)
library(survminer)
load(file = "./Rdata/tcga.alk.signature.Rdata")

# 数据为分类变量
df$age = ifelse(df$age >= median(df$age), "Older", "Young")
df$age = factor(df$age, levels = c("Young", "Older"))

df$race = ifelse(df$race == "white", "White", "Others")
df$race = factor(df$race, levels = c("White", "Others"))

df$N = ifelse(df$N == "N0", "N0", "N1/2")
df$N = factor(df$N, levels = c("N0", "N1/2"))

df$T = ifelse(df$T %in% c("T1","T2"), "T1/2", "T3/4")
df$T = factor(df$T, levels = c("T1/2", "T3/4"))

df$stage = ifelse(df$stage == "stage I", "stage I", "stage II/III/IV")
df$stage = factor(df$stage, levels = c("stage I", "stage II/III/IV"))

df$smoking_history = ifelse(df$smoking_history == 1, 0, 1)

df$diagnosis = factor(df$diagnosis, levels = unique(df$diagnosis))


# 数据预处理
rt = df[,c(1:14, ncol(df))]
rt = rt[,c(1:2, 4:6, 8:12, 15)]
data = data.frame()
Nsample = nrow(rt)

# 构建函数
UniCox = function(x, data){
  BaSurv = Surv(time = data$OS.time, event = data$OS)
  FML = as.formula(paste0("BaSurv~", x))
  GCox = coxph(FML, data = data)
  GSum = summary(GCox)
  HR = round(GSum$coefficients[,2],2)
  PValue = round(GSum$coefficients[,5],3)
  Low_CI = round(GSum$conf.int[,3],2)
  High_CI = round(GSum$conf.int[,4],2)
  
  ## 将结果转化为数据框返回
  Unicox = data.frame("characteristics" = x,
                      "Hazard Ratio" = HR,
                      "CI95_low" = Low_CI,
                      "CI95_high" = High_CI,
                      "P Value" = PValue)
  return(Unicox)
}

## 分析总人群中Cox分析结果
tmp_cox = UniCox(colnames(rt)[ncol(rt)], data = rt)

## 将总人群分析结果记录下来,形成一个数据框,命名为tmp
tmp = data.frame(Variable = "All patients",
                 Count = nrow(rt),
                 Percent = c(nrow(rt)/Nsample*100),
                 Point_Estimate = tmp_cox$Hazard.Ratio,
                 Low = tmp_cox$CI95_low,
                 High = tmp_cox$CI95_high,
                 Pvalue = tmp_cox$P.Value)

## 将tmp放入这个数据框种,后续亚组分析的结果也陆续放在这个数据框
data = rbind(data, tmp)

## 通过循环的方式,计算各个亚组中signature对于生存的影响
for(i in 1:c(ncol(rt) -3)){
  tmp_name = data.frame(Variable = colnames(rt)[i],
                        Count = NA,
                        Percent = NA,
                        Point_Estimate = NA,
                        Low = NA,
                        High = NA,
                        Pvalue = NA)
  
  data = rbind(data, tmp_name)
  
  tmp1 = unique(rt[,i]) %>% as.character()
  
  for(j in 1:length(tmp1)){
    
    tmp_df <- rt[rt[,i] == tmp1[j],]
    
    if(tmp_df[,ncol(tmp_df)] %>% unique() %>% length() == 1){
      next
    }else{
      
      tmp_r = UniCox(colnames(tmp_df)[ncol(tmp_df)], data = tmp_df)
      
      if(tmp_r$CI95_low %in% 'Inf' | tmp_r$CI95_high %in% 'Inf'){
        next
      }else{
        temp = data.frame(Variable = tmp1[j],
                          Count = nrow(tmp_df),
                          Percent = c(nrow(tmp_df)/Nsample*100) %>% round(.,2),
                          Point_Estimate = tmp_r$Hazard.Ratio,
                          Low = tmp_r$CI95_low,
                          High = tmp_r$CI95_high,
                          Pvalue = tmp_r$P.Value)
        data = rbind(data, temp)
      }
    }
    
  }
}

##对于有意义的结果,我们在P值后面标注上*
data$Pvalue[which(data$Pvalue < 0.05 & data$Pvalue >= 0.01)] = data$Pvalue[which(data$Pvalue < 0.05 & data$Pvalue >= 0.01)] %>% paste0(., "*")
data$Pvalue[which(data$Pvalue < 0.01 & data$Pvalue >= 0.001)] = data$Pvalue[which(data$Pvalue < 0.01 & data$Pvalue >= 0.001)] %>% paste0(., "**")
data$Pvalue[which(data$Pvalue < 0.001)] = data$Pvalue[which(data$Pvalue < 0.001)] %>% paste0(., "***")

data$Variable = as.character(data$Variable)
np = ifelse(!is.na(data$Count), paste(data$Count, "(", data$Percent, ")", sep = ""), NA)


#将要在图中展示的文本
tabletext  = cbind(c("\nSubgroup", NA, NA, data$Variable, NA),
                   
                   c("No. of\nPatients(%)", NA, NA, np, NA),
                   
                   c("Hazard Ratio\n(95% CI)", NA, NA, 
                     ifelse(!is.na(data$Count), 
                            paste(format(data$Point_Estimate, nsmall = 2), 
                                  "(", 
                                  format(data$Low, nsmall = 2), 
                                  " to ", 
                                  format(data$High, nsmall = 2), 
                                  ")", 
                                  sep=""),NA), 
                     NA),
                   c("P value", NA, NA, data$Pvalue, NA)
)
view(tabletext)

## 标定划线位置
index_N = nrow(data) + 5

## 计算HR可信区间界值
intervals = seq(min(data$Low %>% na.omit()),
                max(data$High %>% na.omit(),
                    by = c(min(data$Low %>% na.omit()) + max(data$High %>% na.omit()))/10))

tmp_hrzl = list(gpar(lwd = 2, col = "black"), gpar(lwd = 2, col = "black"))
names(tmp_hrzl) = c("3", as.character(index_N))

## 画图

pdf(file = "./Rplot/alk.Subgroup.ForestPlot.pdf", width = 24, height = 25, onefile=FALSE)

forestplot(labeltext = tabletext, #图中的文本
           mean = c(NA,NA,1,data$Point_Estimate,NA),#HR
           lower = c(NA,NA,1,data$Low,NA),#95%置信区间下限
           upper = c(NA,NA,1,data$High,NA),#95%置信区间上限
           Pvalue = c(NA,NA,1,data$Pvalue,NA),#Pvalue
           #title="Hazard Ratio",
           graph.pos = 3,#图在表中的列位置
           graphwidth = unit(.4,"npc"),#图在表中的宽度比例
           fn.ci_norm="fpDrawNormalCI",#box类型选择钻石
           col=fpColors(box="steelblue",lines="black",zero ="black"),#box 颜色
           boxsize=0.7, #box大小根据样本量设置
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           grid = structure(c(data[1,]$Point_Estimate),gp = gpar(col ="black",lty=2,lwd=2)),#虚线(可多条)及其横坐标、颜色、线宽
           # xticks = intervals,
           lwd.xaxis=2,#X轴 线宽
           x1ab="<-Favors Low Score        Favors High Score->",  #X轴 标题
           hrzl_1ines = tmp_hrzl,
           txt_gp=fpTxtGp(label=gpar(cex=1.25), #各种字体大小设置
                          ticks=gpar(cex=1.25), 
                          xlab=gpar(cex = 1.25),
                          title=gpar(cex = 1.25)),
           is.Summary = c(T,rep(F,as.numeric(index_N))),#首行字体类型设置
           lineheight = unit(1,"cm"),#固定行高
           #cex=10,
           colgap = unit(0.5,"cm"),#列之间的间隙
           mar=unit(rep(1.25,times= 4),"cm"),#图形页边距
           new_page = F #是否新页
)


dev.off()

