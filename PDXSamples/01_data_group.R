
# 导入数据，并整理
rm(list = ls())
dir()
exp = read.csv(file = "./RawData/ALK_Brigatinib.csv", header = T)
exp = aggregate(.~ID, mean, data = exp)
# exp = exp[!duplicated(exp$ID),]
rownames(exp) = exp$ID
exp = exp[, -1]

# 分组
Group = c("Control", "Control",
          "Brigatinib_15", "Brigatinib_15",
          "Brigatinib_30", "Brigatinib_30")
Group = factor(Group, levels = c("Control", "Brigatinib_15", "Brigatinib_30"))

# 清洗无效数据
dim(exp)
exp = exp[apply(exp, 1, function(x) sum(x>0) > 0.5*ncol(exp)), ]
dim(exp)
exp[1:4,1:4]

# 保存数据
save(exp, Group, file = "./Rdata/exp_Group.RData")




































