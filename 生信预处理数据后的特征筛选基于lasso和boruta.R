##############################################################################
#基于差异基因分析第一步结果提取基因
#加载包
library("glmnet")
library("survival")
# 设置工作目录
base_dir <- "E:/CRC_model"
setwd(base_dir)
# 读取只有一列且有表头的文件
DEG_data <- read.table("intersectGenes_DESeq2.txt", header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
DEG_gene <- DEG_data$ID
#读取TCGA与临床合并后数据
#读取数据
TCGA <- read.table("E:/CRC_model/PC/TCGA_cli.txt", header=T, sep="\t", check.names=F, row.names=1)

# 定义要保留的列，使用 DEG_gene 和临床数据列
keep_columns <- c("time", "status", DEG_gene)

# 检查所有需要的列是否存在
missing_columns <- setdiff(keep_columns, colnames(TCGA))
if (length(missing_columns) > 0) {
  warning(paste("以下列不存在于数据中:", paste(missing_columns, collapse = ", ")))
  keep_columns <- intersect(keep_columns, colnames(TCGA))
}

# 只保留指定的列
final_data <- TCGA[, keep_columns]

# 查看结果
dim(final_data)
head(final_data)
###############################################################################
#lasso回归
#设置随机种子，K折交叉验证，每次的数据都是随机的，随机数种子一致，就结果给固定住。
rt <- final_data
set.seed(123456)
#构建lasso回归模型
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$time, rt$status))
fit=glmnet(x, y, family="cox", nfolds = 10)
#c-index
cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C", nfolds = 10)
#pdf("lasso.c-index.pdf")
png(file="lasso.c-index_limma.png", width=8, height=8, units="in", res=800)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()
#deviance图形
cvfit=cv.glmnet(x, y, family="cox", type.measure = "deviance", nfolds = 10)
#pdf("lasso.cvfit.pdf")
png(file="lasso.cvfit_limma.png", width=8, height=8, units="in", res=800)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()
#pdf("lasso.lambda.pdf")
png(file="lasso.lambda_limma.png", width=8, height=8, units="in", res=800)
plot(fit, xvar="lambda", label=T)
abline(v=log(cvfit$lambda.min), lty="dashed")
dev.off()
#输出lasso显著基因表达量
coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
row.names(coef)[index]
actCoef
index
lassoSigExp=rt[,c("time", "status", lassoGene)]
lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
write.table(lassoSigExpOut,file="lasso_SigExp_limma.txt",sep="\t",row.names=F,quote=F)
###############################################################################
#Boruta 算法进一步筛选变量
# Joe_2024_10_12
rm(list=ls())
# 加载包
library(tidyverse)
library(survival)
library(Boruta)

# 设置工作目录
base_dir <- "E:/CRC_model"
setwd(base_dir)

# 读取数据
sadata <- read.table("lasso_SigExp_limma.txt", header=TRUE, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
#sadata <- readr::read_csv("F:/Mach_learn_data/sadata.csv")
colnames(sadata)

# 修正变量类型
# 分类变量转factor
#for(i in c(3)){
#  sadata[[i]] <- factor(sadata[[i]])
#}

# 剔除变量-无关变量
#sadata$pid <- NULL
# 剔除样本-含有缺失值的样本，时间非正的样本
# sadata <- na.omit(sadata) # 剔除任意变量下有缺失值的样本
sadata <- sadata %>%
  # na.omit() %>%
  # drop_na(age) %>%  # 剔除指定变量下有缺失值的样本
  filter(time > 0) # 剔除时间非正的样本

#sadata删除id那一列
sadata <- sadata[,-1]

# 数据概况
skimr::skim(sadata)

# 感兴趣的时间节点
range(unique(sadata$time))
itps <- c(1 * c(1, 3, 5))
itps
table(cut(sadata$time, c(0, itps, Inf)))


###################################################

# 数据拆分构建任务对象
set.seed(42)
datasplit <- rsample::initial_split(
  sadata, prop = 0.6, strata = time, breaks = 10
)
traindata <- rsample::training(datasplit) %>%
  sample_n(nrow(.))
testdata <- rsample::testing(datasplit) %>%
  sample_n(nrow(.))

# 拆分得到的数据的生存曲线比较
sadata2 <- sadata
sadata2$set <- "Test"
sadata2$set[datasplit$in_id] <- "Train"
sadata2$set <- factor(sadata2$set)
sfit_set <- survfit(Surv(time, status) ~ set, data=sadata2)
survminer::ggsurvplot(
  sfit_set,
  pval=TRUE, 
  pval.coord = c(0.1, 0.5),
  risk.table=TRUE,
  ggtheme = survminer::theme_survminer() +
    theme(text = element_text(family = "serif")),
  font.family = "serif"
)

###################################################

# 数据预处理
datarecipe_coxph <- 
  recipes::recipe(time + status ~ ., traindata) %>%
  recipes::prep()
datarecipe_coxph

# 按方处理训练集和测试集
traindata2 <- 
  recipes::bake(datarecipe_coxph, new_data = NULL) %>%
  dplyr::select(time, status, everything())
testdata2 <- 
  recipes::bake(datarecipe_coxph, new_data = testdata) %>%
  dplyr::select(time, status, everything())

# 特征筛选
set.seed(42)
result_boruta <- Boruta(
  x = traindata2[,-c(1,2)],
  y = Surv(traindata2$time, traindata2$status), 
  doTrace = 1,
  maxRuns = 100
)
# 筛选结果概况
result_boruta

# 对于存疑变量进一步确认
result_boruta <- TentativeRoughFix(result_boruta)
# 具体结果
attStats(result_boruta)
# 图示
par(las = 3, mar = c(8, 4, 4, 2))
plot(result_boruta, xlab = "", main = "Boruta算法筛选结果")
legend("topleft", 
       legend = c("confirmed", "rejected"),
       pch = 22,
       pt.bg = c("green", "red"))

# 筛选得到的变量
getSelectedAttributes(result_boruta)
cgene <- getSelectedAttributes(result_boruta)
# 筛选之后的数据
traindata <- traindata2 %>%
  dplyr::select(time, status, 
                getSelectedAttributes(result_boruta))
testdata <- testdata2 %>%
  dplyr::select(time, status, 
                getSelectedAttributes(result_boruta))