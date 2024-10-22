#加载
library(survival)
library(survminer)
library(xlsx)

#工作路径
setwd("E:/CRC_model")

library(VennDiagram)
#设置工作目录
setwd("E:/CRC_model")
#读取
DEG <- c("Wilcoxon","limma","edgeR","DESeq2")
for (i in DEG) {
  data=read.table(paste("TCGA.diff.",i,".txt",sep=""), header=F, sep="\t", check.names=F)
  DEG = data[,1]
  data=read.table("Mismatch_Repair_gene.txt", header=F, sep="\t", check.names=F)
  ARG = data[,1]
  
  #画图
  venn.diagram(x = list('ARG' = ARG,i = DEG),
               filename = paste("VN_",i,".png",sep = ""),fill = c("darkorange1","green"))
  #交集基因
  intersectGenes = intersect(DEG,ARG)
  
  write.table(file= paste("intersectGenes_",i,".txt",sep = "") , intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)
}

#工作路径
setwd("E:/CRC_model")
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

rt <- final_data
#0.05,0.01,0.001
p.value = 0.05
#查找预后相关的基因
outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
  #cox分析
  cox <- coxph(Surv(time, status) ~ rt[,i], data = rt)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  if(coxP < p.value){
    outTab=rbind(outTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
}

#输出单因素的结果
write.table(outTab, file="uniCox.txt", sep="\t", row.names=F, quote=F)


#画图
#读取输入文件
rt <- read.table("uniCox.txt",header=T,sep="\t",check.names=F,row.names=1)
#仅仅做前20个,实际运行不运行这一行
#rt = rt[1:20,]
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#绘制森林图
pdf(file="uniCoxforest.pdf", width = 12,height = nrow(rt)/13+15)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "red", "blue")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
axis(1)
dev.off()

rt <- read.table("uniCox.txt",header=T,sep="\t",check.names=F,row.names=1)
#rt = rt[1:20,]
gene <- rownames(rt)
#读入
rt=read.table("TCGA_cli.txt",header=T,sep="\t",check.names=F,row.names=1)   
# 定义要保留的列，使用 DEG_gene 和临床数据列
keep_columns <- c("time", "status", gene)

# 检查所有需要的列是否存在
missing_columns <- setdiff(keep_columns, colnames(rt))
if (length(missing_columns) > 0) {
  warning(paste("以下列不存在于数据中:", paste(missing_columns, collapse = ", ")))
  keep_columns <- intersect(keep_columns, colnames(TCGA))
}

# 只保留指定的列
rt <- rt[, keep_columns]

# 查看结果
dim(final_data)
head(final_data)

#COX模型构建
multiCox=coxph(Surv(time, status) ~ ., data = rt)
#参数多，选取逐步回归方式，step，"both", "backward", "forward"
#参数少，不运行此行
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#输出模型参数
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)
write.table(outTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)

#输出病人风险值
riskScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("time","status",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="risk.txt",
            sep="\t",
            quote=F,
            row.names=F)

#画图
#读取输入文件
rt <- read.table("multiCox.txt",header=T,sep="\t",check.names=F,row.names=1)
gene <- rownames(rt)
hr <- sprintf("%.3f",rt$"HR")
hrLow  <- sprintf("%.3f",rt$"HR.95L")
hrHigh <- sprintf("%.3f",rt$"HR.95H")
Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))

#绘制森林图
pdf(file="multiCoxforest.pdf", width = 7,height = nrow(rt)/13+5)
n <- nrow(rt)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))

xlim = c(0,3)
par(mar=c(4,2.5,2,1))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=0.8
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)

par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, "red", "blue")
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
axis(1)
dev.off()


library(tidyverse)
library(survival)
library(Boruta)

gene <- rownames(rt)
#读取数据
rt <- TCGA
keep_columns <- c("time", "status", gene)
rt <- rt[, keep_columns]
traindata <- rt

testdata <- read.table("E:/CRC_model/PC/GEO_cli.txt", header=T, sep="\t", check.names=F, row.names=1)
testdata <- testdata[, keep_columns]
# 修正变量类型
# 分类变量转factor
#for(i in c(3)){
#  sadata[[i]] <- factor(sadata[[i]])
#}

# 剔除变量-无关变量
#sadata$pid <- NULL
# 剔除样本-含有缺失值的样本，时间非正的样本
# sadata <- na.omit(sadata) # 剔除任意变量下有缺失值的样本
traindata <- traindata %>%
  # na.omit() %>%
  # drop_na(age) %>%  # 剔除指定变量下有缺失值的样本
  filter(time > 0) # 剔除时间非正的样本


# 数据概况
skimr::skim(traindata)

# 感兴趣的时间节点
range(unique(traindata$time))
itps <- c(1 * c(1, 3, 5))
itps
table(cut(traindata$time, c(0, itps, Inf)))


###################################################

# 数据拆分构建任务对象

traindata$set <- "Train"
testdata$set <- "Test"

# 拆分得到的数据的生存曲线比较
sadata2 <- rbind(traindata, testdata)
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
traindata <- rt
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

# 筛选之后的数据
traindata <- traindata2 %>%
  dplyr::select(time, status, 
                getSelectedAttributes(result_boruta))
testdata <- testdata2 %>%
  dplyr::select(time, status, 
                getSelectedAttributes(result_boruta))


##########################################################################
library(glmnet)
library(dplyr)

raw <- traindata

num <- length(raw)
x <- raw[,3:num]

y <- raw[,1:2]
names(y) <- c('time','status')

#z-score
# x <- sapply(x, function(x) (x-mean(x))/sd(x))
#

#使用c-index指标进行交叉验证分析；
cvfit <- cv.glmnet(as.matrix(x), as.matrix(y), family = "cox", type.measure = "C", nfolds = 10,intercept= FALSE)
# plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se

#截距 intercept 即系数中的第一个
coef(cvfit)
#获取系数列表；
coefficient <- coef(cvfit,s=cvfit$lambda.min)
Active.Index <- which(as.numeric(coefficient)!=0)
active.coefficients <- as.numeric(coefficient)[Active.Index]
sig_multi_cox <- rownames(coefficient)[Active.Index]
#
active.coefficients
sig_multi_cox

coef_dataframe = data.frame("features" = sig_multi_cox,'coef' = active.coefficients)
coef_dataframe

feature_names <- coef_dataframe[, 1]
coefficients <- coef_dataframe[, 2]

# 计算训练集的风险评分
traindata$risk_score <- rowSums(as.matrix(traindata[, feature_names, drop = FALSE]) * coefficients)

# 计算测试集的风险评分
testdata$risk_score <- rowSums(as.matrix(testdata[, feature_names, drop = FALSE]) * coefficients)


# 计算训练集风险评分的中位数
train_cutoff <- median(traindata$risk_score, na.rm = TRUE)

# 应用 cutoff，将样本分为高风险和低风险
traindata$risk_group <- ifelse(traindata$risk_score > train_cutoff, "High", "Low")
testdata$risk_group <- ifelse(testdata$risk_score > train_cutoff, "High", "Low")

# 绘制生存曲线
library(survminer)

# 训练集生存曲线
surv_train <- survfit(Surv(time, status) ~ risk_group, data = traindata)
ggsurvplot(surv_train, data = traindata, pval = TRUE, title = "Training Set Survival Curve")

# 测试集生存曲线
surv_test <- survfit(Surv(time, status) ~ risk_group, data = testdata)
ggsurvplot(surv_test, data = testdata, pval = TRUE, title = "Test Set Survival Curve")

###############################################################################################################################
####决策树
library(tidyverse)
library(survival)
library(mlr3verse)
library(mlr3proba)
library(mlr3extralearners)
source("H:/Biolab/Biolab/tidyfuncs4sa.R")
# file.choose()
#testdata <- write.table(testdata,"E:/CRC_model/PC/testdata.txt",sep="\t",row.names=F,quote=F)
#traindata <- write.table(traindata,"E:/CRC_model/PC/traindata.txt",sep="\t",row.names=F,quote=F)


#testdata <- read.table("E:/CRC_model/PC/testdata.txt", header = TRUE, sep = "\t", check.names = FALSE)

#traindata <- read.table("E:/CRC_model/PC/traindata.txt", header = TRUE, sep = "\t", check.names = FALSE)



# 感兴趣的时间节点
range(unique(traindata$time))
itps <- c(1 * c(1, 3, 5))
itps
table(cut(traindata$time, c(0, itps, Inf)))

# 性能指标
measure_sa <- msrs("surv.cindex")
measure_sa

###################################################
# 训练任务对象
traindata <- traindata[,1:7]
testdata <- testdata[,1:7]
task_train <- as_task_surv(
  traindata, 
  time = "time",
  event = "status", 
  type = "right"
)
task_train
# 测试任务对象
task_test <- as_task_surv(
  testdata, 
  time = "time",
  event = "status", 
  type = "right"
)
task_test

###################################################

# 决策树
# https://mlr3extralearners.mlr-org.com/reference/mlr_learners_surv.ctree.html

# 设定模型
learner_ctree <- lrn(
  "surv.ctree",
  alpha = to_tune(c(0.01, 0.05, 0.1)),
  minbucket = to_tune(1, 25)
)
learner_ctree$id <- "ctree"
learner_ctree

# 多核并行
future::plan("multisession")

# 超参数调优设定
set.seed(42)
tune_ctree <- tune(
  tuner = tnr(
    "grid_search", 
    param_resolutions = c(minbucket = 5)
  ),
  task = task_train,
  learner = learner_ctree,
  resampling = rsmp("cv", folds = 5),
  measure = msr("surv.cindex"),
  terminator = trm("none")
)
tune_ctree
as.data.table(tune_ctree$archive) %>%
  as.data.frame() %>%
  select(1:3) %>%
  plotly::plot_ly(
    type = 'parcoords',
    line = list(color = ~surv.cindex, 
                colorscale = 'Jet', 
                showscale = T),
    dimensions = list(
      list(label = 'alpha', values = ~alpha),
      list(label = 'minbucket', values = ~minbucket)
    )
  ) %>%
  plotly::layout(title = "DT HPO Guided by C-Index",
                 font = list(family = "serif"))


# 训练最终模型
learner_ctree$param_set$values <-
  tune_ctree$result_learner_param_vals
set.seed(42)
learner_ctree$train(task_train)
learner_ctree

# 模型概况
learner_ctree
learner_ctree$model
plot(learner_ctree$model)

###################################################

# 预测训练集

summary(predprobtrain_ctree)

predtrain_ctree <- learner_ctree$predict(task_train)
predprobtrain_ctree <- predprob(
  pred = predtrain_ctree, 
  preddata = traindata, 
  etime = "time",
  estatus = "status",
  model = "ctree", 
  dataset = "train", 
  timepoints =itps
)

print("Check predprobtrain_ctree:")
print(head(predprobtrain_ctree))
print(summary(predprobtrain_ctree$time))
print(summary(predprobtrain_ctree$status))


# 性能指标
evaltrain_ctree <- eval4sa(
  predprob = predprobtrain_ctree,
  preddata = traindata,
  etime = "time",
  estatus = "status",
  model = "ctree",
  dataset = "train",
  timepoints = itps,
  plotcalimethod = "nne",  # quantile
  bw4nne = 0.01,
  q4quantile = 5,
  cutoff = "median"
)
evaltrain_ctree$auc
evaltrain_ctree$rocplot
evaltrain_ctree$brierscore
evaltrain_ctree$brierscoretest
evaltrain_ctree$calibrationplot
evaltrain_ctree$riskplot

sadca(
  predprob = predprobtrain_ctree,
  preddata = traindata,
  etime = "time",
  estatus = "status",
  model = "ctree",
  dataset = "train",
  timepoints = itps,
  timepoint = 1, 
  xrange = 0:100 / 100
)

# 预测测试集
predtest_ctree <- learner_ctree$predict(task_test)
predprobtest_ctree <- predprob(
  pred = predtest_ctree, 
  preddata = testdata, 
  etime = "time",
  estatus = "status",
  model = "ctree", 
  dataset = "test", 
  timepoints =itps
)
# 性能指标
predtest_ctree$score(measure_sa)
cindex_bootci(learner_ctree, testdata)

evaltest_ctree <- eval4sa(
  predprob = predprobtest_ctree,
  preddata = testdata,
  etime = "time",
  estatus = "status",
  model = "ctree",
  dataset = "test",
  timepoints = itps,
  plotcalimethod = "nne",  # quantile
  bw4nne = 0.01,
  q4quantile = 5,
  cutoff = "median"
)
evaltest_ctree$auc
evaltest_ctree$rocplot
evaltest_ctree$brierscore
evaltest_ctree$brierscoretest
evaltest_ctree$calibrationplot
evaltest_ctree$riskplot

sadca(
  predprob = predprobtest_ctree,
  preddata = testdata,
  etime = "time",
  estatus = "status",
  model = "ctree",
  dataset = "test",
  timepoints = itps,
  timepoint = 1, 
  xrange = 0:100 / 100
)

