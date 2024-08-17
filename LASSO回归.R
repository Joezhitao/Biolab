#加载包
library("glmnet")
library("survival")

setwd("E:/CORD/LASSO")     #设置工作目录

#读入生存数据
cli=read.table("time_CORD.txt", header=T, sep="\t", check.names=F, row.names=1)
#生存时间以年为单位
cli$time=cli$time/365

#读入
data=read.table("TCGA_CORD_TPM.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#仅保留肿瘤样本
data = data[,group == 0]

#读取单因素cox基因的表达矩阵
unicoxgene=read.table("uniCox.txt", header=T, sep="\t", check.names=F,row.names = 1)
data = data[rownames(unicoxgene),]

#转置
data=t(data)
#样本名仅保留前12字符
rownames(data)=substr(rownames(data),1,12)
#将.改为-
rownames(data) = gsub('[.]', '-', rownames(data))
#获取共同样本
sameSample=intersect(row.names(data), row.names(cli))
#提取共同样本
data=data[sameSample,]
cli=cli[sameSample,]
#合并
rt=cbind(cli, data)

#设置随机种子，K折交叉验证，每次的数据都是随机的，随机数种子一致，就结果给固定住。
set.seed(123456)

#构建lasso回归模型
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$time, rt$state))
fit=glmnet(x, y, family="cox", nfolds = 10)

#family="gaussian" 适用于一维连续因变量（univariate）
#family="mgaussian" 适用于多维连续因变量（multivariate）
#family="poisson" 适用于非负次数因变量（count）
#family="binomial" 适用于二元离散因变量（binary）
#family="multinomial" 适用于多元离散因变量（category）

#c-index
cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C", nfolds = 10)
#pdf("lasso.c-index.pdf")
png(file="lasso.c-index.png", width=8, height=8, units="in", res=800)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()

#deviance图形
cvfit=cv.glmnet(x, y, family="cox", type.measure = "deviance", nfolds = 10)
#pdf("lasso.cvfit.pdf")
png(file="lasso.cvfit.png", width=8, height=8, units="in", res=800)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
dev.off()

#Coefficients图形
#pdf("lasso.lambda.pdf")
png(file="lasso.lambda.png", width=8, height=8, units="in", res=800)
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
lassoSigExp=rt[,c("time", "state", lassoGene)]
lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
