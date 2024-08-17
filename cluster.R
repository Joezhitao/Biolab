#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ConsensusClusterPlus")

#引用包
library(limma)
library(survival)
library(survminer)
library(ConsensusClusterPlus)

rm(list=ls())
#设置工作目录
workDir="E:/CORD/cluster"
setwd(workDir)

#读入
data=read.table("TCGA_CORD_TPM.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#根据正常和肿瘤排序
data = data[,group == 0]

#样本名仅保留前12字符
colnames(data)=substr(colnames(data),1,12)
#将.改为-
colnames(data) = gsub('[.]', '-', colnames(data))

#读取多因素cox基因的表达矩阵
multicoxgene=read.table("multiCox.txt", header=T, sep="\t", check.names=F,row.names = 1)
data = data[rownames(multicoxgene),]

#对样品进行分型
#maxK, 最大的K值，形成一系列梯度
#pItem, 选择80%的样本进行重复抽样
#pfeature, 选择80%的基因进行重复抽样
#reps, 重复抽样的数目，先设置为100，结果不错再设置为1000
#clusterAlg, 聚类的算法,hc,pam,km
#distanc, 距离矩阵的算法,pearson,spearman,euclidean
#title, 输出结果的文件夹名字，包含了输出的图片
#seed, 随机种子，用于固定结果

results=ConsensusClusterPlus(data,
                             maxK=9,
                             reps=100,
                             pItem=0.8,
                             pFeature=1,
                             title=workDir,
                             clusterAlg="pam",
                             distance="euclidean",
                             seed=123,
                             plot="png")

#输出最佳K值
Kvec <- 2:3
x1 <- 0.1; x2 <- 0.9
PAC <- rep(NA, length(Kvec))
names(PAC) <- paste("K=", Kvec, sep="")
for (i in Kvec) {
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])  #计算出共识矩阵
  PAC[i] = Fn(x2) - Fn(x1)
}

optK <- Kvec[which.min(PAC)]
optK

#提取分群信息
duplicated_colnames <- colnames(data)[duplicated(colnames(data))]
print(duplicated_colnames)
colnames(data) <- make.unique(colnames(data))


annCol <- data.frame(Cluster = paste0("Cluster",
                                      results[[optK]][['consensusClass']]),
                     row.names = colnames(data))

table(annCol$Cluster)
identical(rownames(annCol), colnames(data))


#PCA分析，看样本的分布情况
library(FactoMineR)
library(factoextra)

exp <- as.data.frame(t(data))
dat <- PCA(exp, graph = FALSE)
pca <- fviz_pca_ind(dat, 
                    geom.ind =  "point", 
                    col.ind = annCol$Cluster,
                    addEllipses = TRUE ,
                    pointsize = 2.5,
                    ellipse.type = "none",  #外边框的类型，凸多边形"convex"
                    axes.linetype = "blank",  #hide zeroline
                    legend.title = "group")

#pdf("PCA.pdf", width=10, height=10)
png("PCA.png", width=10, height=10, units = "in", res = 800)
print(pca)
dev.off()

#输出分型结果
clusterNum=2      #分成几个亚型
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("cluster")
letter=c("A","B","C","D","E","F","G","H","I","J")
uniqClu=levels(factor(cluster$cluster))
cluster$cluster=letter[match(cluster$cluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="cluster.txt", sep="\t", quote=F, col.names=F)

##############################################################################
#聚类后生存分析
#读入生存数据
cli=read.table("time_CORD.txt", header=T, sep="\t", check.names=F, row.names=1)
#生存时间以年为单位
cli$time=cli$time/365

cli <- cli[rownames(annCol),]
identical(rownames(cli), rownames(annCol))
#去掉cli有NA值的行
cli <- cli[complete.cases(cli),]

# 删除没有的行名

# 找到共同的行名
common_rows <- intersect(rownames(cli), rownames(annCol))

# 子集化两个数据框，保留共同的行
cli_filtered <- cli[common_rows, , drop = FALSE]  # 使用 drop = FALSE 以保持数据框的结构
annCol_filtered <- annCol[common_rows, , drop = FALSE]

# 检查结果
identical(rownames(cli_filtered), rownames(annCol_filtered))  # 应该返回 TRUE

#合并
rt <- cbind(cli_filtered, annCol_filtered)
fitd <- survdiff(Surv(time, state) ~ Cluster, data = rt)
pValue <- fitd$pvalue
#查看p值
pValue
#拟合生存曲线
fit <- survfit(Surv(time, state) ~ Cluster, data = rt)
summary(fit)

pround <- round(pValue, 6)

length = length(levels(factor(rt$Cluster)))

bioCol=c("Firebrick3","MediumSeaGreen","#6E568C","#223D6C")
bioCol=bioCol[1:length]
p <- ggsurvplot(fit, 
                   data=rt,
                   conf.int=F,
                   pval=pround,
                   pval.size=6,
                   legend.title="cluster",
                   legend.labs=levels(factor(rt[,"Cluster"])),
                   legend = c(0.8, 0.8),
                   font.legend=10,
                   xlab="Time(years)",
                   break.time.by = 2,
                   palette = bioCol,
                   surv.median.line = "hv",
                   risk.table=T,
                   cumevents=F,
                   risk.table.height=.25)

#pdf(file=paste("survial.pdf",sep = ""), width=6.5, height=6.25, onefile=FALSE)
png("survial.png", width=6.5, height=6.25, units = "in", res = 800)
print(p)
dev.off()
