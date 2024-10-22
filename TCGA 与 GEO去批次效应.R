#引用包
library(limma)
library(sva)

#设置工作目录
setwd("E:/CRC_model/PC")

#输入文件名称
files=c("TCGA_CRC_TPM.txt", "GSE71187.txt")

geneList=list()
for (i in 1:length(files)) {
  #读入，转化为matrix
  rt=read.table(files[i], header=T, sep="\t",check.names=F,row.names = 1)
  dimnames=list(rownames(rt), colnames(rt))
  rt=matrix(as.numeric(as.matrix(rt)), nrow=nrow(rt), dimnames=dimnames)
  #分别对TCGA取log2及GEO进行标准化
  if(substr(colnames(rt)[1], 1,3) == "TCG"){rt = log2(rt + 1)}
  if(substr(colnames(rt)[1], 1,3) == "GSM"){rt=normalizeBetweenArrays(rt)}
  geneList[[i]] = rt
  #获取基因集共同的基因名
  if (i == 1) {gene3 = rownames(geneList[[i]])}
  else{gene3 = intersect(gene3, rownames(geneList[[i]]))}
}

allmerge2=data.frame()
batchType2 = c()
for (i in 1:length(files)){
  #合并
  if (i == 1) {allmerge2 = geneList[[i]][gene3,]}
  else{allmerge2 = cbind(allmerge2, geneList[[i]][gene3,])}
  #获取不同数据集的分组信息
  batchType2 = c(batchType2, rep(i,ncol(geneList[[i]])))
}

#对数据进行批次矫正，输出矫正后的表达数据
outTab2=ComBat(allmerge2, batchType2, par.prior=TRUE)
outTab2=rbind(ID=colnames(outTab2), outTab2)
write.table(outTab2, file="merge2.txt", sep="\t", quote=F, col.names=F)

# 分离TCGA和GEO数据
tcga_cols <- grep("^TCGA", colnames(outTab2))
geo_cols <- grep("^GSM", colnames(outTab2))

# 创建TCGA数据框并删除第一行和第一列
tcga_data <- outTab2[, tcga_cols]
tcga_data <- tcga_data[-1,-1]
tcga_data <- data.frame(ID = rownames(tcga_data), tcga_data)

# 创建GEO数据框并删除第一行和第一列
geo_data <- outTab2[, geo_cols]
geo_data <- geo_data[-1,-1]
geo_data <- data.frame(ID = rownames(geo_data), geo_data)

# 保存TCGA数据
write.table(tcga_data, file="TCGA.txt", sep="\t", quote=F, row.names=F)

# 保存GEO数据
write.table(geo_data, file="GEO.txt", sep="\t", quote=F, row.names=F)

################################################################################
#合并TCGA与临床数据
#读取输入文件
data=read.table("TCGA.txt", header=T, sep="\t", check.names=F,row.names = 1)
#data=read.table("TCGA_CRC_TPM.txt", header=T, sep="\t", check.names=F,row.names = 1)
# 将列名中的 "." 替换为 "-"
colnames(data) = gsub("\\.", "-", colnames(data))
#与临床数据合并
#读入生存数据
cli=read.table("time_CRC.txt", header=T, sep="\t", check.names=F, row.names=1)
#生存时间以年为单位
cli$time=cli$time/365
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#仅保留肿瘤样本
data = data[,group == 0]

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
# 保存合并TCGA数据
write.table(rt, file="TCGA_cli.txt", sep="\t", quote=F, row.names=T)
#读取数据
A <- read.table("TCGA_cli.txt", header=T, sep="\t", check.names=F, row.names=1)
################################################################################
#合并GEO与临床数据
#pd是已经处理好的临床数据，只要time和status
# 转置表达矩阵
dat <- geo_data
dat_t <- t(dat)

# 确保行名（样本名）是字符型
rownames(dat_t) <- as.character(rownames(dat_t))

# 确保pd的行名（样本名）是字符型
rownames(pd) <- as.character(rownames(pd))

# 找出在临床数据中存在的样本
common_samples <- intersect(rownames(pd), rownames(dat_t))

# 只保留这些共同的样本
pd_filtered <- pd[common_samples, ]
dat_t_filtered <- dat_t[common_samples, ]

# 合并数据
merged_data <- cbind(pd_filtered, dat_t_filtered)
# 保存合并TCGA数据
write.table(merged_data, file="GEO_cli.txt", sep="\t", quote=F, row.names=T)
