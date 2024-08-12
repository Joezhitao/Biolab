#加载R包
library(limma)
library(sva)

#设置工作目录
setwd("E:/CORD/merge")

# 设置文件夹路径
folder_path <- "E:/CORD/merge"

# 获取所有TXT文件的名字
files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = FALSE)

# 输出结果
print(files)


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
write.table(outTab2, file="merge.txt", sep="\t", quote=F, col.names=F)

##############################################################################
#少量样本合并

#读取，转化为matrix
rt1=read.table("TCGA_LUSC_TPM.txt", header=T, sep="\t",check.names=F,row.names = 1)
dimnames=list(rownames(rt1), colnames(rt1))
rt1=matrix(as.numeric(as.matrix(rt1)), nrow=nrow(rt1), dimnames=dimnames)

rt2=read.table("GSE30219.txt", header=T, sep="\t",check.names=F,row.names = 1)
dimnames=list(rownames(rt2), colnames(rt2))
rt2=matrix(as.numeric(as.matrix(rt2)), nrow=nrow(rt2), dimnames=dimnames)

rt3=read.table("GSE74777.txt", header=T, sep="\t",check.names=F,row.names = 1)
dimnames=list(rownames(rt3), colnames(rt3))
rt3=matrix(as.numeric(as.matrix(rt3)), nrow=nrow(rt3), dimnames=dimnames)

#对TCGA数据取log2
rt1 = log2(rt1 + 1)

#对GEO数据进行标准化
rt2=normalizeBetweenArrays(rt2)
rt3=normalizeBetweenArrays(rt3)

#获取共同基因
gene1 = intersect(rownames(rt1),rownames(rt2))
gene2 = intersect(gene1,rownames(rt3))
#length(rownames(rt1))

#合并
allmerge1=data.frame()
allmerge1=cbind(rt1[gene2,], rt2[gene2,])
allmerge1=cbind(allmerge1, rt3[gene2,])
#dim(allmerge1)

#获取不同数据集的分组信息
batchType1 = c()
batchType1 = c(batchType1, rep(1,ncol(rt1)))
batchType1 = c(batchType1, rep(2,ncol(rt2)))
batchType1 = c(batchType1, rep(3,ncol(rt3)))

#对数据进行批次矫正，输出矫正后的表达数据
outTab1=ComBat(allmerge1, batchType1, par.prior=TRUE)
outTab1=rbind(ID=colnames(outTab1), outTab1)
write.table(outTab1, file="merge1.txt", sep="\t", quote=F, col.names=F)