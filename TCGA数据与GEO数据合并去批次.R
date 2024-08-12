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

