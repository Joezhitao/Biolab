#加载包
library(limma)
library(ggpubr)

#设置路径
setwd("E:/CORD/Boxplot")

#读取输入文件
data=read.table("TCGA_CORD_TPM.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#去除正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))

#读取临床数据文件
cli=read.table("clinical.txt", header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))

#设置基因
gene = "FEZF1"

#合并数据
samSample=intersect(row.names(data), row.names(cli))
data=data[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(data, cli)
rt=rt[,c(gene,colnames(cli))]

#临床相关性分析，输出图形结果
for(clinical in colnames(rt[,2:ncol(rt)])){
  data=rt[c(gene, clinical)]
  colnames(data)=c(gene, "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  #设置比较组
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  #绘制箱线图
  boxplot=ggboxplot(data, x="clinical", y=gene, fill="clinical",
                    xlab=clinical,
                    ylab=paste(gene, " expression"),
                    legend.title=clinical)+ 
    stat_compare_means(comparisons = my_comparisons)
  #输出图片
  #pdf(file=paste0("clinicalCor_", clinical, ".pdf"), width=5.5, height=5)
  png(file=paste0("clinicalCor_", clinical, ".png"), width=5.5, height=5, units = "in",res=800)
  print(boxplot)
  dev.off()
}
