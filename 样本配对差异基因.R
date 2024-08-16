#加载
library("limma")
library("ggpubr")

#设置工作目录
setwd("E:/CORD/pairDEG")

#读取输入文件
data=read.table("TCGA_CORD_count.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#去除低表达的基因
data=data[rowMeans(data)>1,]

#提取正常和肿瘤样品表达量
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
normal=data[,group!=0]
tumor=data[,group==0]

#设置需要提取的基因名称
gene = "ABCD3"

#对正常数据进行整理
normal=rbind(normal,gene=normal[gene,])
normal=as.matrix(t(normal[c("gene",gene),]))
rownames(normal)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(normal))

#对肿瘤数据进行整理
tumor=rbind(tumor,gene=tumor[gene,])
tumor=as.matrix(t(tumor[c("gene",gene),]))
rownames(tumor)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(tumor))

#合并数据
samSample=intersect(row.names(normal),row.names(tumor))
normal=normal[samSample,]
tumor=tumor[samSample,]
data=cbind(normal,tumor)
data=as.data.frame(data[,c(1,3)])
colnames(data)=c("Normal","Tumor")

#绘制图形
pdf(file="pairDiff.pdf",width=5.5,height=5)
ggpaired(data, cond1 = "Normal", cond2 = "Tumor",fill = c("red","blue"), palette = "jco",
         xlab="",ylab = paste0(gene," expression"))+
  stat_compare_means(paired = TRUE, label = "p.format", label.x = 1.35)
dev.off()



#加载
library("limma")
library("ggpubr")

#设置工作目录
setwd("E:/CORD/pairDEG") 

#读取输入文件
data=read.table("TCGA_CORD_TPM.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#去除低表达的基因
data=data[rowMeans(data)>1,]

#提取正常和肿瘤样品表达量
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
normal=data[,group!=0]
tumor=data[,group==0]

#对所有基因进行差异表达分析
#对正常数据进行整理
normal=as.matrix(t(normal))
rownames(normal)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(normal))

#对肿瘤数据进行整理
tumor=as.matrix(t(tumor))
rownames(tumor)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*","\\1\\-\\2\\-\\3",rownames(tumor))

#提取共同样本，合并数据
samSample=intersect(row.names(normal),row.names(tumor))
normal=normal[samSample,]
tumor=tumor[samSample,]
rownames(normal) = paste(rownames(normal),"normal",sep="_")
rownames(tumor) = paste(rownames(tumor),"tumor",sep="_")
data=rbind(normal,tumor)
data=as.matrix(t(data))

#分组
Type=c(rep(1,length(samSample)), rep(2,length(samSample)))

#差异分析
outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:length(samSample)])
  treatGeneMeans=mean(data[i,(length(samSample)+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  conMed=median(data[i,1:length(samSample)])
  treatMed=median(data[i,(length(samSample)+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))}
}

pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

logFCfilter=1
fdrFilter=0.05
#输出差异基因
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & 
                   as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="TCGA.diff.Wilcoxon.paired.txt",sep="\t",row.names=F,quote=F)

