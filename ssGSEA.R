#加载包
library(limma)
library(GSEABase)
library(GSVA)
library(pheatmap)
library(reshape2)
library(ggpubr)

#设置工作目录
setwd("E:/CORD/ssGSEAGSVA")

#读取输入文件
data=read.table("TCGA_CORD_count.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)

#读入基因集
geneSets=getGmt("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt", geneIdType=SymbolIdentifier())

#分析 method=c("gsva", "ssgsea", "zscore", "plage") kcdf=c("Gaussian", "Poisson", "none")
#gsvaResult=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)

#若出现以下报错，这是由于GSVA包更新到新版本1.52所导致
#Calling gsva(expr=., gset.idx.list=., method=., ...) is defunct; use a method-specific parameter object (see '?gsva').
#使用以下命令代替上述命令
#不管是从github还是bioconductor都已无法下载安装旧版本，要使用旧版本，只能从有旧版本GSVA包的R中复制

#分析
gsvapar <- gsvaParam(data, geneSets,kcdf='Gaussian',absRanking=TRUE)
gsvaResult <- gsva(gsvapar)

#对打分进行标准化
normalize=function(x){return((x-min(x))/(max(x)-min(x)))}
gsvaResult=normalize(gsvaResult)

#导出
gsvaOut=rbind(id=colnames(gsvaResult), gsvaResult)
write.table(gsvaOut, file="ssgseaOut.txt", sep="\t", quote=F, col.names=F)

#获取分组信息
#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#将2变为1
group=gsub("2", "1", group)
#获取正常及肿瘤组样本数目
conNum=length(group[group==1])
treatNum=length(group[group==0])
Type=c(rep(1,conNum), rep(2,treatNum))
#根据正常和肿瘤排序
gsvaResult1 = gsvaResult[,group == 1]
gsvaResult2 = gsvaResult[,group == 0]
gsvaResult = cbind(gsvaResult1,gsvaResult2)
gsvaResult = cbind(t(gsvaResult),Type)

#提取差异显著的通路
sigGene=c()
for(i in colnames(gsvaResult)[1:(ncol(gsvaResult)-1)]){
  test=wilcox.test(gsvaResult[,i] ~ gsvaResult[,"Type"])
  pvalue=test$p.value
  if(pvalue<0.05){
    sigGene=c(sigGene, i)
  }
}


#热图
hmExp=gsvaResult[,sigGene]
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=rownames(gsvaResult)
Type=as.data.frame(Type)
hmExp = t(hmExp)

#pdf(file="heatmap.pdf", width=10, height=12)
heatmap <- pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5,
         fontsize_col=8)
png(file="heatmap.png", width=10, height=12, units="in", res=800)
print(heatmap)
dev.off()

