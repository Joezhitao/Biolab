
#install.packages("VennDiagram")

#加载
library(VennDiagram)

#设置工作目录
setwd("E:/CORD/VEEN")

#读取
data=read.table("TCGA.diff.Wilcoxon.txt", header=F, sep="\t", check.names=F)
Wilcoxon = data[,1]
data=read.table("TCGA.diff.limma.txt", header=F, sep="\t", check.names=F)
limma = data[,1]
data=read.table("TCGA.diff.edgeR.txt", header=F, sep="\t", check.names=F)
edgeR = data[,1]
data=read.table("TCGA.diff.DESeq2.txt", header=F, sep="\t", check.names=F)
DESeq2 = data[,1]

#画图
venn.diagram(x = list('edgeR' = edgeR,'limma' = limma,'DESeq2' = DESeq2,'Wilcoxon' = Wilcoxon),
             filename = 'VN.png',fill = c("dodgerblue", "goldenrod1", "darkorange1","green"))

#交集基因
intersectGenes1 = intersect(Wilcoxon,limma)
intersectGenes2 = intersect(edgeR,DESeq2)
intersectGenes = intersect(intersectGenes1,intersectGenes2)

write.table(file="intersectGenes.txt", intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)
