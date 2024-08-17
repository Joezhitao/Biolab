#加载包
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringi)	
library(GOplot)
R.utils::setOption("clusterProfiler.download.method",'auto')

#工作路径
setwd("E:/CORD/GOKEGG")

#读入
input_diff =read.table("uniCox.txt",sep="\t",header=T,check.names=F) 
input_gene <- input_diff[,1]
#去除重复基因
input_gene=unique(as.vector(input_gene))

#将gene symbol转换为基因id
#org.Hs.eg为人的物种
#https://www.jianshu.com/p/84e70566a6c6
entrezIDs=mget(input_gene, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
#去除基因id为NA的基因
gene=entrezIDs[entrezIDs!="NA"]
#去除多个ID
gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#筛选条件
pvalueFilter=0.05
qvalueFilter=1       

if(qvalueFilter>0.05){colorSel="pvalue"}else{colorSel="qvalue"}

#GO
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

#保存
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#显示的数目
showNum=20

#柱状图
#pdf(file="GObarplot.pdf", width=17, height=17)
bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, 
            split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
png(file="GObarplot.png", width=16, height=10, res=600, units="in")
print(bar)
dev.off()

#气泡图
#pdf(file="GObubble.pdf", width=17, height=17)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", 
            label_format=130, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
png(file="GObubble.png", width=16, height=10, res=600, units="in")
print(bub)
dev.off()

barplot(kk)  #富集柱形图
dotplot(kk)  #富集气泡图

# 计算术语相似性矩阵
kk <- pairwise_termsim(kk)

# 创建网络图
emapplot(kk)

cnetplot <- cnetplot(kk) #网络图展示富集功能和基因的包含关系
png(file="cnetplot.png", width=16, height=16, res=800, units="in")
print(cnetplot)
dev.off()

emapplot <- emapplot(kk)#网络图展示各富集功能之间共有基因关系
png(file="emapplot.png", width=16, height=16, res=800, units="in")
print(emapplot)
dev.off()

heatplot <- heatplot(kk) #热图展示富集功能和基因的包含关系
png(file="heatplot.png", width=16, height=16, res=800, units="in")
print(heatplot)
dev.off()

#筛选条件
pvalueFilter=0.05
qvalueFilter=1       

if(qvalueFilter>0.05){colorSel="pvalue"}else{colorSel="qvalue"}

#KEGG
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(input_gene[match(strsplit(x,"/")[[1]],
                                                                              as.character(entrezIDs))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

#保存
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=20

#柱状图
#pdf(file="KEGGbarplot.pdf", width=17, height=17)
png(file="KEGGbarplot.png", width=17, height=17, res=600, units="in")
barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
dev.off()

#气泡图
#pdf(file="KEGGbubble.pdf", width = 17, height = 17)
png(file="KEGGbubble.png", width=17, height=17, res=600, units="in")
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
dev.off()
