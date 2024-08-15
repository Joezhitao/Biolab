#加载包
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

#工作目录
setwd("E:/CORD/GSEA")

#读取输入文件,并对输入文件进行整理
rt=read.table("TCGA.all.limma.txt", header=T, sep="\t", check.names=F)
# 将 gene 列移到第一列
rt = rt[, c("gene", setdiff(names(rt), "gene"))]
rt=rt[order(rt[,"logFC"],decreasing=T),]
logFC=as.vector(rt[,"logFC"])
names(logFC)=as.vector(rt[,1])
logFC[1:10]


#读入基因集文件
gmt=read.gmt("c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt")
gmt[,1] = gsub("KEGG_", "", gmt[,1])

#切换GO分析
#gmt=read.gmt("c5.go.v2024.1.Hs.symbols.gmt")
#gmt[,1] = gsub("GO", "", gmt[,1])

#对排序好的基因进行GSEA富集分析
kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.KEGG.txt",sep="\t",quote=F,row.names = F)

#输出实验组富集的图形
termNum=10      #展示通路的数目
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
  showTerm=row.names(kkUp)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Tumor", pvalue_table = T)
  pdf(file="GSEA.tumor.KEGG.pdf", width=13, height=10)
  print(gseaplot)
  dev.off()
}

#输出正常组富集的图形
termNum=10      #展示通路的数目
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
  showTerm=row.names(kkDown)[1:termNum]
  gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Normal", pvalue_table = T)
  pdf(file="GSEA.normal.KEGG.pdf", width=13, height=10)
  print(gseaplot)
  dev.off()
}

