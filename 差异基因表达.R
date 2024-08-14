library(limma)
library(pheatmap)
library(stringr)
library(ggplot2)
library(ggVolcano)

#设置工作目录
setwd("E:/CORD/DEG")

#读取输入文件
data=read.table("TCGA_CORD_TPM.txt", header=T, sep="\t", check.names=F,row.names = 1)
#转化为matrix
dimnames=list(rownames(data), colnames(data))
data=matrix(as.numeric(as.matrix(data)), nrow=nrow(data), dimnames=dimnames)
#去除低表达的基因
data=data[rowMeans(data)>1,]

#正常和肿瘤数目,第14,15字符,01-09是癌症，10-19是正常，20-29是癌旁
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
#将2变为1
group=gsub("2", "1", group)
view(group)
#统计group中0和1的数量
table(group)
#获取正常及肿瘤组样本数目
conNum=length(group[group==1])
treatNum=length(group[group==0])
Type=c(rep(1,conNum), rep(2,treatNum))
#根据正常和肿瘤排序
data1 = data[,group == 1]
data2 = data[,group == 0]
data = cbind(data1,data2)

#差异分析
outTab=data.frame()
for(i in row.names(data)){
  rt=data.frame(expression=data[i,], Type=Type)
  wilcoxTest=wilcox.test(expression ~ Type, data=rt)
  pvalue=wilcoxTest$p.value
  conGeneMeans=mean(data[i,1:conNum])
  treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
  logFC=log2(treatGeneMeans)-log2(conGeneMeans)
  conMed=median(data[i,1:conNum])
  treatMed=median(data[i,(conNum+1):ncol(data)])
  diffMed=treatMed-conMed
  if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){
    outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))}
}

pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

write.table(outTab,file="TCGA.all.Wilcoxon.txt",sep="\t",row.names=F,quote=F)

logFCfilter=1
fdrFilter=0.05
#输出差异基因
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & 
                   as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="TCGA.diff.Wilcoxon.txt",sep="\t",row.names=F,quote=F)

outTab$logFC=as.numeric(as.vector(outTab$logFC))

outTab$change = ifelse(outTab$fdr < 0.05 & abs(outTab$logFC) >= 1, 
                       ifelse(outTab$logFC> 1 ,'Up','Down'),
                       'Stable')
#热图
#设置展示基因的数目
geneNum=200     
outDiff=outDiff[order(as.numeric(as.vector(outDiff$logFC))),]
diffGeneName=as.vector(outDiff[,1])
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=log2(data[hmGene,]+0.01)
Type=c(rep("Normal",conNum),rep("Tumor",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf", width=10, height=25)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 8,
         fontsize_row=5,
         fontsize_col=8)
dev.off()

#火山图
pdf(file="vol.pdf", width=5, height=5)
xMax=6
yMax=max(-log10(outTab$fdr))+1
plot(as.numeric(as.vector(outTab$logFC)), -log10(outTab$fdr), xlab="logFC",ylab="-log10(fdr)",
     main="Volcano", ylim=c(0,yMax),xlim=c(-xMax,xMax),yaxs="i",pch=20, cex=1.2)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))>logFCfilter)
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="red",cex=1.5)
diffSub=subset(outTab, fdr<fdrFilter & as.numeric(as.vector(logFC))<(-logFCfilter))
points(as.numeric(as.vector(diffSub$logFC)), -log10(diffSub$fdr), pch=20, col="green",cex=1.5)
abline(v=0,lty=2,lwd=3)
dev.off()


geneList0 <- c('TUBA1C',"ULBP2")
geneList <- outTab[geneList0,]

library('ggplot2')

rownames(outTab) <- outTab$gene

sum(is.na(outTab$logFC))  # 检查 logFC 列中有多少个 NA 值
sum(is.na(outTab$fdr))    # 检查 fdr 列中有多少个 NA 值
sum(is.na(outTab$change)) # 检查 change 列中有多少个 NA 值

outTab$logFC[is.na(outTab$logFC)] <- 0  # 将 logFC 列中的 NA 值替换为 0
outTab$fdr[is.na(outTab$fdr)] <- 1      # 将 fdr 列中的 NA 值替换为 1


outTab <- na.omit(outTab)
p <- ggplot(# 数据、映射、颜色
  outTab, aes(x = logFC, y = -log10(fdr), colour=change)) +
  geom_point(alpha=0.5, size=3.5) +
  scale_color_manual(values=c("#546de5","#d2dae2","#ff4757"))+
  #突出表示差异基因
  geom_point(data=geneList,aes(x = logFC, y = -log10(fdr)),colour="yellow",size=3.5)+
  #辅助线
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=3,col="black",lwd=0.8) +
  labs(x="log2(fold change)",y="-log10 (p-value)")+   # 坐标轴# 坐标轴和图标题title="Volcano plot",
  theme_bw()+    #去除背景色
  theme(panel.grid = element_blank())+  #去除网格线
  #xlim(-2, 2)+   #设置坐标轴范围
  #图例
  theme(plot.title = element_text(hjust = 0.5,size=24), 
        legend.position="bottom", 
        legend.title = element_blank(),
        legend.text=element_text(size=18),
        legend.key.size = unit(1, 'cm'),
        legend.background = element_rect(fill="gray90", linetype="solid",colour ="gray"),
        axis.title.x =element_text(size=18), 
        axis.title.y=element_text(size=18),
        axis.text=element_text(size=14,face = "bold"))
p

#火山
pdf(file="vol.pdf", width=5, height=5)
print(p)
dev.off()

#标记出5个基因的label
geneList1 <- outTab[rownames(outTab) %in% geneList0,]
geneList1 <- subset(geneList1, select = -change)
geneList1$label <- rownames(geneList1)

pdf("hotplot.pdf", width = 8, height = 10)
library(ggrepel)
p + geom_label_repel(data = geneList1, 
                     aes(x = logFC, y = -log10(fdr), label = label),
                     size = 4,color="black",
                     box.padding = unit(0.4, "lines"), 
                     segment.color = "black",   #连线的颜色
                     segment.size = 0.4,  #连线的粗细
)
dev.off()

