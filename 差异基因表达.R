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
geneNum=50
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
pdf(file="heatmap_CRC_Top50.pdf", width=10, height=8)
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


geneList0 <- c('TUBA1C')
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

#####limma####

#加载包
library(limma)
library(pheatmap)
library(stringr)
library(edgeR)
#设置工作目录
setwd("E:/CORD/DEG")

#读取输入文件
data=read.table("TCGA_CORD_count.txt", header=T, sep="\t", check.names=F,row.names = 1)
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
#获取正常及肿瘤组样本数目
conNum=length(group[group==1])
treatNum=length(group[group==0])
Type=c(rep(1,conNum), rep(2,treatNum))
#根据正常和肿瘤排序
data1 = data[,group == 1]
data2 = data[,group == 0]
data = cbind(data1,data2)
#设置分组信息
Type = factor(Type)
design <- model.matrix(~0+Type)
rownames(design) = colnames(data)
colnames(design) <- levels(Type)

#创建对象
DGElist <- DGEList( counts = data, group = Type )
#过滤掉cpm小于等于1的基因
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 # 自定义
table(keep_gene)
#生成DGEList对象
DGElist <- DGElist[keep_gene, , keep.lib.sizes = FALSE ]
#计算比例因子以将原始库大小转换为有效库大小。
DGElist <- calcNormFactors( DGElist )
#转换logCPM
v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
#非线性最小二乘法
fit <- lmFit(v, design)
colnames(design)=c("normal","tumor")
#构建比较矩阵
cont.matrix <- makeContrasts(contrasts = c('tumor-normal'), levels = design)
#线性拟合模型构建
fit2 <- contrasts.fit(fit, cont.matrix)
#用经验贝叶斯调整t-test中方差的部分
fit2 <- eBayes(fit2)
nrDEG_limma_voom = topTable(fit2, coef = 'tumor-normal', n = Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)

#筛选标准
padj = 0.05
logFC = 1
#筛选
outDiff = nrDEG_limma_voom[(nrDEG_limma_voom$adj.P.Val < padj &
                              (nrDEG_limma_voom$logFC>logFC | nrDEG_limma_voom$logFC<(-logFC))),]
outDiff = outDiff[order(outDiff$logFC),]
write.table(data.frame(ID=rownames(outDiff),outDiff),file="TCGA.diff.limma.txt",sep="\t",row.names=F,quote=F)

#####edgeR#####

#加载包
library(limma)
library(pheatmap)
library(stringr)
library(edgeR)
#设置工作目录
setwd("E:/CORD/DEG")

#读取输入文件
data=read.table("TCGA_CORD_count.txt", header=T, sep="\t", check.names=F,row.names = 1)
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
#获取正常及肿瘤组样本数目
conNum=length(group[group==1])
treatNum=length(group[group==0])
Type=c(rep(1,conNum), rep(2,treatNum))
#根据正常和肿瘤排序
data1 = data[,group == 1]
data2 = data[,group == 0]
data = cbind(data1,data2)
#分组矩阵
Type = factor(Type)
design <- model.matrix(~0+Type)
rownames(design) = colnames(data)
colnames(design) <- levels(Type)

#差异表达矩阵
DGElist <- DGEList( counts = data, group = Type)
keep_gene <- rowSums( cpm(DGElist) > 1 ) >= 2 ## 自定义
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]
DGElist <- calcNormFactors( DGElist )
DGElist <- estimateGLMCommonDisp(DGElist, design)
DGElist <- estimateGLMTrendedDisp(DGElist, design)
DGElist <- estimateGLMTagwiseDisp(DGElist, design)
fit <- glmFit(DGElist, design)
results <- glmLRT(fit, contrast = c(-1, 1))
nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
head(nrDEG_edgeR)

#筛选
padj = 0.05
logFC= 1
outdiff  = nrDEG_edgeR[(nrDEG_edgeR$FDR < padj &
                          (nrDEG_edgeR$logFC>logFC | nrDEG_edgeR$logFC<(-logFC))),]
write.table(data.frame(ID=rownames(outdiff),outdiff),file="TCGA.diff.edgeR.txt",sep="\t",row.names=F,quote=F)

#####DESeq2#####

#加载包
library(limma)
library(pheatmap)
library(stringr)
library(DESeq2)
#设置工作目录
setwd("E:/CORD/DEG")

#读取输入文件
data=read.table("TCGA_CORD_count.txt", header=T, sep="\t", check.names=F,row.names = 1)
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
#获取正常及肿瘤组样本数目
conNum=length(group[group==1])
treatNum=length(group[group==0])
Type=c(rep(1,conNum), rep(2,treatNum))
#根据正常和肿瘤排序
data1 = data[,group == 1]
data2 = data[,group == 0]
data = cbind(data1,data2)

#分组矩阵
condition = factor(Type)
coldata <- data.frame(row.names = colnames(data), condition)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~condition)
View(dds)
dds$condition
# 指定哪一组作为对照组
dds$condition<- relevel(dds$condition, ref = "1")
#差异表达矩阵
dds <- DESeq(dds)
allDEG2 <- as.data.frame(results(dds))

#筛选
padj = 0.05
logFC= 1
outdiff = allDEG2[(allDEG2$padj < padj &
                     (allDEG2$log2FoldChange>logFC | allDEG2$log2FoldChange<(-logFC))),]
write.table(data.frame(ID=rownames(outdiff),outdiff),file="TCGA.diff.DESeq2.txt",sep="\t",row.names=F,quote=F)

#####VennDiagram####

#加载
library(VennDiagram)
#设置工作目录
setwd("E:/CORD/DEG")
#读取
data=read.table("TCGA.diff.Wilcoxon.txt", header=T, sep="\t", check.names=F,row.names = 1)
Wilcoxon = rownames(data)
data=read.table("TCGA.diff.limma.txt", header=T, sep="\t", check.names=F,row.names = 1)
limma = rownames(data)
data=read.table("TCGA.diff.edgeR.txt", header=T, sep="\t", check.names=F,row.names = 1)
edgeR = rownames(data)
data=read.table("TCGA.diff.DESeq2.txt", header=T, sep="\t", check.names=F,row.names = 1)
DESeq2 = rownames(data)

#画图
venn.diagram(x = list('edgeR' = edgeR,'limma' = limma,'DESeq2' = DESeq2,'Wilcoxon' = Wilcoxon),
             filename = 'VN.png',fill = c("dodgerblue", "goldenrod1", "darkorange1","green"))