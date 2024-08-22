library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(stringi)	
library(GOplot)
R.utils::setOption("clusterProfiler.download.method",'auto')
# 设置工作目录
base_dir <- "E:/CRC"
setwd(base_dir)
#读入
input_diff =read.table("uniCox.txt",sep="\t",header=T,check.names=F)
input_diff2 <- read.table("TCGA.diff.limma.txt",sep="\t",header=T,check.names=F)
input_gene <- input_diff[,1]
# 使用基因列表过滤第二个数据集，仅保留ID在input_gene中的行
filtered_diff2 <- input_diff2[input_diff2$ID %in% input_gene, ]
# 查看结果
head(filtered_diff2)
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
write.table(kk,file="eG.txt", sep="\t", quote=F, row.names = F)
library(tidyverse)
eGo <- as.data.frame(kk)
eGo <- separate(data=eGo, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
eGo <- separate(data=eGo, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
eGo <- mutate(eGo, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor

eGoBP <- eGo %>% 
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoCC <- eGo %>% 
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoMF <- eGo %>% 
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGo10 <- rbind(eGoBP,eGoMF,eGoCC)

library(ggplot2)
p <- ggplot(eGo10,aes(enrichment_factor,Description)) + 
  geom_point(aes(size=Count,color=qvalue,shape=ONTOLOGY)) +
  scale_color_gradient(low="red",high = "green") + 
  labs(color="pvalue",size="Count", shape="Ontology",
       x="Enrichment Factor",y="GO term",title="GO enrichment") + 
  theme_bw()
png("GO_enrichment.png",width=12,height=8,units="in",res=600)
print(p)
dev.off()

library(ggplot2)
p <- ggplot(eGo10,aes(enrichment_factor,Description)) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
  scale_color_gradient(low="green",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count", shape="Ontology",
       x="Enrichment Factor",y="GO term",title="GO enrichment") + 
  theme_bw()
png("GO_enrichment2.png",width=12,height=8,units="in",res=600)
p
dev.off()

p <- ggplot(eGo10,aes(enrichment_factor,Description)) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
  scale_color_gradient(low="green",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count", shape="Ontology",
       x="Enrichment Factor",y="GO term",title="GO enrichment") + 
  theme_bw()
png("GO_enrichment3.png",width=12,height=8,units="in",res=600)
p
dev.off()

p <- p + facet_wrap( ~ ONTOLOGY)
png("GO_enrichment4.png",width=12,height=8,units="in",res=600)
print(p)
dev.off()

p <- p + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free')
png("GO_enrichment5.png",width=12,height=8,units="in",res=600)
print(p)
dev.off()
#以用fct_reorder(factor(x), y, .fun = median, .desc = FALSE)函数（将x按照y的顺序排序）对绘图排序
p <- ggplot(eGo10,aes(enrichment_factor, fct_reorder(factor(Description), enrichment_factor))) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
  scale_color_gradient(low="green",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count", shape="Ontology",
       x="Enrichment Factor",y="GO term",title="GO enrichment") + 
  theme_bw()
p <- p + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free')
png("GO_enrichment6.png",width=12,height=8,units="in",res=600)
print(p)
dev.off()

#GOBubble图
library(GOplot)
#准备circle_dat数据
GOterms = data.frame(category = eGo10$ONTOLOGY,
                     ID = eGo10$ID,
                     term = eGo10$Description, 
                     genes = gsub("/", ",", eGo10$geneID), 
                     adj_pval = eGo10$p.adjust)
genelist <- data.frame(ID = filtered_diff2$ID, logFC = filtered_diff2$logFC) #从已有“数据”中提取genelist，1列ID，1列logFC。
circ <- circle_dat(GOterms, genelist)
#绘制GOBubble图
png("GOBubble.png",width=14,height=12,units="in",res=600)
GOBubble(circ, labels = 5, # 标注的界值：-log(adjusted p-value) (默认5)
         table.legend =T, #是否显示右侧legend，默认是
         ID=T, # T标注term名称，F标注ID
         display='single') #是否分屏
dev.off()

#GOBubble2
png("GOBubble2.png",width=14,height=12,units="in",res=600)
GOBubble(circ, labels = 5, # 标注的界值：-log(adjusted p-value) (默认5)
         table.legend =F, #不显示右侧legend
         ID=F, # 标注term名称
         display='single') # 不分屏
dev.off()

#GOCircle图
png("GOCircle.png",width=14,height=12,units="in",res=600)
GOCircle(circ)
dev.off()

#GOCircle图参数
png("GOCircle2.png",width=14,height=12,units="in",res=600)
GOCircle(circ,rad1=2, #内环半径
         rad2=3, #外环半径
         label.size= 5, #标签大小
         label.fontface = 'bold', #标签字体
         nsub=10, #显示的terms数，前10个。（画图前需要先把数据整理好，想要哪些term）
         zsc.col = c('red', 'white', "green"), # z-score颜色
         lfc.col = c('red', 'green')) # 基因up或down的颜色
dev.off()

#GOHeat热图
chord <- chord_dat(circ, #前面的circ数据
                   genelist[1:100,], #选择需要的基因
                   GOterms$term[1:10]) #选择要显示的term
png("GOHeat.png",width=14,height=12,units="in",res=600)
GOHeat(chord, nlfc =1, #数据中有logFC填1，没有填0，默认为0
       fill.col = c('red', 'white', 'blue')) #自定义颜色
dev.off()

#GOChord弦图
png("GOChord.png",width=18,height=14,units="in",res=600)
GOChord(chord)
dev.off()

#参数
png("GOChord2.png",width=18,height=14,units="in",res=600)
GOChord(chord, space = 0, #弦的间隔
        gene.order = 'logFC', #基因排序方式
        gene.space = 0.25, #基因名称和图的间隔
        gene.size = 5, #基因名的大小
        nlfc =1, #是否有logFC
        border.size= NULL, #彩虹边框大小
        lfc.col=c('red','black','cyan')) #自定义logFC 颜色
dev.off()

#树状图
eGoBP <- enrichGO(gene = gene,
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff =0.01, 
                  qvalueCutoff = 0.01,
                  ont="BP", #BP\MF\CC
                  readable =T)
png("GOtree.png",width=14,height=12,units="in",res=600)
plotGOgraph(eGoBP)
dev.off()