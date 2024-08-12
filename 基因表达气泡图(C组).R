library(Seurat)
library(ggplot2)
library(cowplot)
library(xlsx)
library(dplyr)
rm(list=ls())

filepath = "E:/B组数据备份(4.29)/单细胞结果/横向比较/组间/"
genepath = 'E:/B组数据备份(4.29)/横向比较基因/'


filepath <- "E:/B组数据备份(4.29)/鉴定后总样本备份/c_group.RDS"
pbmc <- readRDS(filepath)

levels(pbmc@meta.data$group)
#样本重新排序
cluster_order <- c("Ad1R","Cd1L","Cd1R","Cd15L","Cd15R","Cd30L","Cd30R")
pbmc@meta.data$group <- factor(pbmc@meta.data$group, levels = cluster_order)
#提取肝细胞
cluster <- c("Hepatocytes")
pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% cluster]#抽组
pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)

#读取气泡图基因
gene <- read.xlsx(paste0(genepath,"5265.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
gene <- gene$gene
#选取pbmc中表达的基因进行气泡图绘制
gene <- gene[gene %in% rownames(pbmc)]

sce <- pbmc
png("E:/B组数据备份(4.29)/单细胞结果/气泡图/气泡C组左右肝_衰老基因.png", width = 20, 
    height = 4, units = "in", res = 1000)
DotPlot(sce, features = gene,group.by = 'group')+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust = 1, vjust=0.5, angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#FFCC33','#66CC66','#336699','#330066')) +
  scale_y_discrete(limits = rev(levels(sce@meta.data$group)))
dev.off()
