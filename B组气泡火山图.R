library(ClusterGVis)
library(org.Rn.eg.db)
library(xlsx)
library(dplyr)
library(ggplot2)
library(png)
library(Seurat)
library(xlsx)
library(ggpubr)

color.pals = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
               "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
               "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
               "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

filepath = "E:/B组数据备份(4.29)/单细胞结果/横向比较/组间/"
genepath = 'E:/B组数据备份(4.29)/横向比较基因/'

#seurat对象设置
filepath <- paste("E:/B组数据备份(4.29)/去污染后细胞亚群备份/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)
sce <- pbmc
#抽组
group = c("Ad1R","Bd1R","Bd8R","Bd15R","Bd30R")



for (i in group) {
  pbmc <- sce
  pbmc = pbmc[, pbmc@meta.data$group %in% i]#抽组
  pbmc@meta.data$group <- droplevels(pbmc@meta.data$group)
  
  mouse_gene <- read.xlsx(paste0(genepath,"5265.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
  mouse_gene <- c(mouse_gene$gene)
  
  #画图
  plot <- VlnPlot(pbmc, mouse_gene, stack = TRUE,group.by = "group",
                  sort = FALSE, flip = TRUE, cols = color.pals) +
    theme(legend.position = "none")
  
  png(paste0(filepath,"ViolinPlot_","all",".png",sep = ""), width = 8, height = 16, units = "in", res = 800)
  print(plot)
  dev.off()
}


for (i in group) {
  pbmc <- sce
  pbmc = pbmc[, pbmc@meta.data$group %in% i]#抽组
  pbmc@meta.data$group <- droplevels(pbmc@meta.data$group)
  gene <- read.xlsx(paste0(genepath,"5265.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
  gene <- c(gene$gene)
  plotfile <- paste("E:/B组数据备份(4.29)/单细胞结果/气泡图/seurat_clusters_fabp1_",i,".pdf", sep = "")
  plot <- DotPlot(pbmc, features = gene,group.by = 'seurat_clusters',dot.scale = 8)+
    theme_bw()+
    theme(panel.grid = element_blank(), 
          axis.text.x=element_text(hjust = 1, vjust=0.5, angle=90))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33')) +
    scale_y_discrete(limits = rev(levels(pbmc@meta.data$seurat_clusters)))
  pdf(plotfile, width = 8, height = 5)
  print(plot)
  dev.off()
}
#分组
close()
levels(pbmc@meta.data$seurat_clusters)
cluster <- c("Hepatocytes_3")
pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% cluster]#抽组
pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)

gene <- read.xlsx(paste0(genepath,"5265.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
gene <- c(gene$gene)
plotfile <- paste("E:/B组数据备份(4.29)/单细胞结果/气泡图/Hepatocytes_3_group",".pdf", sep = "")
plot <- DotPlot(pbmc, features = gene,group.by = 'group',dot.scale = 16)+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust = 1, vjust=0.5, angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33')) +
  scale_y_discrete(limits = rev(levels(pbmc@meta.data$group)))

pdf(plotfile, width = 8, height = 5)
print(plot)
dev.off()
#分簇
filepath = "E:/B_group/单细胞结果/横向比较/组间/"
genepath = 'E:/B_group/gene/'
gene <- read.xlsx(paste0(genepath,"gene.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
gene <- c(gene$gene)
plotfile <- paste("E:/B_group/气泡图/气泡_seurat_clusters",".pdf", sep = "")
plot <- DotPlot(pbmc, features = gene,group.by = 'seurat_clusters',dot.scale = 12)+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust = 1, vjust=0.5, angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33')) +
  scale_y_discrete(limits = rev(levels(pbmc@meta.data$seurat_clusters)))

pdf(plotfile, width = 18, height = 6)
print(plot)
dev.off()

###############################################################################
#肝细胞空白上调基因
library(Seurat) 
library(tidyverse)
library(scRNAtoolVis)

group = c("Ad1R","Bd1R","Bd8R","Bd15R","Bd30R")
#group减少"Ad1R"
group = c("Bd1R","Bd8R","Bd15R","Bd30R")
for (j in group) {
  combined_table <- data.frame() # 创建一个空的数据框，用于存储所有子集的数据
  cluster <- levels(pbmc@meta.data$seurat_clusters)
  for (i in cluster) {
    pbmc.markers <- FindMarkers(pbmc, ident.1 = j, ident.2 = "Ad1R", logfc.threshold = 0.25, min.pct = 0.1,
                                group.by = 'group', subset.ident = i)
    # 将 rownames 转换为 gene 列
    pbmc.markers <- pbmc.markers %>% 
      rownames_to_column(var = "gene")
    
    pbmc.markers <- pbmc.markers %>% 
      mutate(cluster = rep(i, nrow(pbmc.markers)), group = rep(j, nrow(pbmc.markers)))# %>%
    #group_by(cluster) %>% top_n(n = 20, wt = abs(avg_log2FC) )
    
    combined_table <- bind_rows(combined_table, pbmc.markers) # 将当前子集的数据添加到 combined_table 中
  }
  
  write.xlsx(combined_table, file = paste("E:/B组数据备份(4.29)/火山图基因/", j, "_combined.xlsx", sep = ""), row.names = F)
}

for (j in group) {
  file = paste("E:/B组数据备份(4.29)/火山图基因/", j, "_combined.xlsx", sep = "")
  combined_markers = openxlsx::read.xlsx(file, sheet = 1)
  
  path = paste("E:/B组数据备份(4.29)/单细胞结果/火山图/Volcano_",j,".png", sep = "")
  
  color.pals = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                 "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
                 "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
                 "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
  
  vol <- jjVolcano(diffData = combined_markers,
                   log2FC.cutoff = 0.1, 
                   size  = 3.5, #设置点的大小
                   fontface = 'italic', #设置字体形式
                   aesCol = c('#87CEFA','#EEA2AD'), #设置点的颜色
                   tile.col = color.pals, #设置cluster的颜色
                   #col.type = "adjustP", #设置矫正方式
                   topGeneN = 20 #设置展示topN的基因
  )
  
  png(path,units = "in",width = 12,height = 12,res = 600)
  print(vol)
  dev.off()
}
#seurat对象设置
filepath <- paste("E:/B组数据备份(4.29)/去污染后细胞亚群备份/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)

# find markers
pbmc.markers <- FindAllMarkers(pbmc, only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = 0)
path = paste("E:/B组数据备份(4.29)/单细胞结果/火山图/Volcano",".pdf", sep = "")

pdf(path, width = 14, height = 12)
jjVolcano(diffData = pbmc.markers,topGeneN = 30)
dev.off()

###############################################################################
##横向小提琴图
#抽组
filepath = "E:/B组数据备份(4.29)/单细胞结果/横向比较/组间/"
group = c("Ad1R","Bd1R","Bd8R","Bd15R","Bd30R")
#sce <- pbmc
for (i in group) {
  pbmc <- sce
  pbmc = pbmc[, pbmc@meta.data$group %in% i]#抽组
  pbmc@meta.data$group <- droplevels(pbmc@meta.data$group)
  
  mouse_gene <- read.xlsx(paste0(genepath,"5265.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
  mouse_gene <- c(mouse_gene$gene)
  
  #画图
  plot <- VlnPlot(pbmc, mouse_gene, stack = TRUE,group.by = "seurat_clusters",
                  sort = FALSE, flip = TRUE, cols = color.pals) +
    theme(legend.position = "none")
  
  pdf(paste0(filepath,"ViolinPlot_",i,".pdf",sep = ""), width = 15, height = 40)
  print(plot)
  dev.off()
}
################################################################################
#featureplot
#每一个基因上跑

mouse_gene <- read.xlsx(paste0(genepath,"5265.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
mouse_gene <- c(mouse_gene$gene)
for (i in mouse_gene) {
  filepath = "E:/B组数据备份(4.29)/单细胞结果/横向比较/组间/"
  gene <- i
  gene_name = paste("肝细胞组间","_",i, sep = "")
  pbmc_hep <- PercentageFeatureSet(sce,features = gene,col.name = gene_name)
  path2 = paste(filepath,"UAMP_",gene_name,".pdf", sep = "")
  p2 <- FeaturePlot(pbmc_hep,gene_name)
  pdf(path2, width = 8, height = 8)
  print(p2)
  dev.off()
  
}
################################################################################
gene_name = paste("肝细胞group","_",'Ppara', sep = "")
#选择seurat对象中有的基因
filepath = "E:/B组数据备份(4.29)/单细胞结果/横向比较/组间/"
genepath = 'E:/B组数据备份(4.29)/横向比较基因/'
mouse_gene <- read.xlsx(paste0(genepath,"5265.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
mouse_gene <- c(mouse_gene$gene)
mouse_gene <- mouse_gene[mouse_gene %in% rownames(pbmc)]
pbmc_hep <- PercentageFeatureSet(pbmc,features = mouse_gene,col.name = gene_name)

factors <- c(levels(pbmc_hep@meta.data$group))
combinations <- combn(factors, 2) 
combinations_list <- split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
library(ggpubr)
my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de')
filepath = "E:/B组数据备份(4.29)/单细胞结果/横向比较/组间/"
path = paste(filepath,"组间横向比较_",gene_name,".pdf", sep = "")
p.percentage <- ggviolin(pbmc_hep@meta.data, x = "group", y = gene_name,
                         color = "group",add = 'mean_sd',fill = 'group',
                         add.params = list(color = "black")) + 
  stat_compare_means(comparisons = combinations_list,label = "p.signif") + 
  scale_color_manual(values = my9color) + 
  scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  #ylim(-0.2, 0.1) +
  NoLegend() + labs(x = '')
pdf(path, width = 8, height = 8)
print(p.percentage)
dev.off()

sce <- pbmc

group = c("Ad1R","Bd1R","Bd8R","Bd15R","Bd30R")
for (i in group) {
  pbmc <- sce
  pbmc = pbmc[, pbmc@meta.data$group %in% i]#抽组
  pbmc@meta.data$group <- droplevels(pbmc@meta.data$group)
  gene_name = paste("肝细胞","_",i,"_",'Ppara', sep = "")
  genepath = 'E:/B组数据备份(4.29)/横向比较基因/'
  mouse_gene <- read.xlsx(paste0(genepath,"5265.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
  mouse_gene <- c(mouse_gene$gene)
  mouse_gene <- mouse_gene[mouse_gene %in% rownames(pbmc)]
  pbmc_hep <- PercentageFeatureSet(pbmc,features = mouse_gene,col.name = gene_name)
  
  factors <- c(levels(pbmc_hep@meta.data$seurat_clusters))
  combinations <- combn(factors, 2) 
  combinations_list <- split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  library(ggpubr)
  my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de')
  filepath = "E:/B组数据备份(4.29)/单细胞结果/横向比较/组间/"
  path = paste(filepath,"组间横向比较_",gene_name,".pdf", sep = "")
  p.percentage <- ggviolin(pbmc_hep@meta.data, x = "seurat_clusters", y = gene_name,
                           color = "seurat_clusters",add = 'mean_sd',fill = 'seurat_clusters',
                           add.params = list(color = "black")) + 
    stat_compare_means(comparisons = combinations_list,label = "p.signif") + 
    scale_color_manual(values = my9color) + 
    scale_fill_manual(values = my9color) +
    theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
    #ylim(-0.2, 0.1) +
    NoLegend() + labs(x = '')
  pdf(path, width = 8, height = 8)
  print(p.percentage)
  dev.off()
}

###
#细胞比例
library(cowplot)

filepath <- paste("E:/B组数据备份(4.29)/鉴定后总样本备份/",'all_sample_decontX035_pc30',".RDS", sep = "")
pbmc <- readRDS(filepath)

cluster <- levels(pbmc@meta.data$seurat_clusters)
dif_cluster <- c("Hepatocytes","Endothelial cells","Hepatic stellate cells","Biliary epithelial cells","Erythrocyte cells")
cluster_filtered <- cluster[!cluster %in% dif_cluster]

pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% cluster_filtered]#抽组
pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)

pbmc$celltype <- Idents(pbmc)
group_b = c("Ad1R","Bd1R","Bd8R","Bd15R","Bd30R")
group_c = c("Bd30R","Bd15R","Bd8R","Bd1R","Ad1R")
#随机生成14个色彩鲜明的颜色
library(RColorBrewer)
colors = brewer.pal(14, "Set3")#集合只有12个颜色
colors <- c(colors, "black")
colors <- c(colors, "white")

cluster_id <- levels(pbmc)
theme_set(theme_cowplot())
result <- setNames( colors,cluster_id)
print(result)

pbmc@meta.data$celltype <- ordered(pbmc@meta.data$celltype, 
                                   levels = cluster_id)
cell <- subset(pbmc, subset = celltype %in% cluster_id)
cell <- ScaleData(cell)
cell_counts <- FetchData(cell, vars = c("celltype","group","orig.ident")) %>%
  mutate(group = factor(group, levels = group_c))
cell_counts_tbl <- cell_counts %>%
  dplyr::count(celltype,group)
write_csv(cell_counts_tbl, path = "E:/B组数据备份(4.29)/单细胞结果/细胞比例/cell_counts_910.csv")
#比例图
#png("E:/B组数据备份(4.29)/单细胞结果/细胞比例/B组细胞比例(9.10).tiff",units = "in",width = 20,height = 8,res = 600)
pdf("E:/B组数据备份(4.29)/单细胞结果/细胞比例/B组细胞比例(9.10).pdf",width = 20,height = 8)
ggplot(data = cell_counts, aes(x = group, fill = celltype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = result) +
  coord_flip() +
  scale_y_reverse()+theme(legend.text = element_text(size = 20),text=element_text(size=20),
                          axis.text = element_text(size=15))
dev.off()