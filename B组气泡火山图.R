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
#抽组
group = c("Ad1R","Bd1R","Bd8R","Bd15R","Bd30R")

sce <- pbmc

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
  plotfile <- paste("E:/B组数据备份(4.29)/单细胞结果/气泡图/气泡再生基因_",i,".png", sep = "")
  plot <- DotPlot(pbmc, features = gene,group.by = 'seurat_clusters',dot.scale = 8)+
    theme_bw()+
    theme(panel.grid = element_blank(), 
          axis.text.x=element_text(hjust = 1, vjust=0.5, angle=90))+
    labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
    scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33')) +
    scale_y_discrete(limits = rev(levels(pbmc@meta.data$seurat_clusters)))
  png(plotfile, width = 8, height = 4, units = "in", res = 1000)
  print(plot)
  dev.off()
}

###############################################################################
#肝细胞空白上调基因
library(Seurat) 
library(tidyverse)
library(scRNAtoolVis)
devtools::install_github('junjunlab/jjPlot')
devtools::install_github("sajuukLyu/ggunchull", type = "source")
install.packages('geomtextpath')
install.packages('ggh4x')
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


file = "E:/B组数据备份(4.29)/单细胞结果/去污染/combined_all_markers.csv"
combined_markers = read.csv(file,header = T,stringsAsFactors = F)
path = paste("E:/B组数据备份(4.29)/单细胞结果/火山图/Volcano",".png", sep = "")

color.pals = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
               "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
               "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
               "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")