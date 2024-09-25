library(Seurat) 
library(tidyverse)
library(scRNAtoolVis)

color.pals = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
               "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
               "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
               "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")


help(jjVolcano)


filepath <- paste("E:/B_group/sub_group/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)

group <- levels(pbmc@meta.data$group)
group <- setdiff(group, "Ad1R")  # 移除 "Ad1R" 从 group 中

combined_table <- data.frame() # 创建一个空的数据框，用于存储所有子集的数据


for (i in group) {
  pbmc.markers <- FindMarkers(pbmc, ident.1 = i, ident.2 = "Ad1R", logfc.threshold = 0.25, min.pct = 0.1,
                              group.by = 'group')
  # 将 rownames 转换为 gene 列
  pbmc.markers <- pbmc.markers %>% 
    rownames_to_column(var = "gene")
  
  pbmc.markers <- pbmc.markers %>% 
    mutate(cluster = rep(i, nrow(pbmc.markers)), group = rep(i, nrow(pbmc.markers)))# %>%
  #group_by(cluster) %>% top_n(n = 20, wt = abs(avg_log2FC) )
  
  combined_table <- bind_rows(combined_table, pbmc.markers) # 将当前子集的数据添加到 combined_table 中
}

write.xlsx(combined_table, file = paste("E:/B_group/火山图基因/", "combined.xlsx", sep = ""), row.names = F)

path = paste("E:/B_group/火山图基因/Volcano",".pdf", sep = "")

#png(path,units = "in",width = 30,height = 20,res = 600)
pdf(path, width = 18, height = 20)
jjVolcano(diffData = combined_table,
          log2FC.cutoff = 0.1, 
          size  = 3.5, #设置点的大小
          fontface = 'italic', #设置字体形式
          aesCol = c('#87CEFA','#EEA2AD'), #设置点的颜色
          tile.col = color.pals, #设置cluster的颜色
          #col.type = "adjustP", #设置矫正方式
          topGeneN = 50 #设置展示topN的基因
)
dev.off()

#环状火山图
# make a polar plot
pdf(path, width = 18, height = 20)
jjVolcano(diffData = combined_table,
          tile.col = corrplot::COL2('RdYlBu', 15)[4:12],
          size  = 3.5,
          fontface = 'italic',
          topGeneN = 40,
          polar = T) +
ylim(-8,10)
dev.off()
library(scRNAtoolVis)
featureCornerAxes(object = pbmc,reduction = 'umap',
                  groupFacet = 'orig.ident',
                  features = c("Ppar","Hp"),
                  cornerTextSize = 3,
                  themebg = 'bwCorner')
