rm(list = ls())
library(SeuratData)
library(Seurat)
library(BuenColors)
library(ggstatsplot)
library(patchwork)
library(xlsx)

filepath <- paste("E:/B_group/sub_group/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)

filepath = "E:/B_group/vio/"
genepath = 'E:/B_group/gene/'

mouse_gene <- read.xlsx(paste0(genepath,"gene.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
mouse_gene <- c(mouse_gene$gene)

#画图
color.pals = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
               "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
               "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
               "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
plot <- VlnPlot(pbmc, mouse_gene, stack = TRUE,
                sort = FALSE, flip = TRUE, cols = color.pals) +
  theme(legend.position = "none")

png(paste0(filepath,"ViolinPlot.png",sep = ""), width = 8, height = 16, units = "in", res = 800)
print(plot)
dev.off()

mouse_gene <- read.xlsx(paste0(genepath,"gene.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
mouse_gene <- c(mouse_gene$gene)
for (i in mouse_gene) {
  filepath = "E:/B_group/UMAP//"
  gene <- i
  gene_name = paste("肝细胞组间","_",i, sep = "")
  pbmc_hep <- PercentageFeatureSet(pbmc,features = gene,col.name = gene_name)
  path2 = paste(filepath,"UAMP_",gene_name,".png", sep = "")
  p2 <- FeaturePlot(pbmc_hep,gene_name,pt.size = 1)
  png(path2, width = 4, height = 4, units = "in", res = 800)
  print(p2)
  dev.off()
  
}
help(FeaturePlot)
sce <- pbmc
cluster <- levels(pbmc@meta.data$seurat_clusters)
for (i in cluster) {
  pbmc <- sce
  pbmc <- pbmc[,pbmc@meta.data$seurat_clusters %in% i]
  pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)
  gene_name = paste("肝细胞组间","_","Mki67_",i, sep = "")
  pbmc_hep <- PercentageFeatureSet(pbmc,features = gene,col.name = gene_name)
  path2 = paste(filepath,"UAMP_",gene_name,".png", sep = "")
  p2 <- FeaturePlot(pbmc_hep,gene_name,pt.size = 1)
  png(path2, width = 4, height = 4, units = "in", res = 800)
  print(p2)
  dev.off()
}
