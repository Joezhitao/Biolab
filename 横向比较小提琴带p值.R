#seurat对象设置
filepath <- paste("E:/B_group/sub_group/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)
levels(pbmc@meta.data$seurat_clusters)

cluster <- c("Hepatocytes_1","Hepatocytes_2")
pbmc <- pbmc[,pbmc@meta.data$seurat_clusters %in% cluster]
pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)

#在基因集上跑
genepath = 'E:/B_group/gene/'
gene <- read.xlsx(paste0(genepath,"gene.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
gene <- c(gene$gene)
gene_name = paste("肝细胞组间","_",'再生基因集', sep = "")
pbmc_hep <- PercentageFeatureSet(pbmc,features = gene,col.name = gene_name)

factors <- c(levels(pbmc_hep@meta.data$group))
combinations <- combn(factors, 2) 
combinations_list <- split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
library(ggpubr)
my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de')

p.percentage <- ggviolin(pbmc_hep@meta.data, x = "group", y = gene_name,
                         color = "group",add = 'mean_sd',fill = 'group',
                         add.params = list(color = "black")) + 
  stat_compare_means(comparisons = combinations_list,label = "p.signif") + 
  scale_color_manual(values = my9color) + 
  scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) #+ 
#ylim(-0.2, 0.1) +
NoLegend() + labs(x = '')

path2 = paste(filepath,"UAMP_",gene_name,".png", sep = "")
p2 <- FeaturePlot(pbmc_hep,gene_name,pt.size = 2)
plotfile <- paste("E:/B_group/气泡图/UMAP",".png", sep = "")
png(plotfile, width = 6, height = 6, units = "in", res = 800)
print(p.percentage)
dev.off()


group <- levels(pbmc@meta.data$group)
sce <- pbmc
for (i in group) {
  pbmc <- sce
  pbmc <- pbmc[,pbmc@meta.data$group %in% i]
  pbmc@meta.data$group <- droplevels(pbmc@meta.data$group)
  #在基因集上跑
  
  
  
  genepath = 'E:/B_group/gene/'
  gene <- read.xlsx(paste0(genepath,"gene.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
  gene <- c(gene$gene)
  #选择seurat对象中有的基因
  gene <- gene[gene %in% rownames(pbmc)]
  gene_name = paste("肝细胞组间","_",'再生基因集', sep = "")
  pbmc_hep <- PercentageFeatureSet(pbmc,features = gene,col.name = gene_name)
  
  factors <- c(levels(pbmc_hep@meta.data$seurat_clusters))
  combinations <- combn(factors, 2) 
  combinations_list <- split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
  library(ggpubr)
  my9color <- c('#5470c6','#91cc75','#fac858','#ee6666','#73c0de')
  
  p.percentage <- ggviolin(pbmc_hep@meta.data, x = "seurat_clusters", y = gene_name,
                           color = "seurat_clusters",add = 'mean_sd',fill = 'seurat_clusters',
                           add.params = list(color = "black")) + 
    stat_compare_means(comparisons = combinations_list,label = "p.signif") + 
    scale_color_manual(values = my9color) + 
    scale_fill_manual(values = my9color) +
    theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) #+ 
  #ylim(-0.2, 0.1) +
  NoLegend() + labs(x = '')
  
  plotfile <- paste("E:/B_group/UMAP/UMAP_横向_",i,".png", sep = "")
  png(plotfile, width = 6, height = 6, units = "in", res = 800)
  print(p.percentage)
  dev.off()
  
  p2 <- FeaturePlot(pbmc_hep,gene_name,pt.size = 2)
  plotfile <- paste("E:/B_group/UMAP/UMAP_cdk_",i,".png", sep = "")
  png(plotfile, width = 6, height = 6, units = "in", res = 800)
  print(p2)
  dev.off()
}

