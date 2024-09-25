#seurat对象设置
filepath <- paste("E:/B_group/sub_group/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)
levels(pbmc@meta.data$seurat_clusters)

cluster <- c("Hepatocytes_3")
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
path = paste(filepath,"组间横向比较_",gene_name,".png", sep = "")
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
plotfile <- paste("E:/B_group/气泡图/UMAP_再生基因",".png", sep = "")
png(plotfile, width = 6, height = 6, units = "in", res = 800)
print(p.percentage)
dev.off()

gene
