library(dplyr)
library(ggplot2)
library(png)
library(Seurat)
library(xlsx)
library(ggpubr)
set.seed(123456)

filepath <- paste("E:/B_group/sub_group/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)

sce <- pbmc
pbmc <- sce

levels(pbmc@meta.data$seurat_clusters)
cluster <- c("Hepatocytes_1","Hepatocytes_2")
pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% cluster]#抽组
pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)

filepath = "E:/B_group/vio/"
genepath = 'E:/B_group/gene/'

mouse_gene <- read.xlsx(paste0(genepath,"gene.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
mouse_gene <- c(mouse_gene$gene)
#选择seurat对象中有的基因
mouse_gene <- mouse_gene[mouse_gene %in% rownames(pbmc)]

gene_name = paste("Hep1_2组间","_",'衰老基因集', sep = "")
pbmc_hep <- PercentageFeatureSet(pbmc,features = mouse_gene,col.name = gene_name)

factors <- c(levels(pbmc_hep@meta.data$group))
combinations <- combn(factors, 2) 
combinations_list <- split(combinations, rep(1:ncol(combinations), each = nrow(combinations)))
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
ggsave(path,width = 8,height = 8)
