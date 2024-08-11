library(dplyr)
library(ggplot2)
library(png)
library(Seurat)
library(xlsx)
library(ggpubr)
set.seed(123456)

#文件保存路径
filepath = "E:/B组数据备份(4.29)/单细胞结果/横向比较/组间/"
sce <- pbmc
#每一个基因上跑
mouse_gene <- c("Mki67")
mouse_gene <- mouse_gene[mouse_gene %in% rownames(sce)]
for (i in mouse_gene) {
  gene <- i
  gene_name = paste("肝细胞组间","_",i, sep = "")
  pbmc_hep <- PercentageFeatureSet(sce,features = gene,col.name = gene_name)
  
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
  ggsave(path,width = 8,height = 8)
  path2 = paste(filepath,"UAMP_",gene_name,".png", sep = "")
  p2 <- FeaturePlot(pbmc_hep,gene_name)
  ggsave(path2,width = 6,height = 6)
  
}
