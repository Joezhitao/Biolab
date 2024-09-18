library(Seurat)
library(stringr)
library(harmony)
library(clustree)
library(dplyr)
library(xlsx)
library(tidyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(magrittr)
library(reshape2)
library(readxl)
library(cowplot)
library(scales)
library(readr)
library(progeny)
library(gplots)
library(tibble)
library(grid)
library(rlang)
library(DoubletFinder)
library(patchwork)
library(dplyr)
rm(list = ls())

#11组数据导入
setwd("E:/R语言单细胞数据备份")
pbmc <- readRDS('pbmc_C_group.RDS')
saveRDS(pbmc, "E:/R语言单细胞数据备份/pbmc_C_group.RDS")
group <- levels(pbmc@meta.data$group)
#分簇(默认分49簇的分辨率)
pc.num = 1:40
pbmc <- FindNeighbors(pbmc,reduction = "harmony", dims = pc.num)
pbmc <- FindClusters(pbmc, resolution = 1.4)
#抽组("Ad1R","Bd1R","Bd8R","Bd15R","Bd30R")
group = c("Ad1R","Cd15L","Cd15R")
#group = c("GSM4041174","GSM4041175")
pbmc = pbmc[, pbmc@meta.data$group %in% group]#抽组
pbmc@meta.data$group <- droplevels(pbmc@meta.data$group)
#鉴定
ident <- read.xlsx("E:/R语言单细胞数据备份/ident.xlsx",
                   sheetIndex = 1,header = T,encoding = "UTF-8")
new.cluster.ids <- c(ident$ident)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
pbmc@meta.data$cluster <- pbmc@meta.data$seurat_clusters
levels(pbmc@meta.data$seurat_clusters) <- new.cluster.ids
levels(pbmc@meta.data$seurat_clusters)
#pbmc$celltype <- Idents(pbmc)
pbmc_all <- SplitObject(pbmc, split.by = "seurat_clusters")
cluster_id <- c("Macrophages","Hepatic stellate cells")
for (i in cluster_id) {
  assign(i, pbmc_all[[i]])
}
#kegg上下调
#group去掉"Ad1R"
group <- c("Cd30L","Cd30R","Cd1L","Cd1R","Cd15L","Cd15R")
for (i in cluster_id) {
  sample <- get(i)
  for (j in group) {
    cell_pbmc <- FindMarkers(sample,ident.1 = j,ident.2 = "Ad1R",group.by = 'group',logfc.threshold = 0.25,min.pct = 0.1)
    #上调
    cell_pbmc_up <- subset(cell_pbmc,p_val_adj<0.05&(avg_log2FC > 0.1))
    genelist <- bitr(row.names(cell_pbmc_up), fromType="SYMBOL",
                     toType="ENTREZID", OrgDb='org.Rn.eg.db')
    genelist <- pull(genelist,ENTREZID)
    ekegg <- enrichKEGG(gene = genelist, organism = 'rno')
    data <- setReadable(ekegg, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
    head(data,5)
    write.table(data,paste("E:/B组数据备份(4.29)/单细胞结果/富集化分析/富集/KEGG_",i,"_",j,"上调.csv",sep = ""),sep = "\t",row.names = F,col.names = T,quote = F)
    
    if (!is.null(ekegg) && nrow(ekegg) > 0) {
      p1 <- barplot(ekegg, showCategory = 20)
      p2 <- dotplot(ekegg, showCategory = 20)
      
      # 创建一个列表,只包含非 NULL 的图
      plot_list <- list(p1, p2)
      plot_list <- plot_list[!sapply(plot_list, is.null)]
      
      # 如果列表为空,则不保存图片
      if (length(plot_list) > 0) {
        # 使用 patchwork 包中的 `/` 运算符垂直排列图片
        plotc <- Reduce("/", plot_list)
        
        path4 <- paste("E:/B组数据备份(4.29)/单细胞结果/富集化分析/富集/上调enrichKEGG_",i,"_",j, ".png", sep = "")
        ggsave(path4, plot = plotc, width = 12, height = 12)
      }
    }
    
  }
  
}
#细胞比例
#细胞比例
pbmc$celltype <- Idents(pbmc)
levels(pbmc$celltype)
group_b = c("Ad1R","Cd1L","Cd1R","Cd15L","Cd15R","Cd30L","Cd30R")
group_c = c("Ad1R","Cd30L","Cd30R","Cd1L","Cd1R","Cd15L","Cd15R")
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
write_csv(cell_counts_tbl, path = "E:/B组数据备份(4.29)/单细胞结果/细胞比例/cell_counts_C.csv")
#比例图
png("E:/B组数据备份(4.29)/单细胞结果/细胞比例/B组细胞比例(C组).tiff",units = "in",width = 20,height = 8,res = 600)
ggplot(data = cell_counts, aes(x = group, fill = celltype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = result) +
  coord_flip() +
  scale_y_reverse()+theme(legend.text = element_text(size = 20),text=element_text(size=20),
                          axis.text = element_text(size=15))
dev.off()