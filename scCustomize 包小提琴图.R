library(scCustomize)
library(xlsx)
library(Seurat)

filepath <- paste("E:/B_group/sub_group/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)

levels(pbmc@meta.data$seurat_clusters)
cluster <- c("Hepatocytes_3")
pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% cluster]#抽组
pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)

filepath = "E:/B_group/vio/"
genepath = 'E:/B_group/gene/'

mouse_gene <- read.xlsx(paste0(genepath,"gene.xlsx", sep = ""),sheetIndex = 1,header = T,encoding = "UTF-8")
mouse_gene <- c(mouse_gene$gene)
#选择seurat对象中有的基因
mouse_gene <- mouse_gene[mouse_gene %in% rownames(pbmc)]


pdf(paste0(filepath,"Hepatocytes_3_violinplot.pdf", sep = ""),width=12, height=25)
Stacked_VlnPlot(seurat_object = pbmc, features = mouse_gene,group.by = "group",
                x_lab_rotate = TRUE)

dev.off()
