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
library(ClusterGVis)
library(org.Rn.eg.db)
#c组数据，质控完成，去污染，未重新聚类
filepath <- paste("E:/B组数据备份(4.29)/鉴定后总样本备份/",'c_group_sample_decontX035',".RDS", sep = "")
pbmc <- readRDS(filepath)
#重新聚类，num：1：30，PC：30，reslution:1.2
pc.num=1:30
pbmc <- pbmc %>% 
  SCTransform() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunHarmony("group", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = pc.num) %>% 
  FindNeighbors(reduction = "harmony", dims = pc.num) %>% 
  FindClusters(resolution = 1.2) %>% 
  identity()

p4 <- DimPlot(pbmc, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)

p4

#细胞鉴定
#手工标注
ident <- read.xlsx("E:/B组数据备份(4.29)/细胞鉴定表格/ident_c.xlsx",
                   sheetIndex = 1,header = T,encoding = "UTF-8")
new.cluster.ids <- c(ident$ident)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
pbmc@meta.data$cluster <- pbmc@meta.data$seurat_clusters
levels(pbmc@meta.data$seurat_clusters) <- new.cluster.ids
levels(pbmc@meta.data$seurat_clusters)



saveRDS(pbmc, file = "E:/B组数据备份(4.29)/鉴定后总样本备份/c_group.RDS")
##抽组，对比空白和C15天左右的差异基因，基于ClusterGVis包
levels(pbmc@meta.data$group)
group = c("Ad1R","Cd15L","Cd15R")
pbmc = pbmc[, pbmc@meta.data$group %in% group]#抽组
pbmc@meta.data$group <- droplevels(pbmc@meta.data$group)

#去掉group中的空白组,!不能与=间有空格
group <- group[group !="Ad1R"]

#提取肝细胞
cluster <- c("Hepatocytes")
pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% cluster]#抽组
pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)

#肝细胞分簇,跑ClusterGVis热图时,样本seurat_clusters至少要大于1簇
pbmc <- pbmc %>% 
  FindNeighbors(reduction = "harmony", dims = pc.num) %>% 
  FindClusters(resolution = 0.2)

#循环跑不同样本间的差异基因，基于ClusterGVis包
for (i in group) {
  group1 <- group[group != i]
  for (j in group1) {
    #找样本间差异基因C15L
    pbmc.C15L <- FindMarkers(pbmc, ident.1 = i, ident.2 = j, group.by = 'group',
                             min.pct = 0.25, logfc.threshold = 0.25)
    # get top 10 genes
    Top.markers <- pbmc.C15L %>%
      dplyr::top_n(n = 80, wt = avg_log2FC)
    #Top.markers添加一列基因，以行名为基因名，列名叫gene
    Top.markers$gene <- rownames(Top.markers)
    #Top.markers添加一列cluster,所有值为C15L
    Top.markers$cluster <- i
    # [3] 以每个细胞单独的表达作为输入
    st.data <- prepareDataFromscRNA(object = pbmc,
                                    diffData = Top.markers,
                                    showAverage = FALSE)
    
    # enrich for clusters
    enrich <- enrichCluster(object = st.data,
                            OrgDb = org.Rn.eg.db,
                            type = "BP",
                            organism = "hsa",
                            pvalueCutoff = 0.5,
                            topn = 5,
                            seed = 5201314)
    #从unique(Top.markers$gene)中随机抽取40个不重复的基因
    markGenes = unique(Top.markers$gene)[sample(1:length(unique(Top.markers$gene)),80,
                                                replace = F)]
    #获取差异基因比较样本的数量
    num <- levels(as.factor(Top.markers$cluster))
    num_elements <- length(num)
    #作图
    plotpath <- paste("E:/B组数据备份(4.29)/单细胞结果/GVis热图/", sep = "")
    png(paste0(plotpath,"diff_gene_GO_",i,"_",j,".png",sep = ""),width = 30,height = 14, units = "in", res = 800)
    visCluster(object = st.data,
               plot.type = "both",
               column_names_rot = 45,
               show_row_dend = F,
               markGenes = markGenes,
               markGenes.side = "left",
               annoTerm.data = enrich,
               line.side = "left",
               col=color.pals ,
               cluster.order = c(1:num_elements),
               go.col = rep(jjAnno::useMyCol("stallion",n = num_elements),each = 5),
               add.bar = T)
    dev.off()
    
  }
  
}

#亚群
#鉴定后数据
filepath <- paste("E:/B组数据备份(4.29)/鉴定后总样本备份/",'c_group',".RDS", sep = "")
pbmc <- readRDS(filepath)

#提取肝细胞
cluster <- c("Hepatocytes")
pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% cluster]#抽组
pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)

pc.num=1:30
pbmc <- pbmc %>%
  SCTransform() %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = pc.num, verbose=FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = pc.num) %>% 
  FindClusters(resolution = 0.1)

filepath <- paste("E:/B组数据备份(4.29)/C组亚群样本//Hepatocytes",".RDS", sep = "")
pbmc <- readRDS(filepath)
saveRDS(pbmc, file= filepath)

pbmc <- FindClusters(pbmc,resolution = 0.2)
p4 <- DimPlot(pbmc, reduction = "umap", group.by = "seurat_clusters",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)

plotpath <- paste("E:/B组数据备份(4.29)/单细胞结果/亚群/", sep = "")
png(paste0(plotpath,"UMAP_鉴定",".png",sep = ""),width = 6,height = 6, units = "in", res = 800)
p4
dev.off()

#亚群命名
new.cluster.ids <- c("Hepatocytes_1","Hepatocytes_2","Hepatocytes_3","Hepatocytes_1","Hepatocytes_2")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
pbmc@meta.data$cluster <- pbmc@meta.data$seurat_clusters
levels(pbmc@meta.data$seurat_clusters) <- new.cluster.ids
levels(pbmc@meta.data$seurat_clusters)

#分组排序
cluster_order <- c("Ad1R","Cd1L","Cd1R","Cd15L","Cd15R","Cd30L","Cd30R")
pbmc@meta.data$group <- factor(pbmc@meta.data$group, levels = cluster_order)

#亚群差异基因
combined_markers <- FindAllMarkers(object = pbmc, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25)

all <- combined_markers %>% group_by(cluster)

write.csv(all, 
          file = "E:/B组数据备份(4.29)/单细胞结果/亚群//combined_all_markers(resolution=0.1).csv", 
          quote = FALSE, 
          row.names = FALSE)

#细胞比例
#细胞比例
levels(pbmc@meta.data$group)
pbmc$celltype <- Idents(pbmc)
group_b = levels(pbmc@meta.data$group)

# 翻转group_b的排序
group_c = rev(group_b)

# 查看翻转后的结果
print(group_c)
#随机生成14个色彩鲜明的颜色
library(RColorBrewer)
colors = brewer.pal(14, "Set3")#集合只有12个颜色

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
write_csv(cell_counts_tbl, path = "E:/B组数据备份(4.29)/单细胞结果/细胞比例/cell_counts_tbl.csv")
#比例图
png("E:/B组数据备份(4.29)/单细胞结果/细胞比例/C组细胞比例(7样本).tiff",units = "in",width = 20,height = 8,res = 600)
ggplot(data = cell_counts, aes(x = group, fill = celltype)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = result) +
  coord_flip() +
  scale_y_reverse()+theme(legend.text = element_text(size = 20),text=element_text(size=20),
                          axis.text = element_text(size=15))
dev.off()