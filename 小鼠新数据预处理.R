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
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

setwd("E:/newdata/")
group_by <- c("ccl4_kong", "ccl4_fang")
for(i in 1:length(group_by)){
  assign(group_by[i], Read10X(data.dir = paste0(group_by[i], '/')))
}

for (i in group_by) {
  df <- get(i)
  seurat_object <- df %>% 
    CreateSeuratObject(min.cells = 3, min.features = 200) %>% 
    PercentageFeatureSet("^mt-", col.name = "percent.mt") %>% 
    subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
  seurat_object$group <- i
  assign(i, seurat_object)
}

for (i in group_by) {
  df <- get(i)
  pc.num=1:30
  #设置去doublet比例，根据细胞数据大小调整，参考官方标准
  num_cells <- ncol(df@assays$RNA@counts)
  numb <- function(x) {
    if (x <= 2000) {
      return(0.008)
    } else if (x > 2000 & x <= 4000) {
      return(0.016)
    } else if (x > 4000 & x <= 6000) {
      return(0.024)
    } else if (x > 6000 & x <= 8000) {
      return(0.032)
    } else if (x > 8000 & x <= 10000) {
      return(0.04)
    } else if (x > 10000 & x <= 12000) {
      return(0.048)
    } else if (x > 12000 & x <= 14000) {
      return(0.056)
    } else if (x > 14000 & x <= 16000) {
      return(0.064)
    } else if (x > 16000 & x <= 18000) {
      return(0.072)
    } else if (x > 18000) {
      return(0.08)
    } else {
      stop("Unexpected input")
    }
  }
  numb_value <- numb(num_cells)
  Find_doublet <- function(data){
    # 寻找最优pk值
    sweep.res.list <- paramSweep_v3(data, PCs = pc.num, sct = T)##sct指是否seurat做了SCT处理，做了为T。
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
    bcmvn <- find.pK(sweep.stats)
    pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)])
    bcmvn$pK[pK_bcmvn]
    
    # 期望doublet数量
    DoubletRate = numb_value  
    homotypic.prop <- modelHomotypic(data$seurat_clusters)   
    nExp_poi <- round(DoubletRate*ncol(data))
    nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
    
    # 鉴定doublets
    data <- doubletFinder_v3(data, PCs = pc.num, pN = 0.25, pK = pK_bcmvn,
                             nExp = nExp_poi.adj, reuse.pANN = F, sct = T)
    colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
    data
  }
  
  seurat_object <- df %>% 
    SCTransform() %>% ##批次校正和归一化
    RunPCA(verbose = F) %>% ##主成分分析（PCA）
    RunUMAP(dims=pc.num) %>% ##UMAP降维分析
    FindNeighbors(dims = pc.num) %>% ##细胞距离计算，并归类到群
    FindClusters(resolution = 0.1) %>% 
    Find_doublet()
  
  assign(i, seurat_object)
  
  filepath <- paste("E:/newdata/doublet_", i, ".RDS", sep = "")
  saveRDS(seurat_object, file= filepath)
}

for (i in group_by) {
  df <- get(i)
  df <- subset(df,subset=doublet_info=="Singlet")
  assign(i, df)
  filepath <- paste("E:/newdata/doublet(clear)_", i, ".RDS", sep = "")
  saveRDS(df, file= filepath)
}

objects_list <- lapply(group_by, get)
pbmc <- merge(x = objects_list[[1]], y = objects_list[-1], add.cell.ids = group_by)
filepath <- paste("E:/newdata/integrated data(doublet_merge)",".RDS", sep = "")
saveRDS(pbmc, file= filepath)

pc.num=1:30
pbmc <- pbmc %>% 
  SCTransform() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>%
  RunPCA(npcs = 50, verbose = FALSE) %>%
  RunHarmony("group", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony", dims = pc.num) %>% 
  FindNeighbors(reduction = "harmony", dims = pc.num) %>% 
  FindClusters(resolution = 0.1) %>% 
  identity()

library(decontX)
scobj <- pbmc
# 提取矩阵(genes x cells)
# 需要注意的是，你应该查看你的矩阵是否已经注释了行名和列明，以及如果你的矩阵是scanpy对象中提取的矩阵应该行列转置
counts <- scobj@assays$RNA@counts
# 计算
# 需要把结果储存在新变量
decontX_results <- decontX(counts) 
# 你可以使用str()查看结果的形式
# RNA污染的计算结果储存在：decontX_results$contamination
# 他是一个数值型向量，长度与细胞数量一致，顺序也与矩阵的colnames一致
# 我们可以直接把他写入metadata
scobj$Contamination =decontX_results$contamination
head(scobj@meta.data) # 可以查看一下结果

# 我们对他进行可视化
library(ggplot2)
png("E:/newdata/污染情况(小鼠).png",units = "in",width = 10,height = 10,
    res = 1000)
FeaturePlot(scobj, 
            features = 'Contamination', 
            raster=FALSE       # 细胞过多时候需要加这个参数
) + 
  scale_color_viridis_c()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  xlab('scVI_UMAP_1')+
  ylab('scVI_UMAP_2')
dev.off()
# 在这里看到，对比前文最开始的状态，我们经过了严格的质控之后，散在的细胞已经少了很多，但是仍然有亚群之间不清晰的界限
# 我们可视化Contamination，惊奇地发现这些边缘的毛躁就是Contamination较高地细胞

low_con_scobj = scobj[,scobj$Contamination < 0.25] #保留小于等于0.25的细胞
filepath <- paste("E:/newdata/",'all_sample_decontX',".RDS", sep = "")
saveRDS(low_con_scobj, file= filepath)

pbmc <- readRDS("E:/newdata/all_sample_decontX.RDS")

levels(pbmc@meta.data$seurat_clusters)
pc.num <- 1:20
pbmc <- FindNeighbors(pbmc,reduction = "harmony", dims = pc.num)
pbmc <- FindClusters(pbmc, resolution = 1.2)

DimPlot(pbmc, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)

combined_markers <- FindAllMarkers(object = pbmc, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25)

all <- combined_markers %>% group_by(cluster)

write.csv(all, 
          file = "E:/newdata/combined_all_markers(resolution=1.2).csv", 
          quote = FALSE, 
          row.names = FALSE)
#细胞鉴定
