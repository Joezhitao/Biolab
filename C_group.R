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
