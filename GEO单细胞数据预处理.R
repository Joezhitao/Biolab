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

setwd("E:/Thyroid_GSM")
group_by <- c("GSM6153597", "GSM5493638")
for(i in 1:length(group_by)){
  assign(group_by[i], Read10X(data.dir = paste0(group_by[i], '/')))
}

for (i in group_by) {
  df <- get(i)
  seurat_object <- df %>% 
    CreateSeuratObject(min.cells = 3, min.features = 200) %>% 
    PercentageFeatureSet("^MT-", col.name = "percent.mt") %>% 
    subset(subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
  seurat_object$group <- i
  assign(i, seurat_object)
}


