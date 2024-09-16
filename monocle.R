library(Seurat)
library(tidyverse)
library(magrittr)
library(monocle)

filepath <- paste("E:/B_group/sub_group/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)
#Transfer into cds
cds <- as.CellDataSet(pbmc)
#Monocle2 process
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
options(future.globals.maxSize = 24 * 1024^3)  # 24GB

all_markers <- FindAllMarkers(pbmc, 
                              only.pos = TRUE, 
                              min.pct = 0.3, 
                              logfc.threshold = 0.25)
deg <- all_markers
deg <- deg[which(deg$cluster %in% unique(pbmc$seurat_clusters)), ]
sel.gene <- unique(deg$gene)
cds <- monocle::setOrderingFilter(cds, sel.gene)
cds <- monocle::reduceDimension(cds, method = 'DDRTree')

## ordering cells
cds <- monocle::orderCells(cds)

## ordering cells by assigning root nodes
GM_state <- function(cds){
  if (length(unique(cds$State)) > 1){
    T0_counts <- table(cds$State, cds$seurat_clusters)[,"Hepatocytes_1"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds <- monocle::orderCells(cds, root_state =  GM_state(cds))

monocle::plot_cell_trajectory(cds, color_by = "Pseudotime")
monocle::plot_cell_trajectory(cds, color_by = "seurat_clusters") 
