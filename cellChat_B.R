suppressMessages(if(!require(CellChat))devtools::install_github("sqjin/CellChat"))
suppressMessages(if(!require(ggplot2))install.packages('ggplot2'))
suppressMessages(if(!require(patchwork))install.packages('patchwork') )
suppressMessages(if(!require(ggalluvial))install.packages('ggalluvial'))
suppressMessages(if(!require(igraph))install.packages('igraph'))
suppressMessages(if(!require(dplyr))install.packages('dplyr'))
suppressMessages(options(stringsAsFactors = FALSE))
suppressMessages(options(futrue.globlas.Maxsize=2*1024**3))
suppressWarnings(suppressMessages(future::plan("multisession", workers = 10)))
library(patchwork)

#for循环跑cellchat
sce <- pbmc
group <- levels(pbmc@meta.data$group)
for (i in group) {
  pbmc <- sce
  pbmc <- pbmc[,pbmc@meta.data$group %in% i]
  pbmc@meta.data$group <- droplevels(pbmc@meta.data$group)
  
  seurat_obj <- pbmc
  data.input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
  identity <- subset(seurat_obj@meta.data, select = "seurat_clusters")
  cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "seurat_clusters")
  groupSize <- as.numeric(table(cellchat@idents))
  
  #CellChatDB <- CellChatDB.human
  CellChatDB <- CellChatDB.mouse
  
  unique(CellChatDB$interaction$annotation)
  dplyr::glimpse(CellChatDB$interaction)
  
  type <- c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")
  for (j in type) {
    CellChatDB.use <- subsetDB(CellChatDB, search = j)#选择"Secreted Signaling"这种类型的细胞通讯
    cellchat@DB <- CellChatDB.use
    options(future.globals.maxSize = 1024*1024*2048)#设置future包中的全局变量最大内存大小
    cellchat <- subsetData(cellchat,features = NULL) %>% 
      identifyOverExpressedGenes() %>% 
      identifyOverExpressedInteractions() %>% 
      projectData(PPI.mouse) %>% 
      computeCommunProb(raw.use = TRUE) %>% 
      filterCommunication(min.cells = 3) %>% 
      computeCommunProbPathway() %>% 
      aggregateNet()
    df.net <- subsetCommunication(cellchat)#查看通路信息以及配体受体信息的表格
    groupSize <- as.numeric(table(cellchat@idents))
    par(mfrow = c(1,2), xpd=TRUE)
    filepath <- paste("E:/C组数据/单细胞结果/cellchat/细胞通讯_",i,"_",j,".png", sep = "")
    png(filepath, units = "in", width = 10, height = 10, res = 800)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    dev.off()
  }
}