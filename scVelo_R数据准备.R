library(Matrix)
library(Seurat)
library(SeuratData)
library(SeuratDisk)

#亚群选择
sample <- c("Endothelial_cells","Hepatocytes")

for (i in sample) {
  filepath <- paste("E:/B_group/sub_group/",i,".RDS", sep = "")
  pbmc <- readRDS(filepath)
  
  #设置默认目录，和scVelo一致:
  # Define the base directory
  dir.create(file.path("E:/B_group/velo/", i), recursive = TRUE, showWarnings = FALSE)
  out_data_dir <- file.path("E:/B_group/velo/", i,"/")
  setwd(out_data_dir)
  
  # assuming that you have some Seurat object called seurat_obj:
  # save metadata table:
  #提取UMAP的前两个维度，在二维平面上绘制这些细胞的UMAP图
  seurat_obj <- pbmc
  seurat_obj$barcode <- colnames(seurat_obj)
  seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
  seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
  write.csv(seurat_obj@meta.data, file='metadata.csv', quote=F, row.names=F)
  
  # write expression counts matrix
  counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
  writeMM(counts_matrix, file=paste0(out_data_dir, 'counts.mtx'))
  
  # write dimesnionality reduction matrix, in this example case pca matrix
  write.csv(seurat_obj@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
  
  # write gene names
  write.table(
    data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
    quote=F,row.names=F,col.names=F
  )
  
}

#亚群分组选择
for (i in sample) {
  filepath <- paste("E:/B_group/sub_group/",i,".RDS", sep = "")
  pbmc <- readRDS(filepath)
  
  group <- levels(pbmc@meta.data$group)
  for (j in group) {
    #设置默认目录，和scVelo一致:
    # Define the base directory
    dir.create(file.path("E:/B_group/velo/", i,"/",j), recursive = TRUE, showWarnings = FALSE)
    out_data_dir <- file.path("E:/B_group/velo/", i,"/",j,"/")
    setwd(out_data_dir)
    
    # assuming that you have some Seurat object called seurat_obj:
    # save metadata table:
    #提取UMAP的前两个维度，在二维平面上绘制这些细胞的UMAP图
    seurat_obj <- pbmc
    
    seurat_obj <- seurat_obj[,seurat_obj@meta.data$group %in% j]
    seurat_obj@meta.data$group <- droplevels(seurat_obj@meta.data$group)
    
    seurat_obj$barcode <- colnames(seurat_obj)
    seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
    seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
    write.csv(seurat_obj@meta.data, file='metadata.csv', quote=F, row.names=F)
    
    # write expression counts matrix
    counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
    writeMM(counts_matrix, file=paste0(out_data_dir, 'counts.mtx'))
    
    # write dimesnionality reduction matrix, in this example case pca matrix
    write.csv(seurat_obj@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)
    
    # write gene names
    write.table(
      data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
      quote=F,row.names=F,col.names=F
    )
  }
}
