library(monocle)
library(Seurat)
library(dplyr)
library(magrittr)
library(patchwork)
library(dplyr)

filepath <- paste("E:/B组数据备份(4.29)/去污染后细胞亚群备份/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)

levels(pbmc@meta.data$seurat_clusters)

group = c("Ad1R","Bd1R","Bd8R","Bd15R","Bd30R")
for (i in group) {
  #替换i中的空格为下划线
  k <- gsub(" ","_",i)
  sample <- pbmc
  sample = sample[, sample@meta.data$group %in% i]#抽组
  sample@meta.data$group <- droplevels(sample@meta.data$group)
  #在E:/B组数据备份(4.29)/单细胞结果/亚群拟时序分析/这个文件夹中创建一个叫k的文件夹
  dir.create(paste0("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/", k))
  
  Mono.cds = sample
  #seurat文件转成CellDataSet格式
  sample_ann <- Mono.cds@meta.data
  Mono_matrix <- as(as.matrix(GetAssayData(Mono.cds, slot = "counts")),
                    "sparseMatrix")
  feature_ann <- data.frame(gene_id = rownames(Mono_matrix),
                            gene_short_name = rownames(Mono_matrix))
  rownames(feature_ann) <- rownames(Mono_matrix)
  Mono_fd <- new("AnnotatedDataFrame", data = feature_ann)
  Mono_pd <- new("AnnotatedDataFrame", data = sample_ann)
  #CellDataSet
  Mono.cds <- newCellDataSet(Mono_matrix, phenoData = Mono_pd, featureData = Mono_fd, 
                             expressionFamily = negbinomial.size())
  #查看phenodata、featuredata
  head(pData(Mono.cds))
  head(fData(Mono.cds))
  
  Mono.cds <- Mono.cds %>% 
    estimateSizeFactors() %>% 
    estimateDispersions(cores = 10, relative_expr = TRUE) %>% 
    detectGenes(min_expr = 0.1)
  
  fData(Mono.cds)$use_for_ordering <- fData(Mono.cds)$num_cells_expressed > 0.05 * ncol(Mono.cds)
  expr.genes <- fData(Mono.cds)$gene_short_name[fData(Mono.cds)$use_for_ordering]
  Mono.cds <- Mono.cds %>% 
    setOrderingFilter(ordering_genes = expr.genes) %>% 
    reduceDimension(max_components = 2, method = "DDRTree") %>% 
    orderCells()
  
  #pic
  pic_plot_ordering_genes <- plot_ordering_genes(Mono.cds)
  filepath = paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/",k,"/","平均表达量",".pdf", sep = "")
  pdf(filepath, width = 8, height = 8)
  print(pic_plot_ordering_genes)
  dev.off()
  
  pic_plot_pc_variance_explained <- plot_pc_variance_explained(Mono.cds, return_all = F)
  filepath = paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/",k,"/","pc_variance_explained",".pdf", sep = "")
  pdf(filepath, width = 8, height = 8)
  print(pic_plot_pc_variance_explained)
  dev.off()
  
  
  Mono.cds_expressed_genes <- row.names(subset(fData(Mono.cds),
                                               num_cells_expressed >= 10))
  Mono.cds_filtered <- Mono.cds[Mono.cds_expressed_genes, ]
  clustering_DEG_genes <- differentialGeneTest(Mono.cds[Mono.cds_expressed_genes,], 
                                               fullModelFormulaStr = "~seurat_clusters", cores = 10)
  #选取差异表达前1000的基因
  Mono.cds_ordering_genes <- 
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
  Mono.cds <- Mono.cds %>% 
    setOrderingFilter(ordering_genes = Mono.cds_ordering_genes) %>% 
    reduceDimension(method = "DDRtree")
  suppressWarnings(Mono.cds <- orderCells(Mono.cds, reverse = T))
  #reverse参数根据先验知识更改轨迹图起点
  #pic
  pic_trajectory <- plot_cell_trajectory(Mono.cds, color_by = "seurat_clusters")|
    plot_cell_trajectory(Mono.cds, color_by = "Pseudotime")
  filepath <- paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/",k,"/","celltype_Pseudotime",
                    ".pdf", sep = "")
  pdf(filepath, width = 18, height = 6)
  print(pic_trajectory)
  dev.off()
  #State
  pic_trajectory_State <- plot_cell_trajectory(Mono.cds, color_by = "State")
  filepath <- paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/",k,"/","State",
                    ".pdf", sep = "")
  pdf(filepath, width = 6, height = 6)
  print(pic_trajectory_State)
  dev.off()
  #separate by celltype
  pic_trajectory_celltype <- plot_cell_trajectory(Mono.cds, color_by = "State") +
    facet_wrap(~seurat_clusters, nrow = 1)
  filepath <- paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/",k,"/","sep_by_celltype",
                    ".pdf", sep = "")
  pdf(filepath, width = 12, height = 4)
  print(pic_trajectory_celltype)
  dev.off()
  
  #寻找随着分化时间逐渐升高或降低的基因表达热图和分布图
  cds_subset <- Mono.cds[row.names(clustering_DEG_genes)
                         [order(clustering_DEG_genes$qval)][1:200], ]
  diff_test_res <- differentialGeneTest(cds_subset,
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  diff_test_res[, c("gene_short_name", "pval", "qval")] %>% head()
  sig.gene <- row.names(subset(diff_test_res, qval < 0.05))
  #sig.gene_plot <- plot_genes_in_pseudotime(cds_subset[sig.gene],  
  #color_by = 'seurat_clusters', ncol = 5)
  #filepath <- paste("E:/B组数据备份(4.29)/单细胞结果/亚群拟时序分析/","sig_200_gene",".png", sep = "")
  #ggsave(filename = filepath, sig.gene_plot, width = 6, height = 6, dpi = 600)
  #numbs是levels(sample@meta.data$seurat_clusters)的个数
  numbs <- length(levels(sample@meta.data$seurat_clusters))
  filepath <- paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/",k,"/","sig_200_gene_heatmap",".pdf", sep = "")
  pdf(filepath, width = 10, height = 25)
  plot_pseudotime_heatmap(cds_subset[sig.gene, ], 
                          num_clusters = numbs, 
                          cores = 10, 
                          show_rownames = T)
  dev.off()

}


sample <- pbmc
Mono.cds = sample
#seurat文件转成CellDataSet格式
sample_ann <- Mono.cds@meta.data
Mono_matrix <- as(as.matrix(GetAssayData(Mono.cds, slot = "counts")),
                  "sparseMatrix")
feature_ann <- data.frame(gene_id = rownames(Mono_matrix),
                          gene_short_name = rownames(Mono_matrix))
rownames(feature_ann) <- rownames(Mono_matrix)
Mono_fd <- new("AnnotatedDataFrame", data = feature_ann)
Mono_pd <- new("AnnotatedDataFrame", data = sample_ann)
#CellDataSet
Mono.cds <- newCellDataSet(Mono_matrix, phenoData = Mono_pd, featureData = Mono_fd, 
                           expressionFamily = negbinomial.size())
#查看phenodata、featuredata
head(pData(Mono.cds))
head(fData(Mono.cds))

Mono.cds <- Mono.cds %>% 
  estimateSizeFactors() %>% 
  estimateDispersions(cores = 10, relative_expr = TRUE) %>% 
  detectGenes(min_expr = 0.1)

fData(Mono.cds)$use_for_ordering <- fData(Mono.cds)$num_cells_expressed > 0.05 * ncol(Mono.cds)
expr.genes <- fData(Mono.cds)$gene_short_name[fData(Mono.cds)$use_for_ordering]
Mono.cds <- Mono.cds %>% 
  setOrderingFilter(ordering_genes = expr.genes) %>% 
  reduceDimension(max_components = 2, method = "DDRTree") %>% 
  orderCells(reverse = T)

#pic
pic_plot_ordering_genes <- plot_ordering_genes(Mono.cds)
filepath = paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/","/","平均表达量",".pdf", sep = "")
pdf(filepath, width = 8, height = 8)
print(pic_plot_ordering_genes)
dev.off()

pic_plot_pc_variance_explained <- plot_pc_variance_explained(Mono.cds, return_all = F)
filepath = paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/","/","pc_variance_explained",".pdf", sep = "")
pdf(filepath, width = 8, height = 8)
print(pic_plot_pc_variance_explained)
dev.off()


Mono.cds_expressed_genes <- row.names(subset(fData(Mono.cds),
                                             num_cells_expressed >= 10))
Mono.cds_filtered <- Mono.cds[Mono.cds_expressed_genes, ]
clustering_DEG_genes <- differentialGeneTest(Mono.cds[Mono.cds_expressed_genes,], 
                                             fullModelFormulaStr = "~seurat_clusters", cores = 10)
#选取差异表达前1000的基因
Mono.cds_ordering_genes <- 
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
Mono.cds <- Mono.cds %>% 
  setOrderingFilter(ordering_genes = Mono.cds_ordering_genes) %>% 
  reduceDimension(method = "DDRtree")
suppressWarnings(Mono.cds <- orderCells(Mono.cds, reverse = T))
#reverse参数根据先验知识更改轨迹图起点
#pic
pic_trajectory <- plot_cell_trajectory(Mono.cds, color_by = "seurat_clusters")|
  plot_cell_trajectory(Mono.cds, color_by = "Pseudotime")
filepath <- paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/","/","celltype_Pseudotime",
                  ".pdf", sep = "")
pdf(filepath, width = 18, height = 6)
print(pic_trajectory)
dev.off()
#State
pic_trajectory_State <- plot_cell_trajectory(Mono.cds, color_by = "State")
filepath <- paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/","/","State",
                  ".pdf", sep = "")
pdf(filepath, width = 6, height = 6)
print(pic_trajectory_State)
dev.off()
#separate by celltype
pic_trajectory_celltype <- plot_cell_trajectory(Mono.cds, color_by = "State") +
  facet_wrap(~seurat_clusters, nrow = 1)
filepath <- paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/","/","sep_by_celltype",
                  ".pdf", sep = "")
pdf(filepath, width = 12, height = 4)
print(pic_trajectory_celltype)
dev.off()

#寻找随着分化时间逐渐升高或降低的基因表达热图和分布图
cds_subset <- Mono.cds[row.names(clustering_DEG_genes)
                       [order(clustering_DEG_genes$qval)][1:200], ]
diff_test_res <- differentialGeneTest(cds_subset,
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res[, c("gene_short_name", "pval", "qval")] %>% head()
sig.gene <- row.names(subset(diff_test_res, qval < 0.05))
#sig.gene_plot <- plot_genes_in_pseudotime(cds_subset[sig.gene],  
#color_by = 'seurat_clusters', ncol = 5)
#filepath <- paste("E:/B组数据备份(4.29)/单细胞结果/亚群拟时序分析/","sig_200_gene",".png", sep = "")
#ggsave(filename = filepath, sig.gene_plot, width = 6, height = 6, dpi = 600)
#numbs是levels(sample@meta.data$seurat_clusters)的个数
numbs <- length(levels(sample@meta.data$seurat_clusters))
filepath <- paste("E:/B组数据备份(4.29)/单细胞结果/拟时序分析/","/","sig_200_gene_heatmap",".pdf", sep = "")
pdf(filepath, width = 10, height = 25)
plot_pseudotime_heatmap(cds_subset[sig.gene, ], 
                        num_clusters = numbs, 
                        cores = 10, 
                        show_rownames = T)
dev.off()