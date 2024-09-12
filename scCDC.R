library(devtools)
install_github("ZJU-UoE-CCW-LAB/scCDC")
install.packages("ddpcr")

library(ddpcr)
library(scCDC)
filepath <- paste("E:/B组数据备份(4.29)/鉴定后总样本备份/",'all_sample_decontX035_pc30',".RDS", sep = "")
pbmc <- readRDS(filepath)
levels(pbmc@meta.data$seurat_clusters)

GCGs = ContaminationDetection(pbmc)
contamination_ratio = ContaminationQuantification(pbmc,rownames(GCGs))
seuratobj_corrected = ContaminationCorrection(pbmc,rownames(GCGs))
DefaultAssay(seuratobj_corrected) = "Corrected"

corrected_count_matrix = data.frame(seuratobj_corrected@assays[["Corrected"]]@counts)

gene_name = paste("肝组织", sep = "")
pbmc_hep <- PercentageFeatureSet(sce,features = "Hp",col.name = gene_name)
p2 <- FeaturePlot(pbmc_hep,gene_name)
seuratobj_corrected@meta.data$cluster
levels(seuratobj_corrected@meta.data$cluster)

sce <- seuratobj_corrected

pc.num = 1:30
sce <- FindNeighbors(sce,reduction = "harmony", dims = pc.num)
sce <- FindClusters(sce, resolution = 1.2, graph.name = "SCT_snn")

sce <- FindClusters(sce, resolution = 1.2)
print(names(sce@graphs))
p4 <- DimPlot(sce, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
pdf("E:/B组数据备份(4.29)/鉴定后总样本备份/umap_去污染前.pdf", width = 18, height = 18)
FeaturePlot(pbmc,features = "Hp")
dev.off()
head(sce@assays$RNA@data["Hp", ])

dev.off()
