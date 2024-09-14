
library(ddpcr)
library(scCDC)
filepath <- paste("E:/B组数据备份(4.29)/鉴定后总样本备份/",'all_sample_scCDC',".RDS", sep = "")
pbmc <- readRDS(filepath)
levels(pbmc@meta.data$seurat_clusters)

GCGs = ContaminationDetection(pbmc)
contamination_ratio = ContaminationQuantification(pbmc,rownames(GCGs))
seuratobj_corrected = ContaminationCorrection(pbmc,rownames(GCGs))
DefaultAssay(seuratobj_corrected) = "Corrected"

corrected_count_matrix = data.frame(seuratobj_corrected@assays[["Corrected"]]@counts)

gene_name = paste("肝组织", sep = "")
pbmc_hep <- PercentageFeatureSet(pbmc,features = "Hp",col.name = gene_name)
FeaturePlot(pbmc_hep,gene_name,pt.size = 2)
gene <- read.xlsx("E:/B_group/gene/gene.xlsx",
                  sheetIndex = 1,header = T,encoding = "UTF-8")
gene <- gene$gene
for (i in gene) {
  gene_name <- paste("肝细胞_", i, sep = "")
  pbmc_hep <- PercentageFeatureSet(pbmc, features = i, col.name = gene_name)
  
  # 打开PDF设备
  pdf_file <- paste("E:/B_group/UMAP/", gene_name, ".pdf", sep = "")
  pdf(pdf_file)  # Start writing to the PDF file
  
  # 绘制特征图
  p <- FeaturePlot(pbmc_hep, gene_name, pt.size = 2)
  print(p)  # 输出图形到PDF
  
  # 关闭PDF设备
  dev.off()
}


pc.num = 1:30
pbmc <- FindNeighbors(pbmc,reduction = "harmony", dims = pc.num)
pbmc <- FindClusters(pbmc, resolution = 1.2, graph.name = "SCT_snn")

sce <- FindClusters(sce, resolution = 1.2)
print(names(sce@graphs))
p4 <- DimPlot(sce, reduction = "umap", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
pdf("E:/B组数据备份(4.29)/鉴定后总样本备份/umap_cluster.pdf", width = 18, height = 18)
print(p4)
dev.off()
head(sce@assays$RNA@data["Hp", ])
FeaturePlot(pbmc,features = "Hp")
dev.off()
###鉴定
#手工标注
library(xlsx)
ident <- read.xlsx("E:/B_group/分簇/ident.xlsx",
                   sheetIndex = 1,header = T,encoding = "UTF-8")
new.cluster.ids <- c(ident$ident)
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
pbmc@meta.data$cluster <- pbmc@meta.data$seurat_clusters
levels(pbmc@meta.data$seurat_clusters) <- new.cluster.ids
levels(pbmc@meta.data$seurat_clusters)
#备份
filepath <- paste("E:/B_group/data_backup/B_group_sample",".RDS", sep = "")
saveRDS(pbmc, file= filepath)
pbmc <- readRDS('E:/B_group/data_backup/B_group_sample.RDS')

filepath <- paste("E:/B_group/sub_group/",'Hepatocytes',".RDS", sep = "")
filepath <- paste("E:/B_group/sub_group/",'Endothelial_cells',".RDS", sep = "")
#saveRDS(pbmc, file= filepath)
pbmc <- readRDS(filepath)

pc.num = 1:30
pbmc <- FindNeighbors(pbmc,reduction = "harmony", dims = pc.num)
pbmc <- FindClusters(pbmc, resolution = 1)

#UMAP
p4 <- DimPlot(pbmc, reduction = "umap", group.by = "ident",   pt.size=2, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
pdf("E:/B_group/UMAP/Hepatocytes.pdf", width = 6, height = 6)
p4
dev.off()

#亚群命名
new.cluster.ids <- c("ECs_1","ECs_2","ECs_3","ECs_1")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
pbmc@meta.data$cluster <- pbmc@meta.data$seurat_clusters
levels(pbmc@meta.data$seurat_clusters) <- new.cluster.ids
levels(pbmc@meta.data$seurat_clusters)
#差异基因
FindAllMarkers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(FindAllMarkers, file = "E:/B_group/高变基因/makers_HEP.csv")

#ClusterGVis
library(ClusterGVis)
library(org.Rn.eg.db)

# 筛选出符合条件的基因
filteredMarkers <- FindAllMarkers %>%
  dplyr::filter(p_val_adj < 0.05, avg_log2FC > 0.1)

# 从筛选后的基因中选择每个簇的前 30 个基因
topMarkers <- filteredMarkers %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 30, wt = avg_log2FC)

# 提取基因列表
topGenes <- topMarkers %>%
  dplyr::pull(gene)

# 查看结果
print(topGenes)

# [3] 以每个细胞单独的表达作为输入
st.data <- prepareDataFromscRNA(object = pbmc,
                                diffData = filteredMarkers,
                                showAverage = FALSE)

# enrich for clusters
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Rn.eg.db,
                        type = "BP",
                        organism = "rno",
                        pvalueCutoff = 0.05,
                        topn = 5,
                        seed = 5201314)

# 在所有簇的前 30 个基因中选择 40 个最特异的基因
markGenes <- topMarkers %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::slice(1:30) %>%
  dplyr::pull(gene)

plotpath <- paste("E:/B_group/高变基因/", sep = "")
pdf(paste0(plotpath,"diff_gene_GO_","HEP",".pdf",sep = ""),width = 14,height = 14)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           markGenes = topGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           col=color.pals ,
           cluster.order = c(1:3),
           go.col = rep(jjAnno::useMyCol("stallion",n = 3),each = 5),
           add.bar = T)
dev.off()

#Go
library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(ggplot2)
library(stringi)	
library(GOplot)
library(tidyverse)
R.utils::setOption("clusterProfiler.download.method",'auto')

cell_pbmc_all <- subset(FindAllMarkers, p_val_adj < 0.05 & avg_log2FC > 0.1 & cluster == "ECs_1")
cell_pbmc_all <- subset(topMarkers, cluster == "ECs_1")

kk=enrichGO(gene=row.names(cell_pbmc_all), 
            OrgDb=org.Rn.eg.db, 
            pvalueCutoff=0.01, 
            qvalueCutoff=0.05, 
            ont="all", 
            keyType = 'SYMBOL',
            readable=T)


eGo <- as.data.frame(kk)
eGo <- separate(data=eGo, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
eGo <- separate(data=eGo, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
eGo <- mutate(eGo, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor
#只用排名前10的term画图，先把BP、MF、CC三种类型各取top10，然后再组合在一起
eGoBP <- eGo %>% 
  filter(ONTOLOGY=="BP") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoCC <- eGo %>% 
  filter(ONTOLOGY=="CC") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGoMF <- eGo %>% 
  filter(ONTOLOGY=="MF") %>%
  filter(row_number() >= 1,row_number() <= 10)
eGo10 <- rbind(eGoBP,eGoMF,eGoCC)

p <- ggplot(eGo10,aes(enrichment_factor,Description)) + 
  geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
  scale_color_gradient(low="green",high ="red") + 
  labs(color=expression(-log[10](p_value)),size="Count", shape="Ontology",
       x="Enrichment Factor",y="GO term",title="GO enrichment") + 
  theme_bw()
p <- p + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free')

pdf("E:/B_group/Go/GO_enrichment5.pdf",width=12,height=8)
print(p)
dev.off()

cluster <- levels(pbmc@meta.data$seurat_clusters)
for (i in cluster) {
  cell_pbmc_all <- subset(topMarkers, cluster == i)
  
  kk=enrichGO(gene=cell_pbmc_all$gene, 
              OrgDb=org.Rn.eg.db, 
              pvalueCutoff=0.05, 
              qvalueCutoff=0.05, 
              ont="all", 
              keyType = 'SYMBOL',
              readable=T)
  
  
  eGo <- as.data.frame(kk)
  eGo <- separate(data=eGo, col=GeneRatio,into = c("GR1", "GR2"), sep = "/") #劈分GeneRatio为2列（GR1、GR2）
  eGo <- separate(data=eGo, col=BgRatio, into = c("BR1", "BR2"), sep = "/") #劈分BgRatio为2列（BR1、BR2）
  eGo <- mutate(eGo, enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2))) #计算Enrichment Factor
  #只用排名前10的term画图，先把BP、MF、CC三种类型各取top10，然后再组合在一起
  eGoBP <- eGo %>% 
    filter(ONTOLOGY=="BP") %>%
    filter(row_number() >= 1,row_number() <= 5)
  eGoCC <- eGo %>% 
    filter(ONTOLOGY=="CC") %>%
    filter(row_number() >= 1,row_number() <= 5)
  eGoMF <- eGo %>% 
    filter(ONTOLOGY=="MF") %>%
    filter(row_number() >= 1,row_number() <= 5)
  eGo10 <- rbind(eGoBP,eGoMF,eGoCC)
  
  p <- ggplot(eGo10,aes(enrichment_factor,Description)) + 
    geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) +
    scale_color_gradient(low="green",high ="red") + 
    labs(color=expression(-log[10](p_value)),size="Count", shape="Ontology",
         x="Enrichment Factor",y="GO term",title="GO enrichment") + 
    theme_bw()
  p <- p + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free')
  
  pdf(paste("E:/B_group/Go/GO_enrichment5_",i,".pdf",sep = ""),width=12,height=8)
  print(p)
  dev.off()
}
