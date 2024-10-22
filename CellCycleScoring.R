library(Seurat)
library(RCurl)
library(AnnotationHub)
library(ensembldb)
library(dplyr)

filepath <- paste("E:/B_group/data_backup/B_group_sample",".RDS", sep = "")
pbmc <- readRDS(filepath)
pbmc@meta.data$seurat_clusters
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv")
cell_cycle_genes <- read.csv("H:/google下载内容/Mus_musculus.csv")

# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah, 
              pattern = c("Mus musculus", "EnsDb"), 
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "S") %>%
  pull("gene_name")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_markers %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("gene_name")

# Perform cell cycle scoring
seurat_phase <- CellCycleScoring(pbmc,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)


# Identify the most variable genes if it hasn't been run
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA and color by cell cycle phase
seurat_phase <- RunPCA(seurat_phase)

# Visualize the PCA, grouping by cell cycle phase
png("E:/B_group/data_backup/CellCycleScoring.png", width = 8, height = 8, units = "in",res = 600)
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase")
dev.off()

plot <- RidgePlot(seurat_phase, features = c("Mki67", "ccnd1"), ncol = 2)

png("E:/B_group/data_backup/B_group_sample_CellCycleScoring.png", width = 16, height = 8, units = "in",res = 600)
print(plot)
dev.off()
