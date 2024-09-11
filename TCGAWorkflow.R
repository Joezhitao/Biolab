library(TCGAWorkflowData)
library(DT)
#从 TCGA 数据库下载转录组分析
query_exp_lgg <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
# Get only first 20 samples to make example faster
query_exp_lgg$results[[1]] <- query_exp_lgg$results[[1]][1:20,]
GDCdownload(query_exp_lgg)
exp_lgg <- GDCprepare(
  query = query_exp_lgg
)

query_exp_gbm <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)
# Get only first 20 samples to make example faster
query_exp_gbm$results[[1]] <- query_exp_gbm$results[[1]][1:20,]
GDCdownload(query_exp_gbm)
exp_gbm <- GDCprepare(
  query = query_exp_gbm
)

###############################################################################
#临床数据下载
library(maftools)
# recovering data from TCGAWorkflowData package.
data(maf_lgg_gbm)

# To prepare for maftools we will also include clinical data
# For a mutant vs WT survival analysis 
# get indexed clinical patient data for GBM samples
gbm_clin <- GDCquery_clinic(project = "TCGA-GBM", type = "Clinical")
# get indexed clinical patient data for LGG samples
lgg_clin <- GDCquery_clinic(project = "TCGA-LGG", type = "Clinical")

###############################################################################
data("TCGA_LGG_Transcriptome_20_samples")
data("TCGA_GBM_Transcriptome_20_samples")

BiocManager::install("DESeq2")

exp_lgg_preprocessed <- TCGAanalyze_Preprocessing(
  object = exp_lgg,
  cor.cut = 0.6,    
  datatype = "unstranded",
  filename = "LGG_IlluminaHiSeq_RNASeqV2.png"
)

exp_gbm_preprocessed <- TCGAanalyze_Preprocessing(
  object = exp_gbm,
  cor.cut = 0.6, 
  datatype = "unstranded",
  filename = "GBM_IlluminaHiSeq_RNASeqV2.png"
)
exp_preprocessed <- cbind(
  exp_lgg_preprocessed, 
  exp_gbm_preprocessed
)

exp_normalized <- TCGAanalyze_Normalization(
  tabDF = cbind(exp_lgg_preprocessed, exp_gbm_preprocessed),
  geneInfo = TCGAbiolinks::geneInfoHT,
  method = "gcContent"
) # 60513   40

exp_filtered <- TCGAanalyze_Filtering(
  tabDF = exp_normalized,
  method = "quantile",
  qnt.cut =  0.25
)  # 44630   40

exp_filtered_lgg <- exp_filtered[
  ,substr(colnames(exp_filtered),1,12) %in% lgg_clin$bcr_patient_barcode
]

exp_filtered_gbm <-   exp_filtered[
  ,substr(colnames(exp_filtered),1,12) %in% gbm_clin$bcr_patient_barcode
]

diff_expressed_genes <- TCGAanalyze_DEA(
  mat1 = exp_filtered_lgg,
  mat2 = exp_filtered_gbm,
  Cond1type = "LGG",
  Cond2type = "GBM",
  fdr.cut = 0.01 ,
  logFC.cut = 1,
  method = "glmLRT"
)
nrow(diff_expressed_genes)
#EA: enrichment analysis
#In order to understand the underlying biological process of DEGs we performed an
#enrichment analysis using TCGAanalyze_EA_complete function.
ansEA <- TCGAanalyze_EAcomplete(
  TFname = "DEA genes LGG Vs GBM", 
  RegulonList = diff_expressed_genes$gene_name
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOBPTab = ansEA$ResBP,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOCCTab = ansEA$ResCC,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  GOMFTab = ansEA$ResMF,
  nRGTab = diff_expressed_genes$gene_name,
  nBar = 20
)

TCGAvisualize_EAbarplot(
  tf = rownames(ansEA$ResBP),
  filename = NULL,
  PathTab = ansEA$ResPat,
  nRGTab = rownames(diff_expressed_genes),
  nBar = 20
)

library(SummarizedExperiment)

# DEGs TopTable
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(
  FC_FDR_table_mRNA = diff_expressed_genes,
  typeCond1 = "LGG",
  typeCond2 = "GBM",
  TableCond1 = exp_filtered[,colnames(exp_filtered_lgg)],
  TableCond2 = exp_filtered[,colnames(exp_filtered_gbm)]
)

dataDEGsFiltLevel$GeneID <- 0

library(clusterProfiler)
# Converting Gene symbol to geneID
eg = as.data.frame(
  bitr(
    dataDEGsFiltLevel$mRNA,
    fromType = "ENSEMBL",
    toType = c("ENTREZID","SYMBOL"),
    OrgDb = "org.Hs.eg.db"
  )
)
eg <- eg[!duplicated(eg$SYMBOL),]
eg <- eg[order(eg$SYMBOL,decreasing=FALSE),]

dataDEGsFiltLevel <- dataDEGsFiltLevel[dataDEGsFiltLevel$mRNA %in% eg$ENSEMBL,]
dataDEGsFiltLevel <- dataDEGsFiltLevel[eg$ENSEMBL,]
rownames(dataDEGsFiltLevel) <- eg$SYMBOL

all(eg$SYMBOL == rownames(dataDEGsFiltLevel))

dataDEGsFiltLevel$GeneID <- eg$ENTREZID

dataDEGsFiltLevel_sub <- subset(dataDEGsFiltLevel, select = c("GeneID", "logFC"))
genelistDEGs <- as.numeric(dataDEGsFiltLevel_sub$logFC)
names(genelistDEGs) <- dataDEGsFiltLevel_sub$GeneID
library(pathview)
# pathway.id: hsa05214 is the glioma pathway
# limit: sets the limit for gene expression legend and color
hsa05214 <- pathview::pathview(
  gene.data  = genelistDEGs,
  pathway.id = "hsa05214",
  species    = "hsa",
  limit = list(gene = as.integer(max(abs(genelistDEGs))))
)

#########################################################
#allDEG2中log2FoldChange的名字改为logFC
# 修改列名
library(clusterProfiler)
library(org.Hs.eg.db)

colnames(allDEG2)[colnames(allDEG2) == "log2FoldChange"] <- "logFC"

eg <- bitr(
  allDEG2$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

df_merged <- merge(allDEG2, eg, by.x = "gene", by.y = "SYMBOL")
genelistDEGs <- df_merged$logFC
names(genelistDEGs) <- df_merged$ENTREZID

library(pathview)
# pathway.id: hsa05214 is the glioma pathway
# limit: sets the limit for gene expression legend and color
	
hsa04014 <- pathview::pathview(
  gene.data  = genelistDEGs,
  pathway.id = "hsa04014",
  species    = "hsa",
  limit = list(gene = as.integer(max(abs(genelistDEGs))))
)


library(clusterProfiler)
library(org.Hs.eg.db)

# 假设你有一组 Entrez ID
gene_list <- c("125058")

# 进行 KEGG 富集分析
kegg_result <- enrichKEGG(
  gene = gene_list,
  organism = "hsa",
  keyType = "kegg", # 因为你使用的是 KEGG 基因 ID
  pvalueCutoff = 0.05
)

# 查看结果
print(kegg_result)