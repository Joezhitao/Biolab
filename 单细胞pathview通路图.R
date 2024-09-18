library(Seurat)
library(clusterProfiler)
library(org.Rn.eg.db)
library(pathview)

filepath <- paste("E:/B_group/sub_group/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)
levels(pbmc@meta.data$group)

DEG <- FindMarkers(pbmc, ident.1 = "Bd1R", ident.2 = "Ad1R", genes.use = NULL,group.by = "group",
                   thresh.use = 0.25, test.use = "bimod", min.pct = 0.1,
                   print.bar = TRUE, only.pos = FALSE, max.cells.per.ident = Inf)
DEG$avg_log2FC
colnames(DEG)[colnames(DEG) == "avg_log2FC"] <- "logFC"
DEG$gene <- rownames(DEG)

eg <- bitr(
  DEG$gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Rn.eg.db
)

df_merged <- merge(DEG, eg, by.x = "gene", by.y = "SYMBOL")
genelistDEGs <- df_merged$logFC
names(genelistDEGs) <- df_merged$ENTREZID


#pathview
rno03320 <- pathview::pathview(
  gene.data  = genelistDEGs,
  pathway.id = "rno03320",
  species    = "rno",
  out.suffix = "my_analysis",
  limit = list(gene = as.integer(max(abs(genelistDEGs))))
)
help("pathview")


filepath <- paste("E:/B_group/sub_group/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)
levels(pbmc@meta.data$group)

group <- levels(pbmc@meta.data$group)
group <- setdiff(group, "Ad1R")
for (i in group) {
  
  DEG <- FindMarkers(pbmc, ident.1 = i, ident.2 = "Ad1R", genes.use = NULL,group.by = "group",
                     thresh.use = 0.25, test.use = "bimod", min.pct = 0.1,
                     print.bar = TRUE, only.pos = FALSE, max.cells.per.ident = Inf)
  DEG$avg_log2FC
  colnames(DEG)[colnames(DEG) == "avg_log2FC"] <- "logFC"
  DEG$gene <- rownames(DEG)
  
  eg <- bitr(
    DEG$gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Rn.eg.db
  )
  
  df_merged <- merge(DEG, eg, by.x = "gene", by.y = "SYMBOL")
  genelistDEGs <- df_merged$logFC
  names(genelistDEGs) <- df_merged$ENTREZID
  
  pathway <- c("rno04979","rno04920","rno04936","rno04910","rno04152")
  for (j in pathway) {
    setwd("E:/B_group/pathway/")
    rno <- pathview::pathview(
      gene.data  = genelistDEGs,
      pathway.id = j,
      species    = "rno",
      out.suffix = paste("_",i,sep = ""),
      limit = list(gene = as.integer(max(abs(genelistDEGs))))
    )
    
  }
}

#查看当前路径
getwd()
dev.off()
