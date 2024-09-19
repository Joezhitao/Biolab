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

#####################
#for循环跑不同组与空白组对比的差异基因的pathview图

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

#进一步提取亚群for循环跑不同组与空白组对比的差异基因的pathview图
filepath <- paste("E:/B_group/sub_group/",'Hepatocytes',".RDS", sep = "")
pbmc <- readRDS(filepath)
sce <- pbmc
group <- levels(pbmc@meta.data$group)
group <- setdiff(group, "Ad1R")


cluster <- levels(pbmc@meta.data$seurat_clusters)

pathway_liver <- c("rno00010","rno00020","rno04010","rno04151","rno04668","rno04621","rno03320","rno04310","rno04064",
                   "rno04210","rno04350","rno04066","rno04630","rno04115","rno04140","rno04110","rno04390","rno04330",
                   "rno04151","rno04150","rno01521","rno04370","rno04012","rno04068","rno04340","rno04611","rno04610",
                   "rno00061","rno00071","rno00062","rno00100","rno00120","rno00121","rno00564","rno00600","rno00561",
                   "rno04910","rno04979","rno04920","rno04215","rno04142","rno03010")

for (k in 1:2) {
  if (k == 1) {
    subset <- cluster[3]  # 仅选择 "Hepatocytes_1"
    cat("Processing:", subset, "\n")
    # 这里添加处理 "Hepatocytes_1" 的代码
    pbmc <- sce
    pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% subset]#抽组
    pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)
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
      
      pathway <- pathway_liver
      for (j in pathway) {
        setwd("E:/B_group/pathway/")
        rno <- pathview::pathview(
          gene.data  = genelistDEGs,
          pathway.id = j,
          species    = "rno",
          out.suffix = paste("_",i,"_hep3",sep = ""),
          limit = list(gene = as.integer(max(abs(genelistDEGs))))
        )
        
      }
    }
  } else if (k == 2) {
    subset <- cluster[1:2]  # 选择 "Hepatocytes_2" 和 "Hepatocytes_3"
    cat("Processing:", paste(subset, collapse = ", "), "\n")
    # 这里添加处理 "Hepatocytes_2" 和 "Hepatocytes_3" 的代码
    pbmc <- sce
    pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% subset]#抽组
    pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)
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
      
      pathway <- pathway_liver
      for (j in pathway) {
        setwd("E:/B_group/pathway/")
        rno <- pathview::pathview(
          gene.data  = genelistDEGs,
          pathway.id = j,
          species    = "rno",
          out.suffix = paste("_",i,"_hep1_2",sep = ""),
          limit = list(gene = as.integer(max(abs(genelistDEGs))))
        )
        
      }
    }
  }
}


##tryCatch 捕获错误
for (k in 1:2) {
  if (k == 1) {
    subset <- cluster[3]  # 仅选择 "Hepatocytes_1"
    cat("Processing:", subset, "\n")
    pbmc <- sce
    pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% subset] # 抽组
    pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)
    for (i in group) {
      DEG <- FindMarkers(pbmc, ident.1 = i, ident.2 = "Ad1R", genes.use = NULL, group.by = "group",
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
      
      pathway <- pathway_liver
      for (j in pathway) {
        setwd("E:/B_group/pathway/")
        # 使用 tryCatch 捕获错误
        tryCatch({
          rno <- pathview::pathview(
            gene.data  = genelistDEGs,
            pathway.id = j,
            species    = "rno",
            out.suffix = paste("_", i, "_hep3", sep = ""),
            limit = list(gene = as.integer(max(abs(genelistDEGs))))
          )
        }, error = function(e) {
          cat("Error in processing pathway", j, ":", conditionMessage(e), "\n")
        })
      }
    }
  } else if (k == 2) {
    subset <- cluster[1:2]  # 选择 "Hepatocytes_2" 和 "Hepatocytes_3"
    cat("Processing:", paste(subset, collapse = ", "), "\n")
    pbmc <- sce
    pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% subset] # 抽组
    pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)
    for (i in group) {
      DEG <- FindMarkers(pbmc, ident.1 = i, ident.2 = "Ad1R", genes.use = NULL, group.by = "group",
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
      
      pathway <- pathway_liver
      for (j in pathway) {
        setwd("E:/B_group/pathway/")
        # 使用 tryCatch 捕获错误
        tryCatch({
          rno <- pathview::pathview(
            gene.data  = genelistDEGs,
            pathway.id = j,
            species    = "rno",
            out.suffix = paste("_", i, "_hep1_2", sep = ""),
            limit = list(gene = as.integer(max(abs(genelistDEGs))))
          )
        }, error = function(e) {
          cat("Error in processing pathway", j, ":", conditionMessage(e), "\n")
        })
      }
    }
  }
}
