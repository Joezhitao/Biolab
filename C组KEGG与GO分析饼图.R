library(dittoSeq)
library(SeuratData)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(org.Rn.eg.db)
library(SingleR)
library(celldex)
library(Seurat)
library(pheatmap)
library(tidyverse)
library(patchwork)
library(Seurat)
library(pheatmap)
library(xlsx)
library(ggrepel)
library(EnhancedVolcano)
library(aPEAR)

#鉴定后数据
filepath <- paste("E:/B组数据备份(4.29)/鉴定后总样本备份/",'c_group',".RDS", sep = "")
pbmc <- readRDS(filepath)
backup <- pbmc
pbmc <- backup

#提取肝细胞
cluster <- c("Hepatocytes")
pbmc = pbmc[, pbmc@meta.data$seurat_clusters %in% cluster]#抽组
pbmc@meta.data$seurat_clusters <- droplevels(pbmc@meta.data$seurat_clusters)

##抽组，每一天左右肝对比
levels(pbmc@meta.data$group)
group = c("Cd15L","Cd15R")
pbmc = pbmc[, pbmc@meta.data$group %in% group]#抽组
pbmc@meta.data$group <- droplevels(pbmc@meta.data$group)

#高变基因分析
group <- levels(pbmc@meta.data$group)# 获取 group 列的所有水平
group <- setdiff(group, c("Cd15R"))
combined_table <- data.frame() # 在循环外部创建一个空的数据框

for (j in group) {
  cell_pbmc <- FindMarkers(pbmc, ident.1 = j, ident.2 = "Cd15R", group.by = "group",
                           logfc.threshold = 0.25, min.pct = 0.1)
  cell_pbmc <- cell_pbmc %>% 
    rownames_to_column(var = "gene")
  
  cell_pbmc$group <- j # 将当前的 j 值赋给 cell_pbmc 的 group 列
  
  combined_table <- bind_rows(combined_table, cell_pbmc) # 将当前子集的数据添加到 combined_table 中
}

write.xlsx(combined_table, file = paste(gofile, "combined.xlsx", sep = "")
           , row.names = F)

#分组
#数据如此格式即可，其他的数据整理成此格式即可
#选择上调基因
df_sig <- combined_table
df_sig <- subset(df_sig,p_val_adj<0.05&(avg_log2FC > 0.1))


group <- data.frame(gene=df_sig$gene,group=df_sig$group)#分组情况

#gene转化为ID
Gene_ID <- bitr(df_sig$gene, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Rn.eg.db")

#构建文件并分析
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')

#KEGG
data_KEGG <- compareCluster(ENTREZID~group,
                            data=data,
                            fun = "enrichKEGG",#函数选择什么定义什么分析
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05,
                            organism= "rno")#物种

data_KEGG <- pairwise_termsim(data_KEGG)#要做分组图，需要先运行这个函数
png("E:/B组数据备份(4.29)/go分析通路信息/KEGG_Cd15L.png", width = 8, height = 8, units = "in", res = 800)
emapplot(data_KEGG, pie="count", cex_category=5, layout="kk")
dev.off()

#for循环跑
sample = c("Cd15L","Cd15R")
for (i in sample) {
  sample_b <- setdiff(sample, i)
  combined_table <- data.frame() # 在循环外部创建一个空的数据框
  
  for (j in sample_b) {
    cell_pbmc <- FindMarkers(pbmc, ident.1 = i, ident.2 = j, group.by = "group",
                             logfc.threshold = 0.25, min.pct = 0.1)
    cell_pbmc <- cell_pbmc %>% 
      rownames_to_column(var = "gene")
    
    cell_pbmc$group <- i # 将当前的 i 值赋给 cell_pbmc 的 group 列
    
    combined_table <- bind_rows(combined_table, cell_pbmc) # 将当前子集的数据添加到 combined_table 中
  }
  
  #write.xlsx(combined_table, file = paste(gofile, "combined.xlsx", sep = "") , row.names = F)
  
  gofile <- paste("E:/B组数据备份(4.29)/单细胞结果/C组Go/",'C_group', sep = "")
  #分组
  #数据如此格式即可，其他的数据整理成此格式即可
  df_sig <- combined_table
  df_sig <- subset(df_sig,p_val_adj<0.05&(avg_log2FC > 0.1))
  
  
  group <- data.frame(gene=df_sig$gene,group=df_sig$group)#分组情况
  
  #gene转化为ID
  Gene_ID <- bitr(df_sig$gene, fromType="SYMBOL", 
                  toType="ENTREZID", 
                  OrgDb="org.Rn.eg.db")
  
  #构建文件并分析
  data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
  
  #KEGG
  data_KEGG <- compareCluster(ENTREZID~group,
                              data=data,
                              fun = "enrichKEGG",#函数选择什么定义什么分析
                              pAdjustMethod = "BH",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05,
                              organism= "rno")#物种
  
  data_KEGG <- pairwise_termsim(data_KEGG)#要做分组图，需要先运行这个函数
  
  keggplot <- emapplot(data_KEGG, pie="count", cex_category=5, layout="kk")
  png(paste0(gofile,"_KEGG_",i,".png",sep = ""), width = 18, height = 18, units = "in", res = 800)
  print(keggplot)
  dev.off()
  
  data_GO <- compareCluster(ENTREZID~group, 
                            data=data, 
                            fun="enrichGO", 
                            OrgDb="org.Rn.eg.db",
                            ont = "ALL",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
  #结果可视化GO
  data_GO <- pairwise_termsim(data_GO)#要做分组图，需要先运行这个函数
  goplot <- emapplot(data_GO, pie="count", cex_category=5, layout="kk",legend_n=2)+
    scale_fill_manual(values = dittoColors())#修改填充颜色
  png(paste0(gofile,"_GO_",i,".png",sep = ""), width = 18, height = 18, units = "in", res = 800)
  print(goplot)
  dev.off()
}
