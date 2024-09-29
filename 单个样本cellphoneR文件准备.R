library(Seurat)
library(biomaRt)#小鼠同源基因转换包
library(dplyr)
options(stringsAsFactors = F)

# 设置大鼠的Mart
rat_mart <- useMart(biomart = "ensembl", dataset = "rnorvegicus_gene_ensembl", host = "nov2020.archive.ensembl.org")
# 设置人类的Mart
human_mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "nov2020.archive.ensembl.org")

#在E:/B组数据备份(4.29)/单细胞结果/亚群拟时序分析/这个文件夹中创建一个叫k的文件夹
path <- paste0("E:/B_group/cellphoneDB准备文件备份/")
dir.create(path)
filepath <- paste("E:/B_group/cellphoneDB准备文件备份/",i,".RDS", sep = "")
exp <- pbmc
setwd(path)
write.table(as.matrix(exp@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)

#准备好注释过的Seurat对象（sce），提取需要的信息
meta_data <- cbind(rownames(exp@meta.data), exp@meta.data[,'seurat_clusters', drop=F])
meta_data <- as.matrix(meta_data) 
meta_data[is.na(meta_data)] = "Unkown"  #去除NA值
meta_data[,1] <- gsub("-",".",meta_data[,1])#后面读取counts.txt文件的时候-被替换成.了，metadata的却没有替换，因此在这里提前替换一下
setwd(path)
sce <- read.table('cellphonedb_count.txt',header = TRUE,row.names = NULL)
names(sce)[1] <-"RGD.symbol"#方便后期的inner join

#转换基因源
rat_genes <- rownames(exp)

gene_conversion <- getLDS(attributes = c("rgd_symbol"), 
                          filters = "rgd_symbol", 
                          values = rat_genes, 
                          mart = rat_mart, 
                          attributesL = c("hgnc_symbol","chromosome_name","start_position"), 
                          martL = human_mart,uniqueRows = T)
gene_conversion <- gene_conversion[,1:2]
gene_conversion <- gene_conversion[!gene_conversion$HGNC.symbol=="",]
hexp <- gene_conversion %>%  inner_join(sce, by = "RGD.symbol")
hexp <- hexp %>% distinct(HGNC.symbol,.keep_all = T)
row.names(hexp) <- hexp[,2]
hexp<-hexp[,-c(1,2)]
#保存文件
setwd(path)
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F) 
write.table(hexp, 'cellphonedb_gene_reverse_count.txt', sep='\t', quote=F)