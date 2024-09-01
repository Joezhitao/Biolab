
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("maftools")

#install.packages("dplyr")

#加载R包
library(maftools)
library(dplyr)
library(mclust)
library(NMF)
library(pheatmap)
library(barplot3d)
library(openxlsx)
library(TRONCO)
BiocManager::install("TRONCO")
BiocManager::install("PoisonAlien/maftools")
# 替换下面的URL和版本号为正确的仓库URL和你想安装的版本号
remotes::install_github("BioinformaticsFMRP/maftools@v2.6.05")

install.packages("barplot3d")
#目录
setwd("E:/CRC/genetic mutation/")
#合并所有数据
unnecessary_files <- list.files(pattern = '*.gz.parcel',recursive = TRUE)
files <- list.files(pattern = '*.gz',recursive = TRUE)
#去掉files中包含unnecessary_files的文件
files <- files[!files %in% unnecessary_files]

all_mut <- data.frame()
for (file in files) {
  mut <- read.delim(file,skip = 7, header = T, fill = TRUE,sep = "\t")
  all_mut <- rbind(all_mut,mut)
}
#仅保留前12个字符
all_mut$Tumor_Sample_Barcode = substr(all_mut$Tumor_Sample_Barcode,1,12)

#数据读入
all_mut <- read.maf(all_mut)

#数据整理
a <- all_mut@data %>%
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>%
  as.data.frame() %>%
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))

#提取基因
gene <- as.character(unique(a$Hugo_Symbol))
#提取样本
sample <- as.character(unique(a$Tumor_Sample_Barcode))

#创建data.frame
mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
#将信息填入
for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}

#创建data.frame
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))
#将信息填入
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

#所有样本突变情况汇总/排序
gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))

colnames(gene_count)[1] = "Gene"
colnames(gene_count)[2] = "Num"
write.table(gene_count,'geneMut.txt', sep="\t", quote=F, row.names = F)
write.mafSummary(maf = all_mut,basename = "input")

#绘制瀑布图oncoplot
#missense_mutation:错义突变 frame_shift_del：移码缺失突变 
#nonsense_mutation：无义突变 frame_shift_ins：移码插入突变 
#splice_site：剪接位点 in_frame_ins：inframe插入 
#in_frame_del：inframe缺失 translation_start_site:转录起始位点 
#nonstop_mutation：终止密码子突变 
plotmafSummary(maf = all_mut, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
pdf(file="maf.pdf", width=6, height=6)
oncoplot(maf = all_mut,
         top = 30, #显示前30个的突变基因信息
         fontSize = 0.6, #设置字体大小
         showTumorSampleBarcodes = F) #不显示病人信息
dev.off()
#通路上的突变
pws = pathways(maf = all_mut, plotType = 'treemap')
pdf(file="pathway.pdf", width=6, height=6)
plotPathways(maf = all_mut, pathlist = pws)
dev.off()
OncogenicPathways(maf = all_mut)
PlotOncogenicPathways(maf = all_mut,pathways="RTK−RAS") 

#计算tmb值
tmb_table = tmb(maf = all_mut,logScale = F)
#保留需要tmb值信息
tmb_table = tmb_table[,c(1,3)]
tmb_table = as.data.frame(tmb_table)
tmb_table[,1]=substr(tmb_table[,1],1,12)
colnames(tmb_table)
tmb_table <- aggregate( . ~ Tumor_Sample_Barcode,data=tmb_table, max)
colnames(tmb_table)[1] = "id"
colnames(tmb_table)[2] = "TMB"
write.table(tmb_table,'TMB.txt', sep="\t", quote=F, row.names = F)

#